# Running the full-data HAPC LaLonde analysis on CESGA FinisTerrae III (FT3)

This file is meant to be read **after cloning the repo onto FT3**. It builds the
environment and submits the full-data version of the HAPC re-analysis described
in [lalonde_hapc_report.md](lalonde_hapc_report.md) (which was a 25 % subsample
dry-run, because degree-2 HAPC on the full CPS-2 pool is memory-bound on a
laptop). On FT3 we run it on the **whole sample** (`HAPC_FRAC=1.0`).

FT3 reference: <https://cesga-docs.gitlab.io/ft3-user-guide/index.html>

---

## 0. What gets computed

`data/lalonde_functionals_hapc.py` estimates the 15 identifying functionals with
`hapc.cv_hapc` base learners (`norm='1'`, `max_degree=2`, `npcs=floor(n^0.6)`,
gaussian outcomes / binomial scores, `approx=False`, CV grid
`log_lambda∈[-12,16]`). Only **10** of the 15 require a base-learner fit; the two
**CPS-2** fits (`OM-CPS-2`, n≈2 629; `IPW-CPS-2`, n≈2 814) dominate the cost.

Outputs (written to `data/`, then copied to `$STORE/shrinkage_results/`):
`lalonde_functionals_hapc.csv`, `lalonde_eb_hapc.csv`, `lalonde_loo_hapc.csv`,
`lalonde_pi_hapc.csv`, `lalonde_functionals_ols_sub.csv` (OLS/logistic baseline
on the same full data).

The original LaTeX (`lalonde_application.tex`) is **never modified**.

---

## 1. Filesystems — where to put things

| What | Filesystem | Why |
|---|---|---|
| The cloned repo + run | **`$LUSTRE`** (1 TB, high-speed, no backup) | FT3 docs: compute jobs run from `$LUSTRE`. |
| Conda env | **`$STORE/.conda`** (500 GB) | `$HOME` is only 10 GB and a full `$HOME` crashes services. |
| Final result CSVs | `$STORE/shrinkage_results/` | `$STORE` has snapshots; the aggregate job copies them there. |

```bash
cd $LUSTRE
git clone <this-repo-url> shrinkage
cd shrinkage
```

Check quotas anytime with `myquota`.

---

## 2. One-time environment setup

```bash
bash hpc/setup_env.sh          # run on a LOGIN node (it compiles hapc); not via sbatch
```

This (see [hpc/setup_env.sh](hpc/setup_env.sh)):
1. redirects conda storage to `$STORE/.conda` (per the FT3 guide);
2. `module load miniconda3` (the **system** miniconda — do **not** install your
   own) plus a compiler + `cmake` to build `hapc`'s C++/pybind11 extension;
3. creates env **`ft3hapc`** with `numpy(<2.3)`, `pandas`, `scikit-learn`,
   `scipy`, `statsmodels`;
4. `git clone https://github.com/meixide/hapc` and `pip install .`;
5. runs a smoke test.

> Verify the exact module name first if the load fails: `module spider miniconda3`
> (and `module spider gcc cmake`). Override the env name or hapc URL with
> `ENV_NAME=... HAPC_REPO=... bash hpc/setup_env.sh`.

---

## 3. Submit the run

### Option A — job array (recommended)

```bash
mkdir -p logs
bash hpc/run_all.sh
```

This submits a **5-task array** ([hpc/submit_array.sh](hpc/submit_array.sh)) —
one task per block (0=RCT, 1=CPS-2, 2=CPS-3, 3=PSID-2, 4=PSID-3) — and an
**aggregation job** ([hpc/submit_aggregate.sh](hpc/submit_aggregate.sh)) chained
with `--dependency=afterok`, so aggregation (EB + LOO + conformal + OLS baseline)
runs automatically once all blocks succeed. Blocks run **concurrently**, so wall
time ≈ the single heaviest block (CPS-2), not the sum.

### Option B — single job (simplest)

```bash
mkdir -p logs
sbatch hpc/submit_single.sh
```

Everything in one process, sequentially ([hpc/submit_single.sh](hpc/submit_single.sh)).
Robust, no orchestration; wall time ≈ sum of all fits.

---

## 4. Why this is the optimal submission

- **Parallelise across the independent unit.** The 5 blocks share no computation,
  so a job array runs them on separate nodes at once. Wall time collapses from
  "sum of all fits" to "the CPS-2 block alone". The aggregation (pooling, EB, LOO,
  conformal) is cheap and depends on all blocks → one dependent job.
- **Parallelise *within* a fit too.** The HAPC cost is the dense degree-2 HAL
  kernel build + eigendecomposition, which is BLAS-heavy. Each task gets
  `-c 16` cores and the scripts export `OMP/OPENBLAS/MKL_NUM_THREADS` to match,
  so numpy threads the linear algebra. (sklearn's LASSO coordinate descent for
  the binomial scores is single-threaded; the kernel/eigendecomp still benefit.)
- **Size memory for the worst case.** On a laptop the CPS-2 binomial fit
  swap-thrashed (the dense HAL design is ~2 GB and peak is a few × that). Array
  tasks request **64 GB** (`-c 16 × --mem-per-cpu=4G`); the single job requests
  **128 GB**. FT3 standard nodes have 180–256 GB, so this never swaps.
- **Partitions / walltime.** `short` = 6 h (default), `medium` = 3 days,
  `long` = 7 days. The CPS-2 block is the long pole (two large fits, ≈1–2 h on
  full data); the array requests `medium` + `-t 12:00:00` and the small blocks
  finish in minutes. The single job requests `medium` + `-t 1-00:00:00`.
- **Grouping note.** Block 1 (CPS-2) carries both heavy fits, so it sets the wall
  time. If you later need it shorter, split CPS-2's `OM` and `IPW` into separate
  array tasks; for a one-off run the 5-block split is the clean sweet spot.

---

## 5. Monitor & inspect

```bash
squeue -u $USER            # or: sqstat            -- queued/running jobs
sacct -j <jobid>           # accounting / exit codes
seff <ARRAYID>_1           # efficiency of the CPS-2 task (the long pole)
tail -f logs/hapc_<ARRAYID>_1.out
```

Each fit prints a CV verdict line, e.g.
`[cv] OM-CPS-2 gaussian best_lambda=... -> INTERIOR`. Watch for
`*** NOT INTERIOR ***`: `BOUNDARY-low` (widen `log_lambda_min`), `BOUNDARY-high` /
`PLATEAU-high(null model)` (genuine null — common only at small n; expected to be
rare on full data).

---

## 6. Configuration knobs (environment variables)

| Var | Default | Meaning |
|---|---|---|
| `HAPC_FRAC` | `0.25` | Subsample fraction; **the HPC scripts set `1.0`** (full data). |
| `HAPC_SEED` | `42` | Subsample RNG seed. |
| `HAPC_BLOCK` | `all` | `0..4` ⇒ compute one block only (set by the array via `SLURM_ARRAY_TASK_ID`). |
| `HAPC_OUTDIR` | script dir | Where CSVs are written. |
| `ENV_NAME` | `ft3hapc` | Conda env name used by the job scripts. |

Locally you can reproduce the dry-run with `HAPC_FRAC=0.25 python data/lalonde_functionals_hapc.py`.

---

## 7. Troubleshooting

- **OOM / `oom-kill`** → raise `--mem-per-cpu` (e.g. `8G`) in the submit script;
  inspect actual use with `seff <jobid>`.
- **`conda: command not found` in the job** → ensure the script has
  `module load miniconda3` then `eval "$(conda shell.bash hook)"` before
  `conda activate` (the provided scripts do this).
- **`hapc` import error** → rebuild: it needs a compiler + `cmake` at build time
  (`module load gcc cmake`); confirm `numpy<2.3`.
- **A `BOUNDARY-low` CV verdict on full data** → lower `GRID['log_lambda_min']`
  in `data/lalonde_functionals_hapc.py` and resubmit just that block
  (`HAPC_BLOCK=<k>`).
- **`logs/` missing** → `mkdir -p logs` before `sbatch` (`run_all.sh` does this).
