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
| Conda env | **`$STORE/conda/envs`** (500 GB) | `$HOME` is only 10 GB and a full `$HOME` crashes services; this account's `~/.condarc` already routes envs/pkgs to `$STORE/conda`. |
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
1. ensures `$STORE/conda/{envs,pkgs}` exist (env placement is governed by this
   account's existing `~/.condarc`, which already routes envs to
   `$STORE/conda/envs`);
2. activates the **system** conda — do **not** install your own. On FT3 plain
   `module load miniconda3` loads `4.11.0` and yields **no working conda**, so
   the script sources the `22.11.1-1` build directly (conda 23.7.4):
   `source /opt/cesga/2020/software/Core/miniconda3/22.11.1-1/etc/profile.d/conda.sh`;
3. creates env **`ft3hapc`** (python 3.11) with `numpy(<2.3)`, `pandas`,
   `scikit-learn`, `scipy`, `statsmodels`;
4. installs the **conda** compiler toolchain and builds `hapc` from
   `git clone https://github.com/meixide/hapc` with it. **This is essential:**
   the gentoo `cesga/2020` gcc links against glibc 2.31, so its `hapc_core.so`
   requires `GLIBC_2.29` — which the system runtime glibc (2.28, RHEL 8.4) the
   conda python uses does **not** provide, giving an `ImportError`. The build
   therefore needs three FT3-specific pieces (all in `setup_env.sh`): pinned
   `gcc_linux-64=12`/`gxx_linux-64=12` (gcc 15's libstdc++ won't compile
   pybind11's `<chrono>` here) on a **`sysroot_linux-64=2.17`** (the default
   2.12 sysroot predates `timespec_get`), with `CXXFLAGS=-D_ISOC11_SOURCE` (the
   2.17 sysroot's pre-C++17 `<features.h>` only declares it under that macro)
   and **`unset LIBRARY_PATH`** (else the gentoo lib64 makes the linker bind
   `exp`/`log`@`GLIBC_2.29`). Result: a `.so` needing only `GLIBC ≤ 2.14`;
5. runs a smoke test (imports `hapc` and runs a tiny `cv_hapc` fit).

> If CESGA retires the `22.11.1-1` build, pick the current one from
> `module avail miniconda3` and update the `source …/conda.sh` path in
> `hpc/*.sh`. Override the env name or hapc URL with
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
- **`conda: command not found` in the job** → on FT3 `module load miniconda3`
  does **not** give a working conda; the scripts instead
  `source /opt/cesga/2020/software/Core/miniconda3/22.11.1-1/etc/profile.d/conda.sh`
  before `conda activate` (the provided scripts do this). If that build is gone,
  pick the current one from `module avail miniconda3` and update the path.
- **`hapc` import error: `/lib64/libm.so.6: version GLIBC_2.29 not found`** →
  the extension was built against the gentoo `cesga/2020` glibc (2.31) but the
  conda python runs on the system glibc (2.28). Rebuild it with the conda
  toolchain exactly as `setup_env.sh` step 4 does — the three things that matter
  are: `sysroot_linux-64=2.17`, `CXXFLAGS=-D_ISOC11_SOURCE`, and
  `unset LIBRARY_PATH` before `pip install --force-reinstall --no-deps .` in the
  hapc clone. Verify with `objdump -T hapc_core*.so | grep -oE 'GLIBC_[0-9.]+' |
  sort -V | tail` showing max ≤ 2.17 (we get 2.14). Also confirm `numpy<2.3`.
  Common cause of a *fresh* failure: the gentoo `gcc_linux-64=12` build hit a
  `'timespec_get' has not been declared` compile error → the 2.17 sysroot/
  `_ISOC11_SOURCE` step was skipped.
- **A `BOUNDARY-low` CV verdict on full data** → lower `GRID['log_lambda_min']`
  in `data/lalonde_functionals_hapc.py` and resubmit just that block
  (`HAPC_BLOCK=<k>`).
- **`logs/` missing** → `mkdir -p logs` before `sbatch` (`run_all.sh` does this).
