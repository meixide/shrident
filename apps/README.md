# NSW application with unit-level influence contributions

The published run (`data/lalonde_functionals_hapc.py`) stored one line per estimator and
nothing per unit, so the covariance between the estimators, and with it a valid standard
error for any pooled estimate, cannot be recovered from it. The scripts here store the
per-row predictions of every nuisance fit and assemble everything else afterwards.

| Script | What it does | Where it runs |
|---|---|---|
| `nsw_prepare.py` | creates the run directory: `units.csv` (rows, persons, folds), `tasks.csv`, `settings.json`, `MANIFEST` | login node, seconds |
| `nsw_fit_nuisances.py` | one HAPC fit per task; writes `nuis/b<block>_<nuisance>_f<fold>.csv` and `.json` | SLURM array |
| `nsw_assemble.py` | estimators, person-level influence matrix `Phi`, `Sigma = Phi'Phi`, pooled estimates, intervals, `report.md` | anywhere, seconds |
| `inference.py` | dispersion estimators, weights, sandwich, thresholding, conformal intervals from `(psi, Phi)` | imported |

## Launch on FinisTerrae III

```bash
cd /mnt/lustre/scratch/nlsas/home/csic/dml/cgm/shrident && git pull && bash hpc/run_nsw_if.sh
```

This prepares `results/nsw_if_hapc_K5_s20261012/`, submits 24 tasks for the CPS-2 block
(2 cores, 24 h limit) and 84 for the other blocks (1 core, 2 h limit), and chains the
assembly. The report is `results/nsw_if_hapc_K5_s20261012/assembled/report.md`.
Commit the run directory the day it finishes.

## What is fitted

Fold 0 is the full-sample fit; folds 1 to K train without the persons of that fold and
predict them (cross-fitting). Training rows are stacked as in the published run, because
the CV folds of `cv_hapc` depend on row order.

| Block | Nuisance | Family | Trained on | Published fit |
|---|---|---|---|---|
| RCT | `e_rct` | binomial | all NSW, label treated | IPW-NB |
| RCT | `m_rct` | gaussian | NSW controls | OR-NB |
| external k | `m_pool` | gaussian | NSW controls and sample k | OM-k |
| external k | `pi_s` | binomial | NSW and sample k, label NSW | IPW-k |
| external k | `m_ec` | gaussian | sample k | new |
| external k | `g_trt` | binomial | NSW and sample k, label treated | new |

The two new fits per block are the second nuisance of each doubly robust estimator.

## What is assembled

| Set | Estimators | Unit-level contributions | Goes to |
|---|---|---|---|
| `onestep_xfit` | one-step estimator of each functional (14), cross-fitted | efficient influence function, out-of-fold nuisances | article |
| `onestep_full` | the same with full-sample nuisances | efficient influence function | supplement |
| `plugin_eif` | the 15 published plug-ins | efficient influence function of each functional | supplement |
| `plugin_known` | the same 15 estimates | fitted nuisances treated as known | supplement |

For every weighted estimator the tables give the effective sample size of the weights
(`ess`) and the largest normalised weight (`max_weight`); look at them before reading
the one-step estimates against the external samples.

The unit is the person: CPS-3 is a subsample of CPS-2 and PSID-3 of PSID-2 (all 429 and
128 rows match), so the data hold 3067 persons in 3624 rows.

## Checks

- `nsw_assemble.py` exits 2 if a full-sample fit is missing and 3 if a plug-in estimate
  differs from `data/lalonde_functionals_hapc.csv` by more than `--tol` (default 0.5).
- On macOS with hapc 2.6.0 the pipeline and the original script give the same three
  RCT-block estimates to all printed digits, and OR-NB equals 1719.70 against the
  published 1764.92: the CV curve of that fit is flat and the two machines choose
  adjacent grid points (lambda 3303 here, 6772 on FT3). Whether FT3 reproduces its own
  published value is what exit status 3 reports.
- Pipeline test on a subsample with a coarse grid (numbers are not meaningful):

```bash
python apps/nsw_prepare.py --run /tmp/nsw_smoke --K 2 --frac 0.12 --smoke
seq 0 53 | xargs -P 6 -I{} python apps/nsw_fit_nuisances.py --run /tmp/nsw_smoke --task {}
python apps/nsw_assemble.py --run /tmp/nsw_smoke
```
