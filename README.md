# shrident

Replication code for the applications of

> García Meixide, C. and Ríos Insua, D. *Shrinkage through multiple identifiability.*

The article pools several estimators of one causal effect, each obtained from a different
identification strategy on the same sample. This repository reproduces its two
applications:

- **National Supported Work (NSW)**, Sections 8.1 and 8.2 and Table 3(a): fourteen
  estimators of the effect of a job-training program, built from the randomized trial and
  from four survey comparison samples, together with a catalogue of 512 doubly robust
  specifications.
- **New Jersey minimum wage (Card and Krueger, 1994)**, Section 8.3 and Table 3(b):
  twenty-five difference-in-differences estimators that use different Pennsylvania control
  groups.

The supplementary material of the article gives the construction of every estimator.

## Layout

| Directory | Contents |
|---|---|
| `apps/` | Current pipeline for both applications. `apps/README.md` documents each script. |
| `data/` | Raw data, the script that builds the clean data sets, and the reference estimates that the pipeline checks against. |
| `results/` | One directory per run, each with a `MANIFEST` that records settings, seeds and software versions. Every number in Section 8 is read from these files. |
| `hpc/` | SLURM scripts for the nuisance fits of the NSW application. They were written for the CESGA FinisTerrae III cluster. |
| `checks/` | Numerical checks and simulations behind the supplementary material. `checks/README.md` lists them. |

## Requirements

- **Python 3.11** with numpy, pandas and scipy. The NSW nuisance fits also need
  [`hapc`](https://github.com/meixide/hapc), highly adaptive principal component regression;
  the reported run used hapc 2.3.0.
- **R 4.3** with haven, dplyr, ggplot2 and sandwich.

## Reproducing the article

Run every command from the repository root.

### 1. Clean data sets

```bash
Rscript data/build_datasets.R
```

This writes `data/lalonde_combined.csv` and `data/ck94.csv` from the raw files.

### 2. New Jersey minimum wage (Section 8.3, Table 3(b))

This step runs on a laptop in seconds.

```bash
Rscript apps/ck_jackknife.R results/ck94_jackknife_stratified
python apps/ck_intervals.py --out results/ck94_intervals_jackknife_B2000_s42 \
    --jackknife results/ck94_jackknife_stratified
Rscript apps/ck_forest.R results/ck94_intervals_jackknife_B2000_s42
```

- `ck_jackknife.R` recomputes the 25 estimates and writes the restaurant-level influence
  contributions, from a delete-one jackknife stratified by state.
- `ck_intervals.py` writes the pooled estimates, sandwich standard errors, contrasts and
  convex weights to `pooled.csv`, `contrasts.csv`, `weights.csv` and `estimates.csv`. It
  uses seed 42.

### 3. NSW catalogue of formulas (Section 8.2)

This step runs on a laptop in seconds.

```bash
python apps/nsw_catalogue.py --out results/nsw_catalogue_dr_parametric_J30_s20261012
```

The script writes:

- `catalogue.csv`: the 512 estimates.
- `draw.csv`: the 31 formulas drawn with seed 20261012.
- `interval.csv`: the split conformal interval and its coverage over 20,000 further draws.

### 4. NSW pooling with cross-fitted nuisances (Section 8.1, Table 3(a))

The nuisance functions are 108 HAPC fits, with five-fold cross-fitting and fold seed
20261012. Together they take about 65 CPU-hours, and the longest fit takes about 8.5 hours,
so they run on a cluster:

```bash
bash hpc/setup_env.sh              # once, on a login node
bash hpc/run_nsw_if.sh 5 20261012  # prepares the run, submits the array, chains the assembly
```

The scripts in `hpc/` assume the FinisTerrae III file systems (`$STORE`, `$LUSTRE`) and its
system conda; adapt those lines for another cluster.

The stored predictions of every fit are in `results/nsw_if_hapc_K5_s20261012/nuis/`. You
can therefore reproduce the article's numbers without refitting:

```bash
python apps/nsw_assemble.py --run results/nsw_if_hapc_K5_s20261012
python apps/nsw_article_numbers.py --run results/nsw_if_hapc_K5_s20261012
Rscript apps/nsw_forest.R results/nsw_if_hapc_K5_s20261012
```

The table rows of the article are the set `onestep_xfit` in `assembled/article/pooled.csv`.
The pairwise contrasts are in `contrasts.csv`, and the figure of the supplement is
`nsw_forest.pdf`.

Influence contributions are indexed by person. CPS-3 and PSID-3 are subsamples of CPS-2 and
PSID-2, so the 3,624 rows of the data correspond to 3,067 distinct persons.

### Numerical reproducibility

HAPC chooses its penalty by cross-validation. When the cross-validation curve is flat, two
machines can select adjacent grid points. For example, the trial outcome regression gives
1719.70 on macOS and 1764.92 on the cluster. Refitting on another machine can therefore
change the last digits of some NSW estimates. Assembling from the stored predictions
reproduces the article exactly.

## Data

The raw data are redistributed here for convenience. Please cite the original sources:

- **NSW and comparison samples**, `data/lalonde/`: Dehejia, R. H. and Wahba, S. (1999),
  *Journal of the American Statistical Association* 94, 1053–1062; and LaLonde, R. J. (1986),
  *American Economic Review* 76, 604–620. Source:
  <https://users.nber.org/~rdehejia/nswdata2.html>.
- **New Jersey and Pennsylvania fast-food survey**, `data/card_krueger/`: Card, D. and
  Krueger, A. B. (1994), *American Economic Review* 84, 772–793. The original codebook and
  `read.me` are included. Source: <https://davidcard.berkeley.edu/data_sets.html>.

## Citation

If you use this code, please cite the article above.
