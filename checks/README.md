# checks/

Numerical checks behind the 2026-09-29 review of the JASA submission
(internal notes: asertii `Multiple/2026-09 JRSSB rejection — review and experiments.md`).

- `limit_experiment.py` — large-sample limit of sqrt(n)(psi_EB - psi*) under exact
  identification: non-Gaussian limit, sandwich coverage, variance vs best single
  estimator vs GLS. Source of Table "Large-sample behaviour" in shrinkage_jasa/main.tex.
- `dispersion_estimators.py` — pairwise vs DerSimonian–Laird tau-hat for NSW and
  Card–Krueger from data/*.csv; also prints the working-independence SE reported
  in arxivv2 for NSW (J=3).

Run with a Python that has numpy/pandas (e.g. ~/miniconda3/bin/python3) from the repo root.
- `iv_environments_sim.py` (+ `.out`) — Example ex:iv_subset: EB vs GLS vs 2SLS on all
  environments vs best single functional, exact and heterogeneous regimes (2026-09-30).
  EB-0 sets the dispersion to zero (Proposition 2(iii) of the article). The threshold rule of
  earlier versions was dropped on 2026-10-12.
- `kolesar_connection.py` (+ `.out`) — Kolesár et al. (JBES 2015) cluster case: the between-
  cluster slope equals the inverse-variance pooled functional of pairwise Walds; the unweighted
  mean of the pairwise discrepancies diverges (Cauchy-like).
- `tsls_if_finite_difference.py` (+ `.out`) — finite-difference check of the 2SLS influence
  function in `iv_environments_sim.tsls_if` (relative error ~2e-7 for gamma = 0 and 2).
  Note: in `iv_environments_sim.py` every environment is observational, including Z=0; this is
  not the randomized-environment design of Example ex:iv_subset.
- `convex_combination.py` (+ `.out`) — the convex (non-negative) minimum-variance combination against
  precision weights, best single estimator and GLS: limit variances, finite samples, IV design.
  Source of the finite-sample table of the supplement.
- `ck_jackknife_crosscheck.R` (+ `.out`) — bootstrap and analytic cross-check of the Card-Krueger
  jackknife standard errors; shows that the published adjusted standard errors were too small.
- `weighted_centring_sim.py` (+ `.out`) — uniform, precision and design weights under four patterns
  of discrepancies (balanced, tied to precision, common sign, centred by design). Source of the
  table that illustrates Theorem 2 and its corollaries in the supplement.
