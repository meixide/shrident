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
  Since 2026-10-12 the thresholded estimator uses the scale-free rule tau2 > log(n) mean(v);
  the earlier rule n tau2 > log n depended on the unit of the outcome and, in this design,
  never set a positive tau2-hat to zero.
- `kolesar_connection.py` (+ `.out`) — Kolesár et al. (JBES 2015) cluster case: the between-
  cluster slope equals the inverse-variance pooled functional of pairwise Walds; the unweighted
  mean of the pairwise discrepancies diverges (Cauchy-like).
- `tsls_if_finite_difference.py` (+ `.out`) — finite-difference check of the 2SLS influence
  function in `iv_environments_sim.tsls_if` (relative error ~2e-7 for gamma = 0 and 2).
  Note: in `iv_environments_sim.py` every environment is observational, including Z=0; this is
  not the randomized-environment design of Example ex:iv_subset.
