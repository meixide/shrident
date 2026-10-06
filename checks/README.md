# checks/

Numerical checks behind the 2026-09-29 review of the JASA submission
(internal notes: asertii `Multiple/2026-09 JASA rejection — review and experiments.md`).

- `limit_experiment.py` — large-sample limit of sqrt(n)(psi_EB - psi*) under exact
  identification: non-Gaussian limit, sandwich coverage, variance vs best single
  estimator vs GLS. Source of Table "Large-sample behaviour" in shrinkage_jasa/main.tex.
- `dispersion_estimators.py` — pairwise vs DerSimonian–Laird tau-hat for NSW and
  Card–Krueger from data/*.csv; also prints the working-independence SE reported
  in arxivv2 for NSW (J=3).

Run with a Python that has numpy/pandas (e.g. ~/miniconda3/bin/python3) from the repo root.
- `iv_environments_sim.py` (+ `.out`) — Example ex:iv_subset: EB vs GLS vs 2SLS on all
  environments vs best single functional, exact and heterogeneous regimes (2026-09-30).
- `kolesar_connection.py` (+ `.out`) — Kolesár et al. (JBES 2015) cluster case: the between-
  cluster slope equals the inverse-variance pooled functional of pairwise Walds; the unweighted
  mean of the pairwise discrepancies diverges (Cauchy-like).
