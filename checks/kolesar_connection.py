"""Kolesár, Chetty, Friedman, Glaeser & Imbens (JBES 2015), Special Case I, in the language of
pooled identifying functionals.

Clusters k = 1..G act as instruments. With direct effects gamma_k and first stages pi_k, the
reduced forms are E[Y | G=k] = delta + beta pi_k + gamma_k and E[X | G=k] = pi_k, so the Wald
functional of the pair (k, l) identifies

    psi_kl = beta + eps_kl,   eps_kl = (gamma_k - gamma_l) / (pi_k - pi_l).

Claims checked here, at the population level (no sampling noise):
 1. The between-cluster slope (the probability limit of their consistent estimators) equals the
    pooled functional sum w_kl psi_kl with w_kl proportional to (pi_k - pi_l)^2, i.e. the
    inverse-variance weights of the pairwise Walds under equal cluster sizes (tau^2 = 0 weights).
 2. Its discrepancy is Cov_G(gamma, pi) / Var_G(pi): their condition Cov(gamma, pi) = 0 is
    exactly zero discrepancy of that pooled functional.
 3. The unweighted mean of eps_kl (what EB weights tend to when tau^2 > 0 and v -> 0) does not
    settle: eps_kl is a ratio of independent normals when (gamma, pi) are independent normal
    (their Assumption 5), so it is Cauchy-like, unbounded, and has no mean.
Run from the repo root:  python checks/kolesar_connection.py
"""
import numpy as np

rng = np.random.default_rng(0)
beta = 1.0
print(f"{'G':>6s} {'pooled IVW':>11s} {'cov/var':>9s} {'unweighted mean eps':>20s} {'median eps':>11s} {'max|eps|':>10s}")
for G in (20, 100, 500, 2000):
    gamma = rng.normal(0, 1, G)             # direct effects, independent of pi
    pi = rng.normal(0, 1, G)                # first-stage effects
    k, l = np.triu_indices(G, 1)
    eps = (gamma[k] - gamma[l]) / (pi[k] - pi[l])
    w = (pi[k] - pi[l]) ** 2
    pooled = beta + (w * eps).sum() / w.sum()
    slope = beta + np.cov(gamma, pi)[0, 1] / np.var(pi, ddof=1)
    print(f"{G:6d} {pooled:11.4f} {slope:9.4f} {eps.mean():20.4f} {np.median(eps):11.4f} {np.abs(eps).max():10.1f}")
