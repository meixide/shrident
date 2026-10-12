"""Finite-difference check of the 2SLS influence function in iv_environments_sim.tsls_if
(aider round 3, F1). For observation i, the empirical influence is
    IF_i ~ n * (beta(w + h e_i) - beta(w)) / h,
with beta(w) the weighted 2SLS slope. Compared with tsls_if() for gamma = 0 (moment
conditions hold) and gamma = 2 (over-identifying restrictions fail, first-stage term matters).
Run from the repo root:  python checks/tsls_if_finite_difference.py
"""
import numpy as np
from iv_environments_sim import simulate, tsls_if


def tsls_weighted(z, a, y, subset, w):
    inS = np.isin(z, subset).astype(float) * w
    Zm = np.column_stack([(z == s).astype(float) for s in subset])
    X = np.column_stack([np.ones(len(y)), a])
    G = (Zm * inS[:, None]).T @ X
    M = (Zm * inS[:, None]).T @ Zm
    m = (Zm * inS[:, None]).T @ y
    W = np.linalg.inv(M)
    return np.linalg.solve(G.T @ W @ G, G.T @ W @ m)[1]


rng = np.random.default_rng(3)
pis = np.linspace(0.2, 0.8, 5)
subset = (0, 1, 2, 3, 4)
for gamma in (0.0, 2.0):
    z, a, y = simulate(4000, pis, 1.0, gamma, rng)
    n = len(y)
    _, infl = tsls_if(z, a, y, subset, len(pis))
    w0 = np.ones(n); b0 = tsls_weighted(z, a, y, subset, w0)
    idx = rng.choice(n, 25, replace=False); h = 1e-4
    fd = np.array([n * (tsls_weighted(z, a, y, subset, w0 + h * (np.arange(n) == i)) - b0) / h for i in idx])
    rel = np.max(np.abs(fd - infl[idx])) / np.max(np.abs(fd))
    print(f"gamma={gamma}: max |FD - IF| / max |FD| over 25 observations = {rel:.2e}")
