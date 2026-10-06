"""Is GLS / GMM worth it? Simulation in the design of Example ex:iv_subset.

q environments Z in {0..q-1} with treatment propensities pi_z; unmeasured confounding.
Potential treatments A(z) = 1{V < pi_z} (monotone), effects Delta = beta + gamma (V - 1/2),
Y(0) = c V + e. For every subset S of environments with |S| >= 2 (J = 2^q - q - 1) the
2SLS estimator of the effect using only units with Z in S is one identifying functional.

  gamma = 0  -> exact identification: every subset targets beta = ATE.
  gamma != 0 -> each subset targets a LATE; pairs target beta + gamma((pi_z+pi_z')/2 - 1/2).
               Symmetric pi around 1/2 -> discrepancies roughly centred; asymmetric -> not.

Estimators of the ATE beta:
  EB      working-independence weights w_j = 1/(v_j + tau2), pairwise tau2   (the paper)
  EB-thr  same, tau2 set to 0 when n tau2 <= log n                           (Prop. 3(iii))
  GLS     weights Omega^+ 1 / 1'Omega^+ 1 (pseudo-inverse: Omega is near singular)
  2SLS    all q environments = the subset S = {0..q-1}, itself one of the J functionals
  best    the functional with the smallest estimated variance
90% intervals: EB and EB-thr with the influence-function sandwich; GLS with 1/(n 1'Omega^+ 1);
2SLS with its robust influence-function SE.

Run from the repo root:  python checks/iv_environments_sim.py [--reps 400]
"""
import argparse
import itertools
import numpy as np

Z90 = 1.6448536269514722


def simulate(n, pis, beta, gamma, rng, c=2.0):
    q = len(pis)
    z = rng.integers(0, q, n)
    v = rng.uniform(size=n)
    a = (v < pis[z]).astype(float)
    delta = beta + gamma * (v - 0.5)
    y = c * v + rng.normal(size=n) + a * delta
    return z, a, y


def tsls_if(z, a, y, subset, q):
    """2SLS on units with Z in subset, instruments = one-hot of Z within subset.
    Returns (slope estimate, n-vector influence function over the full sample).

    The influence function is valid also when the over-identifying restrictions fail
    (heterogeneous effects): beta solves G'W r(beta) = 0 with G = E[1_S Z X'], W = E[1_S Z Z']^{-1},
    r(beta) = E[1_S Z (Y - X'beta)]. Implicit differentiation gives
        dbeta = (G'WG)^{-1} [ dG' W rbar + G' dW rbar + G'W (g_i - rbar) ],
    with dG = 1_S Z X' - G, dW = -W (1_S Z Z' - W^{-1}) W, g_i = 1_S Z (Y - X'beta), rbar = E[g].
    The first two terms (first-stage estimation) vanish only if rbar = 0 (review F1, 2026-09-30)."""
    n = len(y)
    inS = np.isin(z, subset).astype(float)
    Zm = np.zeros((n, len(subset)))
    for k, s in enumerate(subset):
        Zm[:, k] = (z == s)
    Zs = Zm * inS[:, None]
    X = np.column_stack([np.ones(n), a])
    G = Zs.T @ X / n
    M2 = Zs.T @ Zm / n
    W = np.linalg.inv(M2)
    H = G.T @ W @ G
    b = np.linalg.solve(H, G.T @ W @ (Zs.T @ y / n))
    e = y - X @ b
    g = Zs * e[:, None]                                   # n x K moment contributions
    rbar = g.mean(axis=0)
    Wr = W @ rbar
    t1 = (Zs @ Wr)[:, None] * X - (G.T @ Wr)[None, :]     # dG' W rbar, per unit (n x 2)
    t2 = -((Zs @ (W @ G)) * (Zm @ Wr)[:, None] - (G.T @ W @ M2 @ Wr)[None, :])  # G' dW rbar
    t3 = (g - rbar) @ (W @ G)                             # G'W (g_i - rbar)
    infl = np.linalg.solve(H, (t1 + t2 + t3).T).T
    return b[1], infl[:, 1]


def targets(pis, beta, gamma, subsets, n_mc=400000, seed=0):
    """Population target of each subset 2SLS, by one large simulated sample."""
    rng = np.random.default_rng(seed)
    z, a, y = simulate(n_mc, pis, beta, gamma, rng)
    return np.array([tsls_if(z, a, y, S, len(pis))[0] for S in subsets])


def pairwise_tau2(psi, v):
    i, j = np.triu_indices(len(psi), 1)
    return max(0.0, np.mean(((psi[i] - psi[j]) ** 2 - (v[i] + v[j])) / 2))


def one_rep(n, pis, beta, gamma, subsets, rng):
    q = len(pis)
    z, a, y = simulate(n, pis, beta, gamma, rng)
    est, D = zip(*(tsls_if(z, a, y, S, q) for S in subsets))
    psi = np.array(est); D = np.column_stack(D)            # n x J
    Om = D.T @ D / n                                       # Omega-hat
    v = np.diag(Om) / n
    out = {}
    t2 = pairwise_tau2(psi, v)
    for name, tau2 in (("EB", t2), ("EB-thr", t2 if n * t2 > np.log(n) else 0.0)):
        w = 1 / (v + tau2); w /= w.sum()
        out[name] = (w @ psi, np.sqrt(w @ Om @ w / n))
    ev = np.linalg.eigvalsh(Om)
    diag = {"rank": int((ev > 1e-8 * ev.max()).sum()), "fired": float(n * t2 > np.log(n))}
    for rc, name in ((1e-8, "GLS"), (1e-4, "GLS-1e-4")):
        g = np.linalg.pinv(Om, rcond=rc) @ np.ones(len(psi))
        diag["gsum" if rc == 1e-8 else "gsum4"] = g.sum()
        out[name] = (g @ psi / g.sum(), np.sqrt(1 / (n * g.sum())) if g.sum() > 0 else np.nan)
    full = subsets.index(tuple(range(q)))
    out["2SLS"] = (psi[full], np.sqrt(v[full]))
    k = int(np.argmin(v))
    out["best"] = (psi[k], np.sqrt(v[k]))
    return out, t2, diag


def run(n, q, pis, beta, gamma, reps, seed):
    subsets = [S for r in range(2, q + 1) for S in itertools.combinations(range(q), r)]
    tg = targets(pis, beta, gamma, subsets)
    rng = np.random.default_rng(seed)
    res = {m: [] for m in ("EB", "EB-thr", "GLS", "GLS-1e-4", "2SLS", "best")}
    t2s, diags = [], []
    for _ in range(reps):
        o, t2, d = one_rep(n, pis, beta, gamma, subsets, rng)
        t2s.append(t2); diags.append(d)
        for m, (e, se) in o.items():
            res[m].append((e, se))
    print(f"\nq={q} J={len(subsets)} n={n} pis={np.round(pis,2).tolist()} gamma={gamma} "
          f"| targets: mean={tg.mean():.3f} sd={tg.std():.3f} (ATE={beta}) | mean tau2-hat={np.mean(t2s):.4f}")
    D = {k: np.array([d[k] for d in diags]) for k in diags[0]}
    print(f"diagnostics: Omega-hat rank median {np.median(D['rank']):.0f} of {len(subsets)}; "
          f"threshold n*tau2 > log n fired in {100*D['fired'].mean():.1f}% of reps; "
          f"1'Omega^+1: median {np.median(D['gsum']):.3g} (rcond 1e-8) vs {np.median(D['gsum4']):.3g} (1e-4)")
    print(f"coverage below is for the ATE beta; for gamma != 0 the subset targets differ from beta")
    print(f"{'method':8s} {'bias':>8s} {'sd':>8s} {'rmse':>8s} {'mean se':>8s} {'cov90(ATE)':>10s}")
    for m, r in res.items():
        r = np.array(r); err = r[:, 0] - beta
        cov = np.mean(np.abs(err) <= Z90 * r[:, 1])
        print(f"{m:8s} {err.mean():8.4f} {r[:,0].std():8.4f} {np.sqrt(np.mean(err**2)):8.4f} "
              f"{r[:,1].mean():8.4f} {cov:10.3f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=400)
    ap.add_argument("--n", type=int, default=2000)
    a = ap.parse_args()
    sym5 = np.linspace(0.2, 0.8, 5)
    asym5 = np.linspace(0.1, 0.5, 5)
    for pis, gamma in ((sym5, 0.0), (sym5, 2.0), (asym5, 2.0)):
        run(a.n, len(pis), pis, beta=1.0, gamma=gamma, reps=a.reps, seed=1)
    run(a.n, 4, np.linspace(0.2, 0.8, 4), beta=1.0, gamma=0.0, reps=a.reps, seed=2)
