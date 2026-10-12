"""Convex (non-negative) minimum-variance combination: is a stable GLS worth it?

Weights  w = argmin { w' S w : w >= 0, 1'w = 1 },  S = Omega-hat  (convex averaging of
Lavancier & Rochet 2016, Section 2.3). The simplex contains every vertex and the
diagonal precision weights, so the population minimum is no larger than the variance
of the best single estimator or of the working-independence combination. All
estimators share the target, so sqrt(n)(w-hat'psi-hat - psi*) = w-hat' Z_n and the
limit is N(0, min) whenever w-hat converges; no condition on the boundary is needed.

Part 1  large-sample variances in the design of the article's limit table.
Part 2  finite samples in the same design: rmse and coverage of the 90% sandwich
        interval for convex, GLS, working-independence (EB) and best-single weights. The convex
        weights minimise w'Omega-hat w (dispersion set to zero), as in Proposition 3.
Part 3  the IV-environments design of checks/iv_environments_sim.py (Omega-hat singular).

Run from the repo root:  python checks/convex_combination.py
"""
import itertools
import os
import sys

import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "apps"))
import iv_environments_sim as IV
from inference import convex_weights

Z90 = 1.6448536269514722


def pairwise_tau2(psi, v):
    i, j = np.triu_indices(len(psi), 1)
    return max(0.0, np.mean(((psi[i] - psi[j]) ** 2 - (v[i] + v[j])) / 2))


def design(J, rho):
    sig = np.linspace(0.6, 2.0, J)
    return rho * np.outer(sig, sig) + (1 - rho) * np.diag(sig ** 2)


def part1():
    print("PART 1  large-sample variance of sqrt(n)(estimator - psi*), design of the limit table")
    print(f"{'J':>3s} {'rho':>4s} {'diag(0)':>8s} {'best':>8s} {'convex':>8s} {'GLS':>8s} {'support':>8s} {'min GLS w':>10s}")
    for J, rho in [(5, 0.0), (5, 0.5), (5, 0.9), (10, 0.5)]:
        Om = design(J, rho)
        d = 1 / np.diag(Om); d /= d.sum()
        w = convex_weights(Om)
        g = np.linalg.solve(Om, np.ones(J)); g /= g.sum()
        print(f"{J:3d} {rho:4.1f} {d @ Om @ d:8.3f} {np.diag(Om).min():8.3f} {w @ Om @ w:8.3f} "
              f"{g @ Om @ g:8.3f} {int((w > 1e-6).sum()):8d} {g.min():10.3f}")


def part2(reps=2000, seed=3):
    print("\nPART 2  finite samples, same design; psi* = 0; rmse of sqrt(n) estimator and coverage of the 90% sandwich interval")
    print(f"{'J':>3s} {'rho':>4s} {'n':>5s} | " + " | ".join(f"{m:>15s}" for m in ("EB", "best single", "convex", "GLS")))
    rng = np.random.default_rng(seed)
    for (J, rho), n in itertools.product([(5, 0.5), (5, 0.9), (10, 0.5), (10, 0.9)], (100, 400, 2000)):
        Om = design(J, rho)
        C = np.linalg.cholesky(Om)
        res = {m: [] for m in ("EB", "best single", "convex", "GLS")}
        for _ in range(reps):
            D = rng.standard_normal((n, J)) @ C.T
            psi = D.mean(axis=0)
            Dc = D - psi
            Oh = Dc.T @ Dc / n
            v = np.diag(Oh) / n
            t = pairwise_tau2(psi, v)
            W = {"EB": 1 / (v + t), "best single": np.eye(J)[np.argmin(v)],
                 "convex": convex_weights(Oh),   # dispersion set to zero (Proposition 3)
                 "GLS": np.linalg.solve(Oh, np.ones(J))}
            for m, w in W.items():
                w = w / w.sum()
                res[m].append((np.sqrt(n) * w @ psi, np.sqrt(w @ Oh @ w)))
        cells = []
        for m, r in res.items():
            r = np.array(r)
            cells.append(f"{np.sqrt(np.mean(r[:, 0] ** 2)):7.3f} {np.mean(np.abs(r[:, 0]) <= Z90 * r[:, 1]):7.3f}")
        print(f"{J:3d} {rho:4.1f} {n:5d} | " + " | ".join(f"{c:>15s}" for c in cells))


def part3(reps=400, n=2000):
    print("\nPART 3  IV environments (2SLS on every subset of environments), n=2000; error for the ATE beta = 1")
    print(f"{'design':28s} {'method':8s} {'bias':>8s} {'sd':>8s} {'rmse':>8s} {'mean se':>8s} {'cov90':>6s} {'support':>8s}")
    for label, pis, gamma in [("q=5 symmetric, gamma=0", np.linspace(0.2, 0.8, 5), 0.0),
                              ("q=5 symmetric, gamma=2", np.linspace(0.2, 0.8, 5), 2.0)]:
        q = len(pis)
        subsets = [S for r in range(2, q + 1) for S in itertools.combinations(range(q), r)]
        rng = np.random.default_rng(1)
        res = {m: [] for m in ("EB", "convex", "GLS", "best")}
        supp = []
        for _ in range(reps):
            z, a, y = IV.simulate(n, pis, 1.0, gamma, rng)
            est, D = zip(*(IV.tsls_if(z, a, y, S, q) for S in subsets))
            psi, D = np.array(est), np.column_stack(D)
            Om = D.T @ D / n
            v = np.diag(Om) / n
            t = pairwise_tau2(psi, v)
            wc = convex_weights(Om)
            supp.append((wc > 1e-6).sum())
            g = np.linalg.pinv(Om, rcond=1e-8) @ np.ones(len(psi))
            for m, w in (("EB", 1 / (v + t)), ("convex", wc), ("GLS", g), ("best", np.eye(len(psi))[np.argmin(v)])):
                w = w / w.sum()
                res[m].append((w @ psi, np.sqrt(max(w @ Om @ w, 0) / n)))
        for m, r in res.items():
            r = np.array(r); err = r[:, 0] - 1.0
            print(f"{label:28s} {m:8s} {err.mean():8.4f} {r[:, 0].std():8.4f} {np.sqrt(np.mean(err ** 2)):8.4f} "
                  f"{r[:, 1].mean():8.4f} {np.mean(np.abs(err) <= Z90 * r[:, 1]):6.3f} "
                  f"{(np.median(supp) if m == 'convex' else np.nan):8.0f}")


if __name__ == "__main__":
    part1()
    part2()
    part3()
