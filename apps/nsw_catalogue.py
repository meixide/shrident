#!/usr/bin/env python3
"""NSW: a prediction interval over a declared catalogue of formulas, drawn at random.

    python apps/nsw_catalogue.py --out results/nsw_catalogue_dr_parametric_J30_s20261012

Catalogue (512 formulas): the effect on the NSW treated estimated against one survey
sample (CPS-2 or PSID-2) by a doubly robust estimator with a linear outcome model and a
logistic score in a subset S of the 8 covariates, for every subset S; the empty subset is
the difference in means. Variances come from the efficient influence function evaluated
at the fitted nuisances.

The J+1 formulas used are drawn from the catalogue at random without replacement. Given
the data the catalogue's estimates are fixed numbers, so the drawn ones are exchangeable,
and the split conformal interval covers the ESTIMATE of the last drawn formula with
probability at least 1 - alpha over the draw, for every sample size.

Output: catalogue.csv; draw.csv (the recorded draw and the role of each formula);
interval.csv (the recorded interval, and a summary over R further draws: coverage,
endpoints, how often the experimental benchmark is inside); MANIFEST.
"""
import argparse
import hashlib
import itertools
import os
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inference as I

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "lalonde_combined.csv")
COVS = ["age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75"]
SAMPLES = {1: "CPS-2", 3: "PSID-2"}


def logit_fit(X, y, iters=100):
    b = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(X @ b, -30, 30)))
        W = p * (1 - p) + 1e-9
        step = np.linalg.solve(X.T @ (X * W[:, None]) + 1e-8 * np.eye(X.shape[1]), X.T @ (y - p))
        b += step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b


def formula(trt, ec, S):
    """Estimate and variance of one formula: comparison sample ec, covariates S."""
    y1, y0 = trt.Y.to_numpy(float), ec.Y.to_numpy(float)
    n1 = len(y1)
    th1 = y1.mean()
    if not S:
        return th1 - y0.mean(), y1.var(ddof=1) / n1 + y0.var(ddof=1) / len(y0)
    A1, A0 = trt[S].to_numpy(float), ec[S].to_numpy(float)
    both = np.vstack([A1, A0])
    mu, sd = both.mean(0), both.std(0) + 1e-12
    X1 = np.column_stack([np.ones(n1), (A1 - mu) / sd])
    X0 = np.column_stack([np.ones(len(y0)), (A0 - mu) / sd])
    beta = np.linalg.lstsq(X0, y0, rcond=None)[0]
    gam = logit_fit(np.vstack([X1, X0]), np.r_[np.ones(n1), np.zeros(len(y0))])
    g = np.clip(1 / (1 + np.exp(-np.clip(X0 @ gam, -30, 30))), 1e-6, 1 - 1e-6)
    om = g / (1 - g)
    res = y0 - X0 @ beta
    b = (om * res).sum() / om.sum()
    m1 = X1 @ beta
    c1 = ((y1 - th1) - (m1 - m1.mean())) / n1
    c0 = -om * (res - b) / om.sum()
    return th1 - (m1.mean() + b), float((c1 ** 2).sum() + (c0 ** 2).sum())


def interval(psi, v, idx, n_tr, alpha):
    T, C, new = idx[:n_tr], idx[n_tr:-1], idx[-1]
    t = I.tau2_pairwise(psi[T], v[T])
    m = float(I.eb_weights(v[T], t) @ psi[T])
    S = (psi[C] - m) / np.sqrt(t + v[C])
    sc = np.sqrt(t + v[new])
    return m + I._conf_q(S, alpha / 2) * sc, m + I._conf_q(S, 1 - alpha / 2) * sc, t, m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--J", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--R", type=int, default=20000)
    a = ap.parse_args()
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists. Choose a new name; a run directory is never reused.")
    os.makedirs(a.out)

    d = pd.read_csv(DATA)
    trt = d[(d.Z == 0) & (d.A == 1)]
    bench = trt.Y.mean() - d[(d.Z == 0) & (d.A == 0)].Y.mean()
    rows = []
    for zc, name in SAMPLES.items():
        ec = d[d.Z == zc]
        for r in range(len(COVS) + 1):
            for S in itertools.combinations(COVS, r):
                psi, v = formula(trt, ec, list(S))
                rows.append(dict(sample=name, covariates="+".join(S) if S else "(none)", k=len(S), psi=psi, se=np.sqrt(v)))
    cat = pd.DataFrame(rows).rename_axis("formula").reset_index()
    cat.to_csv(os.path.join(a.out, "catalogue.csv"), index=False, float_format="%.4f")
    psi, v = cat.psi.to_numpy(), cat.se.to_numpy() ** 2
    N, J = len(cat), a.J
    n_tr = J // 2
    n_cal = J - n_tr
    alpha = 2 / (n_cal + 1)

    rng = np.random.default_rng(a.seed)
    idx = rng.choice(N, size=J + 1, replace=False)                 # the recorded draw
    lo, hi, t, m = interval(psi, v, idx, n_tr, alpha)
    role = ["training"] * n_tr + ["calibration"] * n_cal + ["new"]
    cat.loc[idx].assign(role=role).to_csv(os.path.join(a.out, "draw.csv"), index=False, float_format="%.4f")
    wT = I.eb_weights(v[idx[:J]], I.tau2_pairwise(psi[idx[:J]], v[idx[:J]]))
    rec = dict(what="recorded draw", nominal=1 - alpha, lo=lo, hi=hi, tau_train=np.sqrt(t), pooled_train=m,
               new_estimate=psi[idx[-1]], new_inside=bool(lo <= psi[idx[-1]] <= hi), benchmark_inside=bool(lo <= bench <= hi),
               tau_drawn=np.sqrt(I.tau2_pairwise(psi[idx[:J]], v[idx[:J]])), pooled_drawn=float(wT @ psi[idx[:J]]))

    cover, los, his, bin_, t0 = [], [], [], [], []
    for _ in range(a.R):
        ix = rng.choice(N, size=J + 1, replace=False)
        l, h, tt, _ = interval(psi, v, ix, n_tr, alpha)
        cover.append(l <= psi[ix[-1]] <= h); los.append(l); his.append(h); bin_.append(l <= bench <= h); t0.append(tt == 0)
    los, his = np.array(los), np.array(his)
    summ = dict(what=f"{a.R} further draws", nominal=1 - alpha, coverage=float(np.mean(cover)),
                lo_median=float(np.median(los)), hi_median=float(np.median(his)),
                lo_q25=float(np.quantile(los, .25)), lo_q75=float(np.quantile(los, .75)),
                hi_q25=float(np.quantile(his, .25)), hi_q75=float(np.quantile(his, .75)),
                benchmark_inside=float(np.mean(bin_)), tau_train_zero=float(np.mean(t0)))
    pd.DataFrame([rec, summ]).to_csv(os.path.join(a.out, "interval.csv"), index=False, float_format="%.4f")

    def sh(*cmd):
        try:
            return subprocess.check_output(cmd, cwd=ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "unavailable"
    with open(DATA, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    lines = [f"run={os.path.basename(os.path.normpath(a.out))}",
             f"created_utc={datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
             f"git_commit={sh('git', 'rev-parse', 'HEAD')}",
             f"git_dirty={'yes' if sh('git', 'status', '--porcelain', '--untracked-files=no') else 'no'}",
             f"data=data/lalonde_combined.csv sha256={sha}",
             f"catalogue={N} formulas: samples {', '.join(SAMPLES.values())} x every subset of {len(COVS)} covariates",
             "estimator=doubly robust for the effect on the treated; linear outcome model, logistic score; empty subset = difference in means",
             "variance=efficient influence function at the fitted nuisances",
             f"experimental_benchmark={bench:.2f}",
             f"J={J} n_train={n_tr} n_cal={n_cal} nominal_level=1-2/(n_cal+1)={1 - alpha:.6f}",
             f"draw_seed={a.seed} (numpy default_rng; the recorded draw is the first, stored in draw.csv)",
             f"further_draws_R={a.R}", "sampling=without replacement from the catalogue",
             f"numpy_version={np.__version__}", f"pandas_version={pd.__version__}"]
    with open(os.path.join(a.out, "MANIFEST"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"catalogue: {N} formulas; estimates {psi.min():.0f} to {psi.max():.0f}, median {np.median(psi):.0f}; "
          f"standard errors {np.sqrt(v.min()):.0f} to {np.sqrt(v.max()):.0f}; pairwise tau {np.sqrt(I.tau2_pairwise(psi, v)):.0f}; "
          f"benchmark {bench:.0f}")
    print(pd.DataFrame([rec]).T.to_string(header=False)); print(pd.DataFrame([summ]).T.to_string(header=False))


if __name__ == "__main__":
    main()
