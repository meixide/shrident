#!/usr/bin/env python3
"""Card-Krueger application: pooled estimates and prediction intervals from the 25 estimates.

    python apps/ck_intervals.py --out results/ck94_intervals_B2000_s42

Reads data/ck94_functionals.csv (written by data/ck94_functionals.R) and refits nothing.
Writes pooled.csv (the four sets of the article's table), intervals.csv (the split
conformal interval of Theorem 3 for one recorded split and its distribution over B
random splits; the leave-one-out interval), split.csv (the recorded split) and MANIFEST.

The leave-one-out scores are checked against the column S_loo of the input; the script
exits 3 if they differ. The leave-one-out interval is centred at the pooled estimate
computed with the pairwise dispersion, the same dispersion that scales it. The R script
centred it at the DerSimonian-Laird pooled estimate (2.54 instead of 2.52), which is why
the article reported [-2.60, 7.85] until 2026-10-12.
"""
import argparse
import hashlib
import os
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inference as I

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "ck94_functionals.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--B", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists. Choose a new name; a run directory is never reused.")
    os.makedirs(a.out)

    d = pd.read_csv(DATA)
    psi, v = d.psi.to_numpy(), d.v.to_numpy()
    J = len(psi)
    v_new = float(np.median(v))

    # ---- pooled estimates, the four sets of the article's table
    rows = []
    for name, keep in [("DIM only", d.method == "DIM"), ("OM only", d.method == "OM"),
                       ("IPW only", d.method == "IPW"), ("All", d.method.notna())]:
        k = keep.to_numpy()
        tp, td = I.tau2_pairwise(psi[k], v[k]), I.tau2_dl(psi[k], v[k])
        rows.append(dict(set=name, J=int(k.sum()), tau_pair=np.sqrt(tp), tau_dl=np.sqrt(td),
                         psi_eb_pair=float(I.eb_weights(v[k], tp) @ psi[k]),
                         psi_eb_dl=float(I.eb_weights(v[k], td) @ psi[k])))
    pooled = pd.DataFrame(rows)
    pooled.to_csv(os.path.join(a.out, "pooled.csv"), index=False, float_format="%.6f")

    # ---- leave-one-out interval, with a check of the scores against the R output
    loo = I.loo_conformal(psi, v, v_new)
    gap = float(np.max(np.abs(loo["scores"] - d.S_loo.to_numpy())))

    # ---- split conformal: one recorded split, and the distribution over B splits
    rng = np.random.default_rng(a.seed)
    n_tr = J // 2
    n_cal = J - n_tr
    alpha = 2 / (n_cal + 1)
    los, his, t2s, first = [], [], [], None
    for b in range(a.B):
        idx = rng.permutation(J)
        T, C = idx[:n_tr], idx[n_tr:]
        t = I.tau2_pairwise(psi[T], v[T])
        m = float(I.eb_weights(v[T], t) @ psi[T])
        S = (psi[C] - m) / np.sqrt(t + v[C])
        sc = np.sqrt(t + v_new)
        los.append(m + I._conf_q(S, alpha / 2) * sc)
        his.append(m + I._conf_q(S, 1 - alpha / 2) * sc)
        t2s.append(t)
        if b == 0:
            first = d[["label"]].assign(role=np.where(np.isin(np.arange(J), T), "train", "calibration"))
    los, his, t2s = np.array(los), np.array(his), np.array(t2s)
    first.to_csv(os.path.join(a.out, "split.csv"), index=False)
    q = lambda x, p: float(np.quantile(x, p))
    inter = pd.DataFrame([
        dict(interval="split, recorded split", level=1 - alpha, lo=los[0], hi=his[0]),
        dict(interval="split, median over splits", level=1 - alpha, lo=q(los, .5), hi=q(his, .5)),
        dict(interval="split, first quartile of each endpoint", level=1 - alpha, lo=q(los, .25), hi=q(his, .25)),
        dict(interval="split, third quartile of each endpoint", level=1 - alpha, lo=q(los, .75), hi=q(his, .75)),
        dict(interval="split, 5% quantile of each endpoint", level=1 - alpha, lo=q(los, .05), hi=q(his, .05)),
        dict(interval="split, 95% quantile of each endpoint", level=1 - alpha, lo=q(los, .95), hi=q(his, .95)),
        dict(interval="leave-one-out", level=loo["level"], lo=loo["lo"], hi=loo["hi"]),
    ])
    inter.to_csv(os.path.join(a.out, "intervals.csv"), index=False, float_format="%.6f")

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
             f"input=data/ck94_functionals.csv sha256={sha}",
             f"J={J} n_train={n_tr} n_cal={n_cal}",
             f"split_level=1-2/(n_cal+1)={1 - alpha:.6f}",
             f"splits_B={a.B}", f"split_seed={a.seed}",
             "recorded_split=first permutation of numpy default_rng(seed); stored in split.csv",
             "v_new=median of the 25 variances", f"v_new_value={v_new:.6f}",
             "dispersion=pairwise estimator, for the training fit and for the leave-one-out fits",
             f"share_of_splits_with_tau2_train_zero={float((t2s == 0).mean()):.4f}",
             f"loo_scores_max_gap_to_R_output={gap:.3g}",
             f"numpy_version={np.__version__}", f"pandas_version={pd.__version__}"]
    with open(os.path.join(a.out, "MANIFEST"), "w") as f:
        f.write("\n".join(lines) + "\n")

    pd.set_option("display.width", 200)
    print(pooled.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print()
    print(inter.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\n" + "\n".join(lines))
    if gap > 1e-6:
        print(f"\nCHECK FAILED: leave-one-out scores differ from the R output by up to {gap:.3g}")
        sys.exit(3)


if __name__ == "__main__":
    main()
