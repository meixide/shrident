#!/usr/bin/env python3
"""Card-Krueger application: pooled estimates and prediction intervals from the 25 estimates.

    python apps/ck_intervals.py --out results/<new directory> --jackknife results/ck94_jackknife_stratified

With --jackknife, estimates, variances and the restaurant-level influence matrix come from
a run of apps/ck_jackknife.R, so every pooled estimate has a sandwich standard error and
the convex combination can be computed. Without it the script uses the variances of
data/ck94_functionals.csv, which omit the Pennsylvania-side variance of the adjusted
estimators (found on 2026-10-12); that mode is kept only to reproduce
results/ck94_intervals_B2000_s42.

Output: pooled.csv (one line per set of estimators), weights.csv and contrasts.csv (with
--jackknife), estimates.csv, intervals.csv (split conformal interval
for one recorded split and its distribution over B random splits; leave-one-out interval),
split.csv (the recorded split), MANIFEST. Nothing is refitted.

The leave-one-out interval is centred at the pooled estimate computed with the pairwise
dispersion, the same dispersion that scales it. In the mode without --jackknife its scores
are checked against the column S_loo of the R output (exit status 3 if they differ).
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
    ap.add_argument("--jackknife", default=None, help="directory written by apps/ck_jackknife.R")
    ap.add_argument("--B", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    a = ap.parse_args()
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists. Choose a new name; a run directory is never reused.")

    d = pd.read_csv(DATA)
    Phi = None
    if a.jackknife:
        jk = pd.read_csv(os.path.join(a.jackknife, "psi.csv"))
        if list(jk.label) != list(d.label) or np.max(np.abs(jk.psi - d.psi)) > 1e-8:
            raise SystemExit("the jackknife run does not match data/ck94_functionals.csv")
        P = pd.read_csv(os.path.join(a.jackknife, "Phi.csv"))
        Phi = P[list(d.label)].to_numpy()
        d = d.assign(v=jk.v.to_numpy(), se=jk.se.to_numpy())
    os.makedirs(a.out)
    psi, v = d.psi.to_numpy(), d.v.to_numpy()
    J = len(psi)
    v_new = float(np.median(v))

    # ---- sets of estimators. Five adjusted pairs coincide (one binary covariate: the outcome
    # model and the weighting estimator are the same number); "distinct" keeps one of each.
    dup = np.zeros(J, bool)
    for i in range(J):
        for j in range(i):
            if d.ctrl[i] == d.ctrl[j] and d.method[j] != "DIM" and abs(psi[i] - psi[j]) < 1e-8:
                dup[i] = True
    sets = [("DIM only", (d.method == "DIM").to_numpy()), ("OM only", (d.method == "OM").to_numpy()),
            ("IPW only", (d.method == "IPW").to_numpy()), ("All", np.ones(J, bool)), ("All distinct", ~dup)]
    rows = []
    for name, k in sets:
        tp, td = I.tau2_pairwise(psi[k], v[k]), I.tau2_dl(psi[k], v[k])
        r = dict(set=name, J=int(k.sum()), tau_pair=np.sqrt(tp), tau_dl=np.sqrt(td),
                 psi_eb_pair=float(I.eb_weights(v[k], tp) @ psi[k]), psi_eb_dl=float(I.eb_weights(v[k], td) @ psi[k]),
                 se_min=float(np.sqrt(v[k].min())), se_max=float(np.sqrt(v[k].max())))
        if Phi is not None:
            q = I.pooled(psi[k], Phi[:, k], n_units=Phi.shape[0])
            r.update({c: q[c] for c in ("tau2_over_mean_v", "dispersion_z",
                                        "se_sandwich_pair", "psi_eb_zero", "se_sandwich_zero", "ci_lo_zero", "ci_hi_zero",
                                        "psi_convex_zero", "se_convex_zero", "ci_lo_convex_zero", "ci_hi_convex_zero",
                                        "support_convex_zero", "psi_convex_pair", "se_convex_pair", "se_working_pair",
                                        "psi_best_single", "se_best_single", "corr_mean", "corr_max", "sigma_rank")})
            r["best_single"] = d.label[np.flatnonzero(k)[q["best_single_index"]]]
        rows.append(r)
    pooled = pd.DataFrame(rows)
    pooled.to_csv(os.path.join(a.out, "pooled.csv"), index=False, float_format="%.6f")

    # ---- with the influence matrix: weights of the full set, and contrasts. Under exact
    # identification every contrast a'psi with 1'a = 0 is centred at zero, with variance a'Sigma a.
    if Phi is not None:
        Sig = Phi.T @ Phi
        w_prec, w_conv, w_unif = I.eb_weights(v, 0.0), I.convex_weights(Sig), np.ones(J) / J
        pd.DataFrame(dict(label=d.label, psi=psi, se=np.sqrt(v), w_precision=w_prec, w_convex=w_conv)).to_csv(
            os.path.join(a.out, "weights.csv"), index=False, float_format="%.6f")
        e = lambda lab: (d.label == lab).to_numpy().astype(float)
        cons = [("precision-weighted minus convex", w_prec - w_conv),
                ("unweighted mean minus convex", w_unif - w_conv),
                ("unweighted mean minus precision-weighted", w_unif - w_prec),
                ("DIM-PA-KFC minus DIM-PA-all", e("DIM-PA-KFC") - e("DIM-PA-all")),
                ("OM-PA-north minus OM-PA-all", e("OM-PA-north") - e("OM-PA-all"))]
        X = np.column_stack([np.ones(J), np.sqrt(v)])
        H = np.linalg.solve(X.T @ X, X.T)          # least-squares slope of the estimates on their standard errors
        rows_c = [dict(contrast=n, estimate=float(c @ psi), se=float(np.sqrt(c @ Sig @ c))) for n, c in cons]
        rows_c.append(dict(contrast="slope of the estimates on their standard errors",
                           estimate=float((H @ psi)[1]), se=float(np.sqrt((H @ Sig @ H.T)[1, 1]))))
        con = pd.DataFrame(rows_c)
        con["z"] = con.estimate / con.se
        con.to_csv(os.path.join(a.out, "contrasts.csv"), index=False, float_format="%.6f")

    # ---- leave-one-out interval
    loo = I.loo_conformal(psi, v, v_new)
    gap = np.nan if a.jackknife else float(np.max(np.abs(loo["scores"] - d.S_loo.to_numpy())))

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
    d.assign(S_loo_pairwise=loo["scores"], duplicate=dup)[
        ["label", "ctrl", "method", "n_ctrl", "psi", "v", "se", "S_loo_pairwise", "duplicate"]].to_csv(
        os.path.join(a.out, "estimates.csv"), index=False, float_format="%.6f")

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
             f"input=data/ck94_functionals.csv sha256={sha}",
             "variances=" + (f"stratified jackknife, from {a.jackknife}" if a.jackknife
                             else "column v of data/ck94_functionals.csv (adjusted estimators: too small)"),
             f"J={J} n_train={n_tr} n_cal={n_cal}",
             f"duplicated_estimators={int(dup.sum())} ({', '.join(d.label[dup])})",
             f"split_level=1-2/(n_cal+1)={1 - alpha:.6f}",
             f"splits_B={a.B}", f"split_seed={a.seed}",
             "recorded_split=first permutation of numpy default_rng(seed); stored in split.csv",
             "v_new=median of the 25 variances", f"v_new_value={v_new:.6f}",
             "dispersion=pairwise estimator, for the training fit and for the leave-one-out fits",
             "pooled estimates: `_zero` sets the dispersion to zero; `_pair` plugs in the pairwise estimate",
             f"share_of_splits_with_tau2_train_zero={float((t2s == 0).mean()):.4f}",
             f"loo_scores_max_gap_to_R_output={gap:.3g}",
             f"numpy_version={np.__version__}", f"pandas_version={pd.__version__}"]
    with open(os.path.join(a.out, "MANIFEST"), "w") as f:
        f.write("\n".join(lines) + "\n")

    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(pooled.T.to_string(float_format=lambda x: f"{x:.3f}"))
    if Phi is not None:
        print()
        print(con.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        print("\nconvex weights:", {l: round(float(x), 3) for l, x in zip(d.label, w_conv) if x > 0})
    print()
    print(inter.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("\n" + "\n".join(lines))
    if not a.jackknife and gap > 1e-6:
        print(f"\nCHECK FAILED: leave-one-out scores differ from the R output by up to {gap:.3g}")
        sys.exit(3)


if __name__ == "__main__":
    main()
