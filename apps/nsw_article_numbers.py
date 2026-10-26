#!/usr/bin/env python3
"""Numbers of the NSW application as they appear in the article, from an assembled run.

    python apps/nsw_article_numbers.py --run results/nsw_if_hapc_K5_s20261012

Reads assembled/ of the run and writes assembled/article/: estimators.csv (the estimators
of the article's set with their weights in the covariate-adjusted analysis), pooled.csv
(the three nested analyses, for every set), contrasts.csv and summary.txt. Nothing is
refitted; this only selects and combines what the assembly stored.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inference as I

SET = "onestep_xfit"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    a = ap.parse_args()
    A = os.path.join(a.run, "assembled")
    out = os.path.join(A, "article")
    os.makedirs(out, exist_ok=True)

    f = pd.read_csv(os.path.join(A, f"functionals_{SET}.csv"))
    lab = list(f.label)
    psi = f.psi.to_numpy()
    Phi = pd.read_csv(os.path.join(A, f"Phi_{SET}.csv"), index_col=0)[lab].to_numpy()
    S = Phi.T @ Phi
    crude = np.array([l.startswith("DIM-") and l != "DIM-RCT" for l in lab])
    adj = ~crude
    Sa = S[np.ix_(adj, adj)]
    w_prec, w_conv = np.zeros(len(lab)), np.zeros(len(lab))
    w_prec[adj], w_conv[adj] = I.eb_weights(np.diag(Sa), 0.0), I.convex_weights(Sa)
    est = f[["label", "block", "psi", "se", "ess", "max_weight"]].assign(
        crude=crude, w_precision_adjusted=np.where(adj, w_prec, np.nan), w_convex_adjusted=np.where(adj, w_conv, np.nan))
    est.to_csv(os.path.join(out, "estimators.csv"), index=False, float_format="%.6f")

    inf = pd.read_csv(os.path.join(A, "inference.csv"))
    keep = ["set", "analysis", "J", "tau_pair", "tau_dl", "dispersion_z", "psi_eb_zero", "se_sandwich_zero", "ci_lo_zero",
            "ci_hi_zero", "se_working_zero", "psi_convex_zero", "se_convex_zero", "ci_lo_convex_zero", "ci_hi_convex_zero",
            "support_convex_zero", "psi_eb_pair", "psi_eb_dl", "psi_convex_pair", "best_single", "se_best_single",
            "corr_mean", "corr_min", "corr_max"]
    inf[keep].to_csv(os.path.join(out, "pooled.csv"), index=False, float_format="%.4f")

    e = lambda l: np.array([x == l for x in lab], float)
    cons = [("precision-weighted minus convex, covariate-adjusted set", w_prec - w_conv)]
    cons += [(f"{l} minus DIM-RCT", e(l) - e("DIM-RCT")) for l in lab if l != "DIM-RCT"]
    rows = [dict(contrast=n, estimate=float(c @ psi), se=float(np.sqrt(c @ S @ c))) for n, c in cons]
    idx = np.flatnonzero(adj)
    best = None
    for i in idx:
        for j in idx:
            if j < i:
                c = e(lab[i]) - e(lab[j])
                z = (c @ psi) / np.sqrt(c @ S @ c)
                if best is None or abs(z) > abs(best[1]):
                    best = (f"{lab[i]} minus {lab[j]} (largest of {len(idx) * (len(idx) - 1) // 2} pairs, covariate-adjusted set)", z, c)
    rows.append(dict(contrast=best[0], estimate=float(best[2] @ psi), se=float(np.sqrt(best[2] @ S @ best[2]))))
    con = pd.DataFrame(rows)
    con["z"] = con.estimate / con.se
    con.to_csv(os.path.join(out, "contrasts.csv"), index=False, float_format="%.4f")

    c = Sa / np.sqrt(np.outer(np.diag(Sa), np.diag(Sa)))
    iu = np.triu_indices(len(Sa), 1)
    fl = pd.read_csv(os.path.join(A, "fit_log.csv"))
    lines = [f"set={SET}", f"persons={Phi.shape[0]}",
             f"adjusted_set_correlations: min {c[iu].min():.3f} median {np.median(c[iu]):.3f} max {c[iu].max():.3f}",
             f"fits={len(fl)} interior={int((fl.cv_verdict == 'interior').sum())} hours={fl.elapsed_s.sum() / 3600:.1f} "
             f"longest_fit_hours={fl.elapsed_s.max() / 3600:.1f} hapc={','.join(sorted(fl.hapc_version.astype(str).unique()))}",
             "not_interior=" + "; ".join(f"{r.nuisance} block {r.block} fold {r.fold}: {r.cv_verdict}"
                                         for r in fl[fl.cv_verdict != "interior"].itertuples())]
    open(os.path.join(out, "summary.txt"), "w").write("\n".join(lines) + "\n")
    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(est.round(3).to_string(index=False)); print(); print(con.round(2).to_string(index=False)); print(); print("\n".join(lines))
    print(inf.loc[inf.set == SET, ["analysis", "J", "tau_pair", "tau_dl", "psi_eb_zero", "se_sandwich_zero", "se_working_zero",
                                   "psi_convex_zero", "se_convex_zero", "psi_eb_pair", "psi_eb_dl", "psi_convex_pair", "se_best_single"]].round(1).to_string(index=False))


if __name__ == "__main__":
    main()
