#!/usr/bin/env python3
"""lalonde_conformal_hapc.py

Conformal prediction intervals + LOO scores for the hapc-based functionals.
Reads lalonde_functionals_hapc.csv (from lalonde_functionals_hapc.py) and applies
the *identical* conformal machinery as data/lalonde_conformal.R. No base learner
is fit here; this only post-processes (psi_j, v_j).
"""
import os
import numpy as np
import pandas as pd

DIR = os.environ.get("HAPC_OUTDIR", os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(42)

def pair_tau2(psi, v):
    i, j = np.triu_indices(len(psi), 1)
    return max(0.0, np.mean(((psi[i] - psi[j]) ** 2 - (v[i] + v[j])) / 2))

def eb_mean(psi, v, tau2):
    w = 1 / (v + tau2); return np.sum(w * psi) / np.sum(w)

def conf_q(scores, p):
    J = len(scores)
    idx = int(np.floor(p * (J + 1))) if p <= 0.5 else int(np.ceil(p * (J + 1)))
    s = np.sort(scores)
    return np.nan if idx < 1 or idx > J else s[idx - 1]

def split_pi(psiT, vT, psiC, vC, v_new, alpha):
    t = pair_tau2(psiT, vT); m = eb_mean(psiT, vT, t)
    S = (psiC - m) / np.sqrt(t + vC)
    qlo, qhi = conf_q(S, alpha / 2), conf_q(S, 1 - alpha / 2)
    sc = np.sqrt(t + v_new)
    lo = np.nan if np.isnan(qlo) else m + qlo * sc
    hi = np.nan if np.isnan(qhi) else m + qhi * sc
    return lo, hi, m

def parametric_pi(psi, v, v_new, alpha=0.10):
    t = pair_tau2(psi, v); m = eb_mean(psi, v, t); z = 1.6448536269514722
    sc = np.sqrt(t + v_new); return m - z * sc, m + z * sc, t, m

def main():
    d = pd.read_csv(f"{DIR}/lalonde_functionals_hapc.csv")
    psi = d.psi.to_numpy(); v = d.v.to_numpy(); v_new = v.mean()
    alpha_split, alpha_loo = 2 / 9, 2 / 16

    pn_lo, pn_hi, t_full, eb_full = parametric_pi(psi, v, v_new, 0.10)

    B = 2000
    los, his = [], []
    for _ in range(B):
        idx = rng.permutation(15); T, C = idx[:7], idx[7:]
        lo, hi, _ = split_pi(psi[T], v[T], psi[C], v[C], v_new, alpha_split)
        los.append(lo); his.append(hi)
    lo_med, hi_med = np.nanmedian(los), np.nanmedian(his)
    sd_lo = np.nanstd(los, ddof=1)

    loo_S = np.array([
        (psi[j] - eb_mean(np.delete(psi, j), np.delete(v, j),
                          pair_tau2(np.delete(psi, j), np.delete(v, j))))
        / np.sqrt(pair_tau2(np.delete(psi, j), np.delete(v, j)) + v[j])
        for j in range(15)])
    qlo, qhi = conf_q(loo_S, alpha_loo / 2), conf_q(loo_S, 1 - alpha_loo / 2)
    lo_loo = eb_full + qlo * np.sqrt(t_full + v_new)
    hi_loo = eb_full + qhi * np.sqrt(t_full + v_new)

    print(f"tau_hat (pairwise, J=15) = ${np.sqrt(t_full):,.0f}   psi_EB(full) = ${eb_full:,.0f}")
    print("\n== Prediction intervals for a new functional (v_new = mean v_j) ==")
    print(f"  Normal 90%            : [{pn_lo:+,.0f}, {pn_hi:+,.0f}]  width {pn_hi-pn_lo:,.0f}")
    print(f"  Split conformal 77.8% : [{lo_med:+,.0f}, {hi_med:+,.0f}]  width {hi_med-lo_med:,.0f}  (SD lo {sd_lo:,.0f}, B={B})")
    print(f"  LOO conformal 87.5%   : [{lo_loo:+,.0f}, {hi_loo:+,.0f}]  width {hi_loo-lo_loo:,.0f}")

    pd.DataFrame({
        "method": ["Normal 90%", "Split conformal 77.8%", "LOO conformal 87.5%"],
        "lower": [pn_lo, lo_med, lo_loo], "upper": [pn_hi, hi_med, hi_loo],
        "width": [pn_hi - pn_lo, hi_med - lo_med, hi_loo - lo_loo],
    }).to_csv(f"{DIR}/lalonde_pi_hapc.csv", index=False)
    print("\nSaved: lalonde_pi_hapc.csv")

if __name__ == "__main__":
    main()
