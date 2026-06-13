#!/usr/bin/env python3
"""lalonde_baseline_subsample.py

Apples-to-apples baseline for the HAPC dry-run: computes the same J=15
functionals with the ORIGINAL base learners (OLS for outcome regressions,
logistic regression for propensity/sampling scores) on the IDENTICAL random
subsample (same SUBSAMPLE_FRAC / SUBSAMPLE_SEED) used by
lalonde_functionals_hapc.py.

This isolates the effect of the base learner (HAPC vs OLS/logistic) from the
effect of subsampling, since the document's published numbers are full-data.

OLS  : statsmodels.OLS                  (matches R lm)
Logit: statsmodels.Logit (unpenalised)  (matches R glm(family=binomial));
       perfect separation on tiny groups is caught and reported.
"""
import os
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.simplefilter("ignore")
DIR = os.environ.get("HAPC_OUTDIR", os.path.dirname(os.path.abspath(__file__)))
COVS = ["age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75"]
SUBSAMPLE_FRAC = float(os.environ.get("HAPC_FRAC", "0.25"))
SUBSAMPLE_SEED = int(os.environ.get("HAPC_SEED", "42"))

def var_mean(y):
    y = np.asarray(y, float); return np.var(y, ddof=1) / len(y)

def var_hajek(y, w):
    y = np.asarray(y, float); w = np.asarray(w, float)
    mu = np.sum(w * y) / np.sum(w); phi = w * (y - mu) / np.mean(w)
    return np.var(phi, ddof=1) / len(phi)

def ols_predict(Xtr, Ytr, Xnew):
    Xtr = sm.add_constant(np.asarray(Xtr, float))
    fit = sm.OLS(np.asarray(Ytr, float), Xtr).fit()
    Xn = sm.add_constant(np.asarray(Xnew, float), has_constant="add")
    return np.asarray(fit.predict(Xn))

def logit_predict(Xtr, Dtr, Xnew):
    Xtr_c = sm.add_constant(np.asarray(Xtr, float))
    try:
        fit = sm.Logit(np.asarray(Dtr, float), Xtr_c).fit(disp=0, maxiter=200)
        Xn = sm.add_constant(np.asarray(Xnew, float), has_constant="add")
        p = np.asarray(fit.predict(Xn))
    except Exception as e:
        print(f"    [logit FAILED -> {type(e).__name__}; falling back to marginal rate]")
        p = np.full(len(Xnew), float(np.mean(Dtr)))
    return np.clip(p, 1e-6, 1 - 1e-6)

def main():
    d = pd.read_csv(f"{DIR}/lalonde_combined.csv")
    d = d.sample(frac=SUBSAMPLE_FRAC, random_state=SUBSAMPLE_SEED).reset_index(drop=True)
    nsw_trt = d[(d.Z == 0) & (d.A == 1)]; nsw_ctrl = d[(d.Z == 0) & (d.A == 0)]
    nsw = d[d.Z == 0]
    ec_list = {"CPS-2": d[d.Z == 1], "CPS-3": d[d.Z == 2], "PSID-2": d[d.Z == 3], "PSID-3": d[d.Z == 4]}
    n1 = len(nsw_trt); theta1 = nsw_trt.Y.mean()
    rows = []
    add = lambda j, l, psi, v, th0: rows.append(dict(j=j, label=l, psi=psi, v=v, se=np.sqrt(v), theta0=th0))

    th0 = nsw_ctrl.Y.mean()
    add(1, "DIM-RCT", theta1 - th0, var_mean(nsw_trt.Y) + var_mean(nsw_ctrl.Y), th0)

    e = logit_predict(nsw[COVS], nsw.A, nsw_ctrl[COVS]); w2 = e / (1 - e)
    th0 = float(np.sum(w2 * nsw_ctrl.Y) / np.sum(w2))
    add(2, "IPW-NB", theta1 - th0, var_mean(nsw_trt.Y) + var_hajek(nsw_ctrl.Y, w2), th0)

    pr = ols_predict(nsw_ctrl[COVS], nsw_ctrl.Y,
                     np.vstack([nsw_trt[COVS].to_numpy(float), nsw_ctrl[COVS].to_numpy(float)]))
    mt, mc = pr[:n1], pr[n1:]; th0 = float(mt.mean())
    v0 = np.var(nsw_ctrl.Y.to_numpy(float) - mc, ddof=1) / len(nsw_ctrl) + np.var(mt, ddof=1) / n1
    add(3, "OR-NB", theta1 - th0, var_mean(nsw_trt.Y) + v0, th0)

    for ki, (k, ec) in enumerate(ec_list.items()):
        base = 3 + 3 * ki
        th0 = ec.Y.mean()
        add(base + 1, f"DIM-{k}", theta1 - th0, var_mean(nsw_trt.Y) + var_mean(ec.Y), th0)
        pool = pd.concat([nsw_ctrl, ec])
        pr = ols_predict(pool[COVS], pool.Y,
                         np.vstack([nsw_trt[COVS].to_numpy(float), pool[COVS].to_numpy(float)]))
        mt, mp = pr[:n1], pr[n1:]; th0 = float(mt.mean())
        v0 = np.var(pool.Y.to_numpy(float) - mp, ddof=1) / len(pool) + np.var(mt, ddof=1) / n1
        add(base + 2, f"OM-{k}", theta1 - th0, var_mean(nsw_trt.Y) + v0, th0)
        pool_s = pd.concat([nsw.assign(S=1), ec.assign(S=0)])
        pi = logit_predict(pool_s[COVS], pool_s.S, ec[COVS]); w = pi / (1 - pi)
        th0 = float(np.sum(w * ec.Y) / np.sum(w))
        add(base + 3, f"IPW-{k}", theta1 - th0, var_mean(nsw_trt.Y) + var_hajek(ec.Y, w), th0)

    res = pd.DataFrame(rows).sort_values("j").reset_index(drop=True)
    res.to_csv(f"{DIR}/lalonde_functionals_ols_sub.csv", index=False)
    pd.set_option("display.width", 200)
    print("==== OLS/logistic functionals on the SAME subsample ====")
    print(res[["j", "label", "psi", "se", "theta0"]].to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print("\nSaved: lalonde_functionals_ols_sub.csv")

if __name__ == "__main__":
    main()
