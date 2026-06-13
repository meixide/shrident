#!/usr/bin/env python3
"""lalonde_functionals_hapc.py

Re-runs the J=15 LaLonde identification functionals using the Python `hapc`
(Highly Adaptive Principal Components) library as the base regressor in place
of every lm()/glm() in data/lalonde_functionals.R.

Base-learner configuration (per user spec):
    * library                 : hapc.cv_hapc   (cross-validated)
    * norm                    : '1'            (LASSO; always)
    * max interaction degree  : max_degree = 2
    * principal components    : npcs = floor(n_train ** 0.6)
    * outcome regressions     : family = 'gaussian'
    * propensity / sampling   : family = 'binomial'
    * approx                  : False          (exact eigendecomposition)

NO linear regression (lm) or logistic GLM (glm) is used. cv_hapc's lambda grid
is in NATURAL-log units (lambda = exp(log_lambda)); the grid spans
log_lambda in [-12, 16] so both the gaussian (earnings-scale, optimum near
lambda~1e3) and binomial (optimum near lambda~0.05) fits turn over at an INTERIOR
minimum. Each fit's CV verdict (interior / boundary / null-plateau) is reported.

Run modes (all controlled by environment variables; nothing to edit):
    HAPC_FRAC   subsample fraction of the whole combined sample (default 0.25
                for a local dry-run; set HAPC_FRAC=1.0 for the full HPC run).
    HAPC_SEED   subsample RNG seed (default 42).
    HAPC_BLOCK  if set to an int 0..4, compute ONLY that block and write a
                partial CSV (for a Slurm job array); if unset/"all", run the
                whole pipeline in one process (functionals + EB + LOO).
                  block 0 = RCT-only      (j=1,2,3)
                  block 1 = CPS-2         (j=4,5,6)
                  block 2 = CPS-3         (j=7,8,9)
                  block 3 = PSID-2        (j=10,11,12)
                  block 4 = PSID-3        (j=13,14,15)
    HAPC_OUTDIR output directory (default: this script's directory).

Influence-function variances, EB shrinkage, the three nested sub-analyses and
the LOO conformal scores use exactly the same formulas as the original scripts;
only the base learner changed. DIM functionals use no base learner and are
therefore identical to the original by construction.
"""
import os
import time
import numpy as np
import pandas as pd
import hapc

DIR = os.environ.get("HAPC_OUTDIR", os.path.dirname(os.path.abspath(__file__)))
COVS = ["age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75"]

SUBSAMPLE_FRAC = float(os.environ.get("HAPC_FRAC", "0.25"))
SUBSAMPLE_SEED = int(os.environ.get("HAPC_SEED", "42"))
_BLOCK_ENV = os.environ.get("HAPC_BLOCK", "all").strip().lower()
BLOCK = None if _BLOCK_ENV in ("", "all", "none") else int(_BLOCK_ENV)

MAX_DEGREE = 2
NORM = "1"
# NB: cv_hapc's grid is in NATURAL-log units, lambda = exp(log_lambda).
# Earnings-scale (gaussian) responses need a large lambda ceiling for the LASSO
# CV curve to turn over at an interior minimum; the binomial scores minimise
# near lambda~0.05. This wide grid brackets both so every fit has an interior
# minimum where the data support one.
GRID = dict(log_lambda_min=-12.0, log_lambda_max=16.0, grid_length=40, nfolds=5)
npcs_rule = lambda n: max(1, int(np.floor(n ** 0.6)))

# block id -> (name, list of j, EC Z-code or None for the RCT block)
BLOCKS = {0: ("RCT", [1, 2, 3], None), 1: ("CPS-2", [4, 5, 6], 1),
          2: ("CPS-3", [7, 8, 9], 2), 3: ("PSID-2", [10, 11, 12], 3),
          4: ("PSID-3", [13, 14, 15], 4)}

_cvlog = []

def _check_interior(label, family, r):
    # CVResult stores the CV error curve in .mses for BOTH families
    # (gaussian: MSE; binomial: logistic deviance).
    err = np.asarray(r.mses)
    imin = int(np.argmin(err))
    rng_err = err.max() - err.min()
    plateau_high = rng_err > 0 and (err[-1] - err.min()) <= 1e-6 * max(abs(err.min()), 1.0)
    if imin == 0:
        pos = "BOUNDARY-low"
    elif imin == len(err) - 1:
        pos = "BOUNDARY-high"
    elif plateau_high:
        pos = "PLATEAU-high(null model)"
    else:
        pos = "INTERIOR"
    _cvlog.append((label, family, r.best_lambda, imin, len(err) - 1, pos))
    flag = "" if pos == "INTERIOR" else "  *** NOT INTERIOR ***"
    print(f"    [cv] {label:<11s} {family:<8s} best_lambda={r.best_lambda:.4g} "
          f"idx {imin}/{len(err)-1} -> {pos}{flag}", flush=True)

def hapc_reg(label, Xtr, Ytr, Xnew):
    """Gaussian outcome regression: fit on (Xtr,Ytr), predict rows of Xnew."""
    Xtr = np.asarray(Xtr, float); Ytr = np.asarray(Ytr, float); Xnew = np.asarray(Xnew, float)
    r = hapc.cv_hapc(Xtr, Ytr, family="gaussian", max_degree=MAX_DEGREE,
                     npcs=npcs_rule(len(Xtr)), norm=NORM, approx=False,
                     predict=Xnew, **GRID)
    _check_interior(label, "gaussian", r)
    return np.asarray(r.predictions, float).ravel()

def hapc_ps(label, Xtr, Dtr, Xnew):
    """Binomial propensity/sampling score: return P(D=1|X) at rows of Xnew."""
    Xtr = np.asarray(Xtr, float); Dtr = np.asarray(Dtr, float); Xnew = np.asarray(Xnew, float)
    r = hapc.cv_hapc(Xtr, Dtr, family="binomial", max_degree=MAX_DEGREE,
                     npcs=npcs_rule(len(Xtr)), norm=NORM, approx=False,
                     predict=Xnew, **GRID)
    _check_interior(label, "binomial", r)
    p = np.asarray(r.predictions, float).ravel()
    return np.clip(p, 1e-6, 1 - 1e-6)   # guard the Hajek ratio p/(1-p)

# ── Influence-function variance helpers (identical to the R script) ──────────
def var_mean(y):
    y = np.asarray(y, float); return np.var(y, ddof=1) / len(y)

def var_hajek(y, w):
    y = np.asarray(y, float); w = np.asarray(w, float)
    mu = np.sum(w * y) / np.sum(w)
    phi = w * (y - mu) / np.mean(w)
    return np.var(phi, ddof=1) / len(phi)

# ── Functional computation, factored so a job array can run one block ────────
def compute_rct(nsw_trt, nsw_ctrl, nsw, theta1, n1):
    rows = []
    add = lambda j, l, psi, v, th0: rows.append(dict(j=j, label=l, psi=psi, v=v, se=np.sqrt(v), theta0=th0))
    th0 = nsw_ctrl.Y.mean()
    add(1, "DIM-RCT", theta1 - th0, var_mean(nsw_trt.Y) + var_mean(nsw_ctrl.Y), th0)
    print("done j=1 DIM-RCT", flush=True)

    e = hapc_ps("IPW-NB", nsw[COVS], nsw.A, nsw_ctrl[COVS]); w2 = e / (1 - e)
    th0 = float(np.sum(w2 * nsw_ctrl.Y) / np.sum(w2))
    add(2, "IPW-NB", theta1 - th0, var_mean(nsw_trt.Y) + var_hajek(nsw_ctrl.Y, w2), th0)
    print("done j=2 IPW-NB", flush=True)

    pr = hapc_reg("OR-NB", nsw_ctrl[COVS], nsw_ctrl.Y,
                  np.vstack([nsw_trt[COVS].to_numpy(float), nsw_ctrl[COVS].to_numpy(float)]))
    mu_trt, mu_ctrl = pr[:n1], pr[n1:]; th0 = float(mu_trt.mean())
    v0 = np.var(nsw_ctrl.Y.to_numpy(float) - mu_ctrl, ddof=1) / len(nsw_ctrl) + np.var(mu_trt, ddof=1) / n1
    add(3, "OR-NB", theta1 - th0, var_mean(nsw_trt.Y) + v0, th0)
    print("done j=3 OR-NB", flush=True)
    return rows

def compute_ec(kname, ec_k, base, nsw_trt, nsw_ctrl, nsw, theta1, n1):
    rows = []
    add = lambda j, l, psi, v, th0: rows.append(dict(j=j, label=l, psi=psi, v=v, se=np.sqrt(v), theta0=th0))
    th0 = ec_k.Y.mean()
    add(base + 1, f"DIM-{kname}", theta1 - th0, var_mean(nsw_trt.Y) + var_mean(ec_k.Y), th0)
    print(f"done j={base+1} DIM-{kname}", flush=True)

    pool = pd.concat([nsw_ctrl, ec_k])
    pr = hapc_reg(f"OM-{kname}", pool[COVS], pool.Y,
                  np.vstack([nsw_trt[COVS].to_numpy(float), pool[COVS].to_numpy(float)]))
    mu_trt, mu_pool = pr[:n1], pr[n1:]; th0 = float(mu_trt.mean())
    v0 = np.var(pool.Y.to_numpy(float) - mu_pool, ddof=1) / len(pool) + np.var(mu_trt, ddof=1) / n1
    add(base + 2, f"OM-{kname}", theta1 - th0, var_mean(nsw_trt.Y) + v0, th0)
    print(f"done j={base+2} OM-{kname}", flush=True)

    pool_s = pd.concat([nsw.assign(S=1), ec_k.assign(S=0)])
    pi = hapc_ps(f"IPW-{kname}", pool_s[COVS], pool_s.S, ec_k[COVS]); w = pi / (1 - pi)
    th0 = float(np.sum(w * ec_k.Y) / np.sum(w))
    add(base + 3, f"IPW-{kname}", theta1 - th0, var_mean(nsw_trt.Y) + var_hajek(ec_k.Y, w), th0)
    print(f"done j={base+3} IPW-{kname}", flush=True)
    return rows

# ── EB shrinkage + LOO (importable by the aggregator) ────────────────────────
def _mom_tau2(psi, v):
    J = len(psi); w = 1 / v; pbar = np.sum(w * psi) / np.sum(w)
    Q = np.sum((psi - pbar) ** 2 / v); C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    return max(0.0, (Q - (J - 1)) / C)

def _pair_tau2(psi, v):
    i, j = np.triu_indices(len(psi), 1)
    return max(0.0, np.mean(((psi[i] - psi[j]) ** 2 - (v[i] + v[j])) / 2))

def _eb(psi, v, tau2):
    w = 1 / (v + tau2); return np.sum(w * psi) / np.sum(w), np.sqrt(1 / np.sum(w))

def compute_eb_loo(res, outdir=DIR):
    """Given the full 15-row functionals frame, compute & save EB + LOO."""
    res = res.sort_values("j").reset_index(drop=True)
    def block(idx, name):
        psi = res.psi.to_numpy()[idx]; v = res.v.to_numpy()[idx]
        tm, tp = _mom_tau2(psi, v), _pair_tau2(psi, v)
        em, sm = _eb(psi, v, tm); ep, sp = _eb(psi, v, tp)
        return dict(analysis=name, J=len(idx), tau_mom=np.sqrt(tm), psi_eb_mom=em,
                    se_eb_mom=sm, tau_pair=np.sqrt(tp), psi_eb_pair=ep, se_eb_pair=sp)
    adj = res.index[~res.label.str.match(r"^DIM-(CPS|PSID)")].to_numpy()
    blocks = pd.DataFrame([
        block(np.array([0, 1, 2]), "RCT-only (J=3)"),
        block(adj, "Covariate-adjusted (J=11)"),
        block(np.arange(15), "Full pool (J=15)")])
    psi15 = res.psi.to_numpy(); v15 = res.v.to_numpy()
    loo = []
    for jj in range(15):
        T = [k for k in range(15) if k != jj]
        tj = _pair_tau2(psi15[T], v15[T]); mj, _ = _eb(psi15[T], v15[T], tj)
        loo.append(dict(j=jj + 1, label=res.label[jj], psi=psi15[jj],
                        S_loo=(psi15[jj] - mj) / np.sqrt(tj + v15[jj])))
    loo = pd.DataFrame(loo)
    res.to_csv(f"{outdir}/lalonde_functionals_hapc.csv", index=False)
    blocks.to_csv(f"{outdir}/lalonde_eb_hapc.csv", index=False)
    loo.to_csv(f"{outdir}/lalonde_loo_hapc.csv", index=False)
    return blocks, loo

def load_data():
    here = os.path.dirname(os.path.abspath(__file__))
    d = pd.read_csv(f"{here}/lalonde_combined.csv")
    if SUBSAMPLE_FRAC < 1.0:
        n_full = len(d)
        d = d.sample(frac=SUBSAMPLE_FRAC, random_state=SUBSAMPLE_SEED).reset_index(drop=True)
        zmap = {0: "NSW", 1: "CPS-2", 2: "CPS-3", 3: "PSID-2", 4: "PSID-3"}
        print(f"*** SUBSAMPLE frac={SUBSAMPLE_FRAC} seed={SUBSAMPLE_SEED}: n={len(d)} of {n_full} ***", flush=True)
        cnt = d.assign(grp=d.Z.map(zmap)).groupby("grp").size()
        nt = ((d.Z == 0) & (d.A == 1)).sum(); nc = ((d.Z == 0) & (d.A == 0)).sum()
        print("    NSW treated=%d NSW control=%d | " % (nt, nc)
              + " ".join(f"{k}={cnt.get(k,0)}" for k in ["CPS-2","CPS-3","PSID-2","PSID-3"]), flush=True)
    else:
        print(f"*** FULL DATA: n={len(d)} ***", flush=True)
    return d

def main():
    t0 = time.time()
    d = load_data()
    nsw_trt = d[(d.Z == 0) & (d.A == 1)]; nsw_ctrl = d[(d.Z == 0) & (d.A == 0)]
    nsw = d[d.Z == 0]; n1 = len(nsw_trt); theta1 = nsw_trt.Y.mean()

    # ── Block mode (one job-array task computes a single block) ──────────────
    if BLOCK is not None:
        name, js, zcode = BLOCKS[BLOCK]
        print(f"=== BLOCK {BLOCK} ({name}): functionals j={js} ===", flush=True)
        if zcode is None:
            rows = compute_rct(nsw_trt, nsw_ctrl, nsw, theta1, n1)
        else:
            rows = compute_ec(name, d[d.Z == zcode], js[0] - 1, nsw_trt, nsw_ctrl, nsw, theta1, n1)
        out = f"{DIR}/lalonde_functionals_hapc_block{BLOCK}.csv"
        pd.DataFrame(rows).to_csv(out, index=False)
        nb = sum(1 for r in _cvlog if r[5] != "INTERIOR")
        print(f"\nCV: {len(_cvlog)-nb}/{len(_cvlog)} interior. Elapsed {time.time()-t0:.1f}s. Saved {out}", flush=True)
        return

    # ── Full single-process mode ─────────────────────────────────────────────
    rows = compute_rct(nsw_trt, nsw_ctrl, nsw, theta1, n1)
    for ki, (kname, zcode) in enumerate([("CPS-2",1),("CPS-3",2),("PSID-2",3),("PSID-3",4)]):
        rows += compute_ec(kname, d[d.Z == zcode], 3 + 3 * ki, nsw_trt, nsw_ctrl, nsw, theta1, n1)
    res = pd.DataFrame(rows).sort_values("j").reset_index(drop=True)
    blocks, loo = compute_eb_loo(res, DIR)

    pd.set_option("display.width", 200, "display.max_columns", 20)
    print("\n==== FUNCTIONALS (hapc base learner: norm='1', max_degree=2, npcs=n^0.6) ====")
    print(res[["j", "label", "psi", "se", "theta0"]].to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print("\n==== EB sub-analyses ====")
    print(blocks.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print("\n==== LOO conformal scores ====")
    print(loo.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nLOO scores: mean={loo.S_loo.mean():.3f}  sd={loo.S_loo.std(ddof=1):.3f}")
    nb = sum(1 for r in _cvlog if r[5] != "INTERIOR")
    print(f"\nCV interior-minimum check: {len(_cvlog)-nb}/{len(_cvlog)} interior"
          + ("" if nb == 0 else f"  ({nb} on boundary/plateau)"))
    print(f"Total elapsed: {time.time()-t0:.1f}s")
    print("Saved: lalonde_functionals_hapc.csv, lalonde_eb_hapc.csv, lalonde_loo_hapc.csv")

if __name__ == "__main__":
    main()
