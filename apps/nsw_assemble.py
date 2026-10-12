#!/usr/bin/env python3
"""Assemble estimators, unit-level influence contributions and inference for the NSW application.

    python apps/nsw_assemble.py --run results/<name>

Reads units.csv and nuis/ of a run, refits nothing, and writes assembled/.

Estimand: effect on the experimental treated, psi = theta1 - theta0, with theta1 their
mean outcome and theta0 their mean outcome without treatment. Each estimator replaces
theta0 by a different functional of the experimental controls or of an external sample.

Estimator sets
  onestep_xfit   one-step (doubly robust) estimator of each functional with out-of-fold
                 nuisances (cross-fitting); 14 estimators. Reported in the article.
  onestep_full   the same with nuisances fitted on the full sample
  plugin_eif     the 15 published plug-in estimates; contributions from the efficient
                 influence function of each functional (right if the plug-in is efficient)
  plugin_known   the same 15 estimates; contributions treat the fitted nuisances as
                 known, which is what the published variances assumed

Unit of observation: the person. CPS-3 and PSID-3 are subsamples of CPS-2 and PSID-2,
so a person can contribute to several estimators; Sigma = Phi'Phi keeps that dependence.

Exit status: 0 complete and consistent; 2 nuisance files missing; 3 the plug-in
estimates differ from the published table by more than --tol.
"""
import argparse
import glob
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import inference as I
import nsw_common as C

PUBLISHED = os.path.join(C.ROOT, "data", "lalonde_functionals_hapc.csv")
CONF_B, CONF_SEED = 2000, 42


# ---------------------------------------------------------------- building blocks
def mean_c(y):
    """Sample mean and the contribution of each row to its error."""
    return y.mean(), (y - y.mean()) / len(y)


def hajek_c(y, w):
    """Weighted mean with normalised weights and the contribution of each row."""
    t = (w * y).sum() / w.sum()
    return t, w * (y - t) / w.sum()


def odds(p):
    p = p.clip(C.CLIP, 1 - C.CLIP)
    return p / (1 - p)


def wdiag(w):
    """Effective sample size of a set of weights and the largest normalised weight."""
    return dict(ess=float(w.sum() ** 2 / (w ** 2).sum()), max_weight=float(w.max() / w.sum()), n_weighted=len(w))


def read_pred(run, block, name, folds):
    files = [C.nuis_path(run, block, name, f) for f in folds]
    missing = [f for f in files if not os.path.exists(f)]
    if missing:
        return None, missing
    s = pd.concat([pd.read_csv(f) for f in files]).set_index("row").pred
    if s.index.has_duplicates:
        raise SystemExit(f"duplicated rows in the predictions of block {block} {name}")
    return s, []


# ---------------------------------------------------------------- estimators
def estimators(d, P, blocks, plugin):
    """All estimators of theta0 for one set of nuisance predictions P[(block, name)].

    Returns {set: [dict(label, block, kind, theta0, c0, ...)]}; c0 holds the row-level
    contributions to the error of theta0. `plugin` adds the published plug-in estimators.
    """
    nsw = d[d.Z == 0]
    T, Cn = nsw[nsw.A == 1], nsw[nsw.A == 0]
    out = {"plugin_known": [], "plugin_eif": [], "onestep": []}

    def add(sets, label, block, kind, theta0, parts, diag=None):
        for s in sets:
            if s == "onestep" or plugin:
                out[s].append(dict(label=label, block=block, kind=kind, theta0=float(theta0),
                                   c0=pd.concat(parts), **(diag or {})))

    def functional(block, tag, m, target, comp, w, plug_label, plug_kind):
        """One identifying functional: mean over `target` of the regression fitted on `comp`."""
        res = comp.Y - m.loc[comp.index]
        t_pl, c_pl = mean_c(m.loc[target.index])
        b, c_res = hajek_c(res, w)
        t_ipw, c_ipw = hajek_c(comp.Y, w)
        dg = wdiag(w)
        if plug_kind == "ipw":
            add(["plugin_known"], plug_label, block, "adj", t_ipw, [c_ipw], dg)
            add(["plugin_eif"], plug_label, block, "adj", t_ipw, [c_pl, c_res], dg)
        else:
            add(["plugin_known"], plug_label, block, "adj", t_pl, [c_pl, mean_c(res)[1]])
            add(["plugin_eif"], plug_label, block, "adj", t_pl, [c_pl, c_res], dg)
        return t_pl + b, [c_pl, c_res], dg

    if 0 in blocks:
        t, c = mean_c(Cn.Y)
        add(out, "DIM-RCT", 0, "dim", t, [c])
        e, m = P[(0, "e_rct")], P[(0, "m_rct")]
        w = odds(e.loc[Cn.index])
        functional(0, "RCT", m, T, Cn, w, "IPW-NB", "ipw")
        t, parts, dg = functional(0, "RCT", m, T, Cn, w, "OR-NB", "or")
        add(["onestep"], "DR-RCT", 0, "adj", t, parts, dg)
    for k in [b for b in blocks if b > 0]:
        name = C.BLOCKS[k]
        ec = d[d.Z == k]
        pool = pd.concat([Cn, ec])
        t, c = mean_c(ec.Y)
        add(out, f"DIM-{name}", k, "dim", t, [c])
        # pooled-control regression averaged over the treated
        t, parts, dg = functional(k, name, P[(k, "m_pool")], T, pool,
                                  odds(P[(k, "g_trt")].loc[pool.index]), f"OM-{name}", "or")
        add(["onestep"], f"DRpool-{name}", k, "adj", t, parts, dg)
        # external-sample regression averaged over the whole experiment
        t, parts, dg = functional(k, name, P[(k, "m_ec")], nsw, ec,
                                  odds(P[(k, "pi_s")].loc[ec.index]), f"IPW-{name}", "ipw")
        add(["onestep"], f"DRext-{name}", k, "adj", t, parts, dg)
    return out, T


def to_matrix(d, T, ests):
    """psi (J) and the person-level influence matrix Phi (persons x J)."""
    th1, c1 = mean_c(T.Y)
    uids = np.sort(d.uid.unique())
    Phi = np.zeros((len(uids), len(ests)))
    for j, e in enumerate(ests):
        if e["c0"].index.has_duplicates:
            raise SystemExit(f"{e['label']}: a row contributes twice")
        phi = c1.reindex(d.index, fill_value=0.0) - e["c0"].reindex(d.index, fill_value=0.0)
        Phi[:, j] = phi.groupby(d.uid).sum().reindex(uids).to_numpy()
    psi = np.array([th1 - e["theta0"] for e in ests])
    return psi, Phi, uids


# ---------------------------------------------------------------- output helpers
def md_table(df, fmt="{:.1f}"):
    cols = list(df.columns)
    def cell(x):
        if isinstance(x, (bool, np.bool_)):
            return "yes" if x else "no"
        if isinstance(x, (int, np.integer)):
            return str(x)
        if isinstance(x, (float, np.floating)):
            return "" if np.isnan(x) else fmt.format(x)
        return str(x)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(cell(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--blocks", default="0,1,2,3,4", help="blocks to assemble (all for a publishable table)")
    ap.add_argument("--tol", type=float, default=0.5, help="largest admissible gap to a published estimate")
    a = ap.parse_args()
    blocks = [int(b) for b in a.blocks.split(",")]
    s = C.settings(a.run)
    K = s["K"]
    d = C.load(a.run).set_index("row", drop=False)
    out_dir = os.path.join(a.run, "assembled")
    os.makedirs(out_dir, exist_ok=True)
    rep = [f"# NSW application: {os.path.basename(os.path.normpath(a.run))}", "",
           f"Blocks assembled: {', '.join(C.BLOCKS[b] for b in blocks)}. Persons: {d.uid.nunique()}; rows: {len(d)}. "
           f"Cross-fitting folds: {K}."]
    publishable = s["frac"] == 1.0 and not s["smoke"] and len(blocks) == len(C.BLOCKS)
    if not publishable:
        rep += ["", "**Not a publishable run**: subsample, coarse grid or a subset of blocks."]

    # ---- nuisance predictions, full-sample and out-of-fold
    versions, missing_all = {}, {}
    for version, folds in (("full", [0]), ("xfit", list(range(1, K + 1)))):
        if not folds:
            continue
        P, missing = {}, []
        for b in blocks:
            for name, (_, _, _, pred) in C.nuisances(d, b).items():
                p, miss = read_pred(a.run, b, name, folds)
                missing += miss
                if p is not None:
                    if set(p.index) != set(pred.row):
                        raise SystemExit(f"block {b} {name} ({version}): predictions do not cover the required rows")
                    P[(b, name)] = p
        if missing:
            missing_all[version] = missing
        else:
            versions[version] = P
    if "full" not in versions:
        print("missing full-sample fits:\n  " + "\n  ".join(missing_all["full"]))
        sys.exit(2)

    # ---- fit log
    metas = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(a.run, "nuis", "*.json")))]
    cv = pd.DataFrame(metas)[["task", "block", "nuisance", "fold", "family", "n_train", "npcs", "best_lambda",
                              "cv_index", "cv_grid_last", "cv_verdict", "elapsed_s", "hapc_version"]].sort_values("task")
    cv.to_csv(os.path.join(out_dir, "fit_log.csv"), index=False)
    bad = cv[cv.cv_verdict != "interior"]
    rep += ["", "## Fits", "",
            f"{len(cv)} fits recorded; {len(cv) - len(bad)} with an interior CV minimum. "
            f"hapc {', '.join(sorted(cv.hapc_version.unique()))}. Total fitting time {cv.elapsed_s.sum() / 3600:.2f} h."]
    if len(bad):
        rep += ["", "Fits whose CV minimum is not interior:", "",
                md_table(bad[["block", "nuisance", "fold", "best_lambda", "cv_index", "cv_verdict"]], "{:.4g}")]
    if "xfit" in missing_all:
        rep += ["", f"Cross-fitted set not assembled: {len(missing_all['xfit'])} out-of-fold files missing."]

    # ---- estimator sets
    # The cross-fitted one-step set is the one reported in the article; the others go to
    # the supplement (decision of 2026-10-12).
    sets = {}
    est_full, T = estimators(d, versions["full"], blocks, plugin=True)
    if "xfit" in versions:
        sets["onestep_xfit"] = estimators(d, versions["xfit"], blocks, plugin=False)[0]["onestep"]
    sets["onestep_full"], sets["plugin_eif"], sets["plugin_known"] = (
        est_full["onestep"], est_full["plugin_eif"], est_full["plugin_known"])
    rep += ["", "Set reported in the article: `onestep_xfit`"
            + ("." if "xfit" in versions else " (**not available in this run**).")
            + " Sets for the supplement: `onestep_full`, `plugin_eif`, `plugin_known`."]

    published = pd.read_csv(PUBLISHED).set_index("label") if os.path.exists(PUBLISHED) else None
    inf_rows, conf_rows, check_fail = [], [], []
    for sname, ests in sets.items():
        psi, Phi, uids = to_matrix(d, T, ests)
        Sigma = Phi.T @ Phi
        v = np.diag(Sigma)
        labels = [e["label"] for e in ests]
        tab = pd.DataFrame(dict(j=np.arange(1, len(ests) + 1), label=labels, block=[C.BLOCKS[e["block"]] for e in ests],
                                psi=psi, se=np.sqrt(v), theta0=[e["theta0"] for e in ests],
                                ess=[e.get("ess", np.nan) for e in ests],
                                max_weight=[e.get("max_weight", np.nan) for e in ests]))
        analyses = {"RCT-only": [j for j, e in enumerate(ests) if e["block"] == 0],
                    "Covariate-adjusted": [j for j, e in enumerate(ests) if e["block"] == 0 or e["kind"] == "adj"],
                    "Full": list(range(len(ests)))}
        for aname, idx in analyses.items():
            if len(idx) < 2 or (aname != "RCT-only" and len(idx) == len(analyses["RCT-only"])):
                continue
            idx = np.array(idx)
            n_units = int((np.abs(Phi[:, idx]).sum(axis=1) > 0).sum())
            r = I.pooled(psi[idx], Phi[:, idx], n_units=n_units)
            r["best_single"] = labels[idx[r.pop("best_single_index")]]
            inf_rows.append(dict(set=sname, analysis=aname, n_persons=n_units, **r))
            if aname == "Full":
                tab["w_eb_pair"] = I.eb_weights(v, I.tau2_pairwise(psi, v))
                tab["w_gls"] = I.gls_weights(Sigma)[0]
            if len(idx) >= 8:
                sc = I.split_conformal(psi[idx], v[idx], v[idx].mean(), B=CONF_B, seed=CONF_SEED)
                lo = I.loo_conformal(psi[idx], v[idx], v[idx].mean())
                conf_rows.append(dict(set=sname, analysis=aname, J=len(idx), split_level=sc["level"],
                                      n_train=sc["n_train"], n_cal=sc["n_cal"], split_lo_first=sc["lo_first"],
                                      split_hi_first=sc["hi_first"], split_lo_median=sc["lo_median"],
                                      split_hi_median=sc["hi_median"], split_lo_sd=sc["lo_sd"], split_hi_sd=sc["hi_sd"],
                                      loo_level=lo["level"], loo_lo=lo["lo"], loo_hi=lo["hi"]))
        if sname.startswith("plugin") and published is not None:
            tab["psi_published"] = published.psi.reindex(labels).to_numpy()
            tab["se_published"] = published.se.reindex(labels).to_numpy()
            tab["gap"] = tab.psi - tab.psi_published
            if sname == "plugin_known":
                rep += ["", "## Agreement with the published table", ""]
                if not publishable:
                    rep += ["Not applicable to this run."]
                else:
                    worst = tab.loc[tab.gap.abs().idxmax()]
                    off = tab[tab.gap.abs() > a.tol]
                    rep += [f"Largest gap {abs(worst.gap):.4g} ({worst.label}); tolerance {a.tol}. "
                            + ("All 15 estimates agree within the tolerance." if off.empty else
                               f"**{len(off)} estimates differ**: "
                               + "; ".join(f"{r.label} {r.psi:.2f} here, {r.psi_published:.2f} published" for r in off.itertuples()))]
                    check_fail = list(off.label)
        tab.to_csv(os.path.join(out_dir, f"functionals_{sname}.csv"), index=False, float_format="%.10g")
        pd.DataFrame(Phi, columns=labels).assign(uid=uids).set_index("uid").to_csv(
            os.path.join(out_dir, f"Phi_{sname}.csv"), float_format="%.8g")
        pd.DataFrame(Sigma, index=labels, columns=labels).to_csv(
            os.path.join(out_dir, f"Sigma_{sname}.csv"), float_format="%.10g")
        rep += ["", f"## Estimators: {sname}", "", md_table(tab.drop(columns=["j"]), "{:.4g}")]

    inf = pd.DataFrame(inf_rows)
    inf.to_csv(os.path.join(out_dir, "inference.csv"), index=False, float_format="%.10g")
    show = ["set", "analysis", "J", "tau_pair", "psi_eb_pair", "se_sandwich_pair", "se_working_pair", "ci_lo_pair",
            "ci_hi_pair", "tau2_over_mean_v", "dispersion_z", "psi_eb_zero", "se_sandwich_zero", "ci_lo_zero", "ci_hi_zero",
            "psi_convex_zero", "se_convex_zero", "ci_lo_convex_zero", "ci_hi_convex_zero", "support_convex_zero", "best_single", "se_best_single",
            "psi_gls", "se_gls", "sigma_rank", "corr_mean"]
    rep += ["", "## Pooled estimates", "",
            "Columns ending in `_zero` set the dispersion to zero: they are the estimates whose sandwich interval is "
            "valid when exact identification is maintained. `psi_convex` is the convex combination (weights "
            "minimising w'(Sigma + t I)w over the simplex); `support_convex` counts its non-zero weights. "
            "`dispersion_z` is the sample variance of the estimates standardised by its mean and standard "
            "deviation under exact identification, computed from Sigma.", "",
            "`se_sandwich` uses the person-level influence matrix; `se_working` is the working-independence "
            "formula of arXiv v2 and is not a valid standard error.", "", md_table(inf[show], "{:.4g}")]
    if conf_rows:
        conf = pd.DataFrame(conf_rows)
        conf.to_csv(os.path.join(out_dir, "conformal.csv"), index=False, float_format="%.10g")
        rep += ["", "## Prediction intervals", "",
                f"Split conformal over {CONF_B} random splits (seed {CONF_SEED}); leave-one-out for comparison.", "",
                md_table(conf, "{:.4g}")]
    with open(os.path.join(out_dir, "report.md"), "w") as f:
        f.write("\n".join(rep) + "\n")
    # The run's MANIFEST describes the fits. This one describes the assembly, which can be
    # redone with later code without refitting; it supersedes any inference setting recorded
    # in the run's MANIFEST (runs prepared before 2026-10-12 record a dispersion threshold
    # that the assembly no longer uses).
    def sh(*cmd):
        try:
            return subprocess.check_output(cmd, cwd=C.ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "unavailable"
    with open(os.path.join(out_dir, "MANIFEST"), "w") as f:
        f.write("\n".join([
            f"assembled_utc={datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
            f"assembly_git_commit={sh('git', 'rev-parse', 'HEAD')}",
            f"assembly_git_dirty={'yes' if sh('git', 'status', '--porcelain', '--untracked-files=no') else 'no'}",
            f"blocks={','.join(C.BLOCKS[b] for b in blocks)}",
            f"sets={','.join(sets)}",
            "article_set=onestep_xfit",
            "inference=pooled estimates with the dispersion set to zero (_zero) and plug-in (_pair, _dl); no threshold",
            "convex=weights minimising w'(Sigma + t I)w over the simplex, SLSQP",
            f"published_table_tolerance={a.tol}",
            "plugin_estimates_differing_from_published="
            + ("not checked (not a publishable run)" if not publishable else ",".join(check_fail) if check_fail else "none"),
            f"conformal_B={CONF_B} conformal_seed={CONF_SEED}",
            f"numpy_version={np.__version__}", f"pandas_version={pd.__version__}"]) + "\n")
    print("\n".join(rep))
    print(f"\nwritten: {out_dir}")
    if check_fail:
        print(f"\nCHECK FAILED: {len(check_fail)} plug-in estimates differ from {PUBLISHED} by more than {a.tol}: "
              + ", ".join(check_fail))
        sys.exit(3)


if __name__ == "__main__":
    main()
