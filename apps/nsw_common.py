"""NSW application with stored nuisance predictions: shared definitions.

An earlier run (whose estimates are kept in data/lalonde_functionals_hapc.csv) kept one
line per estimator and nothing per unit, so the covariance between estimators cannot be recovered from it.
Here every nuisance fit stores its per-row predictions, and the estimators, their
unit-level influence contributions and the inference are assembled afterwards without
refitting.

A *run* is a directory results/<name>/ holding
    settings.json, MANIFEST   what was run (written once by nsw_prepare.py)
    units.csv                 one line per data row: row, uid (person), Z, A, stratum, fold
    tasks.csv                 one line per nuisance fit: task, block, nuisance, fold
    nuis/                     per-task predictions and fit logs   (nsw_fit_nuisances.py)
    assembled/                estimators, influence matrix, inference (nsw_assemble.py)
"""
import json
import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data", "lalonde_combined.csv")
COVS = ["age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75"]

# block id = Z code of the comparison sample (0 is the experimental block)
BLOCKS = {0: "RCT", 1: "CPS-2", 2: "CPS-3", 3: "PSID-2", 4: "PSID-3"}
NESTED = {2: 1, 4: 3}   # CPS-3 is a subsample of CPS-2, PSID-3 of PSID-2: same persons
TAG = {1: "cps2", 2: "cps3", 3: "psid2", 4: "psid3"}

# Base learner of the published run. These values are copied into settings.json when a
# run is prepared; the fits read them from there, never from this dictionary.
HAPC = dict(max_degree=2, norm="1", approx=False, npcs_exponent=0.6,
            log_lambda_min=-12.0, log_lambda_max=16.0, grid_length=40, nfolds=5)
SMOKE = dict(grid_length=8, nfolds=3)   # coarse grid for local pipeline tests only
CLIP = 1e-6                             # propensities are clipped to [CLIP, 1 - CLIP]


# ---------------------------------------------------------------- units
def build_units(K, seed, frac=1.0):
    """Map data rows to persons and give each person a cross-fitting fold.

    The external samples are nested, so one person can sit in two of them. Rows of the
    child sample are matched to rows of the parent on all recorded variables; identical
    records are paired in order of appearance (multiset matching).
    """
    d = pd.read_csv(DATA)
    d["row"] = np.arange(len(d))
    d["uid"] = d["row"]
    key = COVS + ["Y", "A"]
    for child, parent in NESTED.items():
        c, p = d[d.Z == child].copy(), d[d.Z == parent].copy()
        c["occ"], p["occ"] = c.groupby(key).cumcount(), p.groupby(key).cumcount()
        m = c.merge(p[key + ["occ", "row"]], on=key + ["occ"], how="left", suffixes=("", "_parent"))
        if m.row_parent.isna().any():
            raise SystemExit(f"{int(m.row_parent.isna().sum())} rows of {BLOCKS[child]} "
                             f"have no match in {BLOCKS[parent]}: the samples are not nested")
        d.loc[c.index, "uid"] = m.row_parent.to_numpy().astype(int)

    child_of = {parent: child for child, parent in NESTED.items()}
    in_child = {parent: set(d.loc[d.Z == child, "uid"]) for child, parent in NESTED.items()}
    def stratum(r):
        if r.Z == 0:
            return "nsw_treated" if r.A == 1 else "nsw_control"
        parent = NESTED.get(r.Z, r.Z)
        return f"{TAG[parent]}_in_{TAG[child_of[parent]]}" if r.uid in in_child[parent] else f"{TAG[parent]}_only"
    d["stratum"] = d.apply(stratum, axis=1)

    rng = np.random.default_rng(seed)
    persons = d.drop_duplicates("uid")[["uid", "stratum"]].sort_values("uid")
    keep, fold = [], {}
    for _, g in persons.groupby("stratum", sort=True):
        ids = rng.permutation(g.uid.to_numpy())
        if frac < 1.0:
            ids = ids[:max(12, int(round(frac * len(ids))))]
        keep.extend(ids)
        fold.update({u: (1 + i % K if K > 0 else 0) for i, u in enumerate(ids)})
    d = d[d.uid.isin(set(keep))].copy()
    d["fold"] = d.uid.map(fold).astype(int)
    return d[["row", "uid", "Z", "A", "stratum", "fold"]]


def load(run):
    """Units of a run joined with the data, in data order (the order of the published fits)."""
    u = pd.read_csv(os.path.join(run, "units.csv"))
    d = pd.read_csv(DATA)
    d["row"] = np.arange(len(d))
    return u.merge(d[["row", "Y"] + COVS], on="row", how="left", validate="one_to_one").sort_values("row")


def settings(run):
    with open(os.path.join(run, "settings.json")) as f:
        return json.load(f)


# ---------------------------------------------------------------- nuisances
def nuisances(d, block):
    """name -> (family, training rows, label, prediction rows) for one block.

    Training rows are stacked exactly as in the published run, because the CV folds of
    cv_hapc depend on row order. The first two fits of each block are the published
    ones; m_ec and g_trt are the second nuisance of each doubly robust estimator.
    """
    nsw = d[d.Z == 0]
    trt, ctl = nsw[nsw.A == 1], nsw[nsw.A == 0]
    if block == 0:
        return {"e_rct": ("binomial", nsw, nsw.A, nsw),                     # P(A=1 | X) in the experiment
                "m_rct": ("gaussian", ctl, ctl.Y, pd.concat([trt, ctl]))}   # E(Y | X) among experimental controls
    ec = d[d.Z == block]
    pool, allk = pd.concat([ctl, ec]), pd.concat([nsw, ec])
    return {"m_pool": ("gaussian", pool, pool.Y, pd.concat([trt, pool])),   # E(Y | X) among pooled controls
            "pi_s": ("binomial", allk, (allk.Z == 0).astype(float), allk),  # P(experiment | X)
            "m_ec": ("gaussian", ec, ec.Y, allk),                           # E(Y | X) in the external sample
            "g_trt": ("binomial", allk, allk.A, allk)}                      # P(treated | X), treated vs pooled controls


def task_table(K):
    rows = [(b, name, f) for b in BLOCKS for name in nuisances_names(b) for f in range(K + 1)]
    return pd.DataFrame(rows, columns=["block", "nuisance", "fold"]).rename_axis("task").reset_index()


def nuisances_names(block):
    return ["e_rct", "m_rct"] if block == 0 else ["m_pool", "pi_s", "m_ec", "g_trt"]


def nuis_path(run, block, name, fold, ext="csv"):
    return os.path.join(run, "nuis", f"b{block}_{name}_f{fold}.{ext}")
