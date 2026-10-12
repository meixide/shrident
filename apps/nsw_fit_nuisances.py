#!/usr/bin/env python3
"""Fit ONE nuisance (one line of tasks.csv) and store its per-row predictions.

    python apps/nsw_fit_nuisances.py --run results/<name> --task 7

Under SLURM the task index defaults to SLURM_ARRAY_TASK_ID. Fold 0 is the full-sample
fit (the published one for e_rct, m_rct, m_pool, pi_s) and stores predictions on every
row the estimators need. Fold f >= 1 trains on the persons outside fold f and stores
predictions for the persons inside it only (out-of-fold predictions for cross-fitting).

A task whose output already exists is skipped, so a failed array can be resubmitted.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nsw_common as C


def cv_verdict(err):
    """Where the CV curve attains its minimum (same rule as the published run)."""
    err = np.asarray(err, float)
    imin = int(np.argmin(err))
    span = err.max() - err.min()
    plateau = span > 0 and (err[-1] - err.min()) <= 1e-6 * max(abs(err.min()), 1.0)
    pos = ("boundary-low" if imin == 0 else "boundary-high" if imin == len(err) - 1
           else "plateau-high (null model)" if plateau else "interior")
    return imin, pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--task", type=int, default=int(os.environ.get("SLURM_ARRAY_TASK_ID", -1)))
    a = ap.parse_args()
    task = pd.read_csv(os.path.join(a.run, "tasks.csv")).set_index("task").loc[a.task]
    block, name, fold = int(task.block), task.nuisance, int(task.fold)
    out = C.nuis_path(a.run, block, name, fold)
    if os.path.exists(out):
        print(f"task {a.task} already done: {out}")
        return

    import hapc
    s = C.settings(a.run)
    h = s["hapc"]
    d = C.load(a.run)
    family, train, label, pred = C.nuisances(d, block)[name]
    if fold > 0:
        keep = (train.fold != fold).to_numpy()
        train, label = train[keep], label[keep]
        pred = pred[pred.fold == fold]
    n = len(train)
    npcs = max(1, int(np.floor(n ** h["npcs_exponent"])))
    print(f"task {a.task}: block {block} ({C.BLOCKS[block]}) {name} fold {fold} | {family} "
          f"n_train={n} n_pred={len(pred)} npcs={npcs}", flush=True)

    t0 = time.time()
    r = hapc.cv_hapc(train[C.COVS].to_numpy(float), np.asarray(label, float), family=family,
                     max_degree=h["max_degree"], npcs=npcs, norm=h["norm"], approx=h["approx"],
                     predict=pred[C.COVS].to_numpy(float), log_lambda_min=h["log_lambda_min"],
                     log_lambda_max=h["log_lambda_max"], grid_length=h["grid_length"], nfolds=h["nfolds"])
    elapsed = time.time() - t0
    p = np.asarray(r.predictions, float).ravel()
    if len(p) != len(pred) or not np.all(np.isfinite(p)):
        raise SystemExit(f"task {a.task}: {len(p)} predictions for {len(pred)} rows, "
                         f"{int((~np.isfinite(p)).sum())} not finite")
    imin, pos = cv_verdict(r.mses)
    meta = dict(task=int(a.task), block=block, nuisance=name, fold=fold, family=family, n_train=n,
                n_pred=len(pred), npcs=npcs, best_lambda=float(r.best_lambda), cv_index=imin,
                cv_grid_last=len(r.mses) - 1, cv_verdict=pos, elapsed_s=round(elapsed, 1),
                hapc_version=getattr(hapc, "__version__", "unknown"), numpy_version=np.__version__,
                threads=os.environ.get("OMP_NUM_THREADS", "unset"), host=os.uname().nodename)
    with open(C.nuis_path(a.run, block, name, fold, "json"), "w") as f:
        json.dump(meta, f, indent=1)
    tmp = out + ".tmp"
    pd.DataFrame({"row": pred.row.to_numpy(), "pred": p}).to_csv(tmp, index=False, float_format="%.12g")
    os.replace(tmp, out)   # the csv appears only when complete
    print(f"task {a.task} done in {elapsed:.1f}s: lambda={r.best_lambda:.4g} cv {imin}/{len(r.mses) - 1} "
          f"-> {pos}{'' if pos == 'interior' else '  *** NOT INTERIOR ***'}", flush=True)


if __name__ == "__main__":
    main()
