#!/usr/bin/env python3
"""Create a run directory for the NSW application: units, folds, task list, MANIFEST.

    python apps/nsw_prepare.py --run results/nsw_if_hapc_K5_s20261012 --K 5 --seed 20261012

Refuses to write into an existing directory: one run, one directory. Everything a fit
or the assembly needs is read from the files written here, so no published number
depends on a command-line flag or an environment variable that nobody recorded.
"""
import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nsw_common as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--K", type=int, default=5, help="cross-fitting folds (0 = full-sample fits only)")
    ap.add_argument("--seed", type=int, default=20261012, help="seed of the fold assignment")
    ap.add_argument("--frac", type=float, default=1.0, help="person subsample; below 1 only for pipeline tests")
    ap.add_argument("--smoke", action="store_true", help="coarse lambda grid; only for pipeline tests")
    a = ap.parse_args()
    if os.path.exists(a.run):
        raise SystemExit(f"{a.run} exists. Choose a new name; a run directory is never reused.")
    os.makedirs(os.path.join(a.run, "nuis"))

    units = C.build_units(a.K, a.seed, a.frac)
    units.to_csv(os.path.join(a.run, "units.csv"), index=False)
    tasks = C.task_table(a.K)
    tasks.to_csv(os.path.join(a.run, "tasks.csv"), index=False)

    hapc_cfg = dict(C.HAPC, **(C.SMOKE if a.smoke else {}))
    s = dict(K=a.K, fold_seed=a.seed, frac=a.frac, smoke=a.smoke, clip=C.CLIP, hapc=hapc_cfg)
    with open(os.path.join(a.run, "settings.json"), "w") as f:
        json.dump(s, f, indent=1)

    def sh(*cmd):
        try:
            return subprocess.check_output(cmd, cwd=C.ROOT, stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            return "unavailable"
    try:
        import hapc
        hapc_version = getattr(hapc, "__version__", "unknown")
    except Exception as e:                                   # prepared on a machine without hapc
        hapc_version = f"not importable here ({type(e).__name__})"
    with open(C.DATA, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    counts = units.drop_duplicates("uid").groupby("stratum").size()
    lines = [
        f"run={os.path.basename(os.path.normpath(a.run))}",
        f"created_utc={datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        f"host={platform.node()}",
        f"git_commit={sh('git', 'rev-parse', 'HEAD')}",
        f"git_dirty={'yes' if sh('git', 'status', '--porcelain', '--untracked-files=no') else 'no'}",
        f"data=data/lalonde_combined.csv sha256={sha}",
        f"rows={len(units)} persons={units.uid.nunique()}",
        "strata=" + ",".join(f"{k}:{v}" for k, v in counts.items()),
        f"crossfit_folds={a.K}",
        f"fold_seed={a.seed}",
        "fold_assignment=persons permuted within stratum, fold = 1 + position mod K; recorded in units.csv",
        f"person_subsample_frac={a.frac}",
        f"smoke_grid={'yes (NOT a publishable run)' if a.smoke or a.frac < 1 else 'no'}",
        "learner=hapc.cv_hapc",
        *[f"hapc_{k}={v}" for k, v in hapc_cfg.items()],
        "hapc_npcs=floor(n_train ** npcs_exponent)",
        f"propensity_clip={C.CLIP}",
        f"hapc_version={hapc_version}",
        f"numpy_version={np.__version__}",
        f"pandas_version={pd.__version__}",
        f"python_version={platform.python_version()}",
        f"tasks={len(tasks)}",
        "intervals=sandwich at dispersion zero (exact identification maintained) and plug-in versions; apps/inference.py",
        "conformal_B=2000 conformal_seed=42; set in apps/nsw_assemble.py",
    ]
    with open(os.path.join(a.run, "MANIFEST"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\nprepared {a.run}: {len(tasks)} tasks (indices 0-{len(tasks) - 1})")


if __name__ == "__main__":
    main()
