#!/usr/bin/env python3
"""lalonde_aggregate_hapc.py

Aggregates the 5 per-block functional CSVs produced by a job array of
`lalonde_functionals_hapc.py` (HAPC_BLOCK=0..4) into the full 15-functional
table, then computes the EB sub-analyses and LOO conformal scores using the
same code path as the single-process run. Run after all array tasks finish.

    lalonde_functionals_hapc_block0..4.csv  ->  lalonde_functionals_hapc.csv
                                                 lalonde_eb_hapc.csv
                                                 lalonde_loo_hapc.csv
"""
import os
import sys
import pandas as pd
from lalonde_functionals_hapc import compute_eb_loo  # same dir on sys.path

DIR = os.environ.get("HAPC_OUTDIR", os.path.dirname(os.path.abspath(__file__)))

def main():
    parts = []
    for b in range(5):
        f = f"{DIR}/lalonde_functionals_hapc_block{b}.csv"
        if not os.path.exists(f):
            sys.exit(f"ERROR: missing block file {f} -- did array task {b} finish?")
        parts.append(pd.read_csv(f))
    res = pd.concat(parts).sort_values("j").reset_index(drop=True)
    if len(res) != 15:
        sys.exit(f"ERROR: expected 15 functionals, got {len(res)}")
    blocks, loo = compute_eb_loo(res, DIR)
    pd.set_option("display.width", 200, "display.max_columns", 20)
    print("==== FUNCTIONALS (aggregated from job array) ====")
    print(res[["j", "label", "psi", "se", "theta0"]].to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print("\n==== EB sub-analyses ====")
    print(blocks.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print(f"\nLOO scores: mean={loo.S_loo.mean():.3f}  sd={loo.S_loo.std(ddof=1):.3f}")
    print("\nSaved: lalonde_functionals_hapc.csv, lalonde_eb_hapc.csv, lalonde_loo_hapc.csv")

if __name__ == "__main__":
    main()
