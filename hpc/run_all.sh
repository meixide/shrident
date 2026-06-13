#!/bin/bash
# run_all.sh — RECOMMENDED one-liner. Submits the functionals job array, then
# the aggregation job chained with an afterok dependency on the whole array, so
# aggregation starts automatically once (and only if) all 5 blocks succeed.
#
#   cd $LUSTRE/shrinkage && bash hpc/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."          # repo root (so data/ paths resolve)
mkdir -p logs

ARRAY_ID=$(sbatch --parsable hpc/submit_array.sh)
echo "submitted array job:      $ARRAY_ID  (tasks 0-4)"
AGG_ID=$(sbatch --parsable --dependency=afterok:"$ARRAY_ID" hpc/submit_aggregate.sh)
echo "submitted aggregate job:  $AGG_ID  (afterok:$ARRAY_ID)"
echo
echo "monitor with:  squeue -u \$USER        (or: sqstat)"
echo "logs in:       logs/hapc_${ARRAY_ID}_*.{out,err}"
echo "efficiency:    seff ${ARRAY_ID}_1     (CPS-2 task, the long pole)"
