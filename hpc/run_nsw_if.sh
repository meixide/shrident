#!/bin/bash
# run_nsw_if.sh — one-liner for the NSW application with stored nuisance
# predictions. Run on an FT3 login node from the repository:
#
#   bash hpc/run_nsw_if.sh            # K=5 cross-fitting folds, fold seed 20261012
#   bash hpc/run_nsw_if.sh 5 777      # another fold seed -> another run directory
#   bash hpc/run_nsw_if.sh 0          # full-sample fits only (no cross-fitting)
#
# It creates results/nsw_if_hapc_K<K>_s<SEED>/ (refusing to reuse a directory),
# submits one array for the CPS-2 block and one for the other blocks, and chains
# the assembly job on both.
set -euo pipefail
cd "$(dirname "$0")/.."          # repo root
K=${1:-5}
SEED=${2:-20261012}
RUN=results/nsw_if_hapc_K${K}_s${SEED}
mkdir -p logs results

source /opt/cesga/2020/software/Core/miniconda3/22.11.1-1/etc/profile.d/conda.sh
set +u; conda activate "${ENV_NAME:-ft3hapc}"; set -u

python apps/nsw_prepare.py --run "$RUN" --K "$K" --seed "$SEED"

# tasks.csv columns: task,block,nuisance,fold. Block 1 is CPS-2.
BIG=$(awk -F, 'NR>1 && $2==1 {print $1}' "$RUN/tasks.csv" | paste -sd, -)
SMALL=$(awk -F, 'NR>1 && $2!=1 {print $1}' "$RUN/tasks.csv" | paste -sd, -)

# CPS-2: the two published fits took 8.25 h together on one core; each task here is one fit.
A1=$(sbatch --parsable --export=ALL,RUN="$RUN" --array="$BIG" -c 2 -t 24:00:00 -p medium hpc/nsw_if_array.sh)
# Other blocks: the slowest published block (CPS-3, two fits) took 15 minutes.
A2=$(sbatch --parsable --export=ALL,RUN="$RUN" --array="$SMALL" -c 1 -t 02:00:00 -p short hpc/nsw_if_array.sh)
AGG=$(sbatch --parsable --export=ALL,RUN="$RUN" --dependency=afterok:"$A1":"$A2" hpc/nsw_if_assemble.sh)

echo
echo "run directory:   $RUN"
echo "CPS-2 array:     $A1   ($(echo "$BIG" | tr ',' '\n' | wc -l | tr -d ' ') tasks, 24 h limit)"
echo "other blocks:    $A2   ($(echo "$SMALL" | tr ',' '\n' | wc -l | tr -d ' ') tasks, 2 h limit)"
echo "assembly:        $AGG   (afterok:$A1:$A2)"
echo
echo "monitor:   squeue -u \$USER"
echo "progress:  ls $RUN/nuis/*.csv | wc -l      (of $(($(wc -l < "$RUN/tasks.csv") - 1)))"
echo "report:    $RUN/assembled/report.md"
echo "a failed task: resubmit its array index; finished tasks are skipped"
