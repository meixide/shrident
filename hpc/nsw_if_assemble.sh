#!/bin/bash
# nsw_if_assemble.sh — runs after both nuisance arrays finish. Builds the
# estimators, the person-level influence matrix and the inference tables of the
# run, copies the run directory to $STORE, and exits with the status of the
# assembly (3 = the plug-in estimates differ from the published table).
# Submitted by hpc/run_nsw_if.sh; by hand:
#   sbatch --export=ALL,RUN=results/<name> hpc/nsw_if_assemble.sh
#
#SBATCH -J nsw_if_agg
#SBATCH -o logs/nsw_if_agg_%j.out
#SBATCH -e logs/nsw_if_agg_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu=4G
#SBATCH -t 00:30:00
#SBATCH -p short
set -uo pipefail

source /opt/cesga/2020/software/Core/miniconda3/22.11.1-1/etc/profile.d/conda.sh
set +u
conda activate "${ENV_NAME:-ft3hapc}"
set -u

python -u apps/nsw_assemble.py --run "$RUN"
rc=$?

# Permanent copy ($STORE has backups); the scratch filesystem does not.
mkdir -p "$STORE/shrinkage_results"
cp -r "$RUN" "$STORE/shrinkage_results/" && echo "copied $RUN to \$STORE/shrinkage_results/"
echo "assembly exit status: $rc"
exit $rc
