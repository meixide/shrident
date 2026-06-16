#!/bin/bash
# submit_aggregate.sh — runs AFTER the job array finishes. Merges the 5 block
# CSVs into the full functionals table, computes the EB sub-analyses + LOO
# scores, the conformal prediction intervals, and (optionally) the OLS/logistic
# baseline on the same full data for an apples-to-apples comparison.
#
# Submit with a dependency on the array job (see run_all.sh), or by hand:
#   sbatch --dependency=afterok:<ARRAY_JOBID> hpc/submit_aggregate.sh
#
#SBATCH -J hapc_agg
#SBATCH -o logs/hapc_agg_%j.out
#SBATCH -e logs/hapc_agg_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4G
#SBATCH -t 02:00:00
#SBATCH -p short
set -euo pipefail

# FT3: source the system conda directly (`module load miniconda3` gives no
# working conda on this cluster). See hpc/setup_env.sh.
source /opt/cesga/2020/software/Core/miniconda3/22.11.1-1/etc/profile.d/conda.sh
set +u   # env carries the conda build toolchain; its activation script trips set -u
conda activate "${ENV_NAME:-ft3hapc}"
set -u

export HAPC_FRAC=1.0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

python -u data/lalonde_aggregate_hapc.py      # -> functionals + EB + LOO CSVs
python -u data/lalonde_conformal_hapc.py       # -> prediction-interval CSV
python -u data/lalonde_baseline_subsample.py   # OLS/logistic on the SAME full data

# Copy the publishable outputs to permanent storage ($STORE has backups).
mkdir -p "$STORE/shrinkage_results"
cp data/lalonde_functionals_hapc.csv data/lalonde_eb_hapc.csv \
   data/lalonde_loo_hapc.csv data/lalonde_pi_hapc.csv \
   data/lalonde_functionals_ols_sub.csv "$STORE/shrinkage_results/" 2>/dev/null || true
echo "AGGREGATION DONE. Results in data/ and \$STORE/shrinkage_results/"
