#!/bin/bash
# submit_single.sh — SIMPLEST fallback. Runs the whole pipeline (all 15
# functionals + EB + LOO) in ONE job, sequentially, then conformal + baseline.
# Robust and zero-orchestration; wall time ~ sum of all fits (a few hours).
# Prefer hpc/run_all.sh (job array) when wall time matters.
#
#   cd $LUSTRE/shrinkage && sbatch hpc/submit_single.sh
#
#SBATCH -J hapc_all
#SBATCH -o logs/hapc_all_%j.out
#SBATCH -e logs/hapc_all_%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --mem-per-cpu=4G            # 32 x 4G = 128 GB (very safe)
#SBATCH -t 1-00:00:00              # 1 day
#SBATCH -p medium
set -euo pipefail

module purge
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME:-ft3hapc}"

export HAPC_FRAC=1.0               # full data
unset HAPC_BLOCK                   # all blocks in one process
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}

srun python -u data/lalonde_functionals_hapc.py   # functionals + EB + LOO
python -u data/lalonde_conformal_hapc.py           # prediction intervals
python -u data/lalonde_baseline_subsample.py       # OLS/logistic baseline

mkdir -p "$STORE/shrinkage_results"
cp data/lalonde_functionals_hapc.csv data/lalonde_eb_hapc.csv \
   data/lalonde_loo_hapc.csv data/lalonde_pi_hapc.csv \
   data/lalonde_functionals_ols_sub.csv "$STORE/shrinkage_results/" 2>/dev/null || true
echo "SINGLE-JOB RUN DONE."
