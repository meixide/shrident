#!/bin/bash
# nsw_if_array.sh — one nuisance fit of the NSW application per array task.
# Do not sbatch this by hand: hpc/run_nsw_if.sh prepares the run directory and
# submits this script twice, passing --array, --time, -c and the partition for
# the CPS-2 block (long fits) and for the remaining blocks (minutes), with RUN
# exported.
#
# Sizing comes from the seff reports of the published run (logs/hapc_7684731_*):
# every task used one core (6% of 16) and at most 4.4 GB (CPS-2); the other
# blocks stayed under 0.5 GB.
#
#SBATCH -J nsw_if
#SBATCH -o logs/nsw_if_%A_%a.out
#SBATCH -e logs/nsw_if_%A_%a.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=4G
set -euo pipefail

# FT3: source the system conda directly (`module load miniconda3` gives no
# working conda on this cluster). See hpc/setup_env.sh.
source /opt/cesga/2020/software/Core/miniconda3/22.11.1-1/etc/profile.d/conda.sh
set +u   # env carries the conda build toolchain; its activation script trips set -u
conda activate "${ENV_NAME:-ft3hapc}"
set -u

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

echo "host=$(hostname) run=$RUN task=$SLURM_ARRAY_TASK_ID cores=${SLURM_CPUS_PER_TASK:-1}"
srun python -u apps/nsw_fit_nuisances.py --run "$RUN"
