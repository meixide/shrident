#!/bin/bash
# submit_array.sh — RECOMMENDED. Computes the 15 functionals as a 5-task SLURM
# job array (one task per block: 0=RCT, 1=CPS-2, 2=CPS-3, 3=PSID-2, 4=PSID-3).
# Each task is multi-threaded so numpy's BLAS parallelises the HAPC eigendecomp.
# Block 1 (CPS-2) is the long pole (two large fits, n~2.6-2.8k); the others
# finish in minutes. Resources are sized for that worst case.
#
#   cd $LUSTRE/shrinkage && sbatch hpc/submit_array.sh
#
#SBATCH -J hapc_fn
#SBATCH -o logs/hapc_%A_%a.out      # %A=array job id, %a=task id
#SBATCH -e logs/hapc_%A_%a.err
#SBATCH --array=0-4
#SBATCH -N 1                         # one node per task
#SBATCH -n 1                         # one process per task
#SBATCH -c 16                        # cores per task (BLAS threads)
#SBATCH --mem-per-cpu=4G             # 16 x 4G = 64 GB per task (safe for CPS-2)
#SBATCH -t 12:00:00                  # walltime per task (short=6h would also fit)
#SBATCH -p medium
set -euo pipefail

module purge
module load miniconda3
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME:-ft3hapc}"

# Full data, no subsample. Match BLAS thread count to the allocation.
export HAPC_FRAC=1.0
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}
export HAPC_BLOCK=${SLURM_ARRAY_TASK_ID}

echo "host=$(hostname) block=$HAPC_BLOCK cores=$SLURM_CPUS_PER_TASK"
srun python -u data/lalonde_functionals_hapc.py
