#!/bin/bash
# update_hapc.sh — rebuild hapc in the EXISTING conda env from the current head of
# its GitHub repository. Run INTERACTIVELY on an FT3 login node (it compiles the
# C++ extension); do NOT sbatch it.
#
#   bash hpc/update_hapc.sh
#
# It refuses to run while any of your hapc jobs is queued or running: the jobs
# import hapc from this environment, and one run directory must not mix hapc
# versions (each fit records the version it used in its .json).
# After updating, a new run gets a new directory: hpc/run_nsw_if.sh puts the hapc
# version in the directory name.
#
# This is step 4 and 5 of hpc/setup_env.sh; it does not recreate the environment.
set -euo pipefail

ENV_NAME="${ENV_NAME:-ft3hapc}"
if squeue -h -u "$USER" -o %j | grep -q -E '^(nsw_if|hapc_)'; then
  echo "You have hapc jobs queued or running (squeue -u \$USER). Wait until they finish." >&2
  exit 1
fi

source /opt/cesga/2020/software/Core/miniconda3/22.11.1-1/etc/profile.d/conda.sh
set +u; conda activate "$ENV_NAME"; set -u
python -c "import hapc; print('installed before:', getattr(hapc, '__version__', '?'))"

# Same toolchain settings as setup_env.sh (see the comments there).
export CC="${CC:-$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc}"
export CXX="${CXX:-$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++}"
export CXXFLAGS="${CXXFLAGS:-} -D_ISOC11_SOURCE"
unset LIBRARY_PATH

cd "${LUSTRE:-$HOME}/src/hapc"
git pull --ff-only
echo "building hapc at $(git rev-parse --short HEAD)"
rm -rf build _skbuild *.egg-info
pip install --force-reinstall --no-deps .

SO=$(ls "$CONDA_PREFIX"/lib/python*/site-packages/hapc/hapc_core*.so | head -1)
if objdump -T "$SO" | grep -q "GLIBC_2.2[9-9]\|GLIBC_2.3[0-9]"; then
  echo "ERROR: $SO requires GLIBC>=2.29 (gentoo libm leaked into the link)." >&2
  exit 1
fi
cd /tmp
python - <<'PY'
import hapc, numpy as np
print("installed now:", getattr(hapc, "__version__", "?"))
r = hapc.cv_hapc(np.random.rand(60, 3), np.random.rand(60),
                 family="gaussian", max_degree=2, npcs=10, norm="1",
                 log_lambda_min=-12, log_lambda_max=16, grid_length=10, nfolds=3,
                 predict=np.random.rand(2, 3))
print("smoke test best_lambda:", r.best_lambda)
PY
