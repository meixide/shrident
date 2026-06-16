#!/bin/bash
# setup_env.sh — one-time creation of the conda env for the HAPC LaLonde run on
# CESGA FinisTerrae III (FT3). Run this INTERACTIVELY on a login node (it builds
# the hapc C++ extension); do NOT sbatch it.
#
#   bash hpc/setup_env.sh
#
# FT3 notes:
#  * $HOME is only 10 GB and full $HOME crashes services -> conda lives in $STORE.
#  * Use the SYSTEM miniconda; do NOT install your own. `module load miniconda3`
#    is broken here (loads 4.11.0, no working conda) -> we source it directly.
#  * hapc compiles a C++/pybind11 extension. It MUST be built with the conda
#    compiler toolchain, NOT the gentoo cesga/2020 gcc: that gcc links against
#    glibc 2.31, so its .so requires GLIBC_2.29, which the system runtime glibc
#    (2.28, RHEL 8.4) the conda python uses does not provide -> ImportError.
set -euo pipefail

ENV_NAME="${ENV_NAME:-ft3hapc}"
HAPC_REPO="${HAPC_REPO:-https://github.com/meixide/hapc.git}"

# --- 1. conda env storage -----------------------------------------------------
# This account's ~/.condarc already directs envs to $STORE/conda/envs and pkgs
# to $STORE/conda/pkgs (Store has the big quota); no $HOME/.conda symlink needed.
mkdir -p "$STORE/conda/envs" "$STORE/conda/pkgs"

# --- 2. activate the system conda (FT3-specific) -----------------------------
# Plain `module load miniconda3` loads 4.11.0 and yields NO working conda; source
# the system miniconda 22.11.1-1 directly (proven working: conda 23.7.4).
source /opt/cesga/2020/software/Core/miniconda3/22.11.1-1/etc/profile.d/conda.sh

# --- 3. create env and install runtime deps ----------------------------------
conda create -y -n "$ENV_NAME" python=3.11
set +u; conda activate "$ENV_NAME"; set -u
pip install --upgrade pip
# runtime deps for the analysis scripts (hapc itself only needs numpy+sklearn)
pip install "numpy>=1.24,<2.3" pandas scikit-learn scipy statsmodels

# --- 4. build & install hapc with the CONDA compiler toolchain ----------------
# Use conda's gcc 12 + a glibc-2.17 sysroot. (The compilers' default 2.12 sysroot
# predates timespec_get and fails to compile pybind11's <chrono>; gcc 15's
# libstdc++ trips the same header, hence the gcc 12 pin.) One combined install to
# avoid a second slow classic-solver pass.
conda install -y -n "$ENV_NAME" -c conda-forge \
  'gcc_linux-64=12' 'gxx_linux-64=12' 'sysroot_linux-64=2.17'
set +u; conda activate "$ENV_NAME"; set -u     # re-activate to load CC/CXX from the toolchain
export CC="${CC:-$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc}"
export CXX="${CXX:-$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++}"
# Expose C11 timespec_get: the 2.17 sysroot's pre-C++17 <features.h> only
# declares it under _ISOC11_SOURCE.
export CXXFLAGS="${CXXFLAGS:-} -D_ISOC11_SOURCE"
# Drop the gentoo cesga/2020 lib64 (glibc 2.31) from the linker search, else
# exp/log bind to GLIBC_2.29; the conda 2.17 sysroot libm must win instead.
unset LIBRARY_PATH

WORK="${LUSTRE:-$HOME}/src"
mkdir -p "$WORK" && cd "$WORK"
[ -d hapc ] || git clone "$HAPC_REPO" hapc
cd hapc && git pull --ff-only || true
rm -rf build _skbuild *.egg-info                # fresh cmake configure
pip install --force-reinstall --no-deps .       # builds the pybind11/cmake extension

# --- 5. sanity check (glibc baseline + import + tiny fit) ---------------------
SO=$(ls "$CONDA_PREFIX"/lib/python*/site-packages/hapc/hapc_core*.so | head -1)
if objdump -T "$SO" | grep -q "GLIBC_2.2[9-9]\|GLIBC_2.3[0-9]"; then
  echo "ERROR: $SO requires GLIBC>=2.29 (gentoo libm leaked into the link)." >&2
  objdump -T "$SO" | grep -oE "GLIBC_[0-9.]+" | sort -V | uniq | tail -3 >&2
  exit 1
fi
python - <<'PY'
import hapc, numpy as np
print("hapc", getattr(hapc, "__version__", "?"), "OK")
r = hapc.cv_hapc(np.random.rand(60, 3), np.random.rand(60),
                 family="gaussian", max_degree=2, npcs=10, norm="1",
                 log_lambda_min=-12, log_lambda_max=16, grid_length=10, nfolds=3,
                 predict=np.random.rand(2, 3))
print("smoke test best_lambda:", r.best_lambda)
PY
echo "ENV '$ENV_NAME' READY."
