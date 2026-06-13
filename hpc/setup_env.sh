#!/bin/bash
# setup_env.sh — one-time creation of the conda env for the HAPC LaLonde run on
# CESGA FinisTerrae III (FT3). Run this INTERACTIVELY on a login node (it builds
# the hapc C++ extension); do NOT sbatch it.
#
#   bash hpc/setup_env.sh
#
# FT3 notes:
#  * $HOME is only 10 GB and full $HOME crashes services -> conda lives in $STORE.
#  * Use the system miniconda3 module; do NOT install your own miniconda.
#  * hapc compiles a C++ extension (cmake + pybind11) -> needs a compiler.
set -euo pipefail

ENV_NAME="${ENV_NAME:-ft3hapc}"
HAPC_REPO="${HAPC_REPO:-https://github.com/meixide/hapc.git}"

# --- 1. redirect conda storage to $STORE (per FT3 guide) ---------------------
mkdir -p "$STORE/.conda"
if [ -d "$HOME/.conda" ] && [ ! -L "$HOME/.conda" ]; then
  mv "$HOME"/.conda/* "$STORE/.conda/" 2>/dev/null || true
  rm -rf "$HOME/.conda"
fi
[ -L "$HOME/.conda" ] || ln -s "$STORE/.conda" "$HOME/.conda"

# --- 2. load system toolchain (compiler + cmake) and miniconda ---------------
module purge
module load miniconda3
module load gcc cmake 2>/dev/null || module load gcccore cmake 2>/dev/null || true
eval "$(conda shell.bash hook)"

# --- 3. create env and install deps ------------------------------------------
conda create -y -n "$ENV_NAME" python=3.11
conda activate "$ENV_NAME"
pip install --upgrade pip
# runtime deps for the analysis scripts (hapc itself only needs numpy+sklearn)
pip install "numpy>=1.24,<2.3" pandas scikit-learn scipy statsmodels

# --- 4. build & install hapc (C++ extension) ---------------------------------
WORK="${LUSTRE:-$HOME}/src"
mkdir -p "$WORK" && cd "$WORK"
[ -d hapc ] || git clone "$HAPC_REPO" hapc
cd hapc && git pull --ff-only || true
pip install .          # builds the pybind11/cmake extension

# --- 5. sanity check ----------------------------------------------------------
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
