# Source this file: source scripts/activate.sh
if [ -n "${ZSH_VERSION:-}" ]; then
  DFT_ROOT="$(cd -- "$(dirname -- "${(%):-%N}")/.." && pwd)"
else
  DFT_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
fi
export DFT_ROOT
export PATH="$DFT_ROOT/.tools/bin:$DFT_ROOT/.tools/xtb/bin:$PATH"
export UV_CACHE_DIR="$DFT_ROOT/.cache/uv"
export MPLCONFIGDIR="$DFT_ROOT/.cache/matplotlib"
export JUPYTER_RUNTIME_DIR="$DFT_ROOT/.cache/jupyter-runtime"
export IPYTHONDIR="$DFT_ROOT/.cache/ipython"
export OMP_NUM_THREADS="${DFT_THREADS:-4}"
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export MKL_NUM_THREADS=1
if [ -f "$DFT_ROOT/environment/orca.local.sh" ]; then
  source "$DFT_ROOT/environment/orca.local.sh"
fi
source "$DFT_ROOT/.venv/bin/activate"
