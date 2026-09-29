#!/bin/bash
# Isolated MPI environment; the serial ORCA driver launches its MPI modules.
set -euo pipefail
DFT_PARALLEL_ROOT="$(cd -- "$(dirname -- "$0")/.." && pwd)"
export PATH="$DFT_PARALLEL_ROOT/.tools/orca-mpi/bin:$PATH"
export DYLD_LIBRARY_PATH="$DFT_PARALLEL_ROOT/.tools/orca-mpi/lib${DYLD_LIBRARY_PATH:+:$DYLD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
if [ "$#" -ne 1 ]; then
  echo 'Usage: orca_parallel.sh INPUT.inp' >&2
  exit 2
fi
exec "$DFT_PARALLEL_ROOT/.tools/orca-6.1.1/orca" "$1" '--bind-to none'
