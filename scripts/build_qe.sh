#!/bin/bash
set -euo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$root"
version=7.6
src="$root/.tools/src/qe-$version"
build="$root/.tools/build/qe-$version"
prefix="$root/.tools/qe-$version"
mkdir -p "$root/reports" "$root/.tools/src" "$root/.tools/bin"
if [ ! -f "$src/CMakeLists.txt" ]; then
  git clone --depth 1 --branch "qe-$version" --recurse-submodules --shallow-submodules \
    https://github.com/QEF/q-e.git "$src" > reports/qe-download.log 2>&1
fi
cmake -S "$src" -B "$build" \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-16 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX="$prefix" \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DQE_ENABLE_MPI=OFF -DQE_ENABLE_OPENMP=ON -DQE_ENABLE_TEST=ON \
  -DQE_FFTW_VENDOR=FFTW3 -DFFTW3_ROOT=/opt/homebrew/opt/fftw \
  -DBLA_VENDOR=OpenBLAS \
  -DCMAKE_PREFIX_PATH='/opt/homebrew/opt/openblas;/opt/homebrew/opt/fftw' \
  > reports/qe-configure.log 2>&1
cmake --build "$build" --parallel 6 --target \
  qe_pw_exe qe_pp_exe qe_pp_bands_exe qe_pp_dos_exe qe_pp_projwfc_exe \
  qe_phonon_ph_exe qe_phonon_q2r_exe qe_phonon_matdyn_exe qe_phonon_dynmat_exe \
  > reports/qe-build.log 2>&1
mkdir -p "$prefix/bin"
for name in pw ph pp bands dos projwfc q2r matdyn dynmat; do
  cp "$build/bin/$name.x" "$prefix/bin/$name.x"
  ln -sf "../qe-$version/bin/$name.x" "$root/.tools/bin/$name.x"
done
git -C "$src" rev-parse HEAD > environment/qe-commit.txt
git -C "$src" submodule status > environment/qe-submodules.txt
echo "Installed Quantum ESPRESSO $version in $prefix/bin"
