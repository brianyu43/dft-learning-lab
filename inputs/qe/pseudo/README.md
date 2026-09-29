# Silicon pseudopotential for installation verification

`Si_ONCV_PBE-1.1.upf` is copied unchanged from the official Quantum ESPRESSO source tree:

`QEHeat/examples/pseudo/Si_ONCV_PBE-1.1.upf` at the commit recorded in `environment/qe-commit.txt`.

- Source: https://github.com/QEF/q-e/tree/qe-7.6/QEHeat/examples/pseudo
- SHA-256: `d66fa1bb73367f9a2fc34286ecd0876fa1c8c79b25f3b12b30b5b8dc2164148c`
- Header: Si, norm-conserving, PBE, 4 valence electrons. Original PP_INFO and header are retained.
- Purpose: a small executable/SCF installation check. The 30 Ry cutoff and 2×2×2 grid are demonstration settings, not a convergence study.

This is not a selected SSSP production set. Mo/S/Se pseudopotentials and corresponding convergence tests remain part of weeks 3–4. Use fully relativistic potentials when the later SOC calculation requires them.
