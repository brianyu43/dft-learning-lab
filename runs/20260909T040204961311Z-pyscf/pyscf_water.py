"""Real PBE/def2-SVP water SCF, for installation verification only."""
import json
import os
from pathlib import Path
from pyscf import dft, gto, lib

lib.num_threads(int(os.environ.get("OMP_NUM_THREADS", "4")))
mol = gto.M(atom="O 0 0 0; H 0 -0.757 0.587; H 0 0.757 0.587",
            basis="def2-svp", unit="Angstrom", charge=0, spin=0,
            max_memory=2000, verbose=4)
mf = dft.RKS(mol)
mf.xc = "PBE"
mf.grids.level = 3
mf.conv_tol = 1e-9
mf.max_cycle = 100
energy = mf.kernel()
result = {"method": "PBE", "basis": "def2-SVP", "solvent": "vacuum",
          "charge": 0, "multiplicity": 1, "atom_count": 3,
          "scf_converged": bool(mf.converged), "energy_hartree": float(energy),
          "dipole_debye": mf.dip_moment(unit="Debye").tolist(),
          "geometry_optimization": "not_run", "minimum_verification": "not_run"}
Path("result.json").write_text(json.dumps(result, indent=2) + "\n")
print("DFT_LAB_RESULT " + json.dumps(result))
if not mf.converged:
    raise SystemExit(1)
