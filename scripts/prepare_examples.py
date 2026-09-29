"""Generate initial geometries; these are not DFT-optimized structures."""
import hashlib
import json
from pathlib import Path

from ase.build import bulk, mx2
from ase.io import write
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
MOLECULES = {
    "acetaminophen": "CC(=O)Nc1ccc(O)cc1",
    "nicotinamide": "NC(=O)c1cccnc1",
    "isonicotinamide": "NC(=O)c1ccncc1",
    "urea": "NC(N)=O",
    "succinic_acid": "O=C(O)CCC(=O)O",
    "fumaric_acid": "O=C(O)/C=C/C(=O)O",
}


def orca_input(xyz: str, keywords: str) -> str:
    coords = "\n".join(xyz.splitlines()[2:])
    return f"! {keywords}\n%maxcore 2000\n\n* xyz 0 1\n{coords}\n*\n"


def main():
    records = []
    out = ROOT / "structures/molecules"
    out.mkdir(parents=True, exist_ok=True)
    for name, smiles in MOLECULES.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        params = AllChem.ETKDGv3()
        params.randomSeed = 20260909
        params.numThreads = 1
        params.pruneRmsThresh = 0.35
        cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params))
        if not cids or not AllChem.MMFFHasAllMoleculeParams(mol):
            raise RuntimeError(f"Embedding or MMFF parameters missing: {name}")
        energies = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=1000)
        valid = [(i, e) for i, (status, e) in enumerate(energies) if status == 0]
        if not valid:
            raise RuntimeError(f"No converged MMFF conformer: {name}")
        best_index, energy = min(valid, key=lambda pair: pair[1])
        best_id = cids[best_index]
        xyz = Chem.MolToXYZBlock(mol, confId=best_id)
        (out / f"{name}.xyz").write_text(xyz)
        writer = Chem.SDWriter(str(out / f"{name}_initial_conformers.sdf"))
        for cid, (status, en) in zip(cids, energies):
            mol.SetProp("MMFF_status", str(status))
            mol.SetProp("MMFF_energy_kcal_mol", str(en))
            writer.write(mol, confId=cid)
        writer.close()
        records.append({"name": name, "smiles": smiles,
                        "canonical_smiles": Chem.MolToSmiles(Chem.RemoveHs(mol)),
                        "formula": rdMolDescriptors.CalcMolFormula(mol),
                        "atom_count": mol.GetNumAtoms(), "charge": 0, "multiplicity": 1,
                        "seed": 20260909, "requested_conformers": 20,
                        "embedded_conformers": len(cids), "prune_rms_angstrom": 0.35,
                        "selected_id": best_id, "MMFF_energy_kcal_mol": energy,
                        "geometry_status": "MMFF only; DFT optimization not performed",
                        "sha256": hashlib.sha256(xyz.encode()).hexdigest()})
        if name == "acetaminophen":
            for suffix, keywords in (
                ("vacuum_opt", "r2SCAN-3c TightSCF Opt"),
                ("water_opt", "r2SCAN-3c TightSCF Opt CPCM(Water)"),
            ):
                (ROOT / f"inputs/orca/acetaminophen_{suffix}.inp").write_text(orca_input(xyz, keywords))
    (out / "water.xyz").write_text("3\nWater installation-test geometry, angstrom\nO 0 0 0\nH 0 -0.757 0.587\nH 0 0.757 0.587\n")
    water = (out / "water.xyz").read_text()
    (ROOT / "inputs/orca/water_sp.inp").write_text(orca_input(water, "PBE def2-SVP TightSCF SP"))
    (ROOT / "inputs/orca/water_opt_freq.inp").write_text(orca_input(water, "r2SCAN-3c TightSCF Opt Freq"))
    mat = ROOT / "structures/materials"
    mat.mkdir(parents=True, exist_ok=True)
    write(mat / "silicon_initial.cif", bulk("Si", "diamond", a=5.43))
    for name, a, thickness in [("MoS2", 3.18, 3.19), ("MoSe2", 3.32, 3.34)]:
        atoms = mx2(name, kind="2H", a=a, thickness=thickness, vacuum=10.0)
        # ASE vacuum=10 adds 10 A on EACH side: inter-slab vacuum is 20 A.
        write(mat / f"{name}_initial.cif", atoms)
        write(mat / f"{name}_2x2_initial.cif", atoms.repeat((2, 2, 1)))
    (ROOT / "reports/initial-structures.json").write_text(json.dumps(records, indent=2) + "\n")
    print("Generated six molecular starting sets, water, Si, MoS2 and MoSe2 (primitive and 2x2).")


if __name__ == "__main__":
    main()
