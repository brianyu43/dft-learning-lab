"""Audit PubChem identities, enumerate 2D tautomer candidates, and build seed geometries.

Tautomer enumeration is an inventory, not evidence for state populations or an
approved DFT shortlist. Standard InChIKey does not uniquely fix tautomer hydrogens.
"""
import hashlib
import json
from pathlib import Path
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'structures/week1'
OUT.mkdir(parents=True, exist_ok=True)
records = []
for name in ('diflunisal', 'adenine', 'caffeine', 'theobromine', 'xanthine'):
    source = ROOT / f'data/literature/week1/pubchem/{name}.json'
    data = json.loads(source.read_text())['PropertyTable']['Properties'][0]
    molecule = Chem.MolFromSmiles(data['SMILES'])
    assert molecule is not None
    assert rdMolDescriptors.CalcMolFormula(molecule) == data['MolecularFormula']
    assert Chem.MolToInchiKey(molecule) == data['InChIKey']
    assert Chem.GetFormalCharge(molecule) == 0
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in molecule.GetAtoms())
    tautomers = sorted(Chem.MolToSmiles(m) for m in rdMolStandardize.TautomerEnumerator().Enumerate(molecule))
    parent = Chem.MolToSmiles(molecule)
    assert parent in tautomers
    state_records = []
    state_writer = Chem.SDWriter(str(OUT / f'{name}_tautomer_inventory_2d.sdf'))
    for index, smiles in enumerate(tautomers):
        m = Chem.MolFromSmiles(smiles)
        assert rdMolDescriptors.CalcMolFormula(m) == data['MolecularFormula']
        assert Chem.GetFormalCharge(m) == 0
        state_id = f'{name}-t{index:02d}'
        AllChem.Compute2DCoords(m)
        m.SetProp('_Name', state_id)
        m.SetProp('status', 'unreviewed_2D_tautomer_candidate')
        m.SetProp('explicit_state_smiles', smiles)
        state_writer.write(m)
        state_records.append({'state_id': state_id, 'smiles': smiles,
                              'is_pubchem_representation': smiles == parent,
                              'status': 'unreviewed_candidate_not_DFT_selected'})
    state_writer.close()
    geometries = []
    writer = Chem.SDWriter(str(OUT / f'{name}_pubchem_state_initial.sdf'))
    for seed in (17, 43, 89):
        m = Chem.AddHs(molecule)
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.numThreads = 1
        assert AllChem.EmbedMolecule(m, params) == 0
        assert AllChem.MMFFHasAllMoleculeParams(m)
        status = AllChem.MMFFOptimizeMolecule(m, maxIters=1000)
        assert status == 0
        props = AllChem.MMFFGetMoleculeProperties(m)
        energy = AllChem.MMFFGetMoleculeForceField(m, props).CalcEnergy()
        assert m.GetNumAtoms() == sum(a.GetAtomicNum() > 0 for a in m.GetAtoms())
        label = f'{name}-pubchem-seed{seed}'
        m.SetProp('_Name', label)
        m.SetProp('seed', str(seed))
        m.SetProp('MMFF_energy_kcal_mol', str(energy))
        m.SetProp('status', 'MMFF_initial_geometry_only')
        writer.write(m)
        xyz = OUT / f'{label}.xyz'
        xyz.write_text(Chem.MolToXYZBlock(m))
        geometries.append({'seed': seed, 'path': str(xyz.relative_to(ROOT)),
                           'sha256': hashlib.sha256(xyz.read_bytes()).hexdigest(),
                           'MMFF_energy_kcal_mol': energy, 'MMFF_status': status})
    writer.close()
    records.append({'name': name, 'pubchem_cid': data['CID'],
                    'pubchem_url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{data['CID']}",
                    'source_path': str(source.relative_to(ROOT)),
                    'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
                    'formula': data['MolecularFormula'], 'InChIKey': data['InChIKey'],
                    'pubchem_state_smiles': parent,
                    'n_atoms_including_hydrogen': Chem.AddHs(molecule).GetNumAtoms(),
                    'charge': 0, 'multiplicity': 1, 'identity_checks': 'passed',
                    'tautomer_count': len(tautomers), 'tautomer_inventory': state_records,
                    'initial_geometries': geometries,
                    'DFT_state_selection': 'pending_chemical_review',
                    'claim_scope': 'identity and starting geometry only; no DFT energy, conformer population, crystal phase, or tautomer population measured'})
result = {'rdkit_version': rdBase.rdkitVersion, 'protocol_id': 'dft-pilot-20260909-v1',
          'seeds': [17,43,89], 'molecules': records,
          'dimer_generation': 'week2_not_started',
          'state_warning': 'For adenine, include neutral tautomers and plausible proton-transfer complexes in the pre-DFT review. Candidate enumeration is not chemical validation.'}
(ROOT / 'reports/week1-molecule-identities.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({'molecules': len(records), 'tautomer_candidates': sum(r['tautomer_count'] for r in records),
                  'MMFF_initial_geometries': sum(len(r['initial_geometries']) for r in records)}))
