"""Cross-check completed scientific outputs; never launch an engine."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import xml.etree.ElementTree as ET
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.io.espresso import units as espresso_units

# QE 7.6 Modules/constants.f90. ASE espresso's reader retains CODATA 2006.
HARTREE_QE_EV = 4.3597447222071e-18 / 1.602176634e-19

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'structures/optimized/week1'
DEST.mkdir(parents=True, exist_ok=True)
checks = []

def check_recorded_hashes(manifest):
    run = Path(manifest['run_dir'])
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest(run / 'output.log') == manifest['output_sha256']
    assert digest(run / Path(manifest['input']).name) == manifest['source_sha256'][manifest['input']]
    for source, destination in manifest.get('assets', {}).items():
        assert digest(run / destination) == manifest['source_sha256'][source]
    assert digest(run / 'runner_snapshot.py') == manifest['runner_sha256']

for species, count in (('water', 3), ('acetaminophen', 20)):
    manifest = json.loads((ROOT / f'reports/week1-{species}-opt-freq.json').read_text())
    assert manifest['status'] == 'passed'
    check_recorded_hashes(manifest)
    run = Path(manifest['run_dir'])
    stem = Path(manifest['input']).stem
    hess = run / f'{stem}.hess'
    lines = hess.read_text().splitlines()
    start = lines.index('$vibrational_frequencies')
    n = int(lines[start+1])
    frequencies = np.array([float(row.split()[1]) for row in lines[start+2:start+2+n]])
    assert n == 3*count and len(frequencies) == n
    np.testing.assert_allclose(frequencies[:6], 0, atol=1e-8)
    assert np.all(frequencies[6:] > 0)
    np.testing.assert_allclose(frequencies[6:], manifest['result']['frequencies_cm1'], atol=0.0051, rtol=0)
    xyz = run / f'{stem}.xyz'
    atoms = read(xyz)
    assert len(atoms) == count
    comment = xyz.read_text().splitlines()[1]
    xyz_energy = float(re.search(r' E ([-\d.]+)', comment).group(1))
    assert abs(xyz_energy-manifest['result']['energy']) < 1e-10
    assert np.isfinite(atoms.positions).all()
    target = DEST / f'{species}.xyz'
    shutil.copy2(xyz, target)
    result = {'species': species, 'status': 'passed', 'atom_count': len(atoms),
              'formula': atoms.get_chemical_formula(),
              'hessian_sha256': hashlib.sha256(hess.read_bytes()).hexdigest(),
              'optimized_xyz_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
              'hessian_vs_log_frequency_max_difference_cm1': float(np.max(np.abs(frequencies[6:]-manifest['result']['frequencies_cm1']))),
              'hessian_lowest_vibration_cm1': float(min(frequencies[6:])),
              'xyz_comment_energy_matches_log': True,
              'optimized_xyz': str(target.relative_to(ROOT)),
              'energy_change_from_initial_Hartree': manifest['result']['energy']-manifest['result']['energy_history'][0]}
    if species == 'water':
        result['OH_bonds_angstrom'] = [float(atoms.get_distance(0, i)) for i in (1,2)]
        result['HOH_angle_degrees'] = float(atoms.get_angle(1,0,2))
    checks.append(result)

manifest = json.loads((ROOT / 'reports/week1-mos2-relax.json').read_text())
assert manifest['status'] == 'passed'
check_recorded_hashes(manifest)
run = Path(manifest['run_dir'])
images = read(run / 'output.log', index=':', format='espresso-out')
atoms = images[-1]
initial = read(ROOT / 'structures/materials/MoS2_initial.cif')
assert len(atoms) == 3 and atoms.get_chemical_formula() == 'MoS2'
# Console crystal axes have only six decimals; use full-precision XML for cell invariance.
np.testing.assert_allclose(atoms.cell.array, initial.cell.array, atol=2e-6, rtol=0)
xml_path = run / 'scratch/week1_mos2.save/data-file-schema.xml'
xml_output = ET.parse(xml_path).getroot().find('output')
xml_structure = xml_output.find('atomic_structure')
bohr_qe76 = 0.529177210903  # Installed QE Modules/constants.f90.
xml_cell = np.array([np.fromstring(row.text, sep=' ') for row in xml_structure.find('cell')])*bohr_qe76
np.testing.assert_allclose(xml_cell, initial.cell.array, atol=1e-10, rtol=0)
energy = float(atoms.get_potential_energy())
assert abs(energy / espresso_units['Ry'] - manifest['result']['energy']) < 1e-8
forces = atoms.get_forces()
assert np.isfinite(forces).all()
max_force = float(np.max(np.linalg.norm(forces, axis=1)))
assert max_force < 0.01
sites = list(xml_structure.find('atomic_positions'))
xml_positions = np.array([np.fromstring(row.text, sep=' ') for row in sites])*bohr_qe76
np.testing.assert_allclose(atoms.positions, xml_positions, atol=2e-6, rtol=0)
xml_energy_au = float(xml_output.find('total_energy/etot').text)
assert abs(2*xml_energy_au - manifest['result']['energy']) < 1e-8
xml_energy = xml_energy_au*HARTREE_QE_EV
xml_forces = np.fromstring(xml_output.find('forces').text, sep=' ').reshape(3,3)*HARTREE_QE_EV/bohr_qe76
np.testing.assert_allclose(xml_forces, forces, atol=3e-7, rtol=0)
atoms = Atoms([row.attrib['name'] for row in sites], positions=xml_positions, cell=xml_cell, pbc=True)
path = DEST / 'MoS2_fixed_cell_relaxed.cif'
write(path, atoms)
write(DEST / 'MoS2_fixed_cell_relaxed.extxyz', atoms)
checks.append({'species': 'MoS2', 'status': 'passed', 'ASE_ionic_images': len(images),
               'energy_eV': xml_energy, 'ASE_energy_matches_Ry_conversion': True,
               'unit_policy': 'ASE reader checked in its own CODATA 2006 units; final eV/Angstrom from XML and QE 7.6 constants',
               'cell_fixed_verified': True, 'atom_count': len(atoms),
               'cell_verification_source': 'QE full-precision XML, not rounded console axes',
               'QE_XML_sha256': hashlib.sha256(xml_path.read_bytes()).hexdigest(),
               'final_max_atomic_force_eV_angstrom': float(np.max(np.linalg.norm(xml_forces, axis=1))),
               'final_S_to_S_thickness_angstrom': float(np.ptp(atoms.positions[:,2])),
               'final_vacuum_gap_angstrom': float(atoms.cell[2,2]-np.ptp(atoms.positions[:,2])),
               'optimized_cif': str(path.relative_to(ROOT)),
               'optimized_cif_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
               'phonon_stability_verified': False, 'precision_convergence_verified': False})
check_meta = json.loads((ROOT / 'reports/week1-mos2-force-check.json').read_text())
assert check_meta['status'] == 'passed'
check_recorded_hashes(check_meta)
check_run = Path(check_meta['run_dir'])
check_xml_path = check_run / 'scratch/week1_mos2_check.save/data-file-schema.xml'
check_xml = ET.parse(check_xml_path).getroot().find('output')
check_structure = check_xml.find('atomic_structure')
check_cell = np.array([np.fromstring(row.text, sep=' ') for row in check_structure.find('cell')])*bohr_qe76
check_pos = np.array([np.fromstring(row.text, sep=' ') for row in check_structure.find('atomic_positions')])*bohr_qe76
np.testing.assert_allclose(check_cell, xml_cell, atol=1e-10, rtol=0)
np.testing.assert_allclose(check_pos, xml_positions, atol=1e-9, rtol=0)
check_energy = float(check_xml.find('total_energy/etot').text)*HARTREE_QE_EV
check_forces = np.fromstring(check_xml.find('forces').text, sep=' ').reshape(3,3)*HARTREE_QE_EV/bohr_qe76
force_change = float(np.max(np.linalg.norm(check_forces-xml_forces, axis=1)))
check_max_force = float(np.max(np.linalg.norm(check_forces, axis=1)))
assert check_max_force < 0.01
assert abs(check_energy-xml_energy) < 1e-4
checks.append({'species': 'MoS2-strict-SCF', 'status': 'passed',
               'conv_thr_change_Ry': [1e-8, 1e-10], 'geometry_and_cell_unchanged': True,
               'energy_change_eV': check_energy-xml_energy,
               'max_atomic_force_eV_angstrom': check_max_force,
               'max_atomic_force_change_eV_angstrom': force_change,
               'QE_XML_sha256': hashlib.sha256(check_xml_path.read_bytes()).hexdigest(),
               'scope': 'One SCF threshold check only; no complete numerical convergence study'})
report = {'status': 'passed', 'checks': checks,
          'scope': 'Independent file/format checks on real engine outputs, not independent physical methods or experimental validation',
          'cclib_crosscheck': 'unavailable: cclib 1.8.1 raised IndexError while reading ORCA 6.1.1 SCF block; Hessian file and XYZ/log cross-check used instead'}
(ROOT / 'reports/week1-verification.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
