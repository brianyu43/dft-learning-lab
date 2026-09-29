"""Verified QE 7.6 SCF results as ASE energy/forces, without starting an engine."""
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from ase.io.espresso import read_fortran_namelist

from research_run import sha

BOHR_ANGSTROM = 0.529177210903
HARTREE_EV = 4.3597447222071e-18 / 1.602176634e-19
NUMBER = r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, name, atol=1e-12):
    require(np.allclose(actual, expected, rtol=1e-10, atol=atol), f'{name} mismatch')


def scalar(root, path):
    node = root.find(path)
    require(node is not None and node.text is not None, f'Missing XML {path}')
    return float(node.text)


def verify_method(system, root, expected_method):
    name = str(system['input_dft']).lower()
    if expected_method == 'C09':
        require(name == 'vdw-df-c09' and str(system.get('vdw_corr', 'none')).lower() == 'none',
                'Expected C09 without an extra empirical dispersion correction')
    elif expected_method == 'PBE-D3BJ':
        require(name == 'pbe' and str(system.get('vdw_corr', '')).lower() == 'grimme-d3',
                'Expected PBE plus Grimme-D3')
        require(system.get('dftd3_version') == 4 and system.get('dftd3_threebody') is True,
                'PBE-D3BJ requires explicit version 4 and three-body terms')
    else:
        raise ValueError('Unreviewed DFT method')
    for section in ('input/dft', 'output/dft'):
        require(root.findtext(section + '/functional', '').lower() == name, 'XML functional mismatch')
        corr = root.findtext(section + '/vdW/vdw_corr', 'none').lower()
        if expected_method == 'C09':
            require(corr == 'none' and root.findtext(section + '/vdW/non_local_term') == 'VDW1',
                    'C09 nonlocal dispersion missing or changed')
        else:
            require(corr == 'grimme-d3' and scalar(root, section + '/vdW/dftd3_version') == 4
                    and root.findtext(section + '/vdW/dftd3_threebody', '').lower() == 'true',
                    'XML D3(BJ) version or three-body setting mismatch')
            require(root.findtext(section + '/vdW/non_local_term', 'none').lower() in ('none', ''),
                    'Unexpected nonlocal functional combined with D3')


def load_scf(report_path, expected_method='C09', expected_atoms=None):
    report_path = Path(report_path)
    meta = json.loads(report_path.read_text())
    require(meta['engine'] == 'qe' and meta['status'] == 'passed', 'SCF report is not passed')
    require(meta['result']['normal_termination'] and meta['result']['scf_converged'], 'SCF incomplete')
    require(not meta['requested']['optimization'], 'Only fixed-geometry SCF is supported by this bridge')
    run = Path(meta['run_dir'])
    inp = run / Path(meta['input']).name
    require(sha(inp) == meta['source_sha256'][meta['input']], 'Input hash mismatch')
    require(sha(run / 'output.log') == meta['output_sha256'], 'Output hash mismatch')
    require(sha(run / 'runner_snapshot.py') == meta['runner_sha256'], 'Runner snapshot hash mismatch')
    for source, destination in meta['assets'].items():
        require(sha(run / destination) == meta['source_sha256'][source], 'Input asset hash mismatch')
    with inp.open() as stream:
        namelist, _ = read_fortran_namelist(stream)
    require(namelist['control']['calculation'].lower() == 'scf', 'Input is not SCF')
    system, electrons = namelist['system'], namelist['electrons']
    require(system.get('nspin', 1) == 1 and not system.get('noncolin', False)
            and not system.get('lspinorb', False), 'Unexpected spin settings')
    require(system.get('occupations', 'fixed').lower() == 'fixed', 'Unexpected electronic occupations')
    initial = read(inp, format='espresso-in')
    require(len(initial) == meta['n_atoms'], 'Input atom count mismatch')
    paths = list((run / 'scratch').glob('*.save/data-file-schema.xml'))
    require(len(paths) == 1, 'Expected one XML result')
    xml = paths[0]
    root = ET.parse(xml).getroot()
    require(root.findtext('exit_status') == '0' and root.find('closed') is not None, 'XML is not final')
    require(root.find('general_info/creator').get('VERSION') == '7.6', 'Unreviewed QE XML version')
    require(root.findtext('input/control_variables/calculation') == 'scf', 'XML calculation is not SCF')
    require(root.findtext('output/convergence_info/scf_conv/convergence_achieved') == 'true', 'XML SCF not converged')
    verify_method(system, root, expected_method)
    for setting in ('ecutwfc', 'ecutrho'):
        close(2 * scalar(root, 'input/basis/' + setting), float(system[setting]), setting)
        close(2 * scalar(root, 'output/basis_set/' + setting), float(system[setting]), setting)
        if setting + '_Ry' in meta:
            close(float(meta[setting + '_Ry']), float(system[setting]), 'Metadata ' + setting)
    conv_thr = float(electrons['conv_thr'])
    close(2 * scalar(root, 'input/electron_control/conv_thr'), conv_thr, 'SCF threshold', atol=1e-16)
    if 'conv_thr_Ry' in meta:
        close(float(meta['conv_thr_Ry']), conv_thr, 'Metadata SCF threshold', atol=1e-16)
    require(2 * scalar(root, 'output/convergence_info/scf_conv/scf_error') <= conv_thr * 1.00001,
            'Final SCF error exceeds the requested threshold')
    grids = []
    for section in ('input/k_points_IBZ', 'output/band_structure/starting_k_points'):
        grid = root.find(section + '/monkhorst_pack')
        require(grid is not None, 'Automatic k-grid is required')
        grids.append([int(grid.get(k)) for k in ('nk1', 'nk2', 'nk3', 'k1', 'k2', 'k3')])
    require(grids[0] == grids[1] and grids[0][3:] == [0, 0, 0], 'XML k-grid changed or shifted')
    require(grids[0][:3] == meta['k_grid'], 'Metadata k-grid mismatch')
    match = re.search(r'K_POINTS\s+automatic\s*\n\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)', inp.read_text(), re.I)
    require(match is not None and list(map(int, match.groups())) == grids[0], 'Input k-grid mismatch')
    for flag in ('lsda', 'noncolin', 'spinorbit'):
        require(root.findtext('output/band_structure/' + flag) == 'false', 'XML spin setting changed')
    require(root.findtext('output/band_structure/occupations_kind') == 'fixed', 'XML occupations changed')

    structure = root.find('output/atomic_structure')
    species = [item.get('name') for item in structure.find('atomic_positions')]
    positions = np.array([np.fromstring(item.text, sep=' ') for item in structure.find('atomic_positions')]) * BOHR_ANGSTROM
    cell = np.array([np.fromstring(item.text, sep=' ') for item in structure.find('cell')]) * BOHR_ANGSTROM
    forces = np.fromstring(root.findtext('output/forces'), sep=' ').reshape(-1, 3) * HARTREE_EV / BOHR_ANGSTROM
    require(species == initial.get_chemical_symbols() and len(forces) == len(initial), 'Atom ordering or force count changed')
    require(np.isfinite(positions).all() and np.isfinite(cell).all() and np.isfinite(forces).all(), 'Non-finite XML data')
    close(cell, initial.cell.array, 'SCF cell', atol=1e-7)
    close(positions, initial.positions, 'SCF coordinates', atol=1e-7)
    if expected_atoms is not None:
        require(species == expected_atoms.get_chemical_symbols(), 'Requested species differ')
        close(cell, expected_atoms.cell.array, 'Requested cell', atol=1e-7)
        close(positions, expected_atoms.positions, 'Requested coordinates', atol=1e-7)
    valences = {}
    for item in root.find('input/atomic_species'):
        symbol, pseudo_file = item.get('name'), item.findtext('pseudo_file')
        pseudo = run / 'pseudo' / pseudo_file
        require(pseudo.resolve().is_relative_to((run / 'pseudo').resolve()), 'Pseudo path outside run')
        require('pseudo/' + pseudo_file in meta['assets'].values(), 'Unverified pseudopotential')
        valences[symbol] = float(re.search(r'z_valence\s*=\s*"\s*([\d.]+)', pseudo.read_text()).group(1))
    charge = float(system.get('tot_charge', 0))
    close(charge, float(meta['charge']), 'Charge')
    close(scalar(root, 'output/band_structure/nelec'), sum(valences[s] for s in species) - charge, 'Electron count')
    energy_ha = scalar(root, 'output/total_energy/etot')
    require(meta['result']['energy_unit'] == 'Ry', 'Wrong report energy unit')
    close(2 * energy_ha, meta['result']['energy'], 'Report energy', atol=1e-7)
    log = (run / 'output.log').read_text()
    text_forces = re.findall(r'atom\s+\d+\s+type\s+\d+\s+force\s*=\s*(' + NUMBER + r')\s+(' + NUMBER + r')\s+(' + NUMBER + ')', log)
    require(len(text_forces) == len(initial) and 'Ry/au' in log, 'Expected one Ry/au force block')
    converted = np.array([[float(v.replace('D', 'E').replace('d', 'e')) for v in row] for row in text_forces]) * HARTREE_EV / (2 * BOHR_ANGSTROM)
    force_difference = float(np.abs(forces - converted).max())
    require(force_difference <= 0.51e-8 * HARTREE_EV / (2 * BOHR_ANGSTROM), 'XML/text force unit or value mismatch')
    threads = re.search(r'Serial multi-threaded version, running on\s+(\d+) processor cores', log)
    require(threads is not None and int(threads.group(1)) == meta['omp_threads'], 'Actual OpenMP header differs from request')
    atoms = expected_atoms.copy() if expected_atoms is not None else initial.copy()
    atoms.set_positions(positions, apply_constraint=False)
    atoms.set_cell(cell, apply_constraint=False)
    atoms.calc = SinglePointCalculator(atoms, energy=energy_ha * HARTREE_EV, forces=forces.copy())
    audit = dict(job_id=meta['job_id'], method=expected_method, fixed_geometry_verified=True,
                 energy_eV=energy_ha * HARTREE_EV, max_force_eV_angstrom=float(np.linalg.norm(forces, axis=1).max()),
                 energy_conversion_eV_per_Hartree=HARTREE_EV, length_conversion_angstrom_per_Bohr=BOHR_ANGSTROM,
                 xml_text_force_max_difference_eV_angstrom=force_difference,
                 ecutwfc_Ry=system['ecutwfc'], ecutrho_Ry=system['ecutrho'], k_grid=grids[0][:3],
                 k_offset=grids[0][3:], conv_thr_Ry=conv_thr,
                 valence_electrons=valences, openmp_threads_in_log=int(threads.group(1)),
                 xml_nthreads=int(root.findtext('parallel_info/nthreads')),
                 xml_thread_field_is_not_used_as_runtime_proof=True,
                 report_sha256=sha(report_path), xml_sha256=sha(xml), output_sha256=sha(run / 'output.log'),
                 input_sha256=sha(inp), bridge_sha256=sha(Path(__file__)),
                 scope='Verified existing SCF and ASE transfer; no optimization, finite-displacement or D3 engine validation implied')
    return dict(atoms=atoms, energy_eV=energy_ha * HARTREE_EV, forces=forces, audit=audit)
