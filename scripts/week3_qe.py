"""Week-three QE preparation and XML verification, including Se and cell relaxation."""
import hashlib
import io
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
from ase import Atoms
from ase.io import read, write

from research_run import sha
from prepare_week3_drug import immutable

ROOT = Path(__file__).resolve().parents[1]
PSEUDOS = 'data/pseudopotentials/pseudo-dojo-nc-sr-04-pbe-standard'
BOHR = 0.529177210903
HARTREE_EV = 4.3597447222071e-18 / 1.602176634e-19


def json_file(path, value):
    immutable(path, json.dumps(value, indent=2) + '\n')


def job(job_id, atoms, calculation='scf', cut=120, k=12, **extra):
    elements = sorted(set(atoms.get_chemical_symbols()))
    assert all(e in ['Mo', 'S', 'Se'] for e in elements)
    system = {'ecutwfc': cut, 'ecutrho': 4*cut, 'input_dft': 'vdw-df-c09', 'occupations': 'fixed', 'nspin': 1}
    data = {'control': {'calculation': calculation, 'prefix': job_id.replace('-', '_'),
                        'pseudo_dir': './pseudo', 'outdir': './scratch', 'tprnfor': True, 'tstress': True,
                        'nstep': 100, 'etot_conv_thr': 1e-6, 'forc_conv_thr': 1e-4},
            'system': system,
            'electrons': {'conv_thr': 1e-10, 'electron_maxstep': 200, 'mixing_beta': 0.4}}
    if calculation in ['relax', 'vc-relax']:
        data['ions'] = {'ion_dynamics': 'bfgs'}
    if calculation == 'vc-relax':
        data['cell'] = {'cell_dynamics': 'bfgs', 'cell_dofree': '2Dxy', 'press': 0.0, 'press_conv_thr': 0.1}
    inp = ROOT / f'inputs/qe/week3/{job_id}.in'
    text = io.StringIO()
    write(text, atoms, format='espresso-in', pseudopotentials={e: e+'.upf' for e in elements},
          kpts=(k, k, 1), koffset=(0, 0, 0), input_data=data)
    immutable(inp, text.getvalue())
    structure = ROOT / f'structures/week3/materials/{job_id}.extxyz'
    text = io.StringIO(); write(text, atoms, format='extxyz'); immutable(structure, text.getvalue())
    cfg = dict(protocol_id='dft-pilot-20260909-v1-week3', week=3, task_id='SR03', job_id=job_id,
               engine='qe', executable='.tools/qe-7.6/bin/pw.x', input=str(inp.relative_to(ROOT)),
               assets={f'{PSEUDOS}/{e}.upf': f'pseudo/{e}.upf' for e in elements},
               method='vdW-DF-C09; PseudoDojo nc-sr-04 standard PBE scalar NC',
               requested={'optimization': calculation != 'scf', 'frequencies': False}, calculation=calculation,
               n_atoms=len(atoms), charge=0, SOC=False, ecutwfc_Ry=cut, ecutrho_Ry=4*cut, k_grid=[k,k,1],
               initial_vacuum_gap_angstrom=float(atoms.cell[2,2]-np.ptp(atoms.positions[:,2])),
               initial_cell_angstrom=atoms.cell.array.tolist(), geometry_source=str(structure.relative_to(ROOT)),
               conv_thr_Ry=1e-10, timeout_seconds=14400, rss_limit_gib=40, omp_threads=4,
               max_force_eV_angstrom=0.005, accuracy_scope='Composition-specific checks required; no phonon or experimental-cell validation',
               **extra)
    path = ROOT / f'experiments/week3/materials/{job_id}.json'; json_file(path, cfg)
    return str(path.relative_to(ROOT))


def load(job_id):
    meta = json.loads((ROOT / f'reports/{job_id}.json').read_text())
    assert meta['status'] == 'passed', f'Engine requires review: {job_id}'
    run = Path(meta['run_dir'])
    assert sha(run / 'output.log') == meta['output_sha256']
    assert sha(run / Path(meta['input']).name) == meta['source_sha256'][meta['input']]
    for source, dest in meta['assets'].items():
        assert sha(run / dest) == meta['source_sha256'][source]
    files = list((run/'scratch').glob('*.save/data-file-schema.xml')); assert len(files) == 1
    root = ET.parse(files[0]).getroot(); output = root.find('output')
    structure = output.find('atomic_structure')
    symbols = [a.attrib['name'] for a in structure.find('atomic_positions')]
    positions = np.array([np.fromstring(a.text, sep=' ') for a in structure.find('atomic_positions')])*BOHR
    cell = np.array([np.fromstring(a.text, sep=' ') for a in structure.find('cell')])*BOHR
    atoms = Atoms(symbols, positions=positions, cell=cell, pbc=True)
    energy = float(output.find('total_energy/etot').text)
    forces = np.fromstring(output.find('forces').text, sep=' ').reshape(-1,3)*HARTREE_EV/BOHR
    stress_node = output.find('stress')
    stress = (np.fromstring(stress_node.text, sep=' ').reshape(3,3)*HARTREE_EV/BOHR**3
              if stress_node is not None else None)
    nelec = float(output.find('band_structure/nelec').text)
    valence = {e: float(re.search(r'z_valence\s*=\s*"\s*([\d.]+)', (run / f'pseudo/{e}.upf').read_text()).group(1)) for e in set(symbols)}
    assert abs(nelec - sum(valence[e] for e in symbols) + meta['charge']) < 1e-8
    assert abs(energy*2 - meta['result']['energy']) < 1e-7
    initial = read(run / Path(meta['input']).name, format='espresso-in')
    assert symbols == initial.get_chemical_symbols() and np.isfinite(positions).all() and np.isfinite(forces).all()
    assert np.allclose(cell[2], initial.cell[2], rtol=0, atol=1e-7), 'Vacuum cell vector changed'
    if meta['calculation'] == 'scf':
        assert np.allclose(cell, initial.cell, rtol=0, atol=1e-7)
        assert np.allclose(positions, initial.positions, rtol=0, atol=1e-7)
    elif meta['calculation'] == 'relax':
        assert np.allclose(cell, initial.cell, rtol=0, atol=1e-7)
    max_force = float(np.linalg.norm(forces, axis=1).max())
    if meta['requested']['optimization']:
        assert max_force <= meta['max_force_eV_angstrom'], f'Final force {max_force} exceeds review limit'
    if meta['calculation'] == 'vc-relax':
        assert stress is not None, 'Final cell stress is missing'
        assert float(np.abs(stress[:2,:2]).max())*1602.176634 <= 0.5, 'Final-basis in-plane stress requires review'
        assert abs(atoms.cell.lengths()[0]-atoms.cell.lengths()[1]) < 1e-5
        assert abs(atoms.cell.angles()[2]-120) < 1e-4, 'Hexagonal symmetry changed'
    return dict(job_id=job_id, meta=meta, atoms=atoms, energy_eV=energy*HARTREE_EV,
                forces=forces, stress_eV_angstrom3=stress, max_force_eV_angstrom=max_force,
                area_angstrom2=float(np.linalg.norm(np.cross(cell[0],cell[1]))),
                vacuum_gap_angstrom=float(cell[2,2]-np.ptp(positions[:,2])), nelec=nelec,
                xml_sha256=sha(files[0]))


def monolayer(element='S', a=None):
    base = read(ROOT / 'structures/optimized/week1/MoS2_fixed_cell_relaxed.cif')
    base.positions[:,2] -= base.positions[0,2]
    if element == 'Se':
        base.set_chemical_symbols(['Mo','Se','Se'])
        base.positions[:,2] *= 1.08
    if a is not None:
        cell = base.cell.array.copy(); cell[:2] *= a/base.cell.lengths()[0]
        base.set_cell(cell, scale_atoms=True)
    base.center(vacuum=15, axis=2)
    return base


if __name__ == '__main__':
    configs = [job(f'w3-{name}-lattice-vc', monolayer(element, a), calculation='vc-relax',
                   role='in_plane_equilibrium_lattice', initial_lattice_is_only_starting_guess=True)
               for name, element, a in [('mos2','S',3.18), ('mose2','Se',3.30)]]
    json_file(ROOT / 'experiments/week3/materials/lattice-batch.json',
              dict(batch_id='week3-lattice-equilibria', budget_seconds=14400, jobs=configs,
                   execution_enabled=False, review='Run after first drug full-cost review; verify Se sensitivity and final cell before heterobilayers'))
    print(json.dumps(configs))
