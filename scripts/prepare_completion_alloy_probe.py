"""Prepare a three-SCF cost probe from a verified reference, without running QE.

The 27-comparison catalog is an in-memory blueprint. Only CI/r0/D0 and its
two extracted layers become executable inputs. A full batch still requires
measured 24-atom cost and a separate budget review.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
from ase.io import read, write
from ase.io.espresso import read_fortran_namelist

from completion_alloy_geometry import build_family, decoration_equivalence, validate_geometry
from completion_qe_bridge import load_scf, require, scalar, verify_method
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from week3_qe import load as load_reference
from week3_materials_pipeline import bilayer

PSEUDOS = 'data/pseudopotentials/pseudo-dojo-nc-sr-04-pbe-standard'
REFERENCE = 'reports/week3-materials.json'
PLAN = 'experiments/completion/alloy-cost-probe-r1.json'


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def settings_from_reference(selection):
    cut, primitive_k, vacuum = (selection[key] for key in ('ecutwfc_Ry', 'k', 'vacuum_angstrom'))
    require((cut, primitive_k, vacuum) in ((80, 6, 20), (100, 9, 25), (120, 12, 30)),
            'Unreviewed reference settings')
    k = math.ceil(primitive_k / 2)
    return dict(method='C09', ecutwfc_Ry=cut, ecutrho_Ry=4 * cut,
                k_grid=[k, k, 1], k_offset=[0, 0, 0], conv_thr_Ry=1e-10,
                primitive_reference_k=primitive_k, primitive_equivalent_k=2 * k,
                k_rule='ceil(primitive k / 2); no coarser spacing, not identical quadrature',
                vacuum_angstrom=vacuum, charge=0, nspin=1, SOC=False, occupations='fixed')


def calculation_identity(atoms, settings, pseudo_hashes, executable_sha256):
    """Exact computational identity: no centering, wrapping or approximate merge."""
    require(np.isfinite(atoms.positions).all() and np.isfinite(atoms.cell).all(), 'Nonfinite structure')
    return dict(numbers=atoms.numbers.tolist(), positions_angstrom=atoms.positions.tolist(),
                cell_angstrom=atoms.cell.array.tolist(), pbc=atoms.pbc.tolist(), settings=settings,
                canonical_input_sha256=hashlib.sha256(render_scf('alloy_identity', atoms, settings).encode()).hexdigest(),
                pseudopotential_sha256={e: pseudo_hashes[e] for e in sorted(set(atoms.get_chemical_symbols()))},
                executable_sha256=executable_sha256)


def build_catalog(bottom, top, d0, settings, pseudo_hashes, executable_sha256):
    """Build the fixed 27 comparisons and deduplicate only exact identities."""
    require(np.isfinite(d0) and 6.0 <= d0 <= 8.1, 'D0 outside the reviewed bracket')
    require(bottom.pbc.all() and top.pbc.all(), 'Periodic layers required')
    calculations, rows = {}, []
    for rid in range(3):
        for di, distance in enumerate((d0 - .25, d0, d0 + .25)):
            family = build_family(bottom, top, rid, distance, settings['vacuum_angstrom'])
            symmetry = decoration_equivalence(family['CI'], family['CO'])
            require(not symmetry['equivalent'], 'CI and CO are symmetry-equivalent')
            for label, pair in family.items():
                validate_geometry(pair)
                parts = {'pair': pair, 'bottom': pair[pair.arrays['layer_id'] == 0],
                         'top': pair[pair.arrays['layer_id'] == 1]}
                links = {}
                for part, atoms in parts.items():
                    atoms.calc = None
                    atoms.set_constraint()
                    require(np.array_equal(atoms.cell.array, pair.cell.array), 'Reference cell changed')
                    identity = calculation_identity(atoms, settings, pseudo_hashes, executable_sha256)
                    key = digest(identity)
                    if key in calculations:
                        require(calculations[key]['identity'] == identity, 'Identity hash collision')
                    else:
                        calculations[key] = dict(identity=identity, atoms=atoms.copy(), uses=[])
                    calculations[key]['uses'].append(dict(configuration=label, registry_id=rid,
                                                          distance_index=di, part=part))
                    links[part] = key
                rows.append(dict(configuration=label, registry_id=rid, distance_index=di,
                                 Mo_plane_distance_angstrom=distance,
                                 area_angstrom2=float(np.linalg.norm(np.cross(pair.cell[0], pair.cell[1]))),
                                 calculations=links, ci_co_symmetry=symmetry))
    require(len(rows) == 27 and len({r['calculations']['pair'] for r in rows}) == 27,
            'Incomplete pair matrix')
    return rows, calculations


def _raw_paths(meta, report):
    run = Path(meta['run_dir'])
    paths = [report, run / 'output.log', run / Path(meta['input']).name, run / 'runner_snapshot.py']
    paths += list((run / 'scratch').glob('*.save/data-file-schema.xml'))
    paths += [run / dest for dest in meta['assets'].values()]
    for element in set(read(run / Path(meta['input']).name, format='espresso-in').get_chemical_symbols()):
        require(sha(run / 'pseudo' / (element + '.upf')) == sha(ROOT / PSEUDOS / (element + '.upf')),
                'Current pseudopotential differs from the verified reference')
    return paths


def verify_optimized(job_id, summary):
    result = load_reference(job_id)
    meta = result['meta']
    require(meta['requested']['optimization'], 'Expected an optimized reference')
    run = Path(meta['run_dir'])
    require(sha(run / 'runner_snapshot.py') == meta['runner_sha256'], 'Runner changed')
    xml = next((run / 'scratch').glob('*.save/data-file-schema.xml'))
    root = ET.parse(xml).getroot()
    require(root.findtext('exit_status') == '0' and root.find('closed') is not None, 'Unfinished XML')
    for kind in ('scf_conv', 'opt_conv'):
        require(root.findtext(f'output/convergence_info/{kind}/convergence_achieved') == 'true',
                'Reference did not converge')
    with (run / Path(meta['input']).name).open() as handle:
        data, _ = read_fortran_namelist(handle)
    verify_method(data['system'], root, 'C09')
    require(data['system'].get('nspin', 1) == 1 and data['system'].get('occupations', 'fixed') == 'fixed'
            and data['system'].get('tot_charge', 0) == 0 and not data['system'].get('lspinorb', False)
            and not data['system'].get('noncolin', False), 'Unexpected reference electronic state')
    for setting in ('ecutwfc', 'ecutrho'):
        require(abs(2 * scalar(root, 'output/basis_set/' + setting) - meta[setting + '_Ry']) <= 1e-8
                and data['system'][setting] == meta[setting + '_Ry'], 'Reference cutoff mismatch')
    grid = root.find('output/band_structure/starting_k_points/monkhorst_pack')
    require(grid is not None and [int(grid.get(k)) for k in ('nk1', 'nk2', 'nk3')] == meta['k_grid']
            and [int(grid.get(k)) for k in ('k1', 'k2', 'k3')] == [0, 0, 0], 'Reference k-grid mismatch')
    require(2 * scalar(root, 'output/convergence_info/scf_conv/scf_error') <= 1e-10 * 1.00001,
            'Final reference SCF precision failed')
    require(result['xml_sha256'] == summary['xml_sha256'], 'Reference XML provenance mismatch')
    require(abs(result['energy_eV'] - summary['energy_eV']) <= 1e-8, 'Reference energy mismatch')
    require(not result['atoms'].info.get('fixture_only', False), 'Fixture is not a production reference')
    return result


def validate_reference(reference_path=None):
    path = Path(reference_path) if reference_path is not None else ROOT / REFERENCE
    ref = json.loads(path.read_text())
    require(ref.get('status') == 'passed', 'Reference baseline is not passed; no production input written')
    settings = settings_from_reference(ref['setting_selection'])
    require(ref['heterobilayer_numerical_sensitivity']['passed'] is True, 'Reference precision did not pass')
    files = {path}
    layers, equilibria = {}, {}
    for name in ('mos2', 'mose2'):
        for section, destination in [('equilibria', equilibria), ('common_lattice_layers', layers)]:
            summary = ref[section][name]
            result = verify_optimized(summary['job_id'], summary)
            expected_cut, expected_k = ((120, 12) if section == 'equilibria'
                                       else (settings['ecutwfc_Ry'], settings['primitive_reference_k']))
            require(result['meta']['ecutwfc_Ry'] == expected_cut
                    and result['meta']['k_grid'] == [expected_k, expected_k, 1], 'Wrong optimized reference settings')
            destination[name] = result['atoms']
            files.update(_raw_paths(result['meta'], ROOT / f"reports/{summary['job_id']}.json"))
    a_common = float(np.mean([a.cell.lengths()[0] for a in equilibria.values()]))
    require(abs(a_common - ref['common_lattice_angstrom']) <= 1e-7, 'Common lattice is not the mean')
    for name, atoms in layers.items():
        require(abs(atoms.cell.lengths()[0] - a_common) <= 1e-7, 'Wrong common-layer lattice')
        require(abs(a_common / equilibria[name].cell.lengths()[0] - 1 - ref['strain_fraction'][name]) <= 1e-8,
                'Strain provenance mismatch')
    d0 = ref['D0']['Mo_plane_distance_angstrom']
    require(abs(ref['D0']['force_residual_eV_angstrom']) <= .01, 'D0 residual exceeds criterion')
    bracket = ref['D0']['force_bracket_angstrom']
    require(len(bracket) == 2 and bracket[0] <= d0 <= bracket[1], 'D0 outside force bracket')
    cache = {}

    def check_row(row):
        values = {}
        for part, job_id in row['jobs'].items():
            if job_id not in cache:
                report = ROOT / f'reports/{job_id}.json'
                cache[job_id] = load_scf(report)
                files.update(_raw_paths(json.loads(report.read_text()), report))
            values[part] = cache[job_id]
        require(set(values) == {'pair', 'bottom', 'top'}, 'Missing energy component')
        pair = values['pair']['atoms']
        require(pair.get_chemical_symbols() == ['Mo', 'S', 'S', 'Mo', 'Se', 'Se'], 'Wrong reference composition')
        expected_pair = bilayer(layers['mos2'], layers['mose2'], row['registry_id'],
                                row['Mo_plane_distance_angstrom'], row['setting'][2])
        require(np.allclose(pair.positions, expected_pair.positions, atol=1e-7, rtol=0)
                and np.allclose(pair.cell, expected_pair.cell, atol=1e-7, rtol=0),
                'Reference registry or source-layer geometry mismatch')
        for part, indices in [('bottom', slice(0, 3)), ('top', slice(3, 6))]:
            require(values[part]['atoms'].get_chemical_symbols() == pair[indices].get_chemical_symbols(),
                    'Separated species changed')
            require(np.allclose(values[part]['atoms'].positions, pair.positions[indices], atol=1e-7, rtol=0)
                    and np.allclose(values[part]['atoms'].cell, pair.cell, atol=1e-7, rtol=0),
                    'Separated coordinates or cell changed')
        c, k, v = row['setting']
        for value in values.values():
            audit = value['audit']
            require(audit['ecutwfc_Ry'] == c and audit['ecutrho_Ry'] == 4 * c
                    and audit['k_grid'] == [k, k, 1] and audit['conv_thr_Ry'] == 1e-10,
                    'Reference component settings changed')
        require(abs(pair.cell[2, 2] - np.ptp(pair.positions[:, 2]) - v) <= 1e-7, 'Vacuum mismatch')
        require(abs(pair.positions[3, 2] - pair.positions[0, 2] - row['Mo_plane_distance_angstrom']) <= 1e-7,
                'Distance mismatch')
        area = float(np.linalg.norm(np.cross(pair.cell[0], pair.cell[1])))
        gamma = (values['pair']['energy_eV'] - values['bottom']['energy_eV'] - values['top']['energy_eV']) * 1000 / area
        require(abs(gamma - row['gamma_meV_angstrom2']) <= 1e-6, 'Gamma mismatch')
        require(abs(area - row['area_angstrom2']) <= 1e-7, 'Area mismatch')
        fz = float(values['pair']['forces'][3:, 2].sum())
        require(abs(fz - row['upper_layer_force_z_eV_angstrom']) <= 1e-8, 'Force mismatch')
        require(abs(float(values['pair']['forces'][:, 2].sum())) < .002, 'Net force mismatch')
        return gamma, fz

    grid = ref['nine_point_grid']
    require(len(grid) == 9, 'Incomplete reference grid')
    grid_keys, center = set(), None
    selected = [settings['ecutwfc_Ry'], settings['primitive_reference_k'], settings['vacuum_angstrom']]
    for row in grid:
        di = int(np.argmin(np.abs(np.array([d0 - .25, d0, d0 + .25]) - row['Mo_plane_distance_angstrom'])))
        require(abs(row['Mo_plane_distance_angstrom'] - [d0 - .25, d0, d0 + .25][di]) <= 1e-6,
                'Unexpected grid distance')
        key = (row['registry_id'], di)
        require(key not in grid_keys and key[0] in range(3) and row['setting'] == selected, 'Grid coverage/settings mismatch')
        grid_keys.add(key)
        gamma, fz = check_row(row)
        if key == (0, 1):
            center = gamma
            require(abs(fz) <= .01, 'Raw center force failed')
    require(len(grid_keys) == 9 and center is not None, 'Missing reference point')
    for row in grid:
        if row['registry_id'] == 0 and abs(row['Mo_plane_distance_angstrom'] - d0) > .1:
            require(row['gamma_meV_angstrom2'] > center, 'Distance minimum remains unresolved')
    strict = [r for r in ref['comparisons'] if r['tag'] == 'hetero-strict'
              and r['registry_id'] == 0 and abs(r['Mo_plane_distance_angstrom'] - d0) <= 1e-6]
    require(len(strict) == 1, 'Missing unique strict reference')
    require(strict[0]['setting'] == ([120, 12, 30] if selected != [120, 12, 30] else [140, 15, 35]),
            'Strict reference settings were weakened')
    strict_gamma, _ = check_row(strict[0])
    require(abs(strict_gamma - center) <= .1, 'Raw reference precision exceeds 0.1 meV/Angstrom2')
    # Verify the original force sign change, not just the stored D0 number.
    for distance, sign in zip(bracket, (1, -1)):
        rows = [r for r in ref['comparisons'] if r['tag'] == 'hetero' and r['registry_id'] == 0
                and abs(r['Mo_plane_distance_angstrom'] - distance) <= 1e-6 and r['setting'] == selected]
        require(len(rows) == 1 and sign * check_row(rows[0])[1] > 0, 'Raw force bracket failed')
    return dict(layers=layers, d0=d0, settings=settings,
                evidence_sha256={str(p.relative_to(ROOT)): sha(p) for p in sorted(files)},
                common_lattice_angstrom=a_common)


def render_scf(job_id, atoms, settings):
    atoms = atoms.copy()
    atoms.calc = None
    atoms.set_constraint()
    data = dict(control=dict(calculation='scf', prefix=job_id.replace('-', '_'), pseudo_dir='./pseudo',
                             outdir='./scratch', tprnfor=True, tstress=True),
                system=dict(ecutwfc=settings['ecutwfc_Ry'], ecutrho=settings['ecutrho_Ry'],
                            input_dft='vdw-df-c09', occupations='fixed', nspin=1),
                electrons=dict(conv_thr=settings['conv_thr_Ry'], electron_maxstep=200, mixing_beta=.4))
    stream = io.StringIO()
    write(stream, atoms, format='espresso-in', input_data=data,
          pseudopotentials={e: e + '.upf' for e in set(atoms.get_chemical_symbols())},
          kpts=tuple(settings['k_grid']), koffset=(0, 0, 0))
    result = stream.getvalue()
    recovered = read(io.StringIO(result), format='espresso-in')
    require(recovered.get_chemical_symbols() == atoms.get_chemical_symbols(), 'Species serialization failed')
    require(np.allclose(recovered.positions, atoms.positions, atol=1e-9, rtol=0)
            and np.allclose(recovered.cell, atoms.cell, atol=1e-9, rtol=0), 'Geometry serialization failed')
    return result


def prepare():
    target = ROOT / PLAN
    if target.exists():
        saved = json.loads(target.read_text())
        for rel, expected in saved['input_hashes'].items():
            require(sha(ROOT / rel) == expected, f'Changed frozen dependency: {rel}')
        return target
    source = validate_reference()
    settings = source['settings']
    executable = '.tools/qe-7.6/bin/pw.x'
    pseudos = {e: sha(ROOT / PSEUDOS / (e + '.upf')) for e in ('Mo', 'S', 'Se')}
    rows, calculations = build_catalog(source['layers']['mos2'], source['layers']['mose2'], source['d0'],
                                       settings, pseudos, sha(ROOT / executable))
    probe = next(row for row in rows if (row['configuration'], row['registry_id'], row['distance_index']) == ('CI', 0, 1))
    jobs = []
    hashes = dict(source['evidence_sha256'])
    hashes[executable] = sha(ROOT / executable)
    for script in ['prepare_completion_alloy_probe.py', 'completion_alloy_geometry.py', 'completion_qe_bridge.py',
                   'week3_qe.py', 'week3_materials_pipeline.py', 'prepare_week3_drug.py', 'research_run.py']:
        rel = 'scripts/' + script
        hashes[rel] = sha(ROOT / rel)
    for part, key in probe['calculations'].items():
        atoms = calculations[key]['atoms']
        job_id = f'w3-alloy-probe-r1-ci-r0-d0-{part}'
        structure = f'structures/completion/alloy-probe-r1/{job_id}.extxyz'
        stream = io.StringIO(); write(stream, atoms, format='extxyz')
        immutable(ROOT / structure, stream.getvalue())
        inp = f'inputs/qe/completion/{job_id}.in'
        immutable(ROOT / inp, render_scf(job_id, atoms, settings))
        assets = {f'{PSEUDOS}/{e}.upf': f'pseudo/{e}.upf' for e in sorted(set(atoms.get_chemical_symbols()))}
        assets.update({structure: 'initial.extxyz', 'scripts/completion_qe_bridge.py': 'qe_bridge_snapshot.py'})
        cfg = dict(protocol_id='completion-alloy-cost-probe-r1', week=3, task_id='S04-cost',
                   job_id=job_id, engine='qe', executable=executable, input=inp, assets=assets,
                   method='vdW-DF-C09; PseudoDojo nc-sr-04 standard PBE scalar NC', expected_bridge_method='C09',
                   requested=dict(optimization=False, frequencies=False), calculation='scf', n_atoms=len(atoms),
                   charge=0, SOC=False, ecutwfc_Ry=settings['ecutwfc_Ry'], ecutrho_Ry=settings['ecutrho_Ry'],
                   k_grid=settings['k_grid'], conv_thr_Ry=settings['conv_thr_Ry'], omp_threads=4,
                   timeout_seconds=3600, rss_limit_gib=40, geometry_source=structure,
                   initial_cell_angstrom=atoms.cell.array.tolist(),
                   initial_vacuum_gap_angstrom=float(atoms.cell[2, 2] - np.ptp(atoms.positions[:, 2])),
                   role=part, computational_identity_sha256=key,
                   configuration_id='CI', registry_id=0, Mo_plane_distance_angstrom=source['d0'],
                   accuracy_scope='24-atom cost measurement; alloy precision and full matrix remain unverified')
        path = f'experiments/completion/alloy-probe-r1/{job_id}.json'
        immutable(ROOT / path, json.dumps(cfg, indent=2) + '\n')
        jobs.append(path)
        for rel in [path, inp, *assets]: hashes[rel] = sha(ROOT / rel)
    plan = dict(kind='alloy-cost-probe-r1', jobs=jobs, input_hashes=hashes, settings=settings,
                automatic_execution=False, automatic_retries=0, maximum_total_seconds=10800,
                full_matrix_execution_enabled=False, scientific_completion=False,
                probe=probe, matrix_blueprint=rows, candidate_unique_calculations=len(calculations),
                candidate_pair_calculations=27, candidate_separated_layer_calculations=len(calculations) - 27,
                scope='Only one 24-atom CI structure and its two exact 12-atom references are executable',
                follow_up='Verify all three raw SCFs, measure elapsed time/RSS, and review remaining pilot budget before the full matrix')
    immutable(target, json.dumps(plan, indent=2) + '\n')
    return target


if __name__ == '__main__':
    print(prepare())
