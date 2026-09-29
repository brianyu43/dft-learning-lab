"""Prepare three bounded checks; never launch or interrupt a research engine."""
import io
import json
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.io.espresso import read_fortran_namelist

from completion_qe_bridge import load_scf, require
from prepare_week3_drug import immutable
from research_run import ROOT, sha

PSEUDOS = 'data/pseudopotentials/pseudo-dojo-nc-sr-04-pbe-standard'


def prepare_job(job_id, atoms, method, source_report, purpose, cut=80, k=6, conv_thr=1e-12, **details):
    atoms = atoms.copy()
    atoms.calc = None
    atoms.set_constraint()
    system = dict(ecutwfc=cut, ecutrho=4 * cut, occupations='fixed', nspin=1,
                  input_dft='vdw-df-c09' if method == 'C09' else 'PBE')
    require(method in ('C09', 'PBE-D3BJ'), 'Unreviewed method')
    if method == 'PBE-D3BJ':
        system.update(vdw_corr='grimme-d3', dftd3_version=4, dftd3_threebody=True)
    data = dict(control=dict(calculation='scf', prefix=job_id.replace('-', '_'),
                             pseudo_dir='./pseudo', outdir='./scratch', tprnfor=True, tstress=True),
                system=system, electrons=dict(conv_thr=conv_thr, electron_maxstep=200, mixing_beta=0.4))
    elements = sorted(set(atoms.get_chemical_symbols()))
    require(set(elements) <= {'Mo', 'S', 'Se'}, 'Unsupported element')
    structure = ROOT / f'structures/completion/qe-checks/{job_id}.extxyz'
    stream = io.StringIO()
    write(stream, atoms, format='extxyz')
    immutable(structure, stream.getvalue())
    inp = ROOT / f'inputs/qe/completion/{job_id}.in'
    stream = io.StringIO()
    write(stream, atoms, format='espresso-in', pseudopotentials={e: e + '.upf' for e in elements},
          kpts=(k, k, 1), koffset=(0, 0, 0), input_data=data)
    immutable(inp, stream.getvalue())
    recovered = read(inp, format='espresso-in')
    require(recovered.get_chemical_symbols() == atoms.get_chemical_symbols(), 'Input species round trip failed')
    require(np.allclose(recovered.positions, atoms.positions, atol=1e-9, rtol=0), 'Input coordinates round trip failed')
    require(np.allclose(recovered.cell, atoms.cell, atol=1e-9, rtol=0), 'Input cell round trip failed')
    with inp.open() as handle:
        parsed, _ = read_fortran_namelist(handle)
    if method == 'PBE-D3BJ':
        require(parsed['system']['dftd3_version'] == 4 and parsed['system']['dftd3_threebody'] is True,
                'D3(BJ) settings did not survive serialization')
    config = dict(protocol_id='dft-pilot-20260909-v1-completion-qe-checks-r1', week=3,
                  task_id='SR-force-and-method-compatibility', job_id=job_id, engine='qe',
                  executable='.tools/qe-7.6/bin/pw.x', input=str(inp.relative_to(ROOT)),
                  assets={**{f'{PSEUDOS}/{e}.upf': f'pseudo/{e}.upf' for e in elements},
                          str(structure.relative_to(ROOT)): 'initial.extxyz',
                          'scripts/completion_qe_bridge.py': 'qe_bridge_snapshot.py'},
                  method=('vdW-DF-C09' if method == 'C09' else 'PBE-D3(BJ), three-body enabled')
                         + '; PseudoDojo nc-sr-04 standard PBE scalar NC',
                  expected_bridge_method=method, requested=dict(optimization=False, frequencies=False),
                  calculation='scf', n_atoms=len(atoms), charge=0, SOC=False,
                  ecutwfc_Ry=cut, ecutrho_Ry=4 * cut, k_grid=[k, k, 1], conv_thr_Ry=conv_thr,
                  initial_cell_angstrom=atoms.cell.array.tolist(),
                  initial_vacuum_gap_angstrom=float(atoms.cell[2, 2] - np.ptp(atoms.positions[:, 2])),
                  geometry_source=str(structure.relative_to(ROOT)), parent_report=source_report,
                  parent_report_sha256=sha(ROOT / source_report), purpose=purpose,
                  timeout_seconds=1800, rss_limit_gib=40, omp_threads=4,
                  accuracy_scope='Force transfer and method compatibility only; not CI/CO method sensitivity', **details)
    path = ROOT / f'experiments/completion/qe-checks/{job_id}.json'
    immutable(path, json.dumps(config, indent=2) + '\n')
    return path


def main():
    dimer_report = 'reports/w2-bilayer-selected-r0-d6p3-pair.json'
    mono_report = 'reports/w2-mono-combined-c80-k6-v20.json'
    pair = load_scf(ROOT / dimer_report)
    mono = load_scf(ROOT / mono_report)
    jobs = []
    displacement = 0.002
    index, axis = 1, 2
    for sign, suffix in [(-1, 'minus'), (1, 'plus')]:
        atoms = pair['atoms'].copy()
        atoms.positions[index, axis] += sign * displacement
        jobs.append(prepare_job(f'w3-qe-bridge-s-{suffix}', atoms, 'C09', dimer_report,
                                purpose='Central energy derivative for an unconstrained chalcogen direction',
                                displaced_atom_index=index, displaced_cartesian_axis=axis,
                                displacement_angstrom=sign * displacement))
    jobs.append(prepare_job('w3-qe-pbe-d3bj-compatibility', mono['atoms'], 'PBE-D3BJ', mono_report,
                            purpose='Verify real QE 7.6 PBE-D3(BJ) with the existing PBE pseudopotentials'))
    sources = ['scripts/completion_qe_bridge.py', 'scripts/prepare_completion_qe_checks.py',
               '.tools/src/qe-7.6/PW/Doc/INPUT_PW.def', '.tools/src/qe-7.6/Modules/constants.f90',
               'environment/qe-commit.txt', dimer_report, mono_report]
    for report in (dimer_report, mono_report):
        parent = Path(json.loads((ROOT / report).read_text())['run_dir'])
        xml = list((parent / 'scratch').glob('*.save/data-file-schema.xml'))[0]
        sources.append(str(xml.relative_to(ROOT)))
    for path in jobs:
        cfg = json.loads(path.read_text())
        sources.extend([str(path.relative_to(ROOT)), cfg['input'], *cfg['assets'].keys()])
    plan = dict(plan_id='completion-qe-force-and-method-checks-r1', kind='qe_force_and_method_compatibility',
                execution_enabled=False, automatic_execution=False,
                execution_precondition='Only after the active controller releases its engine lock and budget is checked',
                jobs=[str(path.relative_to(ROOT)) for path in jobs],
                input_hashes={path: sha(ROOT / path) for path in sorted(set(sources))},
                maximum_total_seconds=5400,
                measured_parent_seconds=dict(pair=133.097, monolayer=76.467),
                force_check=dict(parent_report=dimer_report, atom_index=index, cartesian_axis=axis,
                                 step_angstrom=displacement, center_force_eV_angstrom=float(pair['forces'][index, axis]),
                                 absolute_force_difference_limit_eV_angstrom=1e-4,
                                 compare='- [E(+h)-E(-h)]/(2h) to original XML force',
                                 on_failure='Review SCF, step-size and force transfer; do not loosen the criterion automatically'),
                method_check=dict(method='PBE-D3BJ', dftd3_version=4, dftd3_threebody=True,
                                  expected_xml_and_force_verification=True,
                                  scope='Runtime option compatibility, not method accuracy or interlayer transfer'),
                scientific_completion=False)
    path = ROOT / 'experiments/completion/qe-force-and-method-checks-r1.json'
    immutable(path, json.dumps(plan, indent=2) + '\n')
    print(json.dumps(dict(plan=str(path.relative_to(ROOT)), jobs=len(jobs),
                         input_hashes=len(plan['input_hashes']), execution_started=False), indent=2))


if __name__ == '__main__':
    main()
