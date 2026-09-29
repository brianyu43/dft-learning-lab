"""Freeze two small method/solvent compatibility jobs; never launch an engine."""
import json
from pathlib import Path

import numpy as np
from ase.io import read

from prepare_completion_drug import GRIDS, LIBS
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from run_week2_frequency_resume import inspect_hessian
from run_week2_xtb_screen import geometry_audit

PLAN = 'experiments/completion/orca-method-solvent-compatibility-r1.json'
SOURCES = {
    'composite_methods': 'https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/3cmethods.html',
    'solvation': 'https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html',
}


def water_source():
    report = ROOT / 'reports/week1-water-opt-freq.json'
    meta = json.loads(report.read_text())
    assert meta['status'] == 'passed' and meta['result']['minimum_verified']
    run = Path(meta['run_dir'])
    assert sha(run / 'output.log') == meta['output_sha256']
    assert sha(run / Path(meta['input']).name) == meta['source_sha256'][meta['input']]
    assert sha(run / 'runner_snapshot.py') == meta['runner_sha256']
    stem = Path(meta['input']).stem
    xyz, hessian = run / (stem + '.xyz'), run / (stem + '.hess')
    atoms = read(xyz)
    assert atoms.numbers.tolist() == [8, 1, 1] and np.isfinite(atoms.positions).all()
    frequencies = inspect_hessian(hessian, atoms, meta['result']['frequencies_cm1'])
    assert len(frequencies) == 3 and np.all(frequencies > 0)
    prior_path = ROOT / 'reports/week1-verification.json'
    prior = json.loads(prior_path.read_text())
    proof = next(row for row in prior['checks'] if row['species'] == 'water')
    assert sha(xyz) == proof['optimized_xyz_sha256']
    assert sha(hessian) == proof['hessian_sha256']
    return xyz, atoms, [report, prior_path, xyz, hessian, run / 'output.log', run / Path(meta['input']).name]


def prepare():
    target = ROOT / PLAN
    if target.exists():
        saved = json.loads(target.read_text())
        for rel, expected in saved['input_hashes'].items():
            assert sha(ROOT / rel) == expected, rel
        return target
    xyz, atoms, source_files = water_source()
    candidate = dict(candidate_id='water-compatibility-r1', species='water', kind='monomer',
                     charge=0, multiplicity=1,
                     graph=dict(atomic_numbers=[8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]],
                                hydrogen_owners={'1': 0, '2': 0}),
                     purpose='Runtime compatibility only; excluded from the four-case drug panel')
    assert geometry_audit(candidate, atoms)['chemical_state_preserved_by_geometry']
    cp = 'experiments/completion/orca-compatibility/water-candidate.json'
    immutable(ROOT / cp, json.dumps(candidate, indent=2) + '\n')
    jobs = []
    g = GRIDS['dense']
    definitions = [
        ('wb97x3c-water-engrad', 'wB97X-3c', 'vacuum', False, 'electronic_energy'),
        ('r2scan3c-smdethanol-water-optfreq', 'r2SCAN-3c', 'SMD(ethanol)', True,
         'orca_smd_total_energy_including_cds'),
    ]
    for label, method, solvent, optimize, quantity in definitions:
        job_id = 'w3-compat-r1-' + label
        geometry = f'structures/completion/orca-compatibility/{job_id}.xyz'
        immutable(ROOT / geometry, xyz.read_text())
        solvation_keyword = ' SMD(ethanol)' if solvent != 'vacuum' else ''
        task = 'VeryTightOpt Opt Freq' if optimize else 'EnGrad'
        solvent_block = '%cpcm\n  SMD18 false\n  draco false\nend\n' if solvent != 'vacuum' else ''
        text = (f"! {method} VeryTightSCF DefGrid3{solvation_keyword} {task}\n"
                f"%maxcore 2000\n%pal nprocs 12 end\n{g['block']}{solvent_block}\n* xyzfile 0 1 initial.xyz\n")
        inp = f'inputs/orca/completion/{job_id}.inp'
        immutable(ROOT / inp, text)
        assets = {f'.tools/orca-mpi/lib/{lib}': lib for lib in LIBS}
        assets.update({geometry: 'initial.xyz', cp: 'candidate.json',
                       'scripts/completion_orca.py': 'completion_parser_snapshot.py',
                       'scripts/week3_recovery.py': 'provenance/gradient_check_snapshot.py',
                       'scripts/run_week2_frequency_resume.py': 'provenance/hessian_check_snapshot.py',
                       'scripts/run_week2_xtb_screen.py': 'provenance/geometry_audit_snapshot.py',
                       'scripts/orca_parallel.sh': 'provenance/launcher_snapshot.sh',
                       'scripts/prepare_completion_orca_compatibility.py': 'provenance/preparer_snapshot.py',
                       'environment/orca-mpi-explicit.txt': 'provenance/orca-mpi-explicit.txt'})
        cfg = dict(protocol_id='completion-orca-method-solvent-compatibility-r1', week=3, task_id='D03-D04-runtime',
                   job_id=job_id, engine='orca', executable='scripts/orca_parallel.sh', input=inp, assets=assets,
                   method=method, solvent=solvent, energy_quantity=quantity, charge=0, multiplicity=1,
                   n_atoms=3, linear=False, requested=dict(optimization=optimize, frequencies=optimize, gradient=True),
                   mpi_processes=12, omp_threads=1, timeout_seconds=1800, rss_limit_gib=40,
                   candidate_id=candidate['candidate_id'], species='water', role='runtime_compatibility_not_pilot',
                   grid_setting='dense', expected_scf_grid=dict(IntAcc=6., AngularGrid=7),
                   initial_structure_sha256=sha(ROOT / geometry), candidate_sha256=sha(ROOT / cp),
                   completion_parser_sha256=sha(ROOT / 'scripts/completion_orca.py'), expected_input_sha256=sha(ROOT / inp),
                   initial_geometry_audit=geometry_audit(candidate, atoms),
                   source_geometry=str(xyz.relative_to(ROOT)), source_geometry_sha256=sha(xyz),
                   main_orca_sha256=sha(ROOT / '.tools/orca-6.1.1/orca'),
                   new_independent_start=False, contributes_to_pilot_scores=False,
                   purpose='Verify runtime options, gradient/geometry transfer and SMD frequency support; no accuracy or full-molecule cost claim')
        config = f'experiments/completion/orca-compatibility/{job_id}.json'
        immutable(ROOT / config, json.dumps(cfg, indent=2) + '\n'); jobs.append(config)
    paths = {str(p.relative_to(ROOT)) for p in source_files}
    paths.update(['scripts/prepare_completion_orca_compatibility.py', 'scripts/prepare_completion_drug.py',
                  'scripts/completion_batch.py', '.tools/orca-6.1.1/orca'])
    for job in jobs:
        cfg = json.loads((ROOT / job).read_text())
        paths.update([job, cfg['input'], *cfg['assets']])
    plan = dict(kind='orca-method-solvent-compatibility-r1', jobs=jobs,
                input_hashes={rel: sha(ROOT / rel) for rel in sorted(paths)},
                automatic_execution=False, automatic_retries=0, maximum_total_seconds=3600,
                execution_precondition='Only after the active controller releases its engine lock, within the existing drug budget',
                sources=SOURCES, documentation_checked_on='2026-09-29',
                required_post_run_review=[
                    'Raw normal termination, charge/spin, MPI count, same input coordinate mapping and engrad energy/gradient',
                    'wB97X-3c composite setup including vDZP/ECP and adapted D4; no manually added basis/gCP/dispersion',
                    'Ethanol solvent header, SMD-CDS descriptors and energy components; CDS must not be added twice',
                    'SMD optimized geometry, analytic Hessian, all three internal frequencies and coordinate agreement',
                    'Observed runtime/peak RSS are for water and cannot estimate the full drug-panel cost'],
                runtime_compatibility_verified=False, full_panel_execution_enabled=False,
                scope='Two small runtime checks only; no full-panel method/solvent sensitivity, numerical accuracy, or formation conclusion')
    immutable(target, json.dumps(plan, indent=2) + '\n')
    return target


if __name__ == '__main__':
    print(prepare())
