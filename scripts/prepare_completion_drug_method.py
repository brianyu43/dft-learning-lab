"""Freeze wB97X-3c EnGrad jobs at validated r2SCAN-3c panel coordinates.

No engine is launched. Full-panel precision and actual method compatibility
must be established before any production input is written. Both grids retain
all original 17 start identities; numerical repeats are not new minima.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from ase.io import read

import completion_drug_panel as panel
import evaluate_completion_orca_compatibility as compatibility
from prepare_completion_drug import GRIDS, LIBS
from prepare_completion_orca_compatibility import PLAN as COMPAT_PLAN
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from run_week2_xtb_screen import geometry_audit

PLAN = 'experiments/completion/wb97x3c-at-uniform-r1.json'
SCORES = 'reports/completion/uniform-vacuum-r4-scores.json'
REVIEW = 'reports/completion/orca-wb97x3c-compatibility-review-r1.json'
COMPAT_CONFIG = 'experiments/completion/orca-compatibility/w3-compat-r1-wb97x3c-water-engrad.json'
TOPICS = ('composite_method', 'basis_ecp', 'dispersion_parameters', 'energy_components')
CONTEXTS = {grid: dict(method='wB97X-3c', environment='vacuum',
                      numerical_protocol_id=f'wb97x3c-{grid}-at-uniform-r1',
                      energy_quantity='electronic_energy') for grid in ('dense', 'strict')}
DEPENDENCIES = [
    'scripts/prepare_completion_drug_method.py', 'scripts/prepare_completion_drug.py',
    'scripts/prepare_week3_drug.py', 'scripts/completion_drug_panel.py',
    'scripts/completion_drug_contrasts.py', 'scripts/evaluate_completion_orca_compatibility.py',
    'scripts/completion_batch.py', 'scripts/research_run.py', '.tools/orca-6.1.1/orca',
]


def require(value, message):
    if not value:
        raise ValueError(message)


def verify_hashes(hashes):
    for rel, expected in hashes.items():
        path = ROOT / rel
        require(path.is_file() and sha(path) == expected, f'Changed or missing evidence: {rel}')


def validate_runtime_review(review, actual, runtime_log, plan_digest, config_digest):
    """Require raw revalidation plus an explicit, traceable parameter review.

    Quoted evidence makes the review inspectable; it is not an automated proof
    of the scientific interpretation. Water does not exercise all drug atoms.
    """
    require(review.get('status') == 'wb97x3c_runtime_reviewed', 'Detailed method review is incomplete')
    require(review.get('compatibility_plan_sha256') == plan_digest
            and review.get('config_path') == COMPAT_CONFIG
            and review.get('config_sha256') == config_digest, 'Compatibility review provenance differs')
    require(review.get('source_sha256') == actual['source_sha256'], 'Compatibility raw evidence differs')
    require(actual['method_evidence']['method'] == 'wB97X-3c'
            and actual['method_evidence']['orbital_basis'].lower() == 'vdzp'
            and actual['solvent'] == 'vacuum' and actual['energy_quantity'] == 'electronic_energy',
            'Wrong compatibility method or energy quantity')
    require('****END OF INPUT****' in runtime_log, 'Missing runtime section')
    runtime = runtime_log.split('****END OF INPUT****', 1)[1]
    for topic in TOPICS:
        row = review.get('checks', {}).get(topic, {})
        require(row.get('accepted') is True and isinstance(row.get('note'), str) and row['note'].strip(),
                f'Missing detailed review: {topic}')
        quotes = row.get('runtime_quotes')
        require(isinstance(quotes, list) and quotes and all(
            isinstance(q, str) and q.strip() and q in runtime for q in quotes),
            f'Missing actual runtime quotation: {topic}')
    require(review.get('water_only_not_drug_accuracy_or_cost') is True,
            'Water compatibility scope must remain explicit')


def prerequisites():
    # Check both prerequisite files before the raw full-panel evaluator can
    # write its immutable scientific report. No partial production preparation.
    require((ROOT / SCORES).is_file() and (ROOT / REVIEW).is_file(),
            'Uniform numerical panel and detailed wB97X-3c review are required; no inputs written')
    compat_plan = json.loads((ROOT / COMPAT_PLAN).read_text())
    verify_hashes(compat_plan['input_hashes'])
    require(COMPAT_CONFIG in compat_plan['jobs'], 'Compatibility job is outside the frozen plan')
    actual = compatibility.audit_job(COMPAT_CONFIG)
    meta = json.loads((ROOT / 'reports/w3-compat-r1-wb97x3c-water-engrad.json').read_text())
    log = (Path(meta['run_dir']) / 'output.log').read_text()
    review = json.loads((ROOT / REVIEW).read_text())
    validate_runtime_review(review, actual, log, sha(ROOT / COMPAT_PLAN), sha(ROOT / COMPAT_CONFIG))
    # The status string alone is insufficient: reread all minima, gradients,
    # Hessians and the strict-grid contrasts through the existing evaluator.
    evaluated = panel.evaluate()
    require(evaluated.resolve() == (ROOT / SCORES).resolve(), 'Unexpected uniform score report')
    scores = json.loads(evaluated.read_text())
    require(scores['status'] == 'uniform_vacuum_numerical_panel_passed'
            and scores['numerical_precision_verified'] is True, 'Full uniform numerical panel did not pass')
    parent, rows = panel.load_uniform_minima()
    require(scores['parent_plan_sha256'] == sha(ROOT / panel.PANEL), 'Uniform parent digest differs')
    hashes = {p: sha(ROOT / p) for p in (SCORES, REVIEW, COMPAT_PLAN, COMPAT_CONFIG, panel.PANEL, panel.CHECKS)}
    hashes.update(compat_plan['input_hashes'])
    hashes.update(actual['source_sha256'])
    hashes.update(parent['input_hashes'])
    checks = json.loads((ROOT / panel.CHECKS).read_text())
    hashes.update(checks['input_hashes'])
    for row in scores['rows']:
        hashes.update(row['source_sha256'])
    return rows, hashes


def validate_rows(rows):
    require(len(rows) == 17, 'The original 17 starts must remain represented')
    lineage = json.loads((ROOT / 'experiments/completion/drug-start-lineage-r1.json').read_text())
    expected = {r['candidate_id']: r for r in lineage['rows']}
    seen = set()
    for row in rows:
        mapping, result = row['mapping'], row['result']
        cid = mapping['candidate_id']
        require(cid in expected and cid not in seen, 'Duplicate or unknown start')
        seen.add(cid)
        original = expected[cid]
        require(mapping['species'] == original['species']
                and mapping['selection_role'] == original['selection_role']
                and result['candidate'] == original['candidate'], 'Start identity or selection role changed')
        source = Path(result['geometry_path']).resolve()
        require(source.is_relative_to(ROOT.resolve()) and sha(source) == result['geometry_sha256'],
                'Source coordinate identity changed')
        atoms = read(source)
        require(np.isfinite(atoms.positions).all()
                and atoms.numbers.tolist() == result['candidate']['graph']['atomic_numbers'],
                'Source atom order or coordinates changed')
        require(result['candidate']['charge'] == 0 and result['candidate']['multiplicity'] == 1,
                'Unexpected charge or multiplicity')
        audit = geometry_audit(result['candidate'], atoms)
        require(audit['heavy_connectivity_preserved'] and not audit['unassigned_hydrogens']
                and not audit.get('separated_geometry', False), 'Source chemical state requires review')
    require(seen == set(expected), 'Incomplete frozen start coverage')


def input_text(grid):
    require(grid in CONTEXTS, 'Only the two frozen method grids are supported')
    return (f"! wB97X-3c VeryTightSCF DefGrid3 EnGrad\n"
            f"%maxcore 2000\n%pal nprocs 12 end\n{GRIDS[grid]['block']}\n* xyzfile 0 1 initial.xyz\n")


def write_job(row, grid):
    mapping, result = row['mapping'], row['result']
    cid = mapping['candidate_id']; source = Path(result['geometry_path'])
    require(sha(source) == result['geometry_sha256'], 'Source changed during preparation')
    job_id = f'w3-method-r1-{cid}-{grid}-engrad'
    base = 'experiments/completion/wb97x3c-r1'
    xyz = f'structures/completion/wb97x3c-r1/{cid}.xyz'
    cp = f'{base}/candidates/{cid}.json'
    inp = f'inputs/orca/completion/wb97x3c-r1/{job_id}.inp'
    immutable(ROOT / xyz, source.read_text())
    require(sha(ROOT / xyz) == result['geometry_sha256'], 'Coordinate copy is not byte-identical')
    immutable(ROOT / cp, json.dumps(result['candidate'], indent=2) + '\n')
    immutable(ROOT / inp, input_text(grid))
    assets = {f'.tools/orca-mpi/lib/{lib}': lib for lib in LIBS}
    assets.update({xyz: 'initial.xyz', cp: 'candidate.json',
        SCORES: 'provenance/uniform-vacuum-scores.json', REVIEW: 'provenance/wb97x3c-compatibility-review.json',
        'scripts/completion_orca.py': 'completion_parser_snapshot.py',
        'scripts/week3_recovery.py': 'provenance/gradient_check_snapshot.py',
        'scripts/run_week2_frequency_resume.py': 'provenance/hessian_check_snapshot.py',
        'scripts/run_week2_xtb_screen.py': 'provenance/geometry_audit_snapshot.py',
        'scripts/prepare_completion_drug_method.py': 'provenance/method_preparer_snapshot.py',
        'scripts/orca_parallel.sh': 'provenance/launcher_snapshot.sh',
        'environment/orca-mpi-explicit.txt': 'provenance/orca-mpi-explicit.txt'})
    cfg = dict(protocol_id='completion-wb97x3c-at-uniform-r1', week=3, task_id='D03',
        job_id=job_id, engine='orca', executable='scripts/orca_parallel.sh', input=inp, assets=assets,
        method='wB97X-3c', solvent='vacuum', energy_quantity='electronic_energy', charge=0, multiplicity=1,
        n_atoms=len(result['candidate']['graph']['atomic_numbers']), linear=False,
        requested=dict(optimization=False, frequencies=False, gradient=True),
        mpi_processes=12, omp_threads=1, timeout_seconds=1800, rss_limit_gib=40,
        candidate_id=cid, species=mapping['species'], comparison_role=mapping['selection_role'],
        role='method_sensitivity_at_fixed_primary_method_minimum', grid_setting=grid,
        comparison_context=CONTEXTS[grid], geometry_origin_context=panel.CONTEXT,
        expected_scf_grid={k: GRIDS[grid][k] for k in ('IntAcc', 'AngularGrid')},
        source_geometry=str(source.relative_to(ROOT)), source_geometry_sha256=result['geometry_sha256'],
        initial_structure_sha256=sha(ROOT / xyz), candidate_sha256=sha(ROOT / cp),
        parent_minimum_job_id=result['meta']['job_id'], parent_minimum_config=mapping['config'],
        completion_parser_sha256=sha(ROOT / 'scripts/completion_orca.py'), expected_input_sha256=sha(ROOT / inp),
        main_orca_sha256=sha(ROOT / '.tools/orca-6.1.1/orca'),
        initial_geometry_audit=geometry_audit(result['candidate'], read(source)),
        new_independent_start=False, selection_uses_formation_labels=False,
        verifies_minimum_under_this_method=False,
        purpose='Fixed r2SCAN-3c coordinates; composite-method and numerical sensitivity, not wB97X-3c optimization')
    config = f'{base}/{job_id}.json'
    immutable(ROOT / config, json.dumps(cfg, indent=2) + '\n')
    return config


def prepare():
    target = ROOT / PLAN
    if target.exists():
        saved = json.loads(target.read_text())
        verify_hashes(saved['input_hashes'])
        require(saved['kind'] == 'wb97x3c-at-uniform-r1' and saved['comparison_contexts'] == CONTEXTS
                and saved['geometry_origin_context'] == panel.CONTEXT
                and saved['numerical_acceptance'] == panel.POLICY, 'Saved method policy changed')
        mapping = saved['comparison_mapping']
        require(len(saved['jobs']) == len(set(saved['jobs'])) == len(mapping) == 34
                and set(saved['jobs']) == {row['config'] for row in mapping}, 'Saved job coverage changed')
        lineage = json.loads((ROOT / 'experiments/completion/drug-start-lineage-r1.json').read_text())
        identities = {r['candidate_id']: r for r in lineage['rows']}
        require({(r['candidate_id'], r['grid']) for r in mapping}
                == {(cid, grid) for cid in identities for grid in CONTEXTS}, 'Saved start/grid coverage changed')
        for row in mapping:
            cfg = json.loads((ROOT / row['config']).read_text())
            original = identities[row['candidate_id']]
            require(all(row[k] == value for k, value in dict(candidate_id=cfg['candidate_id'],
                species=cfg['species'], selection_role=cfg['comparison_role'], grid=cfg['grid_setting'],
                geometry_sha256=cfg['source_geometry_sha256'],
                parent_minimum_config=cfg['parent_minimum_config']).items())
                and row['species'] == original['species'] and row['selection_role'] == original['selection_role'],
                'Saved comparison mapping changed')
        return target
    rows, hashes = prerequisites()
    validate_rows(rows)
    # Resolve every dependency before writing inputs, including the full
    # original lineage; incomplete setup cannot leave a usable production plan.
    dependencies = [*DEPENDENCIES, 'experiments/completion/drug-start-lineage-r1.json']
    hashes.update({p: sha(ROOT / p) for p in dependencies})
    jobs, mapping = [], []
    for grid in ('dense', 'strict'):
        for row in rows:
            config = write_job(row, grid); jobs.append(config)
            cfg = json.loads((ROOT / config).read_text())
            hashes.update({p: sha(ROOT / p) for p in [config, cfg['input'], *cfg['assets']]})
            mapping.append(dict(candidate_id=cfg['candidate_id'], species=cfg['species'],
                selection_role=cfg['comparison_role'], grid=grid, config=config,
                geometry_sha256=cfg['source_geometry_sha256'], parent_minimum_config=cfg['parent_minimum_config']))
    require(len(jobs) == len(set(jobs)) == 34, 'Changed method/numerical comparison scope')
    plan = dict(kind='wb97x3c-at-uniform-r1', jobs=jobs, comparison_mapping=mapping,
        input_hashes=hashes, automatic_execution=False, automatic_retries=0,
        original_start_count=17, chemical_species_count=14, new_independent_start_count=0,
        dense_single_points=17, strict_single_points=17, comparison_contexts=CONTEXTS,
        geometry_origin_context=panel.CONTEXT, full_method_minima_verified=False,
        max_per_job_seconds=1800, sum_of_requested_job_timeouts_seconds=34 * 1800,
        execution_precondition='Existing drug pilot budget and single-engine lock; no full-batch budget reservation implied',
        numerical_acceptance=panel.POLICY,
        required_post_run_review=[
            'Raw runtime composite/basis/ECP/D4 and vacuum evidence for every drug calculation',
            'Byte-identical parent coordinates and valid engrad coordinate/energy mapping',
            'All B/X grid changes and gradient-component differences within the frozen targets',
            'Primary-only and all-pose scopes, baseline-fixed and method-reselected pose scores',
            'All four shared-case contrasts at tie widths 0, 0.5, 1 and 2 kcal/mol'],
        method_sensitivity_complete=False, numerical_precision_verified=False,
        scope='Input freeze only; wB97X-3c SP at r2SCAN-3c minima, no crystal or formation-free-energy conclusion')
    verify_hashes(hashes)
    immutable(target, json.dumps(plan, indent=2) + '\n')
    return target


if __name__ == '__main__':
    print(prepare())
