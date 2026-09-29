"""Read-only raw ORCA validation and full-panel numerical score assembly.

Never reparses through research_run.parse_output: that function overwrites
historical result.json. This reader only writes a new requested report.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
from ase.io import read

from completion_drug_contrasts import analyze_uniform_records, panel_scores
from prepare_completion_drug import GRIDS, prepare_job
from prepare_week3_drug import immutable
from research_run import ROOT, sha
from run_week2_frequency_resume import inspect_hessian
from run_week2_xtb_screen import geometry_audit
from week3_recovery import gradient_in_frame

PANEL = 'experiments/completion/uniform-vacuum-r4.json'
CHECKS = 'experiments/completion/uniform-vacuum-r4-checks.json'
CONTEXT = dict(method='r2SCAN-3c', environment='vacuum',
               numerical_protocol_id='uniform-dense-r4', energy_quantity='electronic_energy')
STRICT_CONTEXT = {**CONTEXT, 'numerical_protocol_id': 'uniform-strict-at-dense-r4'}
POLICY = dict(max_score_change_kcal_mol=.1, max_gradient_component_change_Hartree_Bohr=1e-4)
NUMBER = r'[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_hashes(plan):
    for path, expected in plan['input_hashes'].items():
        require(sha(ROOT / path) == expected, f'Changed frozen input: {path}')


def expected_input(grid, minimum):
    g = GRIDS[grid]
    task = 'VeryTightOpt Opt Freq' if minimum else 'EnGrad'
    return (f"! r2SCAN-3c VeryTightSCF {g['keywords']} {task}\n"
            f"%maxcore 2000\n%pal nprocs 12 end\n{g['block']}\n* xyzfile 0 1 initial.xyz\n")


def audit_raw_result(report_path, *, input_text, grid, minimum):
    """Re-read a saved result, without modifying any source or result file.

    input_text is supplied by the frozen protocol caller. The lower-level
    reader can also inspect historical protocols without approving them for
    the new uniform panel.
    """
    report_path = Path(report_path)
    meta = json.loads(report_path.read_text())
    require(meta['engine'] == 'orca' and meta['status'] == 'passed'
            and meta['exit_code'] == 0 and not meta.get('stop_reason'), 'Run is not normally completed')
    require(meta['method'] == 'r2SCAN-3c' and meta['solvent'] == 'vacuum'
            and meta['charge'] == 0 and meta['multiplicity'] == 1, 'Wrong method/environment/electronic state')
    require(meta['requested'].get('optimization', False) == minimum
            and meta['requested'].get('frequencies', False) == minimum, 'Wrong result type')
    run = Path(meta['run_dir'])
    cfg_path = run / 'job_config.json'
    require(sha(cfg_path) == meta['config_sha256'], 'Archived config hash mismatch')
    config = json.loads(cfg_path.read_text())
    require(all(meta.get(k) == value for k, value in config.items()), 'Config/report settings differ')
    require(sha(run / 'runner_snapshot.py') == meta['runner_sha256'], 'Runner snapshot hash mismatch')
    require(sha(run / 'output.log') == meta['output_sha256'], 'Output hash mismatch')
    paths = [report_path, cfg_path, run / 'runner_snapshot.py', run / 'output.log']
    for source, expected in meta['source_sha256'].items():
        destination = Path(meta['input']).name if source == meta['input'] else meta['assets'][source]
        path = (run / destination).resolve()
        require(path.is_relative_to(run.resolve()) and sha(path) == expected, 'Input/asset hash mismatch')
        paths.append(path)
    inp = run / Path(meta['input']).name
    require(inp.read_text() == input_text, 'Input differs from the frozen numerical protocol')
    log = (run / 'output.log').read_text()
    require('ORCA TERMINATED NORMALLY' in log and 'SCF CONVERGED' in log
            and 'SCF NOT CONVERGED' not in log, 'Raw termination or SCF failure')
    version = re.search(r'Program Version\s+(\S+)', log)
    require(version is not None and version.group(1) == '6.1.1', 'Unreviewed ORCA version')
    ranks = [int(x) for x in re.findall(r'running with\s+(\d+)\s+parallel MPI-processes', log, re.I)]
    require(ranks and set(ranks) == {12} and meta['mpi_processes'] == 12, 'Actual MPI count differs')
    accuracies = re.findall(r'General Integration Accuracy\s+IntAcc\s+\.\.\.\s+([\d.]+)', log)
    angular = re.findall(r'Angular Grid \(max\. ang\.\)\s+AngularGrid\s+\.\.\.\s+(\d+)', log)
    require(accuracies and abs(float(accuracies[0]) - GRIDS[grid]['IntAcc']) < .0006
            and angular and int(angular[0]) == GRIDS[grid]['AngularGrid'], 'Actual SCF grid differs')
    energies = [float(x.replace('D', 'E').replace('d', 'e'))
                for x in re.findall(r'FINAL SINGLE POINT ENERGY\s+(' + NUMBER + ')', log)]
    require(energies and np.isfinite(energies).all() and meta['result']['energy_unit'] == 'Hartree'
            and abs(energies[-1] - meta['result']['energy']) < 1e-9, 'Raw energy mismatch')
    geometry = run / (Path(meta['input']).stem + '.xyz') if minimum else run / 'initial.xyz'
    atoms = read(geometry)
    require(np.isfinite(atoms.positions).all() and len(atoms) == meta['n_atoms'], 'Invalid coordinates')
    initial = read(run / 'initial.xyz')
    require(np.array_equal(initial.numbers, atoms.numbers), 'Atom order/composition changed')
    require(sha(run / 'initial.xyz') == meta['initial_structure_sha256'], 'Initial geometry hash mismatch')
    candidate_path = run / 'candidate.json'
    require(sha(candidate_path) == meta['candidate_sha256'], 'Candidate hash mismatch')
    candidate = json.loads(candidate_path.read_text())
    audit = geometry_audit(candidate, atoms)
    require(audit.get('heavy_connectivity_preserved') and not audit['unassigned_hydrogens']
            and not audit.get('separated_geometry', False), 'Chemical geometry requires review')
    # Proton ownership changes are retained as observed geometric states;
    # they are not silently removed, called electronic charges, or relabelled.
    energy, gradient, gradient_path, gradient_hash = gradient_in_frame(meta, atoms)
    if meta['result'].get('engrad_sha256'):
        require(gradient_hash == meta['result']['engrad_sha256'], 'Saved gradient hash mismatch')
    paths += [geometry, Path(gradient_path)]
    frequencies = None
    if minimum:
        require('THE OPTIMIZATION HAS CONVERGED' in log, 'Raw optimization did not converge')
        require(sha(geometry) == meta['result']['optimized_xyz_sha256'], 'Optimized geometry hash mismatch')
        block = log.rsplit('VIBRATIONAL FREQUENCIES', 1)[-1].split('NORMAL MODES', 1)[0]
        modes = [(int(i), float(f)) for i, f in re.findall(r'^\s*(\d+):\s+([-+]?\d+\.\d+)\s+cm', block, re.M)]
        require([i for i, _ in modes] == list(range(3 * len(atoms)))
                and all(abs(f) <= .001 for _, f in modes[:6]), 'Incomplete internal modes')
        hessian = run / (Path(meta['input']).stem + '.hess')
        require(sha(hessian) == meta['result']['hessian_sha256'], 'Hessian hash mismatch')
        frequencies = inspect_hessian(hessian, atoms, [f for _, f in modes[6:]])
        require(np.all(frequencies > 0), 'Nonpositive internal frequency')
        paths.append(hessian)
    return dict(meta=meta, atoms=atoms, candidate=candidate, energy_Hartree=energy, gradient=gradient,
                geometry_path=geometry, geometry_sha256=sha(geometry), chemical_state=audit,
                minimum_frequency_cm1=float(min(frequencies)) if frequencies is not None else None,
                internal_mode_count=len(frequencies) if frequencies is not None else None,
                evidence_sha256={str(path.relative_to(ROOT)): sha(path) for path in paths
                                 if path.is_relative_to(ROOT)})


def load_uniform_minima():
    path = ROOT / PANEL
    require(path.exists(), 'Uniform panel has not been prepared; no checks or scores written')
    plan = json.loads(path.read_text()); verify_hashes(plan)
    require(plan['comparison_context'] == CONTEXT, 'Unexpected comparison context')
    lineage_path = 'experiments/completion/drug-start-lineage-r1.json'
    require(plan['prerequisite_sha256'][lineage_path] == sha(ROOT / lineage_path), 'Lineage digest mismatch')
    lineage = json.loads((ROOT / lineage_path).read_text())
    frozen = {row['candidate_id']: row for row in lineage['rows']}
    require(len(frozen) == 17 and len(plan['comparison_mapping']) == 17 and len(plan['jobs']) == 17,
            'Changed panel size')
    require(set(plan['jobs']) == {r['config'] for r in plan['comparison_mapping']}, 'Job mapping differs')
    records, seen = [], set()
    for mapping in plan['comparison_mapping']:
        cid = mapping['candidate_id']
        require(cid in frozen and cid not in seen, 'Duplicate or unknown start')
        seen.add(cid)
        require(mapping['species'] == frozen[cid]['species']
                and mapping['selection_role'] == frozen[cid]['selection_role'], 'Frozen selection policy changed')
        cfg = json.loads((ROOT / mapping['config']).read_text())
        require(cfg['candidate_id'] == cid and cfg['grid_setting'] == 'dense', 'Wrong candidate/grid mapping')
        report = ROOT / f"reports/{cfg['job_id']}.json"
        result = audit_raw_result(report, input_text=expected_input('dense', True), grid='dense', minimum=True)
        require(result['meta']['config_sha256'] == sha(ROOT / mapping['config']), 'Current config mismatch')
        require(result['meta']['candidate_id'] == cid and result['meta']['species'] == mapping['species'],
                'Raw result identity mismatch')
        require(result['candidate'] == frozen[cid]['candidate'], 'Raw candidate differs from frozen lineage')
        records.append(dict(mapping=mapping, result=result))
    return plan, records


def prepare_checks():
    target = ROOT / CHECKS
    if target.exists():
        verify_hashes(json.loads(target.read_text()))
        return target
    plan, minima = load_uniform_minima()
    seed_path = 'experiments/completion/geometry-precision-v3-checks.json'
    seed = json.loads((ROOT / seed_path).read_text()); verify_hashes(seed)
    seed_jobs = {json.loads((ROOT / p).read_text())['candidate_id']: p for p in seed['jobs']}
    jobs, parents = [], []
    hashes = {PANEL: sha(ROOT / PANEL)}
    hashes[seed_path] = sha(ROOT / seed_path)
    for row in minima:
        mapping, result = row['mapping'], row['result']
        cid = mapping['candidate_id']
        geometry = result['geometry_path']
        parent_report = f"reports/{result['meta']['job_id']}.json"
        if mapping['reused_representative']:
            cfg_path = seed_jobs[cid]
            strict_cfg = json.loads((ROOT / cfg_path).read_text())
            source = next(p for p, dest in strict_cfg['assets'].items() if dest == 'initial.xyz')
            require(sha(ROOT / source) == result['geometry_sha256'], 'Representative strict geometry differs')
        else:
            cfg_path = prepare_job(f'w3-uniform-r4-{cid}-strict-engrad', cid, geometry,
                                   grid='strict', optimization=False, role='full_panel_strict_at_dense_minimum',
                                   source_evidence=[ROOT / PANEL, ROOT / parent_report])
        jobs.append(cfg_path)
        parents.append(dict(candidate_id=cid, species=mapping['species'], selection_role=mapping['selection_role'],
                            config=mapping['config'], strict_config=cfg_path, report=parent_report,
                            geometry=str(geometry.relative_to(ROOT)), geometry_sha256=result['geometry_sha256'],
                            reused_representative=mapping['reused_representative']))
        hashes.update(result['evidence_sha256'])
        cfg = json.loads((ROOT / cfg_path).read_text())
        for rel in [cfg_path, cfg['input'], *cfg['assets']]: hashes[rel] = sha(ROOT / rel)
    for rel in ['scripts/completion_drug_panel.py', 'scripts/completion_drug_contrasts.py',
                'scripts/prepare_completion_drug.py']:
        hashes[rel] = sha(ROOT / rel)
    require(len(jobs) == 17 and sum(p['reused_representative'] for p in parents) == 2, 'Changed strict panel')
    frozen = dict(kind='uniform-vacuum-r4-strict-checks', jobs=jobs, parents=parents, input_hashes=hashes,
                  comparison_context=CONTEXT, acceptance=POLICY, automatic_retries=0, automatic_execution=False,
                  parent_plan_sha256=sha(ROOT / PANEL), reused_strict_jobs=2, new_strict_jobs=15,
                  full_panel_numerical_precision_verified=False, new_independent_start_count=0)
    immutable(target, json.dumps(frozen, indent=2) + '\n')
    return target


def score_sensitivity(dense_records, strict_records, policy=POLICY):
    require({r['record_id'] for r in dense_records} == {r['record_id'] for r in strict_records},
            'Dense/strict start sets differ')
    for d, s in zip(sorted(dense_records, key=lambda r: r['record_id']),
                    sorted(strict_records, key=lambda r: r['record_id']), strict=True):
        require(all(d[k] == s[k] for k in ('record_id', 'species', 'selection_role')),
                'Dense/strict correspondence changed')
    dense = analyze_uniform_records(dense_records, CONTEXT)
    strict = analyze_uniform_records(strict_records, STRICT_CONTEXT)
    strict_by_id = {r['record_id']: r for r in strict_records}
    changes = []
    complete = True
    for scope in ('primary_only', 'all_selected_poses'):
        fixed = panel_scores({species: strict_by_id[cid]['energy_Hartree']
                              for species, cid in dense['scopes'][scope]['selected_record_ids'].items()})
        for selection, strict_scores in [('dense_selected_poses_fixed', fixed),
                                          ('reselected_on_strict_energy', strict['scopes'][scope]['scores'])]:
            for d, s in zip(dense['scopes'][scope]['scores'], strict_scores, strict=True):
                for score in ('B', 'X'):
                    a, b = d[score]['kcal_mol'], s[score]['kcal_mol']
                    delta = None if a is None or b is None else b - a
                    complete = complete and delta is not None
                    changes.append(dict(scope=scope, selection=selection, coformer=d['coformer'], score=score,
                                        strict_minus_dense_kcal_mol=delta,
                                        within_target=delta is not None and abs(delta) <= policy['max_score_change_kcal_mol']))
    return dict(dense=dense, strict_at_dense_geometries=strict, score_changes=changes,
                all_scores_available=complete, all_score_changes_within_target=all(r['within_target'] for r in changes),
                selected_pose_changes={scope: [species for species, cid in dense['scopes'][scope]['selected_record_ids'].items()
                    if strict['scopes'][scope]['selected_record_ids'].get(species) != cid]
                    for scope in ('primary_only', 'all_selected_poses')})


def evaluate():
    check_path = ROOT / CHECKS
    require(check_path.exists(), 'Full-panel strict checks are missing; no scientific score report written')
    checks = json.loads(check_path.read_text()); verify_hashes(checks)
    require(checks['parent_plan_sha256'] == sha(ROOT / PANEL) and checks['acceptance'] == POLICY,
            'Panel provenance or numerical policy changed')
    require(checks['comparison_context'] == CONTEXT, 'Strict-check comparison context changed')
    _, minima = load_uniform_minima()
    parents = {p['candidate_id']: p for p in checks['parents']}
    require(len(parents) == len(checks['parents']) == len(minima) == 17
            and len(checks['jobs']) == 17
            and set(checks['jobs']) == {p['strict_config'] for p in parents.values()},
            'Strict parent coverage mismatch')
    dense_records, strict_records, raw_rows = [], [], []
    for row in minima:
        mapping, result = row['mapping'], row['result']
        cid = mapping['candidate_id']; parent = parents[cid]
        require(parent['config'] == mapping['config'] and parent['species'] == mapping['species']
                and parent['selection_role'] == mapping['selection_role'], 'Strict lineage mapping changed')
        cfg = json.loads((ROOT / parent['strict_config']).read_text())
        strict = audit_raw_result(ROOT / f"reports/{cfg['job_id']}.json",
                                  input_text=expected_input('strict', False), grid='strict', minimum=False)
        require(strict['meta']['config_sha256'] == sha(ROOT / parent['strict_config'])
                and strict['meta']['candidate_id'] == cid, 'Strict result mapping mismatch')
        require(result['geometry_sha256'] == parent['geometry_sha256'] == strict['geometry_sha256'],
                'Strict coordinates differ from the dense minimum')
        delta_gradient = float(np.abs(strict['gradient'] - result['gradient']).max())
        # Both loaders align to the same saved XYZ frame.
        for output, value, context in [(dense_records, result, CONTEXT), (strict_records, strict, STRICT_CONTEXT)]:
            output.append(dict(record_id=cid, species=mapping['species'], selection_role=mapping['selection_role'],
                               energy_Hartree=value['energy_Hartree'], comparison_context=context))
        raw_rows.append(dict(candidate_id=cid, species=mapping['species'], selection_role=mapping['selection_role'],
                             minimum_frequency_cm1=result['minimum_frequency_cm1'], internal_mode_count=result['internal_mode_count'],
                             chemical_state=result['chemical_state'], geometry_sha256=result['geometry_sha256'],
                             max_gradient_component_change_Hartree_Bohr=delta_gradient,
                             gradient_within_target=delta_gradient <= POLICY['max_gradient_component_change_Hartree_Bohr'],
                             source_sha256={**result['evidence_sha256'], **strict['evidence_sha256']}))
    analysis = score_sensitivity(dense_records, strict_records)
    passed = analysis['all_scores_available'] and analysis['all_score_changes_within_target'] and all(r['gradient_within_target'] for r in raw_rows)
    record = dict(status='uniform_vacuum_numerical_panel_passed' if passed else 'stopped_requires_review',
                  parent_plan_sha256=sha(ROOT / PANEL), check_plan_sha256=sha(check_path), acceptance=POLICY,
                  rows=raw_rows, analysis=analysis, numerical_precision_verified=passed,
                  scientific_project_complete=False, method_sensitivity_complete=False, solvent_sensitivity_complete=False,
                  independent_validation_complete=False, scope='Electronic B/X for the selected vacuum poses only; no crystal structure or formation-probability prediction')
    target = ROOT / 'reports/completion/uniform-vacuum-r4-scores.json'
    immutable(target, json.dumps(record, indent=2) + '\n')
    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('stage', choices=('prepare-checks', 'evaluate'))
    args = parser.parse_args()
    print(prepare_checks() if args.stage == 'prepare-checks' else evaluate())
