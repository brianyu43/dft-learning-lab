"""Separate electronic and harmonic-ZPE B/X for the verified uniform panel.

The primary electronic selection is retained. Reselection by E+ZPE is an
explicit supplementary sensitivity, never a replacement for that analysis.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import completion_drug_panel as panel
from completion_drug_contrasts import analyze_uniform_records, compare_observed_groups, finite_number, panel_scores, TIE_WIDTHS
from completion_drug_zpe import extract, WAVENUMBER_CM1_PER_HARTREE
from prepare_week3_drug import immutable
from research_run import ROOT, sha

PARENT = 'reports/completion/uniform-vacuum-r4-scores.json'
OUTPUT = 'reports/completion/uniform-vacuum-r4-zpe.json'
CONTEXT = {**panel.CONTEXT, 'energy_quantity': 'electronic_plus_unscaled_harmonic_zpe'}
MASS_POLICY_AMU = {1: 1.008, 6: 12.011, 7: 14.007, 8: 15.999, 9: 18.998}
CODE_FILES = ['scripts/completion_drug_zpe_panel.py', 'scripts/completion_drug_zpe.py',
              'scripts/completion_drug_panel.py', 'scripts/completion_drug_contrasts.py',
              'scripts/run_week2_frequency_resume.py', 'scripts/run_week2_xtb_screen.py',
              'scripts/week3_recovery.py', 'scripts/prepare_completion_drug.py',
              'scripts/prepare_week3_drug.py', 'scripts/research_run.py', 'uv.lock']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def mass_policy(atomic_numbers, masses):
    require(len(atomic_numbers) == len(masses) and len(masses) > 0, 'Incomplete atomic mass mapping')
    for number, value in zip(atomic_numbers, masses, strict=True):
        require(number in MASS_POLICY_AMU and abs(finite_number(value) - MASS_POLICY_AMU[number]) <= 1e-7,
                'Isotope or atomic mass policy changed')


def arithmetic(records):
    """Arithmetic on already raw-verified matching electronic/ZPE records.

    Missing controls stay unavailable. This function cannot establish source
    provenance, positive modes, grid convergence or scientific completion.
    """
    electronic = []
    for row in records:
        require(row['comparison_context'] == panel.CONTEXT, 'Wrong electronic comparison context')
        require(finite_number(row['zpe_Hartree']) > 0, 'Nonpositive harmonic ZPE')
        electronic.append(dict(record_id=row['record_id'], species=row['species'], selection_role=row['selection_role'],
            energy_Hartree=finite_number(row['electronic_energy_Hartree']), comparison_context=panel.CONTEXT))
    baseline = analyze_uniform_records(electronic, panel.CONTEXT)
    by_id = {r['record_id']: r for r in records}
    scopes = {}; all_available = True
    for scope, original in baseline['scopes'].items():
        pool = [r for r in records if scope != 'primary_only' or r['selection_role'] == 'primary']
        species = sorted({r['species'] for r in pool})
        corrected_selection = {s: min((r for r in pool if r['species'] == s),
            key=lambda r: (math.fsum([r['electronic_energy_Hartree'], r['zpe_Hartree']]), r['record_id']))['record_id']
            for s in species}
        results = {}
        baseline_scores = {r['coformer']: r for r in original['scores']}
        for selection, identities in [('electronic_selected_poses_fixed', original['selected_record_ids']),
                                      ('reselected_on_electronic_plus_zpe', corrected_selection)]:
            chosen = {s: by_id[cid] for s, cid in identities.items()}
            e = panel_scores({s: r['electronic_energy_Hartree'] for s, r in chosen.items()})
            z = panel_scores({s: r['zpe_Hartree'] for s, r in chosen.items()})
            total = panel_scores({s: math.fsum([r['electronic_energy_Hartree'], r['zpe_Hartree']]) for s, r in chosen.items()})
            differences = []
            for erow, zrow, trow in zip(e, z, total, strict=True):
                coformer = erow['coformer']
                require(zrow['coformer'] == trow['coformer'] == coformer, 'Score order changed')
                for score in ('B', 'X'):
                    ev, zv, tv = (r[score]['kcal_mol'] for r in (erow, zrow, trow))
                    base = baseline_scores[coformer][score]['kcal_mol']
                    available = all(v is not None for v in (ev, zv, tv, base))
                    all_available = all_available and available
                    if available:
                        require(abs(math.fsum([ev, zv, -tv])) < 1e-8, 'Electronic/ZPE decomposition failed')
                    differences.append(dict(coformer=coformer, score=score,
                        electronic_baseline_kcal_mol=base, electronic_at_selected_poses_kcal_mol=ev,
                        harmonic_zpe_correction_kcal_mol=zv, electronic_plus_zpe_kcal_mol=tv,
                        electronic_pose_selection_change_kcal_mol=math.fsum([ev, -base]) if available else None,
                        total_change_from_electronic_baseline_kcal_mol=math.fsum([tv, -base]) if available else None,
                        available=available, missing_species=trow[score]['missing_species']))
            results[selection] = dict(selected_record_ids=identities,
                pose_changes_from_electronic_selection=[s for s in identities if identities[s] != original['selected_record_ids'][s]],
                electronic_scores_at_selected_poses=e, zpe_correction_scores=z, electronic_plus_zpe_scores=total,
                decomposition=differences,
                comparisons={score: [compare_observed_groups({r['coformer']: r[score]['kcal_mol'] for r in total}, width)
                                     for width in TIE_WIDTHS] for score in ('B', 'X')})
        scopes[scope] = results
    return dict(status='arithmetic_only_requires_external_provenance_verification',
                electronic_baseline=baseline, supplementary_comparison_context=CONTEXT, scopes=scopes,
                all_scores_available=all_available, electronic_primary_selection_replaced=False,
                thermal_free_energy_computed=False, independent_observation_count=None)


def evaluate():
    require((ROOT / PARENT).is_file(), 'Uniform electronic panel is incomplete; no ZPE score report written')
    parent_path = panel.evaluate()  # Revalidates every dense and strict source.
    require(parent_path.resolve() == (ROOT / PARENT).resolve(), 'Unexpected parent score report')
    parent = json.loads(parent_path.read_text())
    require(parent['status'] == 'uniform_vacuum_numerical_panel_passed'
            and parent['numerical_precision_verified'], 'Uniform electronic precision did not pass')
    plan, minima = panel.load_uniform_minima()
    require(len(minima) == len(parent['rows']) == 17, 'Incomplete common panel')
    evidence = {p: sha(ROOT / p) for p in [*CODE_FILES, PARENT, panel.PANEL, panel.CHECKS]}
    evidence.update(plan['input_hashes'])
    strict_plan = json.loads((ROOT / panel.CHECKS).read_text())
    evidence.update(strict_plan['input_hashes'])
    for row in parent['rows']:
        evidence.update(row['source_sha256'])
    records, raw_rows = [], []
    for row in minima:
        mapping, result = row['mapping'], row['result']
        meta = result['meta']; run = Path(meta['run_dir'])
        hessian = run / (Path(meta['input']).stem + '.hess')
        require(not meta['linear'] and meta['method'] == 'r2SCAN-3c' and meta['solvent'] == 'vacuum',
                'Unsupported harmonic panel context')
        value = extract(run / 'output.log', hessian, result['atoms'], result['energy_Hartree'])
        mass_policy(value['hessian_atomic_numbers'], value['hessian_masses_amu'])
        value['source_sha256'] = {str(Path(p).relative_to(ROOT)): h for p, h in value['source_sha256'].items()}
        evidence.update(value['source_sha256'])
        records.append(dict(record_id=mapping['candidate_id'], species=mapping['species'],
            selection_role=mapping['selection_role'], electronic_energy_Hartree=result['energy_Hartree'],
            zpe_Hartree=value['zpe_Hartree'], comparison_context=panel.CONTEXT))
        raw_rows.append(dict(candidate_id=mapping['candidate_id'], species=mapping['species'],
            selection_role=mapping['selection_role'], geometry_sha256=result['geometry_sha256'],
            minimum_report=f"reports/{meta['job_id']}.json", **value))
    analysis = arithmetic(records)
    require(analysis['electronic_baseline'] == parent['analysis']['dense'], 'Electronic baseline changed in the supplement')
    require(analysis['all_scores_available'], 'Incomplete ZPE controls; no complete score report written')
    for rel, expected in evidence.items():
        require(sha(ROOT / rel) == expected, f'Source changed during ZPE assembly: {rel}')
    record = dict(status='uniform_harmonic_zpe_supplement_extracted', rows=raw_rows, analysis=analysis,
        source_sha256=evidence, comparison_context=CONTEXT, electronic_comparison_context=panel.CONTEXT,
        mass_policy_amu=MASS_POLICY_AMU, wavenumber_cm1_per_Hartree=WAVENUMBER_CM1_PER_HARTREE,
        frequency_scale_factor=1., electronic_numerical_precision_verified=True,
        frequency_grid_convergence_verified=False, crystal_formation_free_energy_computed=False,
        independent_validation_complete=False, scientific_project_complete=False,
        scope='Harmonic supplement at the common r2SCAN-3c minima; electronic analysis remains primary, no Hessian grid-convergence claim')
    target = ROOT / OUTPUT
    immutable(target, json.dumps(record, indent=2) + '\n')
    return target


if __name__ == '__main__':
    print(evaluate())
