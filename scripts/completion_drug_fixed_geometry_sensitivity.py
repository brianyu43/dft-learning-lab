"""Compare two already verified drug panels at exactly matched coordinates.

This is score arithmetic, not a raw ORCA reader or a scientific gate. It has
no CLI, writes no files, and cannot certify either numerical precision or
method/solvent completion. Production callers must first verify both panels.
"""
from __future__ import annotations

import math
import re

from completion_drug_contrasts import analyze_uniform_records, compare_observed_groups, panel_scores, TIE_WIDTHS


def require(value, message):
    if not value:
        raise ValueError(message)


def validate_contexts(baseline, target, kind):
    require(baseline['method'] == 'r2SCAN-3c' and baseline['environment'] == 'vacuum'
            and baseline['energy_quantity'] == 'electronic_energy', 'Wrong primary comparison context')
    if kind == 'method':
        require(target['method'] == 'wB97X-3c' and target['environment'] == 'vacuum'
                and target['energy_quantity'] == 'electronic_energy',
                'Method comparison must keep the vacuum electronic-energy quantity')
    elif kind == 'solvent':
        require(target['method'] == 'r2SCAN-3c' and target['environment'] == 'SMD(ethanol)'
                and target['energy_quantity'] == 'orca_smd_total_energy_including_cds',
                'Solvent comparison must keep the method and account for CDS exactly once')
    else:
        raise ValueError('Only the frozen method or solvent comparisons are supported')


def compare(baseline_records, target_records, *, baseline_context, target_context, kind):
    """Keep baseline-selected and target-reselected pose scores separate.

    Every record needs record_id, species, selection_role, energy_Hartree,
    comparison_context and geometry_sha256. Hash equality is a correspondence
    check only; this function does not read coordinates or validate the hash.
    Partial but matched panels retain unavailable scores for missing controls.
    """
    require(baseline_records and target_records, 'Empty comparison panel')
    baseline = analyze_uniform_records(baseline_records, baseline_context)
    target = analyze_uniform_records(target_records, target_context)
    validate_contexts(baseline_context, target_context, kind)
    before = {r['record_id']: r for r in baseline_records}
    after = {r['record_id']: r for r in target_records}
    require(set(before) == set(after), 'Panel start identities differ')
    for cid, b in before.items():
        t = after[cid]
        require(all(b[k] == t[k] for k in ('species', 'selection_role')), 'Species or frozen role changed')
        digest = b.get('geometry_sha256')
        require(isinstance(digest, str) and re.fullmatch(r'[0-9a-f]{64}', digest)
                and t.get('geometry_sha256') == digest, 'Coordinates are not matched by SHA-256')

    scopes = {}
    available = True
    for scope, original in baseline['scopes'].items():
        original_ids = original['selected_record_ids']
        fixed_scores = panel_scores({s: after[cid]['energy_Hartree'] for s, cid in original_ids.items()})
        base_by_coformer = {r['coformer']: r for r in original['scores']}
        fixed_by_coformer = {r['coformer']: r for r in fixed_scores}
        choices = [('baseline_selected_poses_fixed', original_ids, fixed_scores),
                   ('reselected_on_target_energy', target['scopes'][scope]['selected_record_ids'],
                    target['scopes'][scope]['scores'])]
        selections = {}
        for selection, selected_ids, scores in choices:
            decomposition = []
            for row in scores:
                coformer = row['coformer']
                for score in ('B', 'X'):
                    terms = (base_by_coformer[coformer][score], fixed_by_coformer[coformer][score], row[score])
                    base, fixed, selected = [r['kcal_mol'] for r in terms]
                    present = all(v is not None for v in (base, fixed, selected))
                    available = available and present
                    effect = math.fsum([fixed, -base]) if present else None
                    reselection = math.fsum([selected, -fixed]) if present else None
                    total = math.fsum([selected, -base]) if present else None
                    if present:
                        require(abs(math.fsum([effect, reselection, -total])) <= 1e-8,
                                'Score change decomposition failed')
                    decomposition.append(dict(coformer=coformer, score=score,
                        baseline_kcal_mol=base, target_at_baseline_selected_poses_kcal_mol=fixed,
                        target_at_selected_poses_kcal_mol=selected,
                        fixed_geometry_context_effect_kcal_mol=effect,
                        target_pose_reselection_effect_kcal_mol=reselection,
                        total_change_from_baseline_kcal_mol=total, available=present,
                        missing_species=sorted({s for r in terms for s in r['missing_species']})))
            comparisons = {score: [compare_observed_groups({r['coformer']: r[score]['kcal_mol'] for r in scores}, width)
                                   for width in TIE_WIDTHS] for score in ('B', 'X')}
            transitions = {}
            for score in ('B', 'X'):
                transitions[score] = []
                for b, t in zip(original['comparisons'][score], comparisons[score], strict=True):
                    rows = []
                    for br, tr in zip(b['rows'], t['rows'], strict=True):
                        rows.append(dict(formed_pair_coformer=br['formed_pair_coformer'],
                            not_formed_pair_coformer=br['not_formed_pair_coformer'],
                            baseline_assessment=br['assessment'], target_assessment=tr['assessment'],
                            assessment_changed=br['assessment'] != tr['assessment']))
                    transitions[score].append(dict(tie_width_kcal_mol=b['tie_width_kcal_mol'], rows=rows))
            selections[selection] = dict(selected_record_ids=selected_ids,
                pose_changes_from_baseline=[s for s in selected_ids if selected_ids[s] != original_ids[s]],
                scores=scores, decomposition=decomposition, comparisons=comparisons,
                assessment_transitions_from_baseline=transitions)
        scopes[scope] = selections
    return dict(status='arithmetic_only_requires_external_provenance_verification', comparison_kind=kind,
        baseline_analysis=baseline, target_analysis=target, scopes=scopes, all_scores_available=available,
        matched_geometry_sha256={cid: r['geometry_sha256'] for cid, r in sorted(before.items())},
        raw_provenance_verified=False, numerical_precision_verified=False,
        target_method_or_solvent_minima_verified=False, scientific_project_complete=False,
        independent_observation_count=None, new_independent_start_count=0,
        scope='Matched-coordinate B/X differences and pose reselection within frozen starts; no formation free energy')
