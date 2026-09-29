"""Analytic artificial panels; these are not method/solvent DFT results."""
from copy import deepcopy
import hashlib
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import completion_drug_fixed_geometry_sensitivity as sensitivity
from completion_drug_contrasts import COFORMERS, HARTREE_KCAL

BASE = dict(method='r2SCAN-3c', environment='vacuum', numerical_protocol_id='test-baseline',
            energy_quantity='electronic_energy')
METHOD = dict(method='wB97X-3c', environment='vacuum', numerical_protocol_id='test-method',
              energy_quantity='electronic_energy')
SOLVENT = dict(method='r2SCAN-3c', environment='SMD(ethanol)', numerical_protocol_id='test-solvent',
               energy_quantity='orca_smd_total_energy_including_cds')


class FixedGeometrySensitivityTests(unittest.TestCase):
    def records(self, context=METHOD):
        energies = {'diflunisal': -100., 'diflunisal--diflunisal': -200.01}
        for i, b in enumerate(COFORMERS):
            e = -20. - i
            energies[b] = e
            energies[f'{b}--{b}'] = 2 * e - .01
            energies[f'diflunisal--{b}'] = -100. + e - .02
        before = [dict(record_id=s + '-primary', species=s, selection_role='primary', energy_Hartree=e,
            comparison_context=dict(BASE), geometry_sha256=hashlib.sha256(s.encode()).hexdigest()) for s, e in energies.items()]
        after = [{**r, 'comparison_context': dict(context)} for r in before]
        return before, after

    def run_comparison(self, b, t, kind='method', context=METHOD):
        return sensitivity.compare(b, t, baseline_context=BASE, target_context=context, kind=kind)

    def test_binding_and_self_competition_changes_have_analytic_values(self):
        b, t = self.records()
        for r in t:
            if r['species'] == 'diflunisal--adenine': r['energy_Hartree'] -= 2. / HARTREE_KCAL
            if r['species'] == 'adenine--adenine': r['energy_Hartree'] -= 6. / HARTREE_KCAL
        result = self.run_comparison(b, t)
        self.assertTrue(result['all_scores_available'])
        rows = result['scopes']['primary_only']['baseline_selected_poses_fixed']['decomposition']
        values = {r['score']: r for r in rows if r['coformer'] == 'adenine'}
        self.assertAlmostEqual(values['B']['total_change_from_baseline_kcal_mol'], -2., places=8)
        self.assertAlmostEqual(values['X']['total_change_from_baseline_kcal_mol'], 1., places=8)
        self.assertFalse(result['raw_provenance_verified'])
        self.assertFalse(result['numerical_precision_verified'])

    def test_method_reference_energy_offsets_cancel_in_complete_contrasts(self):
        b, t = self.records()
        offsets = {'diflunisal': 3., **{s: 10. + i for i, s in enumerate(COFORMERS)}}
        for r in t: r['energy_Hartree'] += sum(offsets[s] for s in r['species'].split('--'))
        result = self.run_comparison(b, t)
        for scope in result['scopes'].values():
            for selection in scope.values():
                for row in selection['decomposition']:
                    self.assertAlmostEqual(row['total_change_from_baseline_kcal_mol'], 0., places=8)

    def test_pose_reselection_is_separate_from_fixed_geometry_effect_and_primary_scope(self):
        b, t = self.records()
        original = next(r for r in b if r['species'] == 'diflunisal--adenine')
        alt = {**original, 'record_id': 'adenine-alt', 'selection_role': 'alternative',
               'energy_Hartree': original['energy_Hartree'] + 1. / HARTREE_KCAL,
               'geometry_sha256': 'a' * 64}
        b.append(alt)
        t.append({**alt, 'comparison_context': dict(METHOD),
                  'energy_Hartree': original['energy_Hartree'] - 1. / HARTREE_KCAL})
        next(r for r in t if r['record_id'] == original['record_id'])['energy_Hartree'] -= .2 / HARTREE_KCAL
        result = self.run_comparison(b, t)
        for scope, expected in [('primary_only', (0., -.2)), ('all_selected_poses', (-.8, -1.))]:
            selected = result['scopes'][scope]['reselected_on_target_energy']
            for row in selected['decomposition']:
                if row['coformer'] == 'adenine':
                    self.assertAlmostEqual(row['fixed_geometry_context_effect_kcal_mol'], -.2, places=8)
                    self.assertAlmostEqual(row['target_pose_reselection_effect_kcal_mol'], expected[0], places=8)
                    self.assertAlmostEqual(row['total_change_from_baseline_kcal_mol'], expected[1], places=8)
            self.assertEqual(selected['pose_changes_from_baseline'], [] if scope == 'primary_only' else ['diflunisal--adenine'])

    def test_solvent_cds_quantity_is_used_once_without_adding_a_second_correction(self):
        b, t = self.records(SOLVENT)
        next(r for r in t if r['species'] == 'diflunisal--caffeine')['energy_Hartree'] += .3 / HARTREE_KCAL
        result = self.run_comparison(b, t, 'solvent', SOLVENT)
        rows = result['scopes']['primary_only']['baseline_selected_poses_fixed']['decomposition']
        for r in rows:
            if r['coformer'] == 'caffeine':
                self.assertAlmostEqual(r['fixed_geometry_context_effect_kcal_mol'], .3, places=8)
        self.assertFalse(result['target_method_or_solvent_minima_verified'])

    def test_all_tie_widths_and_rank_reversal_are_retained(self):
        b, t = self.records()
        next(r for r in b if r['species'] == 'diflunisal--adenine')['energy_Hartree'] -= .5 / HARTREE_KCAL
        next(r for r in t if r['species'] == 'diflunisal--adenine')['energy_Hartree'] += .5 / HARTREE_KCAL
        result = self.run_comparison(b, t)
        rows = result['scopes']['primary_only']['baseline_selected_poses_fixed']['assessment_transitions_from_baseline']['B']
        self.assertEqual([r['tie_width_kcal_mol'] for r in rows], [0., .5, 1., 2.])
        self.assertEqual((rows[0]['rows'][0]['baseline_assessment'], rows[0]['rows'][0]['target_assessment']), ('consistent', 'reversed'))
        self.assertEqual(rows[1]['rows'][0]['target_assessment'], 'tie')
        self.assertIsNone(result['independent_observation_count'])

    def test_missing_control_remains_unavailable_and_unmatched_panels_are_rejected(self):
        b, t = self.records()
        b = [r for r in b if r['species'] != 'adenine--adenine']
        with self.assertRaisesRegex(ValueError, 'identities'): self.run_comparison(b, t)
        t = [r for r in t if r['species'] != 'adenine--adenine']
        result = self.run_comparison(b, t)
        self.assertFalse(result['all_scores_available'])
        row = next(r for r in result['scopes']['primary_only']['baseline_selected_poses_fixed']['decomposition']
                   if r['coformer'] == 'adenine' and r['score'] == 'X')
        self.assertIsNone(row['total_change_from_baseline_kcal_mol'])
        self.assertEqual(row['missing_species'], ['adenine--adenine'])

    def test_changed_geometry_role_species_duplicate_or_nonfinite_energy_is_rejected(self):
        for mutation in ('geometry', 'missing_hash', 'role', 'species', 'duplicate', 'nan', 'boolean'):
            b, t = deepcopy(self.records())
            if mutation == 'geometry': t[0]['geometry_sha256'] = 'b' * 64
            if mutation == 'missing_hash': del t[0]['geometry_sha256']
            if mutation == 'role': t[0]['selection_role'] = 'alternative'
            if mutation == 'species': t[0]['species'] = 'caffeine'
            if mutation == 'duplicate': t.append(dict(t[0]))
            if mutation == 'nan': t[0]['energy_Hartree'] = float('nan')
            if mutation == 'boolean': t[0]['energy_Hartree'] = True
            with self.subTest(mutation=mutation), self.assertRaises(ValueError): self.run_comparison(b, t)

    def test_mixed_method_solvent_and_wrong_energy_quantity_are_rejected(self):
        for context, kind in [({**METHOD, 'environment': 'SMD(ethanol)'}, 'method'),
                              ({**SOLVENT, 'energy_quantity': 'electronic_energy'}, 'solvent'),
                              ({**SOLVENT, 'method': 'wB97X-3c'}, 'solvent')]:
            b, t = self.records(context)
            with self.subTest(context=context), self.assertRaises(ValueError):
                self.run_comparison(b, t, kind, context)


if __name__ == '__main__':
    unittest.main()
