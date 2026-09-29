"""Synthetic supplementary-score arithmetic and fail-closed production gates."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import completion_drug_zpe_panel as supplement
from completion_drug_contrasts import COFORMERS, HARTREE_KCAL


class ZPESupplementTests(unittest.TestCase):
    def records(self):
        energies = {'diflunisal': (-100., .1), 'diflunisal--diflunisal': (-200.01, .202)}
        for i, b in enumerate(COFORMERS):
            e, z = -20. - i, .02 + .001 * i
            energies[b] = (e, z)
            energies[f'{b}--{b}'] = (2 * e - .01, 2 * z + .004)
            energies[f'diflunisal--{b}'] = (-100. + e - .02, .1 + z + .003)
        return [dict(record_id=s + '-primary', species=s, selection_role='primary',
                     electronic_energy_Hartree=e, zpe_Hartree=z, comparison_context=dict(supplement.panel.CONTEXT))
                for s, (e, z) in energies.items()]

    def test_analytic_B_and_X_corrections_are_separate_from_electronic_scores(self):
        result = supplement.arithmetic(self.records())
        self.assertTrue(result['all_scores_available'])
        self.assertFalse(result['electronic_primary_selection_replaced'])
        self.assertFalse(result['thermal_free_energy_computed'])
        for scope in result['scopes'].values():
            for selection in scope.values():
                for row in selection['decomposition']:
                    correction = .003 * HARTREE_KCAL if row['score'] == 'B' else 0.
                    self.assertAlmostEqual(row['harmonic_zpe_correction_kcal_mol'], correction, places=9)
                    self.assertAlmostEqual(row['total_change_from_electronic_baseline_kcal_mol'], correction, places=8)
                    self.assertAlmostEqual(row['electronic_pose_selection_change_kcal_mol'], 0., places=9)
                self.assertEqual(selection['pose_changes_from_electronic_selection'], [])

    def test_pose_switch_separates_electronic_penalty_and_zpe_change(self):
        rows = self.records(); original = next(r for r in rows if r['species'] == 'diflunisal--adenine')
        rows.append({**original, 'record_id': 'adenine-alternative', 'selection_role': 'alternative',
                     'electronic_energy_Hartree': original['electronic_energy_Hartree'] + .3 / HARTREE_KCAL,
                     'zpe_Hartree': original['zpe_Hartree'] - .6 / HARTREE_KCAL})
        result = supplement.arithmetic(rows)
        self.assertEqual(result['electronic_baseline']['scopes']['all_selected_poses']['selected_record_ids']['diflunisal--adenine'],
                         original['record_id'])
        all_poses = result['scopes']['all_selected_poses']
        reselected = all_poses['reselected_on_electronic_plus_zpe']
        self.assertEqual(reselected['pose_changes_from_electronic_selection'], ['diflunisal--adenine'])
        fixed = {r['score']: r for r in all_poses['electronic_selected_poses_fixed']['decomposition'] if r['coformer'] == 'adenine'}
        for row in reselected['decomposition']:
            if row['coformer'] == 'adenine':
                self.assertAlmostEqual(row['electronic_pose_selection_change_kcal_mol'], .3, places=8)
                self.assertAlmostEqual(row['harmonic_zpe_correction_kcal_mol'] - fixed[row['score']]['harmonic_zpe_correction_kcal_mol'], -.6, places=8)
                self.assertAlmostEqual(row['total_change_from_electronic_baseline_kcal_mol'] - fixed[row['score']]['total_change_from_electronic_baseline_kcal_mol'], -.3, places=8)
        self.assertEqual(result['scopes']['primary_only']['reselected_on_electronic_plus_zpe']['pose_changes_from_electronic_selection'], [])

    def test_rank_reversal_and_all_tie_widths_remain_visible(self):
        rows = self.records(); adenine = next(r for r in rows if r['species'] == 'diflunisal--adenine')
        adenine['electronic_energy_Hartree'] -= .5 / HARTREE_KCAL
        adenine['zpe_Hartree'] += 1. / HARTREE_KCAL
        result = supplement.arithmetic(rows)
        base = result['electronic_baseline']['scopes']['primary_only']['comparisons']['B'][0]['rows'][0]
        corrected = result['scopes']['primary_only']['electronic_selected_poses_fixed']['comparisons']['B']
        self.assertEqual(base['assessment'], 'consistent')
        self.assertEqual([r['tie_width_kcal_mol'] for r in corrected], [0., .5, 1., 2.])
        self.assertEqual(corrected[0]['rows'][0]['assessment'], 'reversed')
        self.assertEqual(corrected[1]['rows'][0]['assessment'], 'tie')
        self.assertIsNone(result['independent_observation_count'])

    def test_missing_self_control_does_not_become_zero_or_a_prediction(self):
        rows = [r for r in self.records() if r['species'] != 'adenine--adenine']
        result = supplement.arithmetic(rows); self.assertFalse(result['all_scores_available'])
        for scope in result['scopes'].values():
            for selection in scope.values():
                missing = next(r for r in selection['decomposition'] if r['coformer'] == 'adenine' and r['score'] == 'X')
                self.assertIsNone(missing['electronic_plus_zpe_kcal_mol'])
                self.assertIsNone(missing['total_change_from_electronic_baseline_kcal_mol'])
                self.assertEqual(missing['missing_species'], ['adenine--adenine'])
                self.assertFalse(missing['available'])

    def test_context_nonfinite_energy_duplicate_identity_or_role_change_is_rejected(self):
        for mutation in ('context', 'nan', 'zero', 'boolean', 'duplicate', 'role'):
            rows = deepcopy(self.records())
            if mutation == 'context': rows[0]['comparison_context']['method'] = 'wB97X-3c'
            if mutation == 'nan': rows[0]['zpe_Hartree'] = float('nan')
            if mutation == 'zero': rows[0]['zpe_Hartree'] = 0.
            if mutation == 'boolean': rows[0]['electronic_energy_Hartree'] = True
            if mutation == 'duplicate': rows.append(dict(rows[0]))
            if mutation == 'role': rows[0]['selection_role'] = 'invented'
            with self.subTest(mutation=mutation), self.assertRaises(ValueError): supplement.arithmetic(rows)

    def test_isotope_or_mass_changes_are_not_combined_silently(self):
        supplement.mass_policy([1, 6, 7, 8, 9], [1.008, 12.011, 14.007, 15.999, 18.998])
        for numbers, masses in [([1], [2.014]), ([8], [float('nan')]), ([1, 8], [1.008]), ([16], [32.06])]:
            with self.subTest(numbers=numbers, masses=masses), self.assertRaises(ValueError):
                supplement.mass_policy(numbers, masses)

    def test_missing_full_panel_cannot_write_zpe_results(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with patch.object(supplement, 'ROOT', root), patch.object(supplement.panel, 'evaluate') as evaluate:
                with self.assertRaisesRegex(ValueError, 'no ZPE score report written'): supplement.evaluate()
                evaluate.assert_not_called()
            self.assertEqual(list(root.iterdir()), [])

    def test_parent_pass_flag_is_not_a_substitute_for_raw_revalidation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); parent = root / supplement.PARENT
            parent.parent.mkdir(parents=True); parent.write_text(json.dumps({'status': 'uniform_vacuum_numerical_panel_passed'}))
            with patch.object(supplement, 'ROOT', root), patch.object(supplement.panel, 'evaluate',
                    side_effect=ValueError('raw fixture rejected')) as evaluate:
                with self.assertRaisesRegex(ValueError, 'raw fixture rejected'): supplement.evaluate()
                evaluate.assert_called_once_with()
            self.assertFalse((root / supplement.OUTPUT).exists())


if __name__ == '__main__':
    unittest.main()
