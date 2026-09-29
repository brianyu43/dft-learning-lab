"""Raw-file read-only checks and numerical comparisons, not new experiments."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import completion_drug_panel as panel
from completion_drug_contrasts import COFORMERS, HARTREE_KCAL
from research_run import sha


class RawPanelTests(unittest.TestCase):
    report = ROOT / 'reports/w3-diflunisal--diflunisal-p05-refine-r1-opt-freq.json'

    def read_old(self, report=None, input_text=None):
        meta = json.loads(self.report.read_text())
        inp = ROOT / meta['input']
        return panel.audit_raw_result(report or self.report,
                 input_text=input_text if input_text is not None else inp.read_text(), grid='g3', minimum=True)

    def test_existing_minimum_is_revalidated_without_rewriting_old_results(self):
        meta = json.loads(self.report.read_text()); run = Path(meta['run_dir'])
        paths = [self.report, run / 'manifest.json', run / 'result.json', run / 'output.log']
        before = {p: sha(p) for p in paths}
        result = self.read_old()
        self.assertAlmostEqual(result['energy_Hartree'], -1850.964363429242, places=9)
        self.assertEqual(result['internal_mode_count'], 150)
        self.assertAlmostEqual(result['minimum_frequency_cm1'], 9.18442963, places=7)
        self.assertTrue(result['chemical_state']['chemical_state_preserved_by_geometry'])
        self.assertEqual(result['gradient'].shape, (52, 3))
        self.assertEqual(before, {p: sha(p) for p in paths})

    def test_historical_minimum_is_not_accepted_as_a_dense_panel_result(self):
        with self.assertRaisesRegex(ValueError, 'frozen numerical protocol'):
            self.read_old(input_text=panel.expected_input('dense', True))

    def test_report_pass_flag_cannot_bypass_output_hash_or_raw_protocol(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'report.json'
            meta = json.loads(self.report.read_text()); meta['output_sha256'] = '0' * 64
            path.write_text(json.dumps(meta))
            with self.assertRaisesRegex(ValueError, 'Output hash'):
                self.read_old(path)

    def test_missing_panel_or_checks_cannot_write_inputs_or_scores(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with patch.object(panel, 'ROOT', root):
                with self.assertRaisesRegex(ValueError, 'has not been prepared'):
                    panel.prepare_checks()
                with self.assertRaisesRegex(ValueError, 'checks are missing'):
                    panel.evaluate()
            self.assertEqual(list(root.iterdir()), [])


class ScoreSensitivityTests(unittest.TestCase):
    def records(self):
        energies = {'diflunisal': -100., 'diflunisal--diflunisal': -200.01}
        for i, b in enumerate(COFORMERS):
            energies[b] = -20. - i
            energies[f'{b}--{b}'] = 2 * energies[b] - .01
            energies[f'diflunisal--{b}'] = energies['diflunisal'] + energies[b] - .02
        return [dict(record_id=s + '-primary', species=s, selection_role='primary',
                     energy_Hartree=e, comparison_context=panel.CONTEXT) for s, e in energies.items()]

    def strict(self, dense):
        return [{**r, 'comparison_context': panel.STRICT_CONTEXT} for r in deepcopy(dense)]

    def test_known_small_contrast_change_and_changed_grid_context(self):
        dense = self.records(); strict = self.strict(dense)
        for row in strict:
            if row['species'] == 'diflunisal--adenine': row['energy_Hartree'] += .08 / HARTREE_KCAL
        result = panel.score_sensitivity(dense, strict)
        self.assertTrue(result['all_scores_available'])
        self.assertTrue(result['all_score_changes_within_target'])
        for row in result['score_changes']:
            self.assertAlmostEqual(row['strict_minus_dense_kcal_mol'], .08 if row['coformer'] == 'adenine' else 0., places=8)
        self.assertNotEqual(result['dense']['comparison_context'], result['strict_at_dense_geometries']['comparison_context'])

    def test_failed_contrast_is_not_hidden_by_an_alternative_pose_switch(self):
        dense = self.records()
        primary = next(r for r in dense if r['species'] == 'diflunisal--adenine')
        dense.append({**primary, 'record_id': 'adenine-alternative', 'selection_role': 'alternative',
                      'energy_Hartree': primary['energy_Hartree'] + .03 / HARTREE_KCAL})
        strict = self.strict(dense)
        next(r for r in strict if r['record_id'] == primary['record_id'])['energy_Hartree'] += .2 / HARTREE_KCAL
        result = panel.score_sensitivity(dense, strict)
        self.assertFalse(result['all_score_changes_within_target'])
        self.assertEqual(result['selected_pose_changes']['all_selected_poses'], ['diflunisal--adenine'])
        rows = [r for r in result['score_changes'] if r['scope'] == 'all_selected_poses'
                and r['coformer'] == 'adenine' and r['score'] == 'X']
        by_selection = {r['selection']: r for r in rows}
        self.assertAlmostEqual(by_selection['dense_selected_poses_fixed']['strict_minus_dense_kcal_mol'], .2, places=8)
        self.assertAlmostEqual(by_selection['reselected_on_strict_energy']['strict_minus_dense_kcal_mol'], .03, places=8)

    def test_missing_control_and_mixed_protocol_do_not_pass(self):
        dense = [r for r in self.records() if r['species'] != 'adenine--adenine']
        result = panel.score_sensitivity(dense, self.strict(dense))
        self.assertFalse(result['all_scores_available'])
        self.assertFalse(result['all_score_changes_within_target'])
        dense = self.records(); strict = self.strict(dense)
        strict[0]['comparison_context'] = panel.CONTEXT
        with self.assertRaisesRegex(ValueError, 'Cannot combine'):
            panel.score_sensitivity(dense, strict)

    def test_different_start_sets_and_roles_are_rejected(self):
        dense = self.records(); strict = self.strict(dense)
        with self.assertRaisesRegex(ValueError, 'sets differ'):
            panel.score_sensitivity(dense, strict[:-1])
        strict[0]['selection_role'] = 'alternative'
        with self.assertRaisesRegex(ValueError, 'correspondence'):
            panel.score_sensitivity(dense, strict)


if __name__ == '__main__':
    unittest.main()
