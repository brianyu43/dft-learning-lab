"""Harmonic-ZPE extraction from real archives and explicit malformed copies."""
import json
from pathlib import Path
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import completion_drug_zpe as zpe
from research_run import sha


class HarmonicZPETests(unittest.TestCase):
    def source(self):
        meta = json.loads((ROOT / 'reports/w3-diflunisal--diflunisal-p05-refine-r1-opt-freq.json').read_text())
        run = Path(meta['run_dir']); stem = Path(meta['input']).stem
        return run / 'output.log', run / (stem + '.hess'), read(run / (stem + '.xyz')), meta['result']['energy']

    def test_known_frequency_sum_and_preapplied_scale_not_applied_twice(self):
        frequencies = [1000., 2000., 3000.]
        result = zpe.harmonic_zpe(frequencies, applied_scale_factor=.95)
        self.assertAlmostEqual(result['zpe_Hartree'], 3000 / zpe.WAVENUMBER_CM1_PER_HARTREE, places=14)
        self.assertFalse(result['scale_applied_again'])
        self.assertEqual(result['applied_frequency_scale_factor'], .95)

    def test_zero_imaginary_nonfinite_or_empty_modes_are_not_silently_dropped(self):
        for frequencies in ([], [0., 100.], [-.01, 100.], [np.nan], [np.inf], [[100.]]):
            with self.subTest(frequencies=frequencies), self.assertRaisesRegex(ValueError, 'strictly positive'):
                zpe.harmonic_zpe(frequencies, applied_scale_factor=1.)
        for factor in (0., -1., np.nan, True):
            with self.subTest(factor=factor), self.assertRaisesRegex(ValueError, 'scaling factor'):
                zpe.harmonic_zpe([100.], applied_scale_factor=factor)

    def test_actual_diflunisal_dimer_zpe_matches_raw_summary_without_source_writes(self):
        log, hess, atoms, energy = self.source(); before = (sha(log), sha(hess))
        result = zpe.extract(log, hess, atoms, energy)
        self.assertEqual(result['internal_frequency_count'], 150)
        self.assertAlmostEqual(result['zpe_Hartree'], .36845932400788, places=13)
        self.assertEqual(result['modes_below_100_cm1'], 12)
        self.assertLess(abs(result['sum_minus_printed_zpe_Hartree']), 5e-9)
        self.assertFalse(result['thermal_or_entropy_terms_added'])
        self.assertFalse(result['frequency_grid_convergence_verified'])
        self.assertEqual(before, (sha(log), sha(hess)))

    def test_all_nine_archived_minima_remain_separate_and_source_preserving(self):
        path = ROOT / 'reports/completion/drug-archived-zpe-r1.json'
        before = sha(path); result_path = zpe.archived(); result = json.loads(result_path.read_text())
        self.assertEqual(result_path, path); self.assertEqual(sha(path), before)
        self.assertEqual(result['archived_minimum_count'], 9)
        self.assertEqual(len({r['candidate_id'] for r in result['rows']}), 9)
        self.assertFalse(result['rankings_computed']); self.assertFalse(result['B_X_scores_computed'])
        self.assertTrue(all(r['common_grid_comparison_eligible'] is False for r in result['rows']))
        self.assertLess(result['maximum_absolute_printed_zpe_difference_Hartree'], 5e-9)
        self.assertTrue(all(sha(ROOT / p) == digest for p, digest in result['source_sha256'].items()))

    def test_hessian_must_describe_the_same_geometry(self):
        log, hess, atoms, energy = self.source(); atoms.positions[0, 0] += .02
        with self.assertRaises(AssertionError): zpe.extract(log, hess, atoms, energy)

    def test_duplicate_hessian_mode_or_changed_scale_is_rejected(self):
        log, hess, atoms, energy = self.source(); lines = hess.read_text().splitlines()
        duplicates = list(lines); index = duplicates.index('$vibrational_frequencies') + 2 + 6
        duplicates[index] = re.sub(r'^\s*6\b', '    5', duplicates[index])
        scale = list(lines); scale[scale.index('$frequency_scale_factor') + 1] = '0.95'
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'modified.hess'
            for content, message in [(duplicates, 'mode indices'), (scale, 'Frequency scale')]:
                path.write_text('\n'.join(content) + '\n')
                with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                    zpe.extract(log, path, atoms, energy)

    def test_wrong_summary_cutoff_or_multiple_thermochemistry_blocks_require_review(self):
        log, hess, atoms, energy = self.source(); original = log.read_text()
        variants = [
            (re.sub(r'(^Zero point energy\s+\.\.\.\s+)[-0-9.]+', r'\g<1>0.3', original, flags=re.M),
             'printed ZPE disagree'),
            (re.sub(r'(^Cut-Off Frequency\s+\.\.\.\s+)[-0-9.]+', r'\g<1>100.00', original, flags=re.M),
             'thermochemistry cutoff'),
            (original + '\nTHERMOCHEMISTRY AT 400K\n', 'one thermochemistry block'),
            (original.replace('ORCA TERMINATED NORMALLY', 'incomplete fixture'), 'normal completion'),
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'modified.log'
            for content, message in variants:
                path.write_text(content)
                with self.subTest(message=message), self.assertRaisesRegex(ValueError, message):
                    zpe.extract(path, hess, atoms, energy)

    def test_changed_archived_lineage_cannot_write_new_report(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); p = root / 'experiments/completion/drug-start-lineage-r1.json'
            p.parent.mkdir(parents=True); p.write_text('{}')
            with patch.object(zpe, 'ROOT', root), patch.object(zpe, 'revalidate_lineage', return_value={'changed': True}):
                with self.assertRaisesRegex(ValueError, 'evidence changed'): zpe.archived()
            self.assertFalse((root / 'reports').exists())


if __name__ == '__main__':
    unittest.main()
