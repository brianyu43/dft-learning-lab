"""Synthetic failure fixtures test gates; they are never presented as engine runs."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('research_run', ROOT / 'scripts/research_run.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)

class ParserGates(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name)
        self.meta = {'engine': 'orca', 'n_atoms': 3, 'linear': False, 'exit_code': 0,
                     'requested': {'optimization': True, 'frequencies': True}}
        self.good = '\n'.join([
            'Program Version 6.1.1', 'SCF CONVERGED', 'FINAL SINGLE POINT ENERGY -76.4',
            'THE OPTIMIZATION HAS CONVERGED', 'VIBRATIONAL FREQUENCIES',
            *[f'{i}: {f:.2f} cm**-1' for i, f in enumerate([0]*6 + [1600,3800,3900])],
            'NORMAL MODES', 'ORCA TERMINATED NORMALLY'])

    def tearDown(self):
        self.temp.cleanup()

    def parse(self, text, meta=None):
        (self.directory / 'output.log').write_text(text)
        return runner.parse_output(self.directory, meta or self.meta)

    def test_complete_minimum(self):
        self.assertEqual(self.parse(self.good)['status'], 'passed')

    def test_small_imaginary_mode_is_not_waived(self):
        result = self.parse(self.good.replace('1600.00', '-1.00'))
        self.assertEqual(result['status'], 'needs_review')
        self.assertFalse(result['minimum_verified'])

    def test_truncated_modes_rejected(self):
        self.assertEqual(self.parse(self.good.replace('8: 3900.00 cm**-1', ''))['status'], 'needs_review')

    def test_duplicate_mode_indices_rejected(self):
        self.assertEqual(self.parse(self.good.replace('8: 3900', '7: 3900'))['status'], 'needs_review')

    def test_normal_exit_does_not_imply_opt_convergence(self):
        self.assertEqual(self.parse(self.good.replace('THE OPTIMIZATION HAS CONVERGED', 'OPT FAILED'))['status'], 'needs_review')

    def test_incomplete_engine_exit_rejected(self):
        self.assertEqual(self.parse(self.good.replace('ORCA TERMINATED NORMALLY', ''))['status'], 'needs_review')

    def test_scf_failure_rejected(self):
        self.assertEqual(self.parse(self.good + '\nSCF NOT CONVERGED')['status'], 'needs_review')

    def test_timeout_or_exit_code_rejected(self):
        for changes in ({'exit_code': 1}, {'stop_reason': 'timeout'}):
            self.assertEqual(self.parse(self.good, {**self.meta, **changes})['status'], 'needs_review')

    def test_qe_scf_is_not_relaxation(self):
        meta = {'engine': 'qe', 'exit_code': 0, 'requested': {'optimization': True}}
        output = 'Program PWSCF v.7.6\n! total energy = -100.2 Ry\nconvergence has been achieved\nJOB DONE.'
        self.assertEqual(self.parse(output, meta)['status'], 'needs_review')
        result = self.parse(output + '\nbfgs converged in 4 scf cycles\n', meta)
        self.assertEqual(result['status'], 'passed')
        self.assertIsNone(result['minimum_verified'])

if __name__ == '__main__':
    unittest.main()
