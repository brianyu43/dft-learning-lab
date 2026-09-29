"""Analytic observables and scheduling safeguards; no model inference."""
from copy import deepcopy
import fcntl
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import completion_mace_benchmark as benchmark


class MaceBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.plan = {'cases': {'n6': {'area_angstrom2': 10., 'parts': {
            'pair': {'numbers': [42, 16, 16, 42, 34, 34]},
            'bottom': {'numbers': [42, 16, 16]}, 'top': {'numbers': [42, 34, 34]}}}}}
        self.reference = {'n6': {p: {'energy_eV': e, 'forces_eV_angstrom': np.zeros((n, 3)).tolist()}
                                for p, e, n in [('pair', -101., 6), ('bottom', -40., 3), ('top', -60., 3)]}}

    def compare(self, target):
        return benchmark.deltas(self.plan, self.reference, target)['n6']

    def test_exact_pair_minus_layers_units(self):
        row = self.compare(deepcopy(self.reference))
        self.assertEqual(row['reference_gamma_meV_angstrom2'], -100.)
        self.assertEqual(row['delta_gamma_meV_angstrom2'], 0.)
        self.assertTrue(row['passed'])

    def test_total_energy_can_pass_while_interaction_precision_fails(self):
        target = deepcopy(self.reference)
        target['n6']['pair']['energy_eV'] += .002
        row = self.compare(target)
        self.assertTrue(all(x['passed'] for x in row['parts'].values()))
        self.assertAlmostEqual(row['delta_gamma_meV_angstrom2'], .2)
        self.assertFalse(row['passed'])

    def test_extensive_offset_cancels_in_interaction(self):
        target = deepcopy(self.reference)
        for part, n in [('pair', 6), ('bottom', 3), ('top', 3)]:
            target['n6'][part]['energy_eV'] += n * .0005
        row = self.compare(target)
        self.assertAlmostEqual(row['delta_gamma_meV_angstrom2'], 0., places=10)
        self.assertTrue(row['passed'])

    def test_force_vector_norm_not_largest_component(self):
        target = deepcopy(self.reference)
        target['n6']['pair']['forces_eV_angstrom'][0] = [.003] * 3
        row = self.compare(target)
        self.assertGreater(row['parts']['pair']['max_delta_force_eV_angstrom'], .005)
        self.assertFalse(row['passed'])

    def test_nonfinite_missing_and_misshaped_values_are_rejected(self):
        bad = deepcopy(self.reference); bad['n6']['top']['energy_eV'] = float('nan')
        missing = deepcopy(self.reference); del missing['n6']['bottom']
        forces = deepcopy(self.reference); forces['n6']['pair']['forces_eV_angstrom'] = [[0., 0., 0.]]
        for target in (bad, missing, forces):
            with self.subTest(target=target), self.assertRaises(ValueError): self.compare(target)

    def test_controller_lock_refuses_execution_before_creating_run(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(benchmark, 'ROOT', Path(directory)):
            exp = Path(directory) / 'experiments'; exp.mkdir()
            with (exp / '.completion-controller.lock').open('a') as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                with self.assertRaisesRegex(RuntimeError, 'Another research engine'):
                    benchmark.run(self.plan)
            self.assertFalse((Path(directory) / 'runs').exists())

    def test_direct_worker_cannot_bypass_supervisor(self):
        with patch.dict(os.environ, {'MACE_BENCHMARK_LOCK_FDS': ''}):
            with self.assertRaisesRegex(ValueError, 'supervisor locks'):
                benchmark.worker(self.plan, 'cpu64_4', Path('/tmp/not-created-mace-test.json'))


if __name__ == '__main__': unittest.main()
