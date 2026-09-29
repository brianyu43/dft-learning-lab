"""Archived reference/input tests and mocked lifecycle tests; no QE launches."""
import fcntl
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import numpy as np
from ase.io import read
from ase.io.espresso import read_fortran_namelist

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import completion_qe_provider as provider


class ProviderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = ROOT / 'reports/w2-bilayer-selected-r0-d6p3-pair.json'
        cls.meta = json.loads(cls.report.read_text())
        cls.atoms = provider.load_scf(cls.report)['atoms']

    def make(self, validator=lambda: None):
        return provider.BoundedQESCFProvider(self.report, self.atoms, 'provider-unit-fixture',
                    validate_eligibility=validator, maximum_new_scfs=2,
                    per_scf_timeout_seconds=100, total_seconds_limit=600)

    def test_reference_conditions_survive_new_coordinate_serialization(self):
        atoms = self.atoms.copy(); atoms.positions[1, 2] += .002
        text = provider.render_from_reference(self.meta, atoms, 'serialized-probe')
        output = read(io.StringIO(text), format='espresso-in')
        np.testing.assert_allclose(output.positions, atoms.positions, rtol=0, atol=1e-9)
        np.testing.assert_allclose(output.cell, atoms.cell, rtol=0, atol=1e-9)
        source = Path(self.meta['run_dir']) / Path(self.meta['input']).name
        with source.open() as handle: original, _ = read_fortran_namelist(handle)
        restored, _ = read_fortran_namelist(io.StringIO(text))
        self.assertEqual(dict(original['system']), dict(restored['system']))
        self.assertEqual(dict(original['electrons']), dict(restored['electrons']))
        self.assertEqual(restored['control']['prefix'], 'serialized_probe')
        self.assertNotEqual(provider.geometry_key(atoms), provider.geometry_key(self.atoms))

    def test_failed_eligibility_or_missing_context_cannot_launch(self):
        def blocked(): raise ValueError('Reference and cost prerequisites not passed')
        instance = self.make(blocked)
        with patch.object(instance, '_launch') as launch:
            with self.assertRaisesRegex(ValueError, 'prerequisites'):
                instance.__enter__()
            with self.assertRaisesRegex(ValueError, 'exclusive context'):
                instance(self.atoms)
            launch.assert_not_called()
        self.assertFalse(instance.directory.exists())

    def test_exclusive_controller_lock_prevents_a_second_provider(self):
        instance = self.make()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / 'experiments').mkdir(); (root / 'runs').mkdir()
            instance.directory = root / 'runs/probe'
            with (root / 'experiments/.completion-controller.lock').open('a') as held:
                fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
                # Isolate locking from the separately tested frozen-file audit.
                with patch.object(provider, 'ROOT', root), patch.object(instance, '_verify_frozen'):
                    with self.assertRaises(BlockingIOError): instance.__enter__()
            self.assertFalse(instance.directory.exists())
            self.assertIsNone(instance.lock)

    def test_same_geometry_reuses_only_revalidated_raw_reference(self):
        instance = self.make(); instance.lock = object()
        with patch.object(instance, '_launch') as launch:
            self.assertEqual(instance(self.atoms), self.report)
            launch.assert_not_called()
        self.assertEqual(instance.state['steps'], [])

    def test_changed_cell_charge_dependency_and_unresolved_attempt_stop_before_launch(self):
        for change in ('cell', 'charge', 'dependency', 'failed_attempt'):
            with self.subTest(change=change):
                instance = self.make(); instance.lock = object(); atoms = self.atoms.copy()
                if change == 'cell': atoms.cell[2, 2] += .1
                if change == 'charge': atoms.set_initial_charges(np.ones(len(atoms)))
                if change == 'dependency': instance.frozen[next(iter(instance.frozen))] = '0' * 64
                if change == 'failed_attempt': instance.state['steps'] = [dict(status='prepared')]
                with patch.object(instance, '_launch') as launch:
                    with self.assertRaises(ValueError): instance(atoms)
                    launch.assert_not_called()

    def test_both_trajectory_and_original_pilot_budget_must_cover_next_scf(self):
        instance = self.make(); instance.lock = object(); instance.started = 0
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / 'reports').mkdir()
            (root / 'reports/week2-summary.json').write_text(json.dumps(
                dict(remaining_budget_hours=dict(semiconductor=1))))
            with patch.object(provider, 'ROOT', root), patch.object(provider.time, 'monotonic', return_value=20):
                with patch.object(provider, 'spent', return_value={'qe': 3501}):
                    with self.assertRaisesRegex(ValueError, 'pilot budget'): instance._budget()
                with patch.object(provider, 'spent', return_value={'qe': 3499}):
                    self.assertEqual(instance._budget()['pilot_remaining_seconds'], 101)
            with patch.object(provider.time, 'monotonic', return_value=501):
                with self.assertRaisesRegex(ValueError, 'Trajectory time'): instance._budget()

    def test_subprocess_deadline_and_cleanup_are_wired_without_running_qe(self):
        instance = self.make(); fake = MagicMock(pid=12345)
        fake.wait.return_value = 0; fake.poll.return_value = 0
        with tempfile.TemporaryDirectory() as directory:
            instance.directory = Path(directory)
            with patch.object(provider.subprocess, 'Popen', return_value=fake) as popen, patch.object(instance, '_save'):
                self.assertEqual(instance._launch(Path('/not-executed/config.json')), 0)
                self.assertTrue(popen.call_args.kwargs['start_new_session'])
                fake.wait.assert_called_with(timeout=160)
                fake.wait.side_effect = subprocess.TimeoutExpired('mock runner', 160)
                with patch.object(instance, '_stop_child') as stop:
                    with self.assertRaises(subprocess.TimeoutExpired):
                        instance._launch(Path('/not-executed/config.json'))
                    stop.assert_called_once()


if __name__ == '__main__':
    unittest.main()
