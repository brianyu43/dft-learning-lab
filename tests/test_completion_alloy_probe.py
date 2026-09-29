"""Reference gating and exact component matching, using no new DFT runs."""
from copy import deepcopy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from ase.io import read
from ase.io.espresso import read_fortran_namelist

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import prepare_completion_alloy_probe as probe
from week3_qe import monolayer


class AlloyProbeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # These values are intentionally synthetic, never production evidence.
        cls.bottom, cls.top = monolayer('S', 3.24), monolayer('Se', 3.24)
        cls.settings = probe.settings_from_reference(dict(ecutwfc_Ry=100, k=9, vacuum_angstrom=25))
        cls.pseudos = dict(Mo='fixture-mo', S='fixture-s', Se='fixture-se')
        cls.rows, cls.calculations = probe.build_catalog(cls.bottom, cls.top, 6.7, cls.settings,
                                                       cls.pseudos, 'fixture-executable')

    def test_unfinished_reference_cannot_write_production_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref = root / probe.REFERENCE
            ref.parent.mkdir(parents=True)
            ref.write_text(json.dumps(dict(status='running')))
            with patch.object(probe, 'ROOT', root):
                with self.assertRaisesRegex(ValueError, 'not passed'):
                    probe.prepare()
            self.assertEqual(sorted(str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()),
                             [probe.REFERENCE])

    def test_unweakened_supercell_k_spacing_and_unreviewed_setting_rejection(self):
        for cut, k, vacuum, expected in [(80, 6, 20, 3), (100, 9, 25, 5), (120, 12, 30, 6)]:
            settings = probe.settings_from_reference(dict(ecutwfc_Ry=cut, k=k, vacuum_angstrom=vacuum))
            self.assertEqual(settings['k_grid'], [expected, expected, 1])
            self.assertGreaterEqual(settings['primitive_equivalent_k'], k)
        with self.assertRaisesRegex(ValueError, 'Unreviewed'):
            probe.settings_from_reference(dict(ecutwfc_Ry=60, k=4, vacuum_angstrom=15))

    def test_all_components_preserve_actual_pair_indices_and_cell(self):
        keys = {(r['configuration'], r['registry_id'], r['distance_index']) for r in self.rows}
        self.assertEqual(keys, {(c, r, d) for c in ('C0', 'CI', 'CO') for r in range(3) for d in range(3)})
        for row in self.rows:
            pair = self.calculations[row['calculations']['pair']]['atoms']
            self.assertEqual(len(pair), 24)
            for part, layer in [('bottom', 0), ('top', 1)]:
                component = self.calculations[row['calculations'][part]]['atoms']
                indices = pair.arrays['layer_id'] == layer
                self.assertEqual(len(component), 12)
                np.testing.assert_array_equal(component.positions, pair.positions[indices])
                np.testing.assert_array_equal(component.numbers, pair.numbers[indices])
                np.testing.assert_array_equal(component.cell.array, pair.cell.array)
            self.assertFalse(row['ci_co_symmetry']['equivalent'])

    def test_reuse_exactly_matches_identity_without_centering_or_tolerance(self):
        row = next(r for r in self.rows if (r['configuration'], r['registry_id'], r['distance_index']) == ('CI', 0, 1))
        key = row['calculations']['bottom']
        reused = next(r for r in self.rows if (r['configuration'], r['registry_id'], r['distance_index']) == ('CI', 2, 1))
        self.assertEqual(key, reused['calculations']['bottom'])
        atoms = self.calculations[key]['atoms'].copy()
        atoms.positions[0, 0] += 1e-12
        changed = probe.calculation_identity(atoms, self.settings, self.pseudos, 'fixture-executable')
        self.assertNotEqual(key, probe.digest(changed))
        for field, value in [('settings', {**self.settings, 'conv_thr_Ry': 1e-12}),
                             ('pseudopotential_sha256', {**self.pseudos, 'S': 'changed'}),
                             ('executable_sha256', 'changed')]:
            identity = deepcopy(self.calculations[key]['identity'])
            identity[field] = value
            self.assertNotEqual(key, probe.digest(identity))

    def test_serialized_qe_pair_and_layers_keep_identical_conditions(self):
        row = next(r for r in self.rows if (r['configuration'], r['registry_id'], r['distance_index']) == ('CI', 0, 1))
        recovered, settings = {}, []
        for part, key in row['calculations'].items():
            atoms = self.calculations[key]['atoms']
            text = probe.render_scf('fixture-' + part, atoms, self.settings)
            recovered[part] = read(io.StringIO(text), format='espresso-in')
            data, _ = read_fortran_namelist(io.StringIO(text))
            settings.append((data['system']['ecutwfc'], data['system']['ecutrho'], data['electrons']['conv_thr']))
            self.assertIn('5 5 1  0 0 0', text)
            self.assertEqual(data['control']['calculation'], 'scf')
        self.assertEqual(settings, [(100, 400, 1e-10)] * 3)
        pair = recovered['pair']
        original = self.calculations[row['calculations']['pair']]['atoms']
        for part, layer in [('bottom', 0), ('top', 1)]:
            np.testing.assert_array_equal(recovered[part].positions, pair.positions[original.arrays['layer_id'] == layer])
            np.testing.assert_array_equal(recovered[part].cell.array, pair.cell.array)

    def test_passed_flag_alone_is_not_accepted_as_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'reference.json'
            path.write_text(json.dumps(dict(status='passed', setting_selection=dict(ecutwfc_Ry=100, k=9,
                            vacuum_angstrom=25), heterobilayer_numerical_sensitivity=dict(passed=False))))
            with self.assertRaisesRegex(ValueError, 'precision did not pass'):
                probe.validate_reference(path)


if __name__ == '__main__':
    unittest.main()
