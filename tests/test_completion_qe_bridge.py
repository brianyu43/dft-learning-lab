"""Real archived QE SCF data, deliberate corruption, and ASE state checks."""
import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np
from ase.calculators.calculator import PropertyNotImplementedError

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from completion_qe_bridge import BOHR_ANGSTROM, HARTREE_EV, load_scf


class VerifiedQEBridgeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.run = self.root / 'run'
        self.run.mkdir()
        self.report = self.root / 'report.json'
        self.meta = json.loads((ROOT / 'reports/w2-bilayer-selected-r0-d6p3-pair.json').read_text())
        source = Path(self.meta['run_dir'])
        for name in ['output.log', 'runner_snapshot.py', Path(self.meta['input']).name, *self.meta['assets'].values()]:
            target = self.run / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / name, target)
        original_xml = list((source / 'scratch').glob('*.save/data-file-schema.xml'))[0]
        self.xml = self.run / original_xml.relative_to(source)
        self.xml.parent.mkdir(parents=True)
        shutil.copy2(original_xml, self.xml)
        self.original_tree = ET.parse(self.xml)
        self.meta['run_dir'] = str(self.run)
        self.write_report()

    def write_report(self):
        self.report.write_text(json.dumps(self.meta))

    def test_real_scf_coordinates_energy_forces_and_actual_method(self):
        result = load_scf(self.report)
        out = self.original_tree.getroot().find('output')
        energy = float(out.findtext('total_energy/etot')) * HARTREE_EV
        forces = np.fromstring(out.findtext('forces'), sep=' ').reshape(-1, 3) * HARTREE_EV / BOHR_ANGSTROM
        self.assertAlmostEqual(result['atoms'].get_potential_energy(), energy, places=10)
        np.testing.assert_array_equal(result['atoms'].get_forces(), forces)
        self.assertEqual(result['audit']['openmp_threads_in_log'], 4)
        self.assertLess(result['audit']['xml_text_force_max_difference_eV_angstrom'], 1.4e-7)
        self.assertEqual(result['audit']['valence_electrons'], {'Mo': 14.0, 'S': 6.0})

    def test_calculator_never_reuses_results_at_changed_coordinates(self):
        result = load_scf(self.report)
        result['atoms'].positions[1, 2] += 0.002
        with self.assertRaises(PropertyNotImplementedError):
            result['atoms'].get_forces()

    def test_expected_geometry_mismatch_rejected(self):
        atoms = load_scf(self.report)['atoms']
        atoms.positions[2, 0] += 0.002
        with self.assertRaisesRegex(ValueError, 'Requested coordinates'):
            load_scf(self.report, expected_atoms=atoms)

    def test_c09_cannot_be_reported_as_pbe_d3bj(self):
        with self.assertRaisesRegex(ValueError, 'Expected PBE'):
            load_scf(self.report, expected_method='PBE-D3BJ')

    def test_xml_changes_are_rejected(self):
        changes = [('output/basis_set/ecutwfc', '45'),
                   ('input/electron_control/conv_thr', '5e-8'),
                   ('output/convergence_info/scf_conv/convergence_achieved', 'false'),
                   ('output/band_structure/nelec', '51'),
                   ('output/dft/functional', 'PBE'),
                   ('output/band_structure/spinorbit', 'true'),
                   ('exit_status', '1')]
        for path, text in changes:
            with self.subTest(path=path):
                tree = copy.deepcopy(self.original_tree)
                tree.getroot().find(path).text = text
                tree.write(self.xml)
                with self.assertRaises(ValueError):
                    load_scf(self.report)

    def test_force_sign_and_units_cannot_silently_change(self):
        tree = copy.deepcopy(self.original_tree)
        node = tree.getroot().find('output/forces')
        forces = np.fromstring(node.text, sep=' ')
        node.text = ' '.join(map(str, -forces))
        tree.write(self.xml)
        with self.assertRaisesRegex(ValueError, 'force unit or value'):
            load_scf(self.report)

    def test_xml_geometry_cannot_differ_from_scf_input(self):
        tree = copy.deepcopy(self.original_tree)
        node = tree.getroot().find('output/atomic_structure/atomic_positions')[1]
        position = np.fromstring(node.text, sep=' ')
        position[2] += 0.01
        node.text = ' '.join(map(str, position))
        tree.write(self.xml)
        with self.assertRaisesRegex(ValueError, 'SCF coordinates'):
            load_scf(self.report)

    def test_wrong_k_grid_or_openmp_metadata_rejected(self):
        old = copy.deepcopy(self.meta)
        for key, value in [('k_grid', [9, 9, 1]), ('omp_threads', 1), ('conv_thr_Ry', 1e-8)]:
            with self.subTest(key=key):
                self.meta = copy.deepcopy(old)
                self.meta[key] = value
                self.write_report()
                with self.assertRaises(ValueError):
                    load_scf(self.report)

    def test_unfinished_xml_rejected(self):
        tree = copy.deepcopy(self.original_tree)
        root = tree.getroot()
        root.remove(root.find('closed'))
        tree.write(self.xml)
        with self.assertRaisesRegex(ValueError, 'not final'):
            load_scf(self.report)


if __name__ == '__main__':
    unittest.main()
