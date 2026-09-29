"""Geometric perturbations of marked fixtures, not DFT relaxation results."""
from pathlib import Path
import sys
import unittest

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from completion_qe_trajectory import RegistryGeometryAudit, audit_qe_trajectory
from completion_alloy_geometry import centroids, validate_geometry


class GeometryContinuityTests(unittest.TestCase):
    def setUp(self):
        self.initial = read(ROOT / 'structures/completion/geometry-fixtures/CI-r0-d6.70.extxyz')
        self.audit = RegistryGeometryAudit(self.initial)

    def test_all_27_marked_reference_fixtures_pass(self):
        paths = sorted((ROOT / 'structures/completion/geometry-fixtures').glob('*.extxyz'))
        self.assertEqual(len(paths), 27)
        for path in paths:
            with self.subTest(path=path.name):
                atoms = read(path)
                result = RegistryGeometryAudit(atoms).inspect(atoms)
                self.assertEqual(result['status'], 'geometry_continuity_passed')
                self.assertEqual(result['lost_distance_graph_edges'], [])

    def test_internal_thickness_change_is_allowed_but_actual_vacuum_is_recorded(self):
        atoms = self.initial.copy()
        centers = centroids(atoms)
        indices = np.flatnonzero(atoms.numbers != 42)
        signs = np.sign(atoms.positions[indices, 2] - centers[atoms.arrays['layer_id'][indices], 2])
        atoms.positions[indices, 2] += .025 * signs
        result = self.audit.inspect(atoms)
        self.assertEqual(result['status'], 'geometry_continuity_passed')
        self.assertAlmostEqual(result['current_vacuum_gap_angstrom'],
                               self.initial.info['bilayer_vacuum_gap_angstrom'] - .05)
        # The old validator intentionally covers only frozen structures.
        with self.assertRaisesRegex(ValueError, 'vacuum gap'): validate_geometry(atoms)

    def test_face_crossing_is_detected_even_with_unchanged_layer_labels(self):
        atoms = self.initial.copy(); index = atoms.info['exchange_indices'][0]
        atoms.positions[index, 2] = 2 * centroids(atoms)[0, 2] - atoms.positions[index, 2]
        result = self.audit.inspect(atoms)
        self.assertIn(index, result['changed_face_atom_indices'])
        self.assertEqual(result['status'], 'requires_review')
        np.testing.assert_array_equal(atoms.arrays['layer_id'], self.initial.arrays['layer_id'])

    def test_layer_migration_and_Mo_centroid_motion_are_detected(self):
        atoms = self.initial.copy(); index = atoms.info['exchange_indices'][0]
        atoms.positions[index, 2] = centroids(atoms)[1, 2] + .2
        self.assertIn('spatial_layer_ranges_overlap', self.audit.inspect(atoms)['reasons'])
        atoms = self.initial.copy(); atoms.positions[0, 2] += .002
        self.assertIn('Mo_centroids_moved', self.audit.inspect(atoms)['reasons'])

    def test_in_plane_site_exchange_changes_distance_graph(self):
        atoms = self.initial.copy()
        atoms.positions[[1, 7]] = atoms.positions[[7, 1]]
        result = self.audit.inspect(atoms)
        self.assertEqual(result['changed_face_atom_indices'], [])
        self.assertIn('reference_distance_graph_changed', result['reasons'])
        self.assertTrue(result['lost_distance_graph_edges'])

    def test_severe_chalcogen_approach_is_not_hidden_by_same_faces(self):
        atoms = self.initial.copy(); z = centroids(atoms)[0, 2]
        atoms.positions[1, 2] = z + .01; atoms.positions[2, 2] = z - .01
        result = self.audit.inspect(atoms)
        self.assertEqual(result['changed_face_atom_indices'], [])
        self.assertIn('severe_atomic_approach_requires_review', result['reasons'])

    def test_provenance_cell_and_nonfinite_coordinates_are_rejected(self):
        for change in ('provenance', 'cell', 'nan'):
            with self.subTest(change=change):
                atoms = self.initial.copy()
                if change == 'provenance': atoms.arrays['tile_ij'][0] = [1, 1]
                if change == 'cell': atoms.cell[2, 2] += .1
                if change == 'nan': atoms.positions[0, 0] = np.nan
                with self.assertRaises(ValueError): self.audit.inspect(atoms)

    def test_fixtures_and_missing_raw_reports_cannot_become_production_evidence(self):
        with self.assertRaisesRegex(ValueError, 'Empty or discontinuous'):
            audit_qe_trajectory(self.initial, '/not-a-report', [])
        with self.assertRaisesRegex(ValueError, 'Fixture'):
            audit_qe_trajectory(self.initial, '/not-a-report', [dict(step=0)])
        copied = self.initial.copy(); copied.info.pop('fixture_only')
        with self.assertRaisesRegex(ValueError, 'raw SCF report'):
            audit_qe_trajectory(copied, '/not-a-report', [dict(step=0)])


if __name__ == '__main__':
    unittest.main()
