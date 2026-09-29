"""Scientific input and constraint checks; fixtures are not DFT predictions."""
from pathlib import Path
import io
import json
import sys
import unittest

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.optimize import BFGS
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from completion_alloy_geometry import (FixMoLayerCentroids, build_family, centroids,
                                       decoration_equivalence, mo_groups, validate_geometry)
from week3_qe import monolayer


class HarmonicTrap(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, targets, stiffness):
        super().__init__()
        self.targets = targets
        self.stiffness = stiffness

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        delta = atoms.positions - self.targets
        self.results = {'energy': float((0.5 * self.stiffness[:, None] * delta**2).sum()),
                        'forces': -self.stiffness[:, None] * delta}


def constraint_optimizer_trial():
    """Compare an actual BFGS trajectory with an independent analytic solution."""
    atoms = build_family(monolayer('S', 3.24), monolayer('Se', 3.24), 1, 6.7, 25)['CI']
    initial = atoms.positions.copy()
    cell = atoms.cell.array.copy()
    reference = centroids(atoms)
    targets = initial + np.random.default_rng(17).normal(0, 0.15, initial.shape)
    stiffness = np.linspace(0.7, 2.1, len(atoms))
    exact = targets.copy()
    for group, center in zip(mo_groups(atoms), reference):
        lagrange = (len(group) * center - targets[group].sum(axis=0)) / (1 / stiffness[group]).sum()
        exact[group] += lagrange / stiffness[group, None]
    atoms.calc = HarmonicTrap(targets, stiffness)
    atoms.set_constraint(FixMoLayerCentroids.from_atoms(atoms))
    history = []
    optimizer = BFGS(atoms, logfile=None, alpha=1.0)
    optimizer.attach(lambda: history.append(centroids(atoms).copy()))
    converged = optimizer.run(fmax=1e-8, steps=100)
    return dict(converged=bool(converged), optimizer_steps=optimizer.nsteps,
                maximum_centroid_drift_angstrom=float(np.abs(np.array(history) - reference).max()),
                maximum_solution_error_angstrom=float(np.abs(atoms.positions - exact).max()),
                maximum_cell_drift_angstrom=float(np.abs(atoms.cell.array - cell).max()),
                maximum_displacement_angstrom=float(np.linalg.norm(atoms.positions - initial, axis=1).max()),
                maximum_projected_force_eV_angstrom=float(np.linalg.norm(atoms.get_forces(), axis=1).max()),
                maximum_raw_force_eV_angstrom=float(np.linalg.norm(atoms.get_forces(apply_constraint=False), axis=1).max()),
                scope='Analytic harmonic test only; not a QE relaxation or material energy')


class AlloyGeometryTests(unittest.TestCase):
    def family(self, rid=0, distance=6.7):
        return build_family(monolayer('S', 3.24), monolayer('Se', 3.24), rid, distance, 25)

    def test_full_matrix_composition_registry_and_inner_outer_faces(self):
        for rid in range(3):
            for distance in (6.45, 6.7, 6.95):
                with self.subTest(rid=rid, distance=distance):
                    family = self.family(rid, distance)
                    for label, atoms in family.items():
                        self.assertEqual(validate_geometry(atoms)['total_formula'], 'Mo8S8Se8')
                        self.assertAlmostEqual(centroids(atoms)[1, 2] - centroids(atoms)[0, 2], distance)
                        np.testing.assert_array_equal(atoms.positions, family['C0'].positions)
                        np.testing.assert_array_equal(atoms.cell, family['C0'].cell)
                    for layer in (0, 1):
                        i = family['CI'].info['exchange_indices'][layer]
                        o = family['CO'].info['exchange_indices'][layer]
                        np.testing.assert_allclose(family['CI'].positions[i, :2], family['CO'].positions[o, :2], atol=1e-10)

    def test_primitive_registry_and_parent_sources_preserved(self):
        bottom, top = monolayer('S', 3.24), monolayer('Se', 3.24)
        bpos, tpos = bottom.positions.copy(), top.positions.copy()
        atoms = build_family(bottom, top, 1, 6.7, 25)['CI']
        shift = centroids(atoms)[1] - centroids(atoms)[0]
        np.testing.assert_allclose(shift[:2], ((bottom.cell[0] + 2 * bottom.cell[1]) / 3)[:2], atol=1e-12)
        np.testing.assert_array_equal(bottom.positions, bpos)
        np.testing.assert_array_equal(top.positions, tpos)

    def test_wrong_lattice_overlapping_layers_and_corrupt_annotations_rejected(self):
        with self.assertRaisesRegex(ValueError, 'in-plane'):
            build_family(monolayer('S', 3.18), monolayer('Se', 3.30), 0, 6.7, 25)
        with self.assertRaisesRegex(ValueError, 'overlap'):
            self.family(distance=1.0)
        atoms = self.family()['CI']
        atoms.info['configuration_id'] = 'CO'
        with self.assertRaisesRegex(ValueError, 'wrong face'):
            validate_geometry(atoms)
        atoms = self.family()['CI']
        atoms.positions[0, 0] += 0.1
        with self.assertRaisesRegex(ValueError, 'registry'):
            validate_geometry(atoms)

    def test_ci_co_inequivalence_all_nine_geometries_and_three_tolerances(self):
        for rid in range(3):
            for distance in (6.45, 6.7, 6.95):
                family = self.family(rid, distance)
                for tolerance in (1e-5, 1e-4, 1e-3):
                    with self.subTest(rid=rid, distance=distance, tolerance=tolerance):
                        result = decoration_equivalence(family['CI'], family['CO'], tolerance)
                        self.assertFalse(result['equivalent'])
                        self.assertEqual(result['operations_checked'], result['parent_operation_count'])
                        self.assertGreater(result['parent_operation_count'], 1)

    def test_equivalence_positive_control_permutation_wrapping_and_origin_shift(self):
        first = self.family(2)['CI']
        second = first[np.random.default_rng(43).permutation(24)]
        second.positions += np.array([0.2, -0.3, 0.4]) @ second.cell.array
        second.wrap()
        result = decoration_equivalence(first, second)
        self.assertTrue(result['equivalent'])
        self.assertLess(result['max_residual_angstrom'], 1e-10)

    def test_different_parent_rejected(self):
        family = self.family()
        family['CO'].positions[0, 0] += 0.02
        with self.assertRaisesRegex(ValueError, 'parent geometries differ'):
            decoration_equivalence(family['CI'], family['CO'])

    def test_nontrivial_rotation_positive_control(self):
        first = build_family(monolayer('S', 3.24), monolayer('Se', 3.24), 0, 6.7, 25,
                             top_exchange_tile=(0, 1))['CI']
        angle = 2 * np.pi / 3
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                             [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        second = first.copy()
        second.positions = second.positions @ rotation.T
        second.wrap()
        result = decoration_equivalence(first, second)
        self.assertTrue(result['equivalent'])
        self.assertFalse(np.array_equal(result['rotation_fractional'], np.eye(3)))

    def test_structure_round_trip_remains_verifiable_and_json_serializable(self):
        for atoms in self.family(1).values():
            stream = io.StringIO()
            write(stream, atoms, format='extxyz')
            recovered = read(io.StringIO(stream.getvalue()), format='extxyz')
            check = validate_geometry(recovered)
            self.assertEqual(json.loads(json.dumps(check))['total_formula'], 'Mo8S8Se8')


class LayerConstraintTests(unittest.TestCase):
    def atoms(self):
        atoms = build_family(monolayer('S', 3.24), monolayer('Se', 3.24), 1, 6.7, 25)['CI']
        atoms.set_constraint(FixMoLayerCentroids.from_atoms(atoms))
        return atoms

    def test_multiple_steps_do_not_move_other_groups_or_freeze_internal_motion(self):
        atoms = self.atoms()
        reference = centroids(atoms)
        free = atoms.numbers != 42
        rng = np.random.default_rng(89)
        for _ in range(100):
            trial = atoms.positions + rng.normal(0, 0.01, (24, 3))
            expected_chalcogens = trial[free].copy()
            atoms.set_positions(trial)
            np.testing.assert_array_equal(atoms.positions[free], expected_chalcogens)
            np.testing.assert_allclose(centroids(atoms), reference, atol=1e-12, rtol=0)
        self.assertEqual(atoms.constraints[0].get_removed_dof(atoms), 6)

    def test_projected_force_matches_tangent_energy_derivative(self):
        atoms = self.atoms()
        rng = np.random.default_rng(17)
        direction = rng.normal(size=(24, 3))
        for group in mo_groups(atoms):
            direction[group] -= direction[group].mean(axis=0)
        target = atoms.positions + rng.normal(0, 0.1, (24, 3))
        atoms.calc = HarmonicTrap(target, np.linspace(0.7, 2.1, 24))
        original = atoms.positions.copy()
        forces = atoms.get_forces()
        raw = atoms.get_forces(apply_constraint=False)
        np.testing.assert_array_equal(forces[atoms.numbers != 42], raw[atoms.numbers != 42])
        for group in mo_groups(atoms):
            np.testing.assert_allclose(forces[group].sum(axis=0), 0, atol=1e-14)
        eps = 1e-5
        atoms.set_positions(original + eps * direction)
        plus = atoms.get_potential_energy()
        atoms.set_positions(original - eps * direction)
        minus = atoms.get_potential_energy()
        self.assertAlmostEqual((plus - minus) / (2 * eps), float(-(forces * direction).sum()), places=8)

    def test_actual_optimizer_matches_constrained_analytic_minimum(self):
        result = constraint_optimizer_trial()
        self.assertTrue(result['converged'])
        self.assertLess(result['maximum_centroid_drift_angstrom'], 1e-12)
        self.assertLess(result['maximum_solution_error_angstrom'], 2e-8)
        self.assertEqual(result['maximum_cell_drift_angstrom'], 0)
        self.assertGreater(result['maximum_displacement_angstrom'], 0.1)
        self.assertGreater(result['maximum_raw_force_eV_angstrom'], 1e-3)

    def test_serialization_permutation_and_slice_safety(self):
        atoms = self.atoms()
        restored = FixMoLayerCentroids.fromdict(atoms.constraints[0].todict())
        np.testing.assert_array_equal(restored.reference_centroids, centroids(atoms))
        permuted = atoms[np.random.default_rng(43).permutation(24)]
        permuted.set_positions(permuted.positions + 0.1)
        np.testing.assert_allclose(centroids(permuted), centroids(atoms), atol=1e-12)
        with self.assertRaisesRegex(ValueError, 'must not be sliced'):
            _ = atoms[:12]


if __name__ == '__main__':
    unittest.main()
