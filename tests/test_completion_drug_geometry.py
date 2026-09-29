"""Symmetry/proton-state guards using archived coordinates; no new DFT."""
from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from completion_drug_geometry import compare, describe, align, inline_orca_geometry


class DrugGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = json.loads((ROOT / 'experiments/completion/drug-start-lineage-r1.json').read_text())
        cls.rows = {r['candidate_id']: r for r in source['rows']}

    def load(self, name):
        row = self.rows[name]
        return deepcopy(row['candidate']), read(ROOT / row['source_geometry'])

    def test_rotation_translation_and_identical_H_permutation(self):
        c, a = self.load('caffeine-t00-s17'); b = a.copy(); before = a.positions.copy()
        b.rotate(37, 'z'); b.translate([3, -2, 7])
        owners = c['graph']['hydrogen_owners']
        for owner in set(owners.values()):
            group = [int(h) for h, x in owners.items() if x == owner]
            b.positions[group] = b.positions[np.roll(group, 1)]
        result = compare(c, a, c, b, context_a={'test': 'same'}, context_b={'test': 'same'}, energy_a=-1., energy_b=-1.)
        self.assertTrue(result['duplicate_minimum_verified'])
        self.assertLess(result['all_atom_RMSD_upper_bound_angstrom'], 1e-8)
        np.testing.assert_array_equal(a.positions, before)

    def test_identical_homodimer_fragments_can_exchange_as_whole_system(self):
        c, a = self.load('diflunisal--diflunisal-p05'); b = a.copy(); n = c['fragment_atom_counts'][0]
        b.positions = np.concatenate([a.positions[n:], a.positions[:n]])
        result = compare(c, a, c, b, context_a={'test': 'same'}, context_b={'test': 'same'}, energy_a=-1., energy_b=-1.)
        self.assertTrue(result['duplicate_minimum_verified'])
        self.assertTrue(result['fragments_aligned_together'])
        self.assertGreaterEqual(result['heavy_graph_mapping_count'], 2)

    def test_energy_context_and_tolerance_are_required_for_duplicate_certificate(self):
        c, a = self.load('caffeine-t00-s17')
        for ca, cb, energy in [({'grid': 1}, {'grid': 2}, -1.), ({'grid': 1}, {'grid': 1}, -1. + .02 / 627.509474)]:
            result = compare(c, a, c, a, context_a=ca, context_b=cb, energy_a=-1., energy_b=energy)
            self.assertFalse(result['duplicate_minimum_verified'])
            if ca != cb: self.assertIsNone(result['energy_gap_kcal_mol'])

    def test_proton_shift_changes_state_even_with_identical_heavy_coordinates(self):
        c, ion = self.load('diflunisal--adenine-p36'); neutral = ion.copy()
        vector = neutral.positions[2] - neutral.positions[28]; vector /= np.linalg.norm(vector)
        neutral.positions[18] = neutral.positions[2] + vector
        self.assertEqual(describe(c, ion)['hydrogen_owners'][18], 28)
        self.assertEqual(describe(c, neutral)['hydrogen_owners'][18], 2)
        result = compare(c, ion, c, neutral)
        self.assertEqual(result['status'], 'different_geometry_assigned_state_or_graph')

    def test_existing_caffeine_poses_have_different_geometries_and_contacts(self):
        c, a = self.load('diflunisal--caffeine-p14'); d, b = self.load('diflunisal--caffeine-p06')
        first, second = describe(c, a), describe(d, b)
        self.assertEqual((first['contacts'][0]['donor'], first['contacts'][0]['acceptor']), (17, 29))
        self.assertEqual((second['contacts'][0]['donor'], second['contacts'][0]['acceptor']), (2, 32))
        result = compare(c, a, d, b)
        self.assertEqual(result['status'], 'geometrically_separated_under_declared_tolerance')
        self.assertFalse(result['energy_comparable'])
        self.assertGreater(result['all_atom_RMSD_lower_bound_angstrom'], .05)

    def test_alignment_does_not_use_reflection(self):
        p = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        q = p.copy(); q[:, 0] *= -1
        rms, rotation, _, _ = align(p, q)
        self.assertAlmostEqual(np.linalg.det(rotation), 1.)
        self.assertGreater(rms, .1)

    def test_original_inline_input_is_read_and_H_overlap_rejected(self):
        atoms = inline_orca_geometry(ROOT / 'inputs/orca/week2/w2-diflunisal-opt-freq.inp')
        self.assertEqual(len(atoms), 26)
        c, a = self.load('caffeine-t00-s17')
        h, owner = next(iter(c['graph']['hydrogen_owners'].items()))
        a.positions[int(h)] = a.positions[owner]
        with self.assertRaisesRegex(ValueError, 'overlaps'): describe(c, a)


if __name__ == '__main__':
    unittest.main()
