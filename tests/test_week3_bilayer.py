"""Reject mismatched cells and verify the physical definition of the reference grid."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from week3_qe import monolayer
from week3_materials_pipeline import bilayer


class BilayerDefinitionTests(unittest.TestCase):
    def test_mismatched_lattices_rejected(self):
        with self.assertRaises(AssertionError):bilayer(monolayer('S',3.18),monolayer('Se',3.30),0,6.6,20)

    def test_registry_uses_primitive_vectors_and_vacuum_uses_outer_atoms(self):
        b=monolayer('S',3.24);t=monolayer('Se',3.24)
        zero=bilayer(b,t,0,6.7,25);shift=bilayer(b,t,1,6.7,25)
        self.assertEqual(shift.get_chemical_formula(),'Mo2S2Se2')
        np.testing.assert_allclose(shift.positions[:3],zero.positions[:3],rtol=0,atol=1e-12)
        expected=np.tile((b.cell[0]+2*b.cell[1])/3,(3,1))
        np.testing.assert_allclose(shift.positions[3:]-zero.positions[3:],expected,rtol=0,atol=1e-12)
        self.assertAlmostEqual(shift.positions[3,2]-shift.positions[0,2],6.7)
        self.assertAlmostEqual(shift.cell[2,2]-np.ptp(shift.positions[:,2]),25)
        np.testing.assert_allclose(shift[:3].cell,shift.cell)


if __name__=='__main__':unittest.main()
