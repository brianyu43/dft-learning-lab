"""Geometry/state bookkeeping failures must not silently become screening evidence."""
import importlib.util
import json
from pathlib import Path
import unittest
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('screen',ROOT/'scripts/run_week2_xtb_screen.py')
screen=importlib.util.module_from_spec(spec);spec.loader.exec_module(screen)

class GeometryAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data=json.loads((ROOT/'experiments/week2/chemistry-candidates.json').read_text())
        cls.neutral=next(c for c in data['dimers'] if c['species']=='diflunisal--caffeine' and c['motif']=='hydrogen_bond')
        cls.ion=next(c for c in data['dimers'] if c['motif']=='proton_transfer_start')

    def test_neutral_state_and_linear_hbond(self):
        c=self.neutral;a=read(ROOT/c['path']);result=screen.geometry_audit(c,a)
        self.assertTrue(result['chemical_state_preserved_by_geometry'])
        detail=c['motif_detail']
        self.assertAlmostEqual(a.get_angle(detail['donor_atom'],detail['hydrogen_atom'],detail['acceptor_atom']),180,places=5)

    def test_intended_proton_transfer_is_reported(self):
        c=self.ion;r=screen.geometry_audit(c,read(ROOT/c['path']))
        self.assertEqual(len(r['cross_fragment_proton_transfer']),1)
        self.assertTrue(r['heavy_connectivity_preserved'])
        self.assertFalse(r['chemical_state_preserved_by_geometry'])

    def test_separated_pair_not_treated_as_bound(self):
        c=self.neutral;a=read(ROOT/c['path']);a.positions[c['fragment_atom_counts'][0]:]+=np.array([100,0,0])
        self.assertTrue(screen.geometry_audit(c,a)['separated_geometry'])

    def test_composition_change_detected(self):
        c=self.neutral;a=read(ROOT/c['path']);a.numbers[0]=8 if a.numbers[0]!=8 else 6
        self.assertEqual(screen.geometry_audit(c,a)['status'],'atom_order_or_composition_changed')

if __name__=='__main__':unittest.main()
