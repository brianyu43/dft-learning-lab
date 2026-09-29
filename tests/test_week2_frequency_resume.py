"""Regression checks against a real completed Hessian plus corrupted copies."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from run_week2_frequency_resume import inspect_hessian

class HessianResumeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        meta=json.loads((ROOT/'reports/week1-water-opt-freq.json').read_text())
        run=Path(meta['run_dir']);stem=Path(meta['input']).stem
        cls.hess=run/f'{stem}.hess';cls.atoms=read(run/f'{stem}.xyz')
        cls.freq=meta['result']['frequencies_cm1']

    def test_real_hessian_matches_translated_geometry(self):
        atoms=self.atoms.copy();atoms.positions+=[1.2,-3.4,2.1]
        frequencies=inspect_hessian(self.hess,atoms,self.freq)
        self.assertTrue(np.all(frequencies>0))

    def test_changed_bond_geometry_rejected(self):
        atoms=self.atoms.copy();atoms.positions[1,0]+=.01
        with self.assertRaises(AssertionError):inspect_hessian(self.hess,atoms,self.freq)

    def test_changed_atom_identity_rejected(self):
        atoms=self.atoms.copy();atoms[0].symbol='N'
        with self.assertRaises(AssertionError):inspect_hessian(self.hess,atoms,self.freq)

    def test_incomplete_frequency_log_rejected(self):
        with self.assertRaises(AssertionError):inspect_hessian(self.hess,self.atoms,self.freq[:-1])

if __name__=='__main__':unittest.main()
