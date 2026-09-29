"""Archived real SCFs for the adapter; exact quadratic test for the optimizer."""
from pathlib import Path
import sys
import unittest

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from completion_alloy_geometry import centroids, mo_groups
from completion_qe_bridge import load_scf
from completion_qe_relaxation import VerifiedQESCFCalculator, minimize_fixed_registry, settings_fingerprint


class Quadratic(Calculator):
    implemented_properties=['energy','forces']
    def __init__(self,target,stiffness):
        super().__init__();self.target=target;self.stiffness=stiffness
    def calculate(self,atoms=None,properties=('energy','forces'),system_changes=all_changes):
        super().calculate(atoms,properties,system_changes)
        d=atoms.positions-self.target
        self.results={'energy':float(.5*np.sum(self.stiffness[:,None]*d*d)),
                      'forces':-self.stiffness[:,None]*d}


class VerifiedSCFOptimizerTests(unittest.TestCase):
    def setUp(self):
        self.r0=ROOT/'reports/w2-bilayer-selected-r0-d6p3-pair.json'
        self.r1=ROOT/'reports/w2-bilayer-selected-r1-d6p3-pair.json'
        self.base=load_scf(self.r0)

    def test_archived_real_SCFs_update_only_at_the_requested_geometry(self):
        calls=[]
        def provider(atoms):calls.append(atoms.copy());return self.r1
        atoms=self.base['atoms'].copy();atoms.calc=VerifiedQESCFCalculator(atoms,self.r0,provider)
        np.testing.assert_array_equal(atoms.get_forces(),self.base['forces'])
        self.assertEqual(atoms.get_potential_energy(),self.base['energy_eV'])
        self.assertEqual(len(calls),0)
        other=load_scf(self.r1)
        atoms.set_positions(other['atoms'].positions)
        np.testing.assert_array_equal(atoms.get_forces(),other['forces'])
        self.assertEqual(atoms.get_potential_energy(),other['energy_eV'])
        self.assertEqual(len(calls),1)
        self.assertEqual(atoms.calc.revalidate_history(),2)

    def test_stale_report_rejected_and_cached_forces_cleared(self):
        atoms=self.base['atoms'].copy();atoms.calc=VerifiedQESCFCalculator(atoms,self.r0,lambda _:self.r0)
        atoms.get_forces();atoms.positions[1,2]+=.002
        with self.assertRaisesRegex(ValueError,'Requested coordinates'):atoms.get_forces()
        self.assertEqual(atoms.calc.results,{})
        self.assertEqual(len(atoms.calc.evaluation_history),1)

    def test_fixed_cell_species_or_charge_change_rejected_before_backend(self):
        for change in ('cell','species','charge'):
            with self.subTest(change=change):
                calls=[];atoms=self.base['atoms'].copy()
                atoms.calc=VerifiedQESCFCalculator(atoms,self.r0,lambda a:calls.append(a))
                if change=='cell':atoms.cell[2,2]+=.1
                if change=='species':atoms.numbers[1]=34
                if change=='charge':atoms.set_initial_charges(np.ones(len(atoms)))
                with self.assertRaises(ValueError):atoms.get_forces()
                self.assertEqual(calls,[])

    def test_numerical_settings_have_a_distinct_fingerprint(self):
        self.assertEqual(settings_fingerprint(self.r0),settings_fingerprint(self.r1))
        strict=ROOT/'reports/w2-bilayer-strict-r0-d6p3-pair.json'
        self.assertNotEqual(settings_fingerprint(self.r0),settings_fingerprint(strict))

    def test_analytic_24_atom_optimization_has_exact_constrained_solution(self):
        atoms=read(ROOT/'structures/completion/geometry-fixtures/CI-r0-d6.70.extxyz')
        origin=atoms.positions.copy();centers=centroids(atoms)
        rng=np.random.default_rng(712)
        target=origin+rng.normal(scale=.11,size=origin.shape)
        stiffness=np.linspace(.8,2.5,len(atoms))
        exact=target.copy()
        for indices in mo_groups(atoms):
            # Lagrange-multiplier solution for nonuniform quadratic springs.
            residual=target[indices].sum(axis=0)-origin[indices].sum(axis=0)
            exact[indices]-=(1/stiffness[indices])[:,None]*residual/np.sum(1/stiffness[indices])
        calc=Quadratic(target,stiffness)
        steps=[]
        result=minimize_fixed_registry(atoms,calc,fmax=1e-7,max_steps=100,maxstep=.05,
                                       initial_curvature=1.,on_step=steps.append,allow_analytic_test=True)
        self.assertTrue(result['converged'])
        self.assertEqual(result['evidence_kind'],'analytic_test_only')
        np.testing.assert_allclose(result['atoms'].positions,exact,atol=1e-6,rtol=0)
        np.testing.assert_allclose(centroids(result['atoms']),centers,atol=1e-12,rtol=0)
        np.testing.assert_array_equal(result['atoms'].cell,atoms.cell)
        self.assertGreater(result['trace'][-1]['raw_max_force_eV_angstrom'],.005)
        self.assertLess(result['trace'][-1]['projected_max_force_eV_angstrom'],1e-7)
        self.assertTrue(all(r['scf_report'] is None for r in steps))
        self.assertEqual(len(steps),result['trace'][-1]['step']+1)

    def test_default_curvature_reaches_planned_force_tolerance_in_analytic_test(self):
        atoms=read(ROOT/'structures/completion/geometry-fixtures/CI-r0-d6.70.extxyz')
        target=atoms.positions+np.random.default_rng(712).normal(scale=.11,size=atoms.positions.shape)
        result=minimize_fixed_registry(atoms,Quadratic(target,np.linspace(.8,2.5,len(atoms))),
                                       fmax=.005,max_steps=100,maxstep=.05,allow_analytic_test=True)
        self.assertTrue(result['converged'])
        self.assertEqual(result['initial_BFGS_curvature_eV_angstrom2'],70.)

    def test_short_optimization_budget_is_not_reported_as_convergence(self):
        atoms=read(ROOT/'structures/completion/geometry-fixtures/CI-r0-d6.70.extxyz')
        target=atoms.positions.copy();target[1,2]+=.3
        result=minimize_fixed_registry(atoms,Quadratic(target,np.ones(len(atoms))),
                                       fmax=1e-7,max_steps=1,maxstep=.03,allow_analytic_test=True)
        self.assertFalse(result['converged'])
        self.assertGreater(result['trace'][-1]['projected_max_force_eV_angstrom'],1e-7)

    def test_analytic_calculator_requires_an_explicit_test_label(self):
        atoms=read(ROOT/'structures/completion/geometry-fixtures/CI-r0-d6.70.extxyz')
        with self.assertRaisesRegex(ValueError,'explicitly marked analytic'):
            minimize_fixed_registry(atoms,Quadratic(atoms.positions,np.ones(len(atoms))),
                                    fmax=.005,max_steps=10,maxstep=.03)


if __name__=='__main__':unittest.main()
