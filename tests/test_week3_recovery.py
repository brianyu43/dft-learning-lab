"""Recovery gates: archived real modes/gradients and isolated corruption fixtures."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
from ase.io import read

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from research_run import sha
from inspect_week3_mode import read_modes, rigid_fraction
from week3_recovery import PLAN, load_recovery, gradient_in_frame, verify_numerical_check, HARTREE_KCAL


class ArchivedModeAndGradient(unittest.TestCase):
    def setUp(self):
        self.meta=json.loads((ROOT/'reports/w3-diflunisal--diflunisal-p05-opt-freq.json').read_text())
        self.run=Path(self.meta['run_dir']);self.stem=Path(self.meta['input']).stem

    def test_archived_mode_is_cartesian_internal_fragment_motion(self):
        atoms,modes=read_modes(self.run/(self.stem+'.hess'))
        vector=modes[:,6].reshape(-1,3)
        self.assertAlmostEqual(float(np.linalg.norm(vector)),1,places=8)
        self.assertLess(np.linalg.norm(np.sum(atoms.get_masses()[:,None]*vector,axis=0)),1e-7)
        for ids in [np.arange(26),np.arange(26,52)]:
            self.assertGreater(rigid_fraction(atoms,vector,ids)['fraction_explained'],.95)

    def test_real_engrad_matches_geometry_up_to_proper_rotation(self):
        atoms=read(self.run/(self.stem+'.xyz'))
        angle=.37
        rot=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
        e,g,_,_=gradient_in_frame(self.meta,atoms)
        atoms.positions=atoms.positions@rot+np.array([2.,-1.,.7])
        e2,g2,_,_=gradient_in_frame(self.meta,atoms)
        self.assertEqual(e,e2)
        np.testing.assert_allclose(g2,g@rot,rtol=0,atol=1e-10)
        atoms.positions[0,0]+=.001
        with self.assertRaises(AssertionError):gradient_in_frame(self.meta,atoms)


class FrozenPlan(unittest.TestCase):
    def test_real_parent_is_retained_and_changed_parent_is_rejected(self):
        plan,items=load_recovery()
        self.assertEqual(len(items),1)
        self.assertEqual(plan['new_independent_start_count'],0)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            paths=[PLAN,'experiments/week3/drug/panel-plan.json',plan['inspection']]
            for item in plan['replacements']:
                paths += [item[k] for k in ['parent_config','parent_report','replacement_config']]
                retry=json.loads((ROOT/item['replacement_config']).read_text())
                paths += [retry['input'],next(p for p,v in retry['assets'].items() if v=='initial.xyz')]
            for p in paths:
                dst=root/p;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(ROOT/p,dst)
            load_recovery(root)
            p=root/plan['replacements'][0]['parent_report']
            p.write_text(p.read_text()+' ')
            with self.assertRaises(AssertionError):load_recovery(root)

    def test_numerical_gate_withholds_inconsistent_or_failed_result(self):
        # Synthetic metadata exercises the promotion gate, not an engine validation.
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            def put(p,data):
                dest=root/p;dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_text(json.dumps(data));return sha(dest)
            policy={'max_absolute_energy_delta_kcal_mol':.1,'max_gradient_component_delta_Hartree_Bohr':1e-4}
            planhash=put(PLAN,{'common_grid_policy':policy})
            rhash=put('refine-config.json',{'job_id':'refine'})
            bhash=put('base-config.json',{})
            def meta(config_hash,energy,minimum):
                return {'status':'passed','config_sha256':config_hash,
                        'result':{'energy':energy,'minimum_verified':minimum,'week3_verification_passed':minimum}}
            rreport=put('refine-report.json',meta(rhash,-100.,True))
            breport=put('base-report.json',meta(bhash,-100.+.05/HARTREE_KCAL,False))
            item={'replacement_config':'refine-config.json','replacement_config_sha256':rhash,'common_grid_job_id':'base'}
            check=dict(status='passed',refinement_report='refine-report.json',refinement_report_sha256=rreport,
                common_grid_report='base-report.json',common_grid_report_sha256=breport,
                common_grid_config='base-config.json',common_grid_config_sha256=bhash,
                refinement_job_id='refine',refinement_config_sha256=rhash,common_grid_job_id='base',
                recovery_plan_sha256=planhash,same_geometry_verified=True,refined_minimum_verified=True,
                refined_energy_Hartree=-100.,common_grid_energy_Hartree=-100.+.05/HARTREE_KCAL,
                delta_energy_kcal_mol=.05,max_gradient_component_delta_Hartree_Bohr=1e-5)
            self.assertTrue(verify_numerical_check(check,item,root))
            check['max_gradient_component_delta_Hartree_Bohr']=2e-4
            with self.assertRaises(AssertionError):verify_numerical_check(check,item,root)
            check['status']='needs_review'
            self.assertFalse(verify_numerical_check(check,item,root))
            failed=meta(rhash,-100.,False);failed['status']='needs_review'
            check['refinement_report_sha256']=put('refine-report.json',failed)
            with self.assertRaises(AssertionError):verify_numerical_check(check,item,root)


if __name__=='__main__':unittest.main()
