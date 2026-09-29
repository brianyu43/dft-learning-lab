"""Real archived outputs, temporary metadata, and deliberate corruption; no engines."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_week3_drug_v2 as drug
import week3_qe as qe
from research_run import sha


class DrugVerification(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.run=Path(self.tmp.name)
        self.meta=json.loads((ROOT/'reports/week1-water-opt-freq.json').read_text())
        old=Path(self.meta['run_dir']);stem=Path(self.meta['input']).stem
        for name in ['output.log',Path(self.meta['input']).name,stem+'.xyz',stem+'.hess']:
            shutil.copy2(old/name,self.run/name)
        shutil.copy2(old/(stem+'.xyz'),self.run/'initial.xyz')
        shutil.copy2(ROOT/'scripts/run_week3_drug_v2.py',self.run/'week3_parser_snapshot.py')
        candidate={'kind':'monomer','graph':{'atomic_numbers':[8,1,1],
                    'bonds':[[0,1,1],[0,2,1]],'hydrogen_owners':{'1':0,'2':0}}}
        (self.run/'candidate.json').write_text(json.dumps(candidate))
        self.meta.update(candidate_sha256=sha(self.run/'candidate.json'),
                         initial_structure_sha256=sha(self.run/'initial.xyz'),
                         week3_parser_sha256=sha(self.run/'week3_parser_snapshot.py'))

    def test_real_water_monomer_has_no_fragment_separation_field(self):
        result=drug.checked_parse(self.run,self.meta)
        self.assertEqual(result['status'],'passed')
        self.assertTrue(result['hessian_geometry_matches'])
        self.assertNotIn('separated_geometry',result['geometry_audit'])

    def test_changed_candidate_metadata_rejected(self):
        with (self.run/'candidate.json').open('a') as f:f.write(' ')
        result=drug.checked_parse(self.run,self.meta)
        self.assertEqual(result['status'],'needs_review')

    def test_mismatched_hessian_geometry_rejected(self):
        from ase.io import read,write
        p=self.run/(Path(self.meta['input']).stem+'.xyz')
        atoms=read(p);atoms.positions[1,0]+=.02;write(p,atoms)
        result=drug.checked_parse(self.run,self.meta)
        self.assertEqual(result['status'],'needs_review')
        self.assertFalse(result['minimum_verified'])


class QEVerification(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);(self.root/'reports').mkdir()
        self.job='w2-mono-combined-c80-k6-v20'
        self.meta=json.loads((ROOT/f'reports/{self.job}.json').read_text())
        self.meta['calculation']='scf'
        self.path=self.root/f'reports/{self.job}.json'
        self.path.write_text(json.dumps(self.meta))

    def test_real_qe_xml_energy_force_and_pseudo_electrons(self):
        with patch.object(qe,'ROOT',self.root):r=qe.load(self.job)
        self.assertEqual(r['nelec'],26)
        self.assertAlmostEqual(r['energy_eV'],self.meta['result']['energy']*qe.HARTREE_EV/2,places=6)
        self.assertEqual(len(r['forces']),3)

    def test_wrong_charge_metadata_rejected(self):
        self.meta['charge']=1;self.path.write_text(json.dumps(self.meta))
        with patch.object(qe,'ROOT',self.root),self.assertRaises(AssertionError):qe.load(self.job)


if __name__=='__main__':unittest.main()
