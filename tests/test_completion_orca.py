"""Scientific parser checks against preserved real ORCA output, in temporary copies."""
import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import completion_orca as runner
from research_run import sha


class CompletionParser(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.run=Path(self.tmp.name)
        self.meta=json.loads((ROOT/'reports/w3-diflunisal--diflunisal-p05-refine-r1-opt-freq.json').read_text())
        old=Path(self.meta['run_dir']);stem=Path(self.meta['input']).stem
        names=['output.log',Path(self.meta['input']).name,stem+'.xyz',stem+'.hess',stem+'.engrad','initial.xyz','candidate.json']
        for name in names:shutil.copy2(old/name,self.run/name)
        shutil.copy2(ROOT/'scripts/completion_orca.py',self.run/'completion_parser_snapshot.py')
        self.meta['expected_input_sha256']=sha(self.run/Path(self.meta['input']).name)
        self.meta['completion_parser_sha256']=sha(self.run/'completion_parser_snapshot.py')
        self.meta['assets']={'fixture/initial.xyz':'initial.xyz','fixture/candidate.json':'candidate.json',
                             'fixture/parser.py':'completion_parser_snapshot.py'}
        self.meta['source_sha256']={self.meta['input']:self.meta['expected_input_sha256'],
            **{p:sha(self.run/d) for p,d in self.meta['assets'].items()}}
        self.stem=stem

    def test_real_refined_minimum_hessian_and_gradient(self):
        r=runner.checked_parse(self.run,self.meta)
        self.assertEqual(r['status'],'passed')
        self.assertTrue(r['scientific_verification_passed'])
        self.assertEqual(len(r['full_precision_frequencies_cm1']),150)
        self.assertAlmostEqual(min(r['full_precision_frequencies_cm1']),9.1844296317,places=7)

    def test_gradient_coordinate_corruption_cannot_pass(self):
        p=self.run/(self.stem+'.engrad')
        rows=p.read_text().splitlines()
        # The final block is atomic number followed by Cartesian coordinates.
        fields=rows[-1].split();fields[1]=str(float(fields[1])+.05);rows[-1]=' '.join(fields)
        p.write_text('\n'.join(rows)+'\n')
        r=runner.checked_parse(self.run,self.meta)
        self.assertEqual(r['status'],'needs_review')
        self.assertFalse(r['scientific_verification_passed'])

    def test_unexpected_input_is_rejected_even_if_engine_result_passed(self):
        p=self.run/Path(self.meta['input']).name;p.write_text(p.read_text()+'# changed\n')
        r=runner.checked_parse(self.run,self.meta)
        self.assertEqual(r['status'],'needs_review')

    def test_recorded_grid_matches_real_ground_state_output(self):
        self.meta['expected_scf_grid']={'IntAcc':4.959,'AngularGrid':6}
        r=runner.checked_parse(self.run,self.meta)
        self.assertEqual(r['status'],'passed')
        self.assertEqual(r['observed_mpi_process_counts'],[12])
        self.assertEqual(r['applied_scf_grid_verified'],self.meta['expected_scf_grid'])

    def test_wrong_grid_metadata_is_rejected(self):
        self.meta['expected_scf_grid']={'IntAcc':6.,'AngularGrid':7}
        self.assertEqual(runner.checked_parse(self.run,self.meta)['status'],'needs_review')

    def test_wrong_process_metadata_is_rejected(self):
        self.meta['mpi_processes']=16
        self.assertEqual(runner.checked_parse(self.run,self.meta)['status'],'needs_review')


if __name__=='__main__':unittest.main()
