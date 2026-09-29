"""Solvent input assembly and prerequisite failures; no new scientific runs.

Temporary fixtures deliberately bypass scientific prerequisites to exercise
serialization. They use historical start coordinates, not completed uniform
minima, and cannot be used as evidence that a production method panel exists.
"""
from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import prepare_completion_drug_solvent as prep
from research_run import sha


class SolventPreparationTests(unittest.TestCase):
    def fixture(self, root):
        lineage = json.loads((ROOT / 'experiments/completion/drug-start-lineage-r1.json').read_text())
        files = [*prep.DEPENDENCIES, prep.SCORES, prep.REVIEW,
            'scripts/completion_orca.py', 'scripts/week3_recovery.py',
            'scripts/run_week2_frequency_resume.py', 'scripts/run_week2_xtb_screen.py',
            'scripts/orca_parallel.sh', 'environment/orca-mpi-explicit.txt',
            *(f'.tools/orca-mpi/lib/{lib}' for lib in prep.LIBS)]
        for rel in files:
            path = root / rel; path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('TEST FIXTURE ONLY: not executable and not scientific evidence\n')
        lp = root / 'experiments/completion/drug-start-lineage-r1.json'
        lp.parent.mkdir(parents=True, exist_ok=True); lp.write_text(json.dumps(lineage))
        rows = []
        for original in lineage['rows']:
            source = root / 'fixture-sources' / (original['candidate_id'] + '.xyz')
            source.parent.mkdir(exist_ok=True)
            source.write_bytes((ROOT / original['source_geometry']).read_bytes())
            rows.append(dict(mapping=dict(candidate_id=original['candidate_id'],
                species=original['species'], selection_role=original['selection_role'],
                config='fixture-parent-' + original['candidate_id']),
                result=dict(candidate=original['candidate'], geometry_path=source,
                            geometry_sha256=sha(source), meta=dict(job_id='fixture-only-no-engine'))))
        return rows

    def review(self):
        actual = dict(source_sha256={'fixture-output.log': 'b' * 64},
            method_evidence=dict(method='r2SCAN-3c', orbital_basis='def2-mTZVPP'),
            solvent='SMD(ethanol)', energy_quantity='orca_smd_total_energy_including_cds',
            smd_energy_evidence=dict(cds_added_again=False), internal_frequencies_cm1=[1000., 2000., 3000.])
        review = dict(status='smdethanol_runtime_reviewed', compatibility_plan_sha256='plan-fixture',
            config_path=prep.COMPAT_CONFIG, config_sha256='config-fixture',
            source_sha256=actual['source_sha256'], water_only_not_drug_accuracy_or_cost=True,
            checks={topic: dict(accepted=True, note='Synthetic review for a unit check, not ORCA evidence',
                               runtime_quotes=['synthetic runtime evidence']) for topic in prep.TOPICS})
        log = 'synthetic credits only\n****END OF INPUT****\nsynthetic runtime evidence\n'
        return review, actual, log

    def validate(self, review, actual, log):
        prep.validate_runtime_review(review, actual, log, 'plan-fixture', 'config-fixture')

    def test_missing_scientific_prerequisites_write_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with patch.object(prep, 'ROOT', root), patch.object(prep.panel, 'evaluate') as evaluate:
                with self.assertRaisesRegex(ValueError, 'no inputs written'):
                    prep.prepare()
                evaluate.assert_not_called()
            self.assertEqual(list(root.iterdir()), [])

    def test_runtime_review_requires_actual_hashes_topics_and_runtime_quotes(self):
        review, actual, log = self.review(); self.validate(review, actual, log)
        for mutate, error in [
            (lambda r: r.update(source_sha256={}), 'raw evidence differs'),
            (lambda r: r['checks'].pop('ethanol_smd_parameters'), 'detailed review'),
            (lambda r: r['checks']['composite_method'].update(runtime_quotes=['synthetic credits only']),
             'actual runtime quotation'),
            (lambda r: r.update(water_only_not_drug_accuracy_or_cost=False), 'scope'),
        ]:
            changed = deepcopy(review); mutate(changed)
            with self.subTest(error=error), self.assertRaisesRegex(ValueError, error):
                self.validate(changed, actual, log)

    def test_pass_flags_do_not_bypass_actual_compatibility_revalidation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for rel, value in [(prep.SCORES, dict(status='uniform_vacuum_numerical_panel_passed')),
                               (prep.REVIEW, dict(status='smdethanol_runtime_reviewed')),
                               (prep.COMPAT_PLAN, dict(input_hashes={}, jobs=[prep.COMPAT_CONFIG]))]:
                path = root / rel; path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(value))
            with patch.object(prep, 'ROOT', root), patch.object(prep.compatibility, 'audit_job',
                    side_effect=ValueError('Raw fixture output was rejected')) as audit:
                with self.assertRaisesRegex(ValueError, 'Raw fixture output was rejected'): prep.prepare()
                audit.assert_called_once_with(prep.COMPAT_CONFIG)
            self.assertFalse((root / 'inputs').exists())

    def test_wrong_environment_energy_quantity_or_hessian_cannot_enable_smd(self):
        review, actual, log = self.review()
        for key, value in [('method_evidence', dict(method='wB97X-3c', orbital_basis='vDZP')),
                           ('solvent', 'vacuum'), ('energy_quantity', 'electronic_energy'),
                           ('smd_energy_evidence', None), ('smd_energy_evidence', dict(cds_added_again=True)),
                           ('internal_frequencies_cm1', [-1., 1000., 2000.]),
                           ('internal_frequencies_cm1', [1000., 2000.])]:
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                self.validate(review, {**actual, key: value}, log)

    def test_all_starts_both_grids_same_coordinates_and_resource_caps(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); rows = self.fixture(root)
            original = {str(r['result']['geometry_path']): sha(r['result']['geometry_path']) for r in rows}
            with patch.object(prep, 'ROOT', root), patch.object(prep, 'prerequisites', return_value=(rows, {})):
                plan_path = prep.prepare(); plan = json.loads(plan_path.read_text())
            self.assertEqual(len(plan['jobs']), 34)
            self.assertEqual(plan['new_independent_start_count'], 0)
            self.assertEqual(len({r['species'] for r in plan['comparison_mapping']}), 14)
            self.assertFalse(plan['solvent_sensitivity_complete'])
            self.assertFalse(plan['automatic_execution'])
            self.assertFalse(plan['solvent_geometry_optimization_enabled'])
            by_id = {}
            for path in plan['jobs']:
                cfg = json.loads((root / path).read_text())
                inp = (root / cfg['input']).read_text()
                self.assertEqual(inp, prep.input_text(cfg['grid_setting']))
                self.assertEqual(cfg['method'], 'r2SCAN-3c')
                self.assertEqual(cfg['solvent'], 'SMD(ethanol)')
                self.assertEqual(cfg['energy_quantity'], 'orca_smd_total_energy_including_cds')
                self.assertEqual(inp.count('SMD(ethanol)'), 1)
                self.assertIn('SMD18 false', inp)
                self.assertIn('draco false', inp)
                self.assertEqual(cfg['requested'], dict(optimization=False, frequencies=False, gradient=True))
                self.assertEqual((cfg['mpi_processes'], cfg['omp_threads'], cfg['timeout_seconds'], cfg['rss_limit_gib']),
                                 (12, 1, 1800, 40))
                self.assertFalse(cfg['verifies_minimum_in_this_solvent'])
                xyz = root / next(p for p, dst in cfg['assets'].items() if dst == 'initial.xyz')
                self.assertEqual(sha(xyz), cfg['source_geometry_sha256'])
                by_id.setdefault(cfg['candidate_id'], set()).add(sha(xyz))
                if cfg['grid_setting'] == 'strict':
                    self.assertIn('GridPruning Unpruned', inp)
                    self.assertIn('HGridReduced false', inp)
            self.assertEqual(len(by_id), 17)
            self.assertTrue(all(len(hashes) == 1 for hashes in by_id.values()))
            self.assertEqual(original, {p: sha(Path(p)) for p in original})

    def test_missing_start_or_changed_role_cannot_write_production_plan(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); rows = self.fixture(root)
            changed = deepcopy(rows); changed[-1]['mapping']['selection_role'] = 'invented'
            for altered, error in [(rows[:-1], '17 starts'), (changed, 'selection role')]:
                with patch.object(prep, 'ROOT', root), patch.object(prep, 'prerequisites', return_value=(altered, {})):
                    with self.assertRaisesRegex(ValueError, error): prep.prepare()
                self.assertFalse((root / prep.PLAN).exists())
                self.assertFalse((root / 'inputs').exists())

    def test_changed_source_geometry_is_rejected_before_writes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); rows = self.fixture(root)
            rows[-1]['result']['geometry_path'].write_text('corrupted XYZ\n')
            with patch.object(prep, 'ROOT', root), patch.object(prep, 'prerequisites', return_value=(rows, {})):
                with self.assertRaisesRegex(ValueError, 'coordinate identity changed'): prep.prepare()
            self.assertFalse((root / 'inputs').exists())

    def test_saved_plan_is_idempotent_and_tampered_dependency_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); rows = self.fixture(root)
            with patch.object(prep, 'ROOT', root), patch.object(prep, 'prerequisites', return_value=(rows, {})) as gate:
                target = prep.prepare(); first = sha(target)
                self.assertEqual(prep.prepare(), target); self.assertEqual(sha(target), first)
                self.assertEqual(gate.call_count, 1)
                saved_text = target.read_text(); changed = json.loads(saved_text)
                changed['comparison_mapping'].pop(); target.write_text(json.dumps(changed))
                with self.assertRaisesRegex(ValueError, 'Saved job coverage changed'): prep.prepare()
                target.write_text(saved_text)
                (root / 'scripts/completion_orca.py').write_text('changed test fixture')
                with self.assertRaisesRegex(ValueError, 'Changed or missing evidence'): prep.prepare()


if __name__ == '__main__':
    unittest.main()
