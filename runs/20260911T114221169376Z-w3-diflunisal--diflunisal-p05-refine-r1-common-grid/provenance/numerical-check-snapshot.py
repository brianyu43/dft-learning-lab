"""Frozen, explicit replacement mapping and a same-geometry numerical check."""
import json
from pathlib import Path

import numpy as np
from ase.io import read

from prepare_week3_drug import immutable
from research_run import ROOT, sha, save, now
from run_week2_frequency_resume import BOHR_ANGSTROM

PLAN = 'experiments/week3/drug/recovery-plan-r1.json'
HARTREE_KCAL = 627.509474


def load_recovery(root=ROOT):
    path = root / PLAN
    if not path.exists():
        return None, {}
    plan = json.loads(path.read_text())
    assert plan['schema'] == 1 and plan['automatic_retries'] == 0
    assert plan['reviewed_refinement_attempts'] == 1 and plan['new_independent_start_count'] == 0
    assert plan['parent_panel_plan_sha256'] == sha(root/'experiments/week3/drug/panel-plan.json')
    assert plan['inspection_sha256'] == sha(root/plan['inspection'])
    inspection = json.loads((root/plan['inspection']).read_text())
    replacements = {}
    assert len(plan['replacements']) == 1
    for item in plan['replacements']:
        for field in ['parent_config', 'parent_report', 'replacement_config']:
            assert sha(root/item[field]) == item[field+'_sha256'], field
        parent = json.loads((root/item['parent_config']).read_text())
        retry = json.loads((root/item['replacement_config']).read_text())
        failed = json.loads((root/item['parent_report']).read_text())
        assert failed['status'] == 'needs_review' and not failed['result']['minimum_verified']
        assert failed['config_sha256'] == item['parent_config_sha256']
        assert inspection['parent_report_sha256'] == item['parent_report_sha256']
        assert retry['parent_job_id'] == parent['job_id'] and retry['refinement_attempt'] == 1
        for key in ['candidate_id','species','role','method','charge','multiplicity','solvent','n_atoms','candidate_sha256']:
            assert parent[key] == retry[key], key
        assert retry['method'] == 'r2SCAN-3c' and retry['mpi_processes'] == 12
        assert retry['omp_threads'] == 1 and retry['rss_limit_gib'] <= 40 and retry['timeout_seconds'] <= 14400
        assert retry['initial_structure_sha256'] == inspection['geometry_sha256']
        assert (root/retry['input']).read_text() == '! r2SCAN-3c VeryTightSCF VeryTightOpt DefGrid3 Opt Freq\n%maxcore 2000\n%pal nprocs 12 end\n\n* xyzfile 0 1 initial.xyz\n'
        initial=next(source for source,target in retry['assets'].items() if target=='initial.xyz')
        assert sha(root/initial)==retry['initial_structure_sha256']
        replacements[item['parent_config']] = item
    return plan, replacements


def checked_report(config_path):
    cfg = json.loads((ROOT/config_path).read_text())
    path = ROOT/f"reports/{cfg['job_id']}.json"
    meta = json.loads(path.read_text())
    assert meta['status'] == 'passed' and meta['config_sha256'] == sha(ROOT/config_path)
    run = Path(meta['run_dir'])
    assert meta['output_sha256'] == sha(run/'output.log')
    for source, digest in meta['source_sha256'].items():
        target = Path(meta['input']).name if source == meta['input'] else meta['assets'][source]
        assert sha(run/target) == digest
    return cfg, meta, path


def gradient_in_frame(meta, reference):
    run = Path(meta['run_dir'])
    path = run/(Path(meta['input']).stem+'.engrad')
    values = [s.strip() for s in path.read_text().splitlines() if s.strip() and not s.lstrip().startswith('#')]
    n = int(values[0]); energy = float(values[1])
    gradient = np.array([float(v) for v in values[2:2+3*n]]).reshape(n,3)
    coords = np.array([list(map(float, v.split())) for v in values[2+3*n:]])
    assert n == len(reference) and coords.shape == (n,4)
    assert coords[:,0].astype(int).tolist() == reference.numbers.tolist()
    assert np.isfinite(gradient).all() and np.isfinite(coords).all()
    assert abs(energy-meta['result']['energy']) < 1e-9
    p = coords[:,1:] * BOHR_ANGSTROM
    p -= p.mean(axis=0)
    q = reference.positions-reference.positions.mean(axis=0)
    u,_,vt = np.linalg.svd(p.T@q)
    rotation = u @ np.diag([1,1,np.linalg.det(u@vt)]) @ vt
    np.testing.assert_allclose(p@rotation,q,rtol=0,atol=1e-7)
    return energy, gradient@rotation, str(path), sha(path)


def verify_numerical_check(check, item, root=ROOT):
    for field in ['refinement_report','common_grid_report','common_grid_config']:
        assert sha(root/check[field]) == check[field+'_sha256']
    retry = json.loads((root/item['replacement_config']).read_text())
    assert check['refinement_job_id'] == retry['job_id']
    assert check['common_grid_job_id'] == item['common_grid_job_id']
    assert check['refinement_config_sha256'] == item['replacement_config_sha256']
    assert check['recovery_plan_sha256'] == sha(root/PLAN)
    assert check['same_geometry_verified'] and check['refined_minimum_verified']
    refined=json.loads((root/check['refinement_report']).read_text())
    base=json.loads((root/check['common_grid_report']).read_text())
    assert refined['status']==base['status']=='passed'
    assert refined['result']['minimum_verified'] and refined['result']['week3_verification_passed']
    assert refined['config_sha256']==item['replacement_config_sha256']
    assert base['config_sha256']==check['common_grid_config_sha256']
    assert abs(check['refined_energy_Hartree']-refined['result']['energy'])<1e-9
    assert abs(check['common_grid_energy_Hartree']-base['result']['energy'])<1e-9
    assert abs(check['delta_energy_kcal_mol']-(check['common_grid_energy_Hartree']-check['refined_energy_Hartree'])*HARTREE_KCAL)<1e-8
    assert check['max_gradient_component_delta_Hartree_Bohr']>=0
    policy = json.loads((root/PLAN).read_text())['common_grid_policy']
    expected = (abs(check['delta_energy_kcal_mol']) <= policy['max_absolute_energy_delta_kcal_mol']
                and check['max_gradient_component_delta_Hartree_Bohr'] <= policy['max_gradient_component_delta_Hartree_Bohr'])
    assert (check['status'] == 'passed') == expected
    return expected


def complete_numerical_check(item, run_config):
    plan, _ = load_recovery()
    cfg, meta, report = checked_report(item['replacement_config'])
    assert meta['result']['week3_verification_passed'] and meta['result']['minimum_verified']
    path = ROOT/item['common_grid_check_report']
    if path.exists():
        check = json.loads(path.read_text())
        assert verify_numerical_check(check,item), 'Existing numerical check requires review'
        return check
    final = Path(meta['run_dir'])/(Path(meta['input']).stem+'.xyz')
    assert sha(final) == meta['result']['optimized_xyz_sha256']
    job = item['common_grid_job_id']
    initial = ROOT/f'structures/week3/drug_initial/{job}.xyz'
    immutable(initial,final.read_text())
    inp = ROOT/f'inputs/orca/week3/{job}.inp'
    immutable(inp,'! r2SCAN-3c TightSCF EnGrad\n%maxcore 2000\n%pal nprocs 12 end\n\n* xyzfile 0 1 initial.xyz\n')
    assets = {k:v for k,v in cfg['assets'].items() if v.endswith('.dylib')}
    assets.update({str(initial.relative_to(ROOT)):'initial.xyz',
                   str(report.relative_to(ROOT)):'provenance/refinement-report.json',
                   'scripts/week3_recovery.py':'provenance/numerical-check-snapshot.py',
                   PLAN:'provenance/recovery-plan-r1.json'})
    single = dict(protocol_id=cfg['protocol_id'],week=3,task_id='DR03',job_id=job,engine='orca',
                  executable=cfg['executable'],input=str(inp.relative_to(ROOT)),
                  method='r2SCAN-3c',charge=0,multiplicity=1,solvent='vacuum',n_atoms=52,linear=False,
                  requested={'optimization':False,'frequencies':False,'gradient':True},
                  mpi_processes=12,omp_threads=1,timeout_seconds=900,rss_limit_gib=40,assets=assets,
                  initial_structure_sha256=sha(initial),refinement_job_id=cfg['job_id'],
                  purpose='Original numerical settings on the verified refined geometry; no new pose or independent minimum')
    config_path = ROOT/f'experiments/week3/drug/{job}.json'
    immutable(config_path,json.dumps(single,indent=2)+'\n')
    run_config(str(config_path.relative_to(ROOT)))
    _, base, base_report = checked_report(str(config_path.relative_to(ROOT)))
    reference = read(final)
    e_ref,g_ref,p_ref,h_ref = gradient_in_frame(meta,reference)
    e_base,g_base,p_base,h_base = gradient_in_frame(base,reference)
    delta = (e_base-e_ref)*HARTREE_KCAL
    grad_delta = float(np.max(np.abs(g_base-g_ref)))
    policy = plan['common_grid_policy']
    passed = abs(delta)<=policy['max_absolute_energy_delta_kcal_mol'] and grad_delta<=policy['max_gradient_component_delta_Hartree_Bohr']
    check = dict(checked_at=now(),status='passed' if passed else 'needs_review',
        refinement_job_id=cfg['job_id'],refinement_config_sha256=sha(ROOT/item['replacement_config']),
        refinement_report=str(report.relative_to(ROOT)),refinement_report_sha256=sha(report),
        common_grid_job_id=job,common_grid_config=str(config_path.relative_to(ROOT)),common_grid_config_sha256=sha(config_path),
        common_grid_report=str(base_report.relative_to(ROOT)),common_grid_report_sha256=sha(base_report),
        recovery_plan_sha256=sha(ROOT/PLAN),same_geometry_verified=True,refined_minimum_verified=True,
        geometry_sha256=sha(final),refined_energy_Hartree=e_ref,common_grid_energy_Hartree=e_base,
        delta_energy_kcal_mol=delta,max_gradient_component_delta_Hartree_Bohr=grad_delta,
        refined_engrad=p_ref,refined_engrad_sha256=h_ref,common_grid_engrad=p_base,common_grid_engrad_sha256=h_base,
        energy_used_for_comparison='common_grid_energy_Hartree',policy=policy)
    save(path,check)
    assert verify_numerical_check(check,item), 'Grid/SCF numerical sensitivity exceeds the declared limit'
    return check
