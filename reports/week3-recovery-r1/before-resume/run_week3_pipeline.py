"""Authorized finite week-three calculation pipeline; one engine, budgeted, no retries."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

import psutil

from research_run import ROOT, now, save, sha
from summarize_week3 import summarize
from week3_materials_pipeline import execute as run_materials
from plot_week3 import main as make_figures


def spent(engine):
    total=0.0
    for p in (ROOT/'runs').glob('*-w3-*/manifest.json'):
        m=json.loads(p.read_text())
        if m.get('engine')==engine:
            total+=m.get('elapsed_seconds',0.0)
    return total


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--execute',action='store_true');args=parser.parse_args()
    if not args.execute:raise SystemExit('Use --execute only within the authorized week-three scope')
    plan_path=ROOT/'experiments/week3/drug/panel-plan.json';plan=json.loads(plan_path.read_text())
    original=json.loads((ROOT/'reports/week2-summary.json').read_text())
    remaining={'orca':original['remaining_budget_hours']['drug']*3600,
               'qe':original['remaining_budget_hours']['semiconductor']*3600}
    for path,digest in plan['input_hashes'].items():assert sha(ROOT/path)==digest
    assert plan['new_start_count_including_first']==14
    out=ROOT/'runs'/('week3-pipeline-'+now().replace(':','').replace('.','').replace('+','_'))
    out.mkdir()
    code_files=['run_week3_pipeline.py','week3_materials_pipeline.py','week3_qe.py','summarize_week3.py',
                'run_week3_drug.py','run_week3_drug_v2.py','prepare_week3_drug.py','research_run.py',
                'run_week2_xtb_screen.py','run_week2_frequency_resume.py','week2_qe_results.py','orca_parallel.sh','plot_week3.py']
    code_hashes={f'scripts/{p}':sha(ROOT/'scripts'/p) for p in code_files}
    for p in code_files:shutil.copy2(ROOT/'scripts'/p,out/p)
    shutil.copy2(plan_path,out/'panel-plan.json')
    state=dict(status='waiting_for_first_full_job',started_at=now(),pid=os.getpid(),controller_dir=str(out.relative_to(ROOT)),
               plan_sha256=sha(plan_path),code_sha256=code_hashes,steps=[],automatic_retries=0,
               concurrent_engines=1,remaining_after_week2_seconds=remaining,
               scope='Week-three bounded 14 new ORCA starts and adaptive MoS2/MoSe2 reference; no week-four jobs')
    proc=None

    def status(stage_name):
        state.update(stage=stage_name,updated_at=now(),spent_week3_seconds={e:spent(e) for e in remaining})
        save(out/'state.json',state);save(ROOT/'reports/week3-run-state.json',state)
        summarize()

    def check_code():
        for p,h in code_hashes.items():assert sha(ROOT/p)==h, f'Frozen pipeline dependency changed: {p}'

    def run_config(config_path):
        nonlocal proc
        check_code();cfg=json.loads((ROOT/config_path).read_text());report=ROOT/f"reports/{cfg['job_id']}.json"
        if report.exists():
            m=json.loads(report.read_text())
            assert m['status']=='passed' and m['config_sha256']==sha(ROOT/config_path),'Existing attempt requires review'
            return m
        budget=remaining[cfg['engine']]-spent(cfg['engine'])
        assert budget>=cfg['timeout_seconds'],'Insufficient remaining budget for the bounded job timeout'
        assert cfg['week']==3 and cfg['rss_limit_gib']<=40
        command=[sys.executable,'scripts/run_week3_drug_v2.py','--config',config_path] if cfg['engine']=='orca' else [sys.executable,'scripts/research_run.py','--config',config_path]
        state.update(current_job=cfg['job_id'],current_config=config_path,current_config_sha256=sha(ROOT/config_path))
        status(state.get('stage','running'))
        with (out/'pipeline.log').open('a') as log:
            log.write(json.dumps({'started_at':now(),'command':command})+'\n');log.flush()
            proc=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            code=proc.wait();proc=None
        state['steps'].append(dict(job_id=cfg['job_id'],config=config_path,exit_code=code,ended_at=now()))
        status(state['stage'])
        assert code==0,f"Job requires review: {cfg['job_id']}"
        m=json.loads(report.read_text());assert m['status']=='passed'
        if cfg['engine']=='orca':assert m['result']['week3_verification_passed']
        return m

    def interrupted(*_):raise KeyboardInterrupt('user termination')
    signal.signal(signal.SIGTERM,interrupted)
    with (ROOT/'experiments/.week3-controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            status('first_drug_full_job_cost')
            first_cfg=json.loads((ROOT/plan['first_job']).read_text());first_report=ROOT/f"reports/{first_cfg['job_id']}.json"
            deadline=time.monotonic()+first_cfg['timeout_seconds']+300
            while not first_report.exists():
                if time.monotonic()>deadline:raise RuntimeError('First-job completion report not available within bounded wait')
                time.sleep(10)
            check_code();first=json.loads(first_report.read_text())
            assert first['config_sha256']==sha(ROOT/plan['first_job'])
            assert first['status']=='passed' and first['result']['week3_verification_passed'] and first['result']['minimum_verified']
            # A planning scenario, not a measured full-panel prediction or guaranteed upper bound.
            scenario=first['elapsed_seconds']*(1+2*len(plan['remaining_jobs']))
            review=dict(reviewed_at=now(),first_job=first_cfg['job_id'],first_elapsed_seconds=first['elapsed_seconds'],
                        first_peak_sampled_process_tree_rss_bytes=first['peak_sampled_process_tree_rss_bytes'],
                        first_complete_opt_freq_verified=True,remaining_starts=len(plan['remaining_jobs']),
                        planning_scenario_seconds=scenario,scenario='First measured full job for pilot plus twice that cost per remaining start; molecular cost can differ',
                        remaining_after_week2_seconds=remaining['orca'],passes=scenario<=remaining['orca'],
                        first_report_sha256=sha(first_report),panel_plan_sha256=sha(plan_path))
            save(out/'first-cost-review.json',review);save(ROOT/'reports/week3-first-cost-review.json',review)
            assert review['passes'],'First measured cost exceeds the precommitted panel budget scenario'
            state['status']='running';status('drug_contact_and_self_association_controls')
            for path in plan['remaining_jobs']:run_config(path)
            drug=summarize();assert drug['drug_complete_new_jobs']==14
            make_figures()
            status('semiconductor_reference');run_materials(run_config,status)
            make_figures()
            state.update(status='completed',ended_at=now(),current_job=None,figure_visual_review='pending human or agent inspection');status('completed')
        except BaseException as exc:
            if proc is not None and proc.poll() is None:
                try:
                    owner=psutil.Process(proc.pid);members=owner.children(recursive=True)+[owner]
                    for p in reversed(members):
                        try:p.terminate()
                        except psutil.Error:pass
                    _,alive=psutil.wait_procs(members,timeout=10)
                    for p in alive:
                        try:p.kill()
                        except psutil.Error:pass
                except psutil.Error:pass
            state.update(status='stopped_requires_review',ended_at=now(),reason=f'{type(exc).__name__}: {exc}')
            material=ROOT/'reports/week3-materials.json'
            if material.exists():
                m=json.loads(material.read_text())
                if m['status']=='running':m.update(status='requires_review',reason=state['reason']);save(material,m)
            status('stopped_requires_review')
            raise
    return 0


if __name__=='__main__':raise SystemExit(main())
