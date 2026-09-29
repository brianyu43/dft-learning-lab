"""One bounded scientific engine at a time, explicit batches, immutable run snapshots."""
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

from research_run import ROOT,sha,save,now


def spent():
    used={'orca':0.,'qe':0.}
    for p in (ROOT/'runs').glob('*/manifest.json'):
        m=json.loads(p.read_text())
        if m.get('week',0)>=3 and m.get('engine') in used:
            used[m['engine']]+=m.get('elapsed_seconds',0.)
    return used


def verify_report(config_path):
    cfg=json.loads((ROOT/config_path).read_text());p=ROOT/f"reports/{cfg['job_id']}.json"
    m=json.loads(p.read_text());run=Path(m['run_dir'])
    assert m['status']=='passed' and m['config_sha256']==sha(ROOT/config_path)
    assert m['output_sha256']==sha(run/'output.log')
    for source,digest in m['source_sha256'].items():
        destination=Path(m['input']).name if source==m['input'] else m['assets'][source]
        assert sha(run/destination)==digest
    if cfg['engine']=='orca':
        assert m['result']['scientific_verification_passed']
        if cfg['requested'].get('frequencies'):assert m['result']['minimum_verified']
    return m


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--plan',type=Path);group.add_argument('--materials',action='store_true')
    parser.add_argument('--execute',action='store_true');args=parser.parse_args()
    assert args.execute, 'Reviewed local execution requires --execute'
    out=ROOT/'runs'/('completion-'+('materials-' if args.materials else 'drug-')+now().replace(':','').replace('.','').replace('+','_'))
    out.mkdir()
    files=['completion_batch.py','research_run.py']
    if args.materials:files+=['week3_materials_pipeline.py','week3_qe.py','week2_qe_results.py','prepare_week3_drug.py']
    else:files+=['completion_orca.py','week3_recovery.py','prepare_week3_drug.py','run_week2_xtb_screen.py','run_week2_frequency_resume.py','orca_parallel.sh']
    code_hashes={f'scripts/{p}':sha(ROOT/'scripts'/p) for p in files}
    for p in files:shutil.copy2(ROOT/'scripts'/p,out/p)
    plan=json.loads(args.plan.read_text()) if args.plan else None
    if plan:
        shutil.copy2(args.plan,out/'plan.json')
        for p,digest in plan['input_hashes'].items():assert sha(ROOT/p)==digest
    previous=json.loads((ROOT/'reports/week2-summary.json').read_text())
    budgets={'orca':previous['remaining_budget_hours']['drug']*3600,'qe':previous['remaining_budget_hours']['semiconductor']*3600}
    state=dict(status='running',started_at=now(),pid=os.getpid(),controller_dir=str(out.relative_to(ROOT)),
               kind='materials' if args.materials else plan['kind'],stage='initializing',code_sha256=code_hashes,
               plan_sha256=sha(args.plan) if plan else None,steps=[],automatic_retries=0,
               remaining_after_week2_seconds=budgets)
    latest=ROOT/'reports/completion'/('materials-run-state.json' if args.materials else 'drug-run-state.json')
    latest.parent.mkdir(exist_ok=True)
    proc=None

    def status(stage=None):
        if stage:state['stage']=stage
        state.update(updated_at=now(),spent_after_week2_seconds=spent())
        save(out/'state.json',state);save(latest,state)

    def frozen():
        for p,digest in code_hashes.items():assert sha(ROOT/p)==digest,f'Changed running dependency: {p}'
        if plan:
            assert sha(args.plan)==state['plan_sha256']
            for p,digest in plan['input_hashes'].items():assert sha(ROOT/p)==digest

    def run_config(config_path):
        nonlocal proc
        frozen();cfg=json.loads((ROOT/config_path).read_text())
        report=ROOT/f"reports/{cfg['job_id']}.json"
        if report.exists():return verify_report(config_path)
        assert not list((ROOT/'runs').glob(f"*-{cfg['job_id']}")), 'Inspect an existing attempt before launching a duplicate'
        assert cfg['rss_limit_gib']<=40 and cfg['timeout_seconds']<=14400
        assert budgets[cfg['engine']]-spent()[cfg['engine']]>=cfg['timeout_seconds'], 'Budget review required'
        command=[sys.executable,'scripts/completion_orca.py' if cfg['engine']=='orca' else 'scripts/research_run.py','--config',config_path]
        state.update(current_job=cfg['job_id'],current_config=config_path,current_config_sha256=sha(ROOT/config_path));status()
        with (out/'controller.log').open('a') as log:
            log.write(json.dumps({'started_at':now(),'command':command})+'\n');log.flush()
            proc=subprocess.Popen(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            state.update(child_pid=proc.pid);status()
            while proc.poll() is None:
                time.sleep(2)
                if time.time()-datetime_timestamp(state['updated_at'])>20:status()
            code=proc.returncode;proc=None
        state['steps'].append(dict(job_id=cfg['job_id'],config=config_path,exit_code=code,ended_at=now()));status()
        assert code==0,f"Review failed job: {cfg['job_id']}"
        return verify_report(config_path)

    signal.signal(signal.SIGTERM,lambda *_:(_ for _ in ()).throw(KeyboardInterrupt('terminated')))
    with (ROOT/'experiments/.completion-controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            status()
            if args.materials:
                from week3_materials_pipeline import execute
                execute(run_config,status)
            else:
                status(plan['kind'])
                for cfg in plan['jobs']:run_config(cfg)
            state.update(status='completed',ended_at=now(),current_job=None,child_pid=None);status('completed')
        except BaseException as exc:
            if proc is not None and proc.poll() is None:
                try:
                    parent=psutil.Process(proc.pid);members=parent.children(recursive=True)+[parent]
                    for p in members:
                        try:p.terminate()
                        except psutil.Error:pass
                    _,alive=psutil.wait_procs(members,timeout=10)
                    for p in alive:
                        try:p.kill()
                        except psutil.Error:pass
                except psutil.Error:pass
            state.update(status='stopped_requires_review',ended_at=now(),reason=f'{type(exc).__name__}: {exc}');status()
            raise


def datetime_timestamp(value):
    from datetime import datetime
    return datetime.fromisoformat(value).timestamp()


if __name__=='__main__':main()
