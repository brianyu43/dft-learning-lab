"""Finite continuation of the already-authorized week-two pipeline; no recurrence.

Wait for the active three-job ORCA batch, run the reviewed 16-job QE batch, then
analyze existing outputs. No retries, no week-three jobs, no parameter adjustment.
"""
from datetime import datetime,timezone
import fcntl
import hashlib
import json
from pathlib import Path
import signal
import subprocess
import sys
import time
import psutil

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/week2-finish-controller';OUT.mkdir(exist_ok=True)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def save(state):
    p=OUT/'state.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(p)
def stop(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
signal.signal(signal.SIGTERM,stop)

with (OUT/'.lock').open('a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    upstream=ROOT/'experiments/week2/drug-first-dft-batch.json'
    expected_hash=sha(upstream);upstream_jobs=json.loads(upstream.read_text())['jobs']
    downstream=ROOT/'experiments/week2/bilayer-batch.json';downstream_hash=sha(downstream)
    assert len(upstream_jobs)==3 and len(json.loads(downstream.read_text())['jobs'])==16
    snapshot=OUT/'runner_snapshot.py'
    if snapshot.exists():assert sha(snapshot)==sha(Path(__file__)), 'Preserve/review the old controller before changing it'
    else:snapshot.write_bytes(Path(__file__).read_bytes())
    state={'status':'waiting_for_active_ORCA_batch','started_at':now(),'scope':'week2 only, finite pipeline',
           'pid':__import__('os').getpid(),'upstream_batch_sha256':expected_hash,'downstream_batch_sha256':downstream_hash,
           'runner_sha256':sha(Path(__file__)),'steps':[],'automatic_retries':0}
    save(state);proc=None;start=time.monotonic()
    try:
        while True:
            try:report=json.loads((ROOT/'reports/week2-drug-first-dft.json').read_text())
            except (FileNotFoundError,json.JSONDecodeError):report={'status':'running'}
            if 'batch_sha256' in report:assert report['batch_sha256']==expected_hash
            if report['status']=='completed':break
            if report['status']!='running':raise RuntimeError('Upstream ORCA batch requires review: '+report['status'])
            if time.monotonic()-start>21600:raise RuntimeError('Six-hour continuation waiting budget reached')
            time.sleep(5)
        for path in upstream_jobs:
            config=json.loads((ROOT/path).read_text());meta=json.loads((ROOT/f"reports/{config['job_id']}.json").read_text())
            assert meta['status']=='passed' and meta['config_sha256']==sha(ROOT/path)
        assert sha(downstream)==downstream_hash
        commands=[
            ['scripts/run_week2_batch.py','experiments/week2/bilayer-batch.json'],
            ['scripts/analyze_week2_monolayer.py'],['scripts/check_week2_combined.py'],
            ['scripts/analyze_week2_xtb.py'],['scripts/verify_week2_xtb.py'],
            ['scripts/verify_week2_drug.py'],['scripts/verify_week2_bilayer_inputs.py'],
            ['scripts/analyze_week2_bilayers.py'],['scripts/plot_week2.py'],
            ['scripts/summarize_week2.py'],['scripts/verify_week2_integrity.py']]
        for args in commands:
            state.update(status='running_finite_pipeline',current_command=args);save(state)
            print(json.dumps({'stage':'week2-continuation','command':args}),flush=True)
            with (OUT/'output.log').open('a') as output:
                proc=subprocess.Popen([sys.executable,*args],cwd=ROOT,stdout=output,stderr=subprocess.STDOUT,start_new_session=True)
                code=proc.wait()
            state['steps'].append({'command':args,'exit_code':code,'finished_at':now()});save(state);proc=None
            if code:raise RuntimeError(f'Command failed with {code}: {args}')
        state.update(status='finite_pipeline_finished',ended_at=now(),
                     scientific_status=json.loads((ROOT/'reports/week2-summary.json').read_text())['status'])
        save(state);print(json.dumps(state),flush=True)
    except BaseException as exc:
        if proc is not None and proc.poll() is None:
            try:
                parent=psutil.Process(proc.pid);owned=parent.children(recursive=True)+[parent]
                for process in reversed(owned):
                    try:process.terminate()
                    except psutil.Error:pass
                _,alive=psutil.wait_procs(owned,timeout=10)
                for process in alive:
                    try:process.kill()
                    except psutil.Error:pass
            except psutil.Error:pass
        state.update(status='stopped_requires_review',ended_at=now(),reason=f'{type(exc).__name__}: {exc}')
        save(state);raise
