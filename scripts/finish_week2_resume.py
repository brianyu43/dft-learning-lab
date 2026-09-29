"""Finite authorized resume: Freq checkpoint, 16 QE controls, verification/report.

No retries or week-three jobs. The original stopped controller remains archived.
"""
from datetime import datetime, timezone
import fcntl
import json
from pathlib import Path
import signal
import subprocess
import sys
import psutil
from research_run import sha

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'runs/week2-resume-20260911'
def now():return datetime.now(timezone.utc).isoformat()
def write(state):
    target=OUT/'state.json';tmp=target.with_suffix('.tmp');tmp.write_text(json.dumps(state,indent=2)+'\n');tmp.replace(target)
def interrupted(*_):raise KeyboardInterrupt('user termination')
signal.signal(signal.SIGTERM,interrupted)

with (OUT/'.controller.lock').open('a') as lock:
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    bench=json.loads((ROOT/'reports/week2-parallel-benchmark.json').read_text());assert bench['status']=='passed'
    batch=ROOT/'experiments/week2/bilayer-batch.json'
    assert len(json.loads(batch.read_text())['jobs'])==16
    snapshot=OUT/'resume_controller_snapshot.py'
    if snapshot.exists():assert sha(snapshot)==sha(Path(__file__)), 'Review an existing controller before retry'
    else:snapshot.write_bytes(Path(__file__).read_bytes())
    state=dict(status='running_finite_resume',started_at=now(),pid=__import__('os').getpid(),
               scope='Week 2 only; one engine at a time; no automatic retries',steps=[],
               runner_sha256=sha(Path(__file__)),bilayer_batch_sha256=sha(batch),
               parallel_benchmark_sha256=sha(ROOT/'reports/week2-parallel-benchmark.json'))
    write(state);proc=None
    commands=[['scripts/run_week2_frequency_resume.py','--execute'],
              ['scripts/verify_week2_drug.py'],
              ['scripts/run_week2_batch.py','experiments/week2/bilayer-batch.json'],
              ['scripts/analyze_week2_monolayer.py'],['scripts/check_week2_combined.py'],
              ['scripts/analyze_week2_xtb.py'],['scripts/verify_week2_xtb.py'],
              ['scripts/verify_week2_bilayer_inputs.py'],['scripts/analyze_week2_bilayers.py'],
              ['scripts/plot_week2.py'],['scripts/summarize_week2.py'],['scripts/verify_week2_integrity.py']]
    try:
        for args in commands:
            assert sha(batch)==state['bilayer_batch_sha256']
            state.update(current_command=args);write(state)
            print(json.dumps({'stage':'week2-resume','command':args}),flush=True)
            with (OUT/'continuation.log').open('a') as log:
                proc=subprocess.Popen([sys.executable,*args],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
                code=proc.wait()
            state['steps'].append(dict(command=args,exit_code=code,ended_at=now()));write(state);proc=None
            if code:raise RuntimeError(f'Command requires review: {args}, exit {code}')
        state.update(status='finite_pipeline_finished',ended_at=now(),
                     scientific_status=json.loads((ROOT/'reports/week2-summary.json').read_text())['status']);write(state)
    except BaseException as exc:
        if proc and proc.poll() is None:
            try:
                parent=psutil.Process(proc.pid);owned=parent.children(recursive=True)+[parent]
                for p in reversed(owned):
                    try:p.terminate()
                    except psutil.Error:pass
                _,alive=psutil.wait_procs(owned,timeout=10)
                for p in alive:
                    try:p.kill()
                    except psutil.Error:pass
            except psutil.Error:pass
        state.update(status='stopped_requires_review',ended_at=now(),reason=f'{type(exc).__name__}: {exc}');write(state);raise
