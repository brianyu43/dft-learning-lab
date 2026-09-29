"""Finite reviewed job list, sequential, fail-stop, no retries; durable per-job outputs."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('batch',type=Path)
args=parser.parse_args()
batch=json.loads(args.batch.read_text())
status={'batch_id':batch['batch_id'],'started_at':datetime.now(timezone.utc).isoformat(),
        'batch_sha256':hashlib.sha256(args.batch.read_bytes()).hexdigest(),'status':'running','jobs':[],'automatic_retries':0}
path=ROOT/'reports'/f"{batch['batch_id']}.json"
def save():
    path.write_text(json.dumps(status,indent=2)+'\n')
save()
start=time.monotonic()
for cfg in batch['jobs']:
    if time.monotonic()-start>batch['budget_seconds']:
        status['status']='budget_review_required';break
    config=json.loads((ROOT/cfg).read_text())
    job=config['job_id']
    report=ROOT/f'reports/{job}.json'
    if report.exists():
        old=json.loads(report.read_text())
        if old['status']=='passed' and old['config_sha256']==hashlib.sha256((ROOT/cfg).read_bytes()).hexdigest():
            status['jobs'].append({'job_id':job,'status':'reused_existing_completed_run','run_dir':old['run_dir']});save();continue
        status.update(status='existing_attempt_requires_review',blocked_job=job);break
    run=subprocess.run([sys.executable,str(ROOT/'scripts/research_run.py'),'--config',str(ROOT/cfg)],cwd=ROOT)
    status['jobs'].append({'job_id':job,'exit_code':run.returncode,'report':str(report.relative_to(ROOT))})
    save()
    if run.returncode:
        status.update(status='job_failed_requires_review',blocked_job=job);break
else:
    status['status']='completed'
status.update(ended_at=datetime.now(timezone.utc).isoformat(),batch_elapsed_seconds=round(time.monotonic()-start,3))
save();print(json.dumps(status),flush=True)
raise SystemExit(0 if status['status']=='completed' else 1)
