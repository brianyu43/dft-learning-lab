"""Read-only status of finite week-two engine jobs and their current stage."""
from datetime import datetime,timezone
import json
from pathlib import Path
import re

ROOT=Path(__file__).resolve().parents[1]
jobs=[]
for path in sorted([*(ROOT/'runs').glob('*-w2-*/manifest.json'),*(ROOT/'runs').glob('*-p2-*/manifest.json')]):
    try:jobs.append(json.loads(path.read_text()))
    except json.JSONDecodeError:continue
active=[]
for job in jobs:
    if job['status']!='running':continue
    run=Path(job['run_dir']);log=(run/'output.log').read_text(errors='replace') if (run/'output.log').exists() else ''
    row={'job_id':job['job_id'],'run_dir':str(run.relative_to(ROOT)),'engine':job['engine']}
    try:progress=json.loads((run/'progress.json').read_text())
    except (FileNotFoundError,json.JSONDecodeError):progress={}
    row.update(elapsed_minutes=round(progress.get('elapsed_seconds',0)/60,2),
               peak_sampled_RSS_GiB=round(progress.get('peak_sampled_process_tree_rss_bytes',0)/1024**3,3))
    if job['engine']=='orca':
        markers=[('VIBRATIONAL FREQUENCIES','frequency_output_and_finalization'),
                 ('SCF HESSIAN','Hessian_assembly'),('ORCA SCF RESPONSE CALCULATION','electronic_response'),
                 ('GEOMETRIC PERTURBATIONS','frequency_derivative_integrals'),
                 ('THE OPTIMIZATION HAS CONVERGED','post_optimization_integrals')]
        row['stage']=next((stage for marker,stage in markers if marker in log),'geometry_optimization')
        if row['stage']=='geometry_optimization' and not job['requested'].get('optimization'):
            row['stage']='frequency_SCF' if job['requested'].get('frequencies') else 'SCF_or_gradient'
        cycles=re.findall(r'GEOMETRY OPTIMIZATION CYCLE\s+(\d+)',log)
        row['geometry_cycle']=int(cycles[-1]) if cycles else None
        if row['stage']=='electronic_response':
            it=re.findall(r'^\s*ITERATION\s+(\d+):',log,re.M);row['response_iteration']=int(it[-1]) if it else None
    else:
        it=re.findall(r'iteration #\s*(\d+)',log)
        row['stage']='SCF_converged_force_and_output' if 'convergence has been achieved' in log else 'SCF_iterations'
        row['SCF_iteration']=int(it[-1]) if it else None
    active.append(row)
core=[j for j in jobs if j['job_id'].startswith('w2-')]
parallel=[j for j in jobs if j['job_id'].startswith('p2-')]
result={'observed_at':datetime.now(timezone.utc).isoformat(),'DFT_jobs_passed':sum(j['status']=='passed' for j in core),
        'passed_by_engine':{e:sum(j['status']=='passed' and j['engine']==e for j in core) for e in ('orca','qe')},
        'DFT_jobs_requiring_review':[j['job_id'] for j in core if j['status'] not in ('running','passed')],
        'parallel_validation_passed':sum(j['status']=='passed' for j in parallel),
        'active_engine_jobs':active}
state=ROOT/'runs/week2-finish-controller/state.json'
resumed=ROOT/'runs/week2-resume-20260911/state.json'
if resumed.exists():state=resumed
if state.exists():
    controller=json.loads(state.read_text());result['continuation_status']=controller['status']
    if controller.get('current_command'):result['continuation_command']=controller['current_command']
pause=ROOT/'experiments/week2/pause-checkpoint.json'
status=ROOT/'experiments/week2/status.json'
if status.exists():result['task_status']=json.loads(status.read_text())['status']
if not resumed.exists() and result.get('task_status')=='resumed_parallel_validation':
    result['continuation_status']='resumed_parallel_validation'
if pause.exists():
    result['automatic_restart']=False
    result['paused_checkpoint']='experiments/week2/pause-checkpoint.json'
    paused=json.loads(pause.read_text())['source_run_dir']
    result['user_paused_jobs']=[j['job_id'] for j in jobs if str(Path(j['run_dir']).relative_to(ROOT))==paused]
    result['DFT_jobs_requiring_review']=[j for j in result['DFT_jobs_requiring_review'] if j not in result['user_paused_jobs']]
print(json.dumps(result,indent=2))
