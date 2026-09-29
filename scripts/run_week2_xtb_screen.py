"""Sequential, bounded GFN2-xTB optimization and graph/geometry audit of starts.

No automatic retries. Geometry-based proton ownership is not an electronic charge.
Reuse checks both candidate and executable hashes, and preserves all failed outputs.
"""
import argparse
from datetime import datetime,timezone
import fcntl
import hashlib
import json
from pathlib import Path
import os
import re
import signal
import subprocess
import time
import numpy as np
import psutil
from ase.io import read
from rdkit import Chem

ROOT=Path(__file__).resolve().parents[1]
P=Chem.GetPeriodicTable()
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(path,value):path.write_text(json.dumps(value,indent=2)+'\n')

def geometry_audit(candidate,atoms):
    z=candidate['graph']['atomic_numbers'];p=atoms.positions
    if len(atoms)!=len(z) or atoms.numbers.tolist()!=z:return {'status':'atom_order_or_composition_changed'}
    heavy=[i for i,a in enumerate(z) if a>1]
    expected={tuple(sorted((i,j))) for i,j,_ in candidate['graph']['bonds'] if z[i]>1 and z[j]>1}
    observed={(i,j) for i in heavy for j in heavy if i<j and np.linalg.norm(p[i]-p[j])<1.25*(P.GetRcovalent(z[i])+P.GetRcovalent(z[j]))}
    changes=[];unassigned=[]
    for h,owner in candidate['graph']['hydrogen_owners'].items():
        h=int(h);nearest=min(heavy,key=lambda j:np.linalg.norm(p[h]-p[j]));distance=float(np.linalg.norm(p[h]-p[nearest]))
        if distance>1.45*(P.GetRcovalent(1)+P.GetRcovalent(z[nearest])):unassigned.append(h)
        if nearest!=owner:changes.append({'hydrogen':h,'from':owner,'to':nearest,'bond_distance_angstrom':distance})
    result={'status':'geometry_checked','heavy_connectivity_preserved':expected==observed,
            'lost_heavy_bonds':[list(x) for x in sorted(expected-observed)],'new_heavy_bonds':[list(x) for x in sorted(observed-expected)],
            'hydrogen_owner_changes':changes,'unassigned_hydrogens':unassigned,
            'chemical_state_preserved_by_geometry':expected==observed and not changes and not unassigned}
    if candidate['kind']=='dimer':
        n=candidate['fragment_atom_counts'][0]
        result['closest_interfragment_distance_angstrom']=float(np.min(np.linalg.norm(p[:n,None,:]-p[None,n:,:],axis=2)))
        result['separated_geometry']=result['closest_interfragment_distance_angstrom']>4.0
        result['cross_fragment_proton_transfer']=[x for x in changes if (x['from']<n)!=(x['to']<n)]
    return result

def parse(run,candidate,meta):
    text=(run/'output.log').read_text(errors='replace')
    error=(run/'stderr.log').read_text(errors='replace')
    values=re.findall(r'TOTAL ENERGY\s+([-+\d.Ee]+)\s+Eh',text)
    normal='normal termination of xtb' in text+error
    opt='GEOMETRY OPTIMIZATION CONVERGED' in text
    scf='convergence criteria satisfied' in text
    energy=float(values[-1]) if values else None
    result={'normal_termination':normal,'geometry_optimization_converged':opt,'scf_converged':scf,
            'energy_Hartree':energy,'minimum_verified':False,'method':'GFN2-xTB, semiempirical, not DFT',
            'geometry_audit':geometry_audit(candidate,read(run/'xtbopt.xyz')) if (run/'xtbopt.xyz').is_file() else None}
    version=re.search(r'xtb version\s+(\S+)',text);result['engine_version']=version.group(1) if version else None
    result['status']='passed' if meta.get('exit_code')==0 and not meta.get('stop_reason') and normal and opt and scf and energy is not None and np.isfinite(energy) else 'needs_review'
    return result

def main():
    def interrupted(signum, frame):
        raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM, interrupted)
    parser=argparse.ArgumentParser();parser.add_argument('--phase',choices=['monomers','dimers'],required=True)
    parser.add_argument('--parse-existing',action='store_true')
    parser.add_argument('--reviewed-continuation',type=Path,help='Explicit saved review of the normal-SCF/geometry-iteration-limit failure class')
    args=parser.parse_args()
    inventory_path=ROOT/'experiments/week2/chemistry-candidates.json';inventory=json.loads(inventory_path.read_text())
    exe=ROOT/'.tools/xtb/bin/xtb';exe_hash=sha(exe)
    out=ROOT/f'runs/week2-xtb-{args.phase}';out.mkdir(parents=True,exist_ok=True)
    review=None
    if args.reviewed_continuation:
        review=json.loads(args.reviewed_continuation.read_text())
        assert review['phase']==args.phase and review['action']=='continue_unattempted_candidates_without_retry'
        assert review['reviewed_failure_class']=='normal_scf_converged_geometry_iteration_limit'
        for entry in review['evidence']:
            m=json.loads((ROOT/entry['run_dir']/'manifest.json').read_text())
            assert sha(ROOT/entry['run_dir']/'output.log')==entry['output_sha256']==m['output_sha256']
            assert m['exit_code']==0 and m['result']['status']=='needs_review'
        previous=ROOT/f'reports/week2-xtb-{args.phase}.json'
        archived=out/'batch_before_reviewed_continuation.json'
        if previous.exists() and not archived.exists():archived.write_bytes(previous.read_bytes())
    frozen_inventory=out/'chemistry-candidates.json'
    if frozen_inventory.exists():
        assert sha(frozen_inventory)==sha(inventory_path), 'Candidate metadata changed; review before resuming'
    else:frozen_inventory.write_bytes(inventory_path.read_bytes())
    snapshot=out/'runner_snapshot.py'
    if snapshot.exists():
        if sha(snapshot)!=sha(Path(__file__)):
            assert review or args.parse_existing, 'Runner changed; review existing attempts before resuming'
            snapshot=out/f'runner_revision_{sha(Path(__file__))[:12]}.py'
            if snapshot.exists():assert sha(snapshot)==sha(Path(__file__))
            else:snapshot.write_bytes(Path(__file__).read_bytes())
    else:snapshot.write_bytes(Path(__file__).read_bytes())
    index={'phase':args.phase,'inventory_sha256':sha(inventory_path),'executable_sha256':exe_hash,'status':'running','results':[],
           'runner_sha256':sha(Path(__file__)),'started_at':datetime.now(timezone.utc).isoformat(),'automatic_retries':0}
    if review:index['reviewed_continuation']={'path':str(args.reviewed_continuation),'sha256':sha(args.reviewed_continuation),'details':review}
    report=ROOT/f'reports/week2-xtb-{args.phase}.json'
    env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
    failures=0;start=time.monotonic()
    with (ROOT/'experiments/.research-run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for candidate in inventory[args.phase]:
            run=out/candidate['candidate_id'];manifest_path=run/'manifest.json'
            if manifest_path.exists():
                meta=json.loads(manifest_path.read_text())
                assert meta['candidate_sha256']==candidate['sha256'] and meta['executable_sha256']==exe_hash
                assert json.loads((run/'candidate.json').read_text())==candidate
                assert sha(run/'input.xyz')==candidate['sha256']
                if 'exit_code' not in meta or ('result' not in meta and not args.parse_existing):
                    index.update(status='incomplete_existing_attempt_requires_review',blocked_candidate=candidate['candidate_id']);break
                if args.parse_existing:
                    result=parse(run,candidate,meta);meta['result']=result;save(manifest_path,meta)
                else:result=meta['result']
            else:
                if args.parse_existing:continue
                if time.monotonic()-start>3600 or sum(r['elapsed_seconds'] or 0 for r in index['results'])>3600:
                    index['status']='screening_budget_review';break
                run.mkdir();source=ROOT/candidate['path'];assert sha(source)==candidate['sha256']
                (run/'input.xyz').write_bytes(source.read_bytes())
                (run/'candidate.json').write_text(json.dumps(candidate,indent=2)+'\n')
                cmd=[str(exe),'input.xyz','--gfn','2','--opt','tight','--chrg','0','--uhf','0']
                meta={'candidate_sha256':candidate['sha256'],'executable_sha256':exe_hash,'command':cmd,'stop_reason':None,
                      'runner_sha256':sha(Path(__file__)),'runner_snapshot':str(snapshot.relative_to(ROOT)),
                      'started_at':datetime.now(timezone.utc).isoformat(),'thread_environment':{k:env[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS')},
                      'timeout_seconds':120,'rss_limit_gib':40}
                save(manifest_path,meta)
                began=time.monotonic();peak=0;samples=0;denied=0
                with (run/'output.log').open('w') as stdout,(run/'stderr.log').open('w') as stderr:
                    proc=subprocess.Popen(cmd,cwd=run,env=env,stdout=stdout,stderr=stderr,start_new_session=True)
                    try:
                        while proc.poll() is None:
                            try:
                                p=psutil.Process(proc.pid);rss=sum(x.memory_info().rss for x in [p,*p.children(recursive=True)]);peak=max(peak,rss);samples+=1
                            except psutil.AccessDenied:denied+=1
                            except (psutil.Error,OSError):pass
                            if time.monotonic()-began>120 or peak>40*1024**3:
                                meta['stop_reason']='time_or_rss_limit';os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=10);break
                            time.sleep(0.1)
                    except BaseException:
                        meta['stop_reason']='interrupted_or_runner_error'
                        raise
                    finally:
                        if proc.poll() is None:
                            os.killpg(proc.pid,signal.SIGKILL);proc.wait()
                        meta.update(exit_code=proc.returncode,elapsed_seconds=round(time.monotonic()-began,4),
                                    peak_sampled_process_tree_rss_bytes=peak or None,memory_samples=samples,memory_access_denied=denied)
                        save(manifest_path,meta)
                meta['output_sha256']=sha(run/'output.log')
                if (run/'xtbopt.xyz').exists():meta['optimized_xyz_sha256']=sha(run/'xtbopt.xyz')
                result=parse(run,candidate,meta);meta['result']=result;save(manifest_path,meta)
            row={'candidate_id':candidate['candidate_id'],'species':candidate['species'],'run_dir':str(run.relative_to(ROOT)),
                 'elapsed_seconds':meta.get('elapsed_seconds'),'candidate':candidate,'result':result}
            index['results'].append(row);save(report,index)
            known_iteration_limit=bool(review and meta.get('exit_code')==0 and not meta.get('stop_reason')
                                      and result['normal_termination'] and result['scf_converged']
                                      and not result['geometry_optimization_converged']
                                      and re.search(r'FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN \d+ ITERATIONS',
                                                    (run/'output.log').read_text(errors='replace')))
            # A reviewed geometry failure stays excluded; only orchestration changes.
            failures=failures+1 if result['status']!='passed' and not known_iteration_limit else 0
            if len(index['results'])%10==0:print(json.dumps({'phase':args.phase,'completed':len(index['results']),'latest_energy_Ha':result['energy_Hartree']}),flush=True)
            if failures>=3:index['status']='three_consecutive_failures_requires_review';break
        else:index['status']='completed'
    index.update(ended_at=datetime.now(timezone.utc).isoformat(),batch_elapsed_seconds=round(time.monotonic()-start,3),
                 calculation_wall_seconds=sum(r['elapsed_seconds'] or 0 for r in index['results']))
    save(report,index)
    print(json.dumps({'phase':args.phase,'status':index['status'],'completed':len(index['results'])}),flush=True)
    return 0 if index['status']=='completed' else 1
if __name__=='__main__':raise SystemExit(main())
