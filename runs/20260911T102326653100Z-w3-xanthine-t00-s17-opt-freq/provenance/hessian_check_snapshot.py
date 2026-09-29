"""A new Freq-only run linked to the preserved, completed 87-cycle optimization.

The original interrupted run is never relabelled as normally completed. Separate
optimization evidence, frozen coordinates and the new Hessian establish the minimum.
"""
import argparse
import json
from pathlib import Path
import signal

import numpy as np
from ase.io import read

import research_run
from week2_resume_provenance import ROOT, sha, validate_checkpoint

CHECKPOINT=ROOT/'experiments/week2/pause-checkpoint.json'
LIBS=['libmpi.40.dylib','libopen-rte.40.dylib','libopen-pal.40.dylib','libopen-orted-mpir.dylib','libz.1.dylib']
BASE_PARSE=research_run.parse_output
BOHR_ANGSTROM=0.529177210903

def input_text(nprocs):
    return f'! r2SCAN-3c TightSCF Freq\n%maxcore 2000\n%pal nprocs {nprocs} end\n\n* xyzfile 0 1 saved-optimization.xyz\n'

def inspect_hessian(path, expected_atoms, logged_frequencies):
    lines=path.read_text().splitlines()
    start=lines.index('$vibrational_frequencies');n=int(lines[start+1])
    frequencies=np.array([float(line.split()[1]) for line in lines[start+2:start+2+n]])
    assert n==3*len(expected_atoms) and len(frequencies)==n and np.isfinite(frequencies).all()
    np.testing.assert_allclose(frequencies[:6],0,rtol=0,atol=1e-8)
    np.testing.assert_allclose(frequencies[6:],logged_frequencies,rtol=0,atol=.0051)
    start=lines.index('$atoms');natoms=int(lines[start+1]);assert natoms==len(expected_atoms)
    records=[line.split() for line in lines[start+2:start+2+natoms]]
    assert [r[0] for r in records]==expected_atoms.get_chemical_symbols()
    positions=np.array([[float(x) for x in r[2:5]] for r in records])*BOHR_ANGSTROM
    # ORCA may translate/rotate to its Hessian frame. Preserve atom mapping and
    # require the same internal geometry under a proper rigid transformation.
    p=positions-positions.mean(axis=0);q=expected_atoms.positions-expected_atoms.positions.mean(axis=0)
    u,_,vt=np.linalg.svd(p.T@q);rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    np.testing.assert_allclose(p@rotation,q,rtol=0,atol=1e-7)
    return frequencies[6:]

def resume_parse(run,meta):
    result=BASE_PARSE(run,meta)
    result.update(optimization_from_checkpoint_verified=False,minimum_verified=False,status='needs_review')
    try:
        proof=validate_checkpoint(CHECKPOINT)
        assert meta['optimization_provenance']==proof
        assert meta['method']=='r2SCAN-3c' and meta['charge']==0 and meta['multiplicity']==1 and meta['n_atoms']==50
        assert meta['requested']=={'optimization':False,'frequencies':True}
        assert (run/Path(meta['input']).name).read_text()==input_text(meta['mpi_processes'])
        assert sha(run/'saved-optimization.xyz')==proof['geometry_sha256']
        assert sha(run/'provenance/optimization-output.log')==proof['optimization_log_sha256']
        assert result['geometry_optimization_converged'] is None
        assert not result['optimization_cycles'], 'Freq-only continuation must not optimize'
        atoms=read(run/'saved-optimization.xyz')
        hess=run/(Path(meta['input']).stem+'.hess')
        frequencies=inspect_hessian(hess,atoms,result['frequencies_cm1'])
        delta=result['energy']-proof['optimization_energy_Hartree']
        assert abs(delta)<=1e-6, 'Changed method/state or numerical consistency needs review'
        complete=(meta.get('exit_code')==0 and not meta.get('stop_reason') and result['normal_termination']
                  and result['scf_converged'] and result['frequency_analysis_complete'])
        minimum=bool(complete and np.all(frequencies>0))
        result.update(optimization_from_checkpoint_verified=True,minimum_verified=minimum,
                      status='passed' if minimum else 'needs_review',
                      optimization_provenance=proof,frequency_geometry_sha256=sha(run/'saved-optimization.xyz'),
                      hessian_sha256=sha(hess),hessian_coordinates_match_checkpoint=True,
                      hessian_geometry_match_scope='Same atom mapping up to translation and proper rigid rotation; 1e-7 Angstrom coordinate tolerance',
                      full_precision_frequencies_cm1=frequencies.tolist(),
                      energy_difference_from_optimization_Hartree=delta,
                      verification_scope='Completed optimization in preserved original log plus independent Freq-only Hessian at identical coordinates')
    except (AssertionError,ValueError,KeyError,FileNotFoundError,IndexError) as exc:
        result['resume_validation_error']=f'{type(exc).__name__}: {exc}'
    result['resume_parser_sha256']=sha(Path(__file__))
    research_run.save(run/'result.json',result)
    return result

def prepare():
    bench=json.loads((ROOT/'reports/week2-parallel-benchmark.json').read_text())
    assert bench['status']=='passed' and bench['water_OptFreq_backend_checked']
    proof=validate_checkpoint(CHECKPOINT)
    assert bench['plan']['geometry_source_sha256']==proof['geometry_sha256']
    assert bench['plan']['main_ORCA_sha256']==sha(ROOT/'.tools/orca-6.1.1/orca')
    assert bench['plan']['libmpi_sha256']==sha(ROOT/'.tools/orca-mpi/lib/libmpi.40.dylib')
    fastest=bench['fastest_observed_parallel_nprocs'];times=bench['median_wall_seconds_by_nprocs']
    n=12 if times['12']<=1.05*times[str(fastest)] else fastest
    assert n in [8,12,16]
    job=f'w2-diflunisal--caffeine-freq-resume-n{n}'
    inp=ROOT/f'inputs/orca/week2/{job}.inp'
    if inp.exists():assert inp.read_text()==input_text(n)
    else:inp.write_text(input_text(n))
    assets={f'.tools/orca-mpi/lib/{lib}':lib for lib in LIBS}
    assets.update({proof['optimized_geometry']:'saved-optimization.xyz',
                   proof['optimization_log']:'provenance/optimization-output.log',
                   'experiments/week2/pause-checkpoint.json':'provenance/pause-checkpoint.json',
                   'environment/orca-mpi-explicit.txt':'provenance/orca-mpi-explicit.txt',
                   'scripts/run_week2_frequency_resume.py':'resume_parser_snapshot.py',
                   'scripts/week2_resume_provenance.py':'resume_provenance_snapshot.py',
                   'scripts/orca_parallel.sh':'parallel_launcher_snapshot.sh'})
    config=dict(protocol_id='week2-resume-20260911',week=2,task_id='DR01',job_id=job,engine='orca',
                executable='scripts/orca_parallel.sh',input=str(inp.relative_to(ROOT)),method='r2SCAN-3c',
                charge=0,multiplicity=1,n_atoms=50,linear=False,solvent='vacuum',
                requested={'optimization':False,'frequencies':True},mpi_processes=n,omp_threads=1,
                timeout_seconds=7200,rss_limit_gib=40,assets=assets,optimization_provenance=proof,
                purpose='Complete interrupted frequency validation at identical saved optimized coordinates',
                process_selection_rule='Use 12 processes when its measured median EnGrad time is within 5 percent of the fastest tested setting; otherwise use the fastest',
                parallel_benchmark_sha256=sha(ROOT/'reports/week2-parallel-benchmark.json'))
    path=ROOT/'experiments/week2/frequency-resume.json';text=json.dumps(config,indent=2)+'\n'
    if path.exists():assert path.read_text()==text,'Existing resume config must be preserved'
    else:path.write_text(text)
    return path

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--execute',action='store_true')
    parser.add_argument('--parse-run',type=Path);args=parser.parse_args()
    if args.parse_run:
        run=args.parse_run.resolve();meta=json.loads((run/'manifest.json').read_text())
        assert meta['status']!='running' and 'exit_code' in meta
        result=resume_parse(run,meta);meta.update(result=result,status=result['status'])
        research_run.save(run/'manifest.json',meta);research_run.save(ROOT/f"reports/{meta['job_id']}.json",meta)
        print(json.dumps(result,indent=2));return 0 if result['status']=='passed' else 1
    path=prepare();print(json.dumps({'config':str(path.relative_to(ROOT))}),flush=True)
    if not args.execute:return 0
    config=json.loads(path.read_text());report=ROOT/f"reports/{config['job_id']}.json"
    assert not report.exists(),'A previous resume attempt exists; review it rather than retry automatically'
    research_run.parse_output=resume_parse
    return research_run.run_job(path)

if __name__=='__main__':
    signal.signal(signal.SIGTERM,lambda *_: (_ for _ in ()).throw(KeyboardInterrupt('terminated')))
    raise SystemExit(main())
