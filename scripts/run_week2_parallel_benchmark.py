"""Bounded fixed-geometry MPI validation after the week-two core engines finish.

A small water Opt/Freq checks the parallel backend. Two reversed orders measure
1/8/12/16-process EnGrad costs on the same 50-atom geometry; no week-three panel.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import signal
import statistics
import time

import numpy as np
from ase.io import read
from research_run import run_job, sha
from week2_resume_provenance import validate_checkpoint

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'experiments/week2/parallel-benchmark'
REVISION=1

def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2)+'\n')

def frozen(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.read_text()==text, f'Changed frozen input: {path}'
    else:
        path.write_text(text)

def prepare(species, xyz, nprocs, repeat):
    atoms = read(xyz)
    job = f'p2-{species}-n{nprocs}-r{repeat}' if REVISION==1 else f'p2-v{REVISION}-{species}-n{nprocs}-r{repeat}'
    inp = ROOT/f'inputs/orca/week2/{job}.inp'
    coordinates = '\n'.join(f'{s} {p[0]:.14f} {p[1]:.14f} {p[2]:.14f}' for s,p in zip(atoms.symbols,atoms.positions))
    water_smoke = species=='water'
    task = 'TightOpt Opt Freq' if water_smoke else 'EnGrad'
    text = f'! r2SCAN-3c TightSCF {task}\n%maxcore 2000\n%pal nprocs {nprocs} end\n\n* xyz 0 1\n{coordinates}\n*\n'
    frozen(inp, text)
    config = dict(protocol_id='week2-parallel-performance-extension', week=2, job_id=job,
                  task_id='PERF2', engine='orca', executable='scripts/orca_parallel.sh',
                  input=str(inp.relative_to(ROOT)), method='r2SCAN-3c', charge=0,
                  multiplicity=1, solvent='vacuum', n_atoms=len(atoms), linear=False,
                  requested={'optimization':water_smoke,'frequencies':water_smoke,'gradient':True},
                  timeout_seconds=900, rss_limit_gib=40, omp_threads=1,
                  mpi_processes=nprocs, benchmark_repeat=repeat,
                  geometry_source=str(xyz.relative_to(ROOT)), geometry_source_sha256=sha(xyz),
                  purpose='Water MPI Opt/Freq backend smoke' if water_smoke else 'Fixed-geometry energy and gradient throughput; not Opt/Freq speed or a new minimum',
                  assets={'environment/orca-mpi-explicit.txt':'provenance/orca-mpi-explicit.txt',
                          **{f'.tools/orca-mpi/lib/{lib}':lib for lib in ['libmpi.40.dylib','libopen-rte.40.dylib','libopen-pal.40.dylib','libopen-orted-mpir.dylib','libz.1.dylib']}})
    path = OUT/f'{job}.json'
    frozen(path, json.dumps(config,indent=2)+'\n')
    return path

def inspect(config_path):
    config=json.loads(config_path.read_text())
    report=ROOT/f"reports/{config['job_id']}.json"
    m=json.loads(report.read_text()); run=Path(m['run_dir'])
    assert m['status']=='passed' and m['config_sha256']==sha(config_path)
    assert m['output_sha256']==sha(run/'output.log')
    assert m['executable_sha256']==sha(ROOT/'scripts/orca_parallel.sh')
    gradient=run/f"{config['job_id']}.engrad"
    data=[line.strip() for line in gradient.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
    n=int(data[0]); energy=float(data[1]); g=np.array([float(x) for x in data[2:2+3*n]]).reshape(n,3)
    assert n==config['n_atoms'] and np.isfinite(g).all()
    assert abs(energy-m['result']['energy'])<1e-9
    log=(run/'output.log').read_text()
    ranks=[int(x) for x in re.findall(r'running with\s+(\d+)\s+parallel MPI-processes',log,re.I)]
    if config['mpi_processes']>1:
        assert ranks and all(x==config['mpi_processes'] for x in ranks), 'MPI rank log requires review'
    return dict(job_id=config['job_id'], nprocs=config['mpi_processes'], repeat=config['benchmark_repeat'],
                elapsed_seconds=m['elapsed_seconds'], energy_Hartree=energy, gradient_Hartree_Bohr=g.tolist(),
                peak_sampled_RSS_GiB=m['peak_sampled_process_tree_rss_bytes']/1024**3,
                run_dir=m['run_dir'], MPI_rank_log=ranks, engrad_sha256=sha(gradient),
                config_sha256=m['config_sha256'], output_sha256=m['output_sha256'],
                vibrational_frequencies_cm1=m['result']['frequencies_cm1'], minimum_verified=m['result']['minimum_verified'])

def main():
    global OUT,REVISION
    parser=argparse.ArgumentParser();parser.add_argument('--execute',action='store_true')
    parser.add_argument('--checkpoint',type=Path);parser.add_argument('--revision',type=int,default=1);args=parser.parse_args()
    REVISION=args.revision
    if REVISION>1:OUT=OUT/f'revision-{REVISION}'
    proof=None
    if args.checkpoint:
        proof=validate_checkpoint(args.checkpoint.resolve());xyz=ROOT/proof['optimized_geometry']
        for manifest in (ROOT/'runs').glob('*-w2-*/manifest.json'):
            assert json.loads(manifest.read_text())['status']!='running', 'Core engine is active'
    else:
        controller=json.loads((ROOT/'runs/week2-finish-controller/state.json').read_text())
        assert controller['status']=='finite_pipeline_finished', 'Keep the already-running core pipeline exclusive'
        verification=json.loads((ROOT/'reports/week2-drug-verification.json').read_text());assert verification['status']=='passed'
        dimer=next(r for r in verification['checks'] if r['species']=='diflunisal--caffeine')
        xyz=ROOT/dimer['optimized_xyz'];assert sha(xyz)==dimer['optimized_xyz_sha256']
    smoke=prepare('water',ROOT/'structures/optimized/week1/water.xyz',2,0)
    configs=[prepare('dimer',xyz,n,rep) for rep,order in [(1,[1,8,12,16]),(2,[16,12,8,1])] for n in order]
    plan=dict(status='prepared', purpose='User-requested core-count review',
              optimization_provenance=proof,
              geometry_source_sha256=sha(xyz), main_ORCA_sha256=sha(ROOT/'.tools/orca-6.1.1/orca'),
              mpirun_sha256=sha(ROOT/'.tools/orca-mpi/bin/mpirun'), libmpi_sha256=sha(ROOT/'.tools/orca-mpi/lib/libmpi.40.dylib'),
              launcher_sha256=sha(ROOT/'scripts/orca_parallel.sh'),
              jobs=[str(p.relative_to(ROOT)) for p in [smoke,*configs]],
              energy_tolerance_Hartree=1e-6, maximum_gradient_difference_tolerance_Hartree_Bohr=1e-5,
              budget_seconds=2700, automatic_retries=0,
              scope='Water Opt/Freq backend check, followed by same dimer geometry, method, SCF tolerance and fresh initial guess. Two EnGrad repeats, reverse order; OS/thermal variance remains. Dimer analytical Hessian scaling is not measured.')
    frozen(OUT/'plan.json',json.dumps(plan,indent=2)+'\n')
    if not args.execute:
        print(json.dumps(plan,indent=2));return
    started=time.monotonic(); rows=[]
    for config in [smoke,*configs]:
        if time.monotonic()-started>plan['budget_seconds']:raise RuntimeError('Benchmark review budget reached')
        job=json.loads(config.read_text())['job_id']
        if not (ROOT/f'reports/{job}.json').exists():
            assert run_job(config)==0, f'Engine needs review: {job}'
        row=inspect(config);rows.append(row)
        save(ROOT/'reports/week2-parallel-benchmark.json',dict(status='running',rows=rows,plan=plan))
        if config==smoke:
            serial=json.loads((ROOT/'reports/week1-water-opt-freq.json').read_text())['result']
            assert row['minimum_verified'] is True
            assert abs(row['energy_Hartree']-serial['energy'])<=1e-6
            np.testing.assert_allclose(row['vibrational_frequencies_cm1'],serial['frequencies_cm1'],rtol=0,atol=.1)
    samples=rows[1:]; ref=samples[0]
    for r in samples:
        r['energy_delta_reference_Hartree']=r['energy_Hartree']-ref['energy_Hartree']
        r['max_gradient_delta_reference_Hartree_Bohr']=float(np.max(np.linalg.norm(np.array(r['gradient_Hartree_Bohr'])-ref['gradient_Hartree_Bohr'],axis=1)))
    good=all(abs(r['energy_delta_reference_Hartree'])<=plan['energy_tolerance_Hartree'] and r['max_gradient_delta_reference_Hartree_Bohr']<=plan['maximum_gradient_difference_tolerance_Hartree_Bohr'] for r in samples)
    timings={n:statistics.median(r['elapsed_seconds'] for r in samples if r['nprocs']==n) for n in [1,8,12,16]}
    fastest=min([8,12,16],key=lambda n:timings[n])
    result=dict(status='passed' if good else 'numerical_agreement_review_required', rows=rows, plan=plan,
                median_wall_seconds_by_nprocs=timings, speedup_vs_one={n:timings[1]/t for n,t in timings.items()},
                fastest_observed_parallel_nprocs=fastest, recommendation_requires_Freq_scaling_check=True,
                total_engine_wall_seconds=sum(r['elapsed_seconds'] for r in rows),
                no_new_target_optimization_or_minimum=True, water_OptFreq_backend_checked=True, no_production_jobs_changed=True)
    save(ROOT/'reports/week2-parallel-benchmark.json',result)
    files=[OUT/'plan.json',*configs,smoke,ROOT/'environment/orca-mpi-explicit.txt',Path(__file__),ROOT/'scripts/orca_parallel.sh']
    for r in rows:
        files.extend(p for p in Path(r['run_dir']).iterdir() if p.is_file() and p.suffix in ['.inp','.engrad','.hess','.xyz','.json','.log','.py','.txt'])
    save(ROOT/'reports/week2-parallel-evidence-sha256.json',dict(algorithm='SHA-256',files={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(files))}))
    print(json.dumps({k:v for k,v in result.items() if k not in ['rows','plan']},indent=2))
    if not good:raise RuntimeError('Parallel energy/gradient comparison needs review')

if __name__=='__main__':
    signal.signal(signal.SIGTERM,lambda *_: (_ for _ in ()).throw(KeyboardInterrupt('terminated')))
    main()
