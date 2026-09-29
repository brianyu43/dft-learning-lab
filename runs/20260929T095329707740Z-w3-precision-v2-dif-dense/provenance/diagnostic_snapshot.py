"""Fixed-geometry DIF/DI2 grid test: six controlled EnGrad jobs, no retries."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from ase.io import read

from prepare_week3_drug import immutable, LIBS
from research_run import ROOT, sha, save, now
from week3_recovery import gradient_in_frame

OUT=ROOT/'reports/completion'
PLAN=ROOT/'experiments/completion/precision-diagnostic.json'
HARTREE_KCAL=627.509474


def prepare():
    OUT.mkdir(exist_ok=True)
    prior=json.loads((ROOT/'reports/week2-drug-verification.json').read_text())
    mono=next(r for r in prior['checks'] if r['species']=='diflunisal')
    assert mono['minimum_verified'] and sha(ROOT/mono['optimized_xyz'])==mono['optimized_xyz_sha256']
    dimer_report=ROOT/'reports/w3-diflunisal--diflunisal-p05-refine-r1-opt-freq.json'
    dimer=json.loads(dimer_report.read_text())
    assert dimer['status']=='passed' and dimer['result']['minimum_verified']
    dpath=Path(dimer['run_dir'])/(Path(dimer['input']).stem+'.xyz')
    assert sha(dpath)==dimer['result']['optimized_xyz_sha256']
    sources={'dif':ROOT/mono['optimized_xyz'],'dif2':dpath}
    settings={'g2':('DefGrid2',''),'g3':('DefGrid3',''),
              'dense':('DefGrid3','%method\n  AngularGrid 7\n  IntAcc 6.0\nend\n')}
    jobs=[]
    for setting,(keyword,extra) in settings.items():
        for species,source in sources.items():
            job=f'w3-precision-v2-{species}-{setting}'
            structure=ROOT/f'structures/completion/{job}.xyz'
            immutable(structure,source.read_text())
            text=f'! r2SCAN-3c VeryTightSCF {keyword} EnGrad\n%maxcore 2000\n%pal nprocs 12 end\n{extra}\n* xyzfile 0 1 initial.xyz\n'
            inp=ROOT/f'inputs/orca/completion/{job}.inp';immutable(inp,text)
            cfg=dict(protocol_id='dft-completion-20260929-precision-v2',week=3,task_id='C02',
                job_id=job,engine='orca',executable='scripts/orca_parallel.sh',input=str(inp.relative_to(ROOT)),
                method='r2SCAN-3c',charge=0,multiplicity=1,solvent='vacuum',n_atoms=len(read(source)),linear=False,
                requested={'optimization':False,'frequencies':False,'gradient':True},
                timeout_seconds=1200,rss_limit_gib=40,mpi_processes=12,omp_threads=1,
                assets={**{f'.tools/orca-mpi/lib/{lib}':lib for lib in LIBS},
                        str(structure.relative_to(ROOT)):'initial.xyz',
                        'scripts/completion_precision_diagnostic.py':'provenance/diagnostic_snapshot.py',
                        'scripts/orca_parallel.sh':'provenance/launcher_snapshot.sh'},
                species=species,grid_setting=setting,source_geometry=str(source.relative_to(ROOT)),
                source_geometry_sha256=sha(source),initial_structure_sha256=sha(structure),
                purpose='Same-geometry grid sensitivity in DIF self-association, not a new minimum or pose')
            cfgpath=ROOT/f'experiments/completion/{job}.json'
            immutable(cfgpath,json.dumps(cfg,indent=2)+'\n');jobs.append(str(cfgpath.relative_to(ROOT)))
    plan=dict(protocol_id='dft-completion-20260929-precision-v2',jobs=jobs,
        grid_definitions={k:{'simple_keyword':v[0],'extra_method_block':v[1]} for k,v in settings.items()},
        original_failed_gate='reports/week3-recovery-r1/numerical-check.json',
        original_failed_gate_sha256=sha(ROOT/'reports/week3-recovery-r1/numerical-check.json'),
        fixed_sources={k:{'path':str(v.relative_to(ROOT)),'sha256':sha(v)} for k,v in sources.items()},
        contrast='E(DIF-DIF at saved refined geometry) - 2 E(DIF at saved optimized geometry)',
        scf='VeryTightSCF in all six jobs',
        target_contrast_change_kcal_mol=.1,
        decision='Measure first. No automatic change of prior failed gate or launch of the pending panel.',
        reuse_scope='Two fixed geometries; numerical convergence of other candidates must be checked separately',
        input_hashes={p:sha(ROOT/p) for p in jobs})
    immutable(PLAN,json.dumps(plan,indent=2)+'\n')
    return plan


def validate(config_path):
    cfg=json.loads((ROOT/config_path).read_text());report=ROOT/f"reports/{cfg['job_id']}.json"
    meta=json.loads(report.read_text());run=Path(meta['run_dir'])
    assert meta['status']=='passed' and meta['config_sha256']==sha(ROOT/config_path)
    assert sha(run/'output.log')==meta['output_sha256']
    for source,digest in meta['source_sha256'].items():
        target=Path(meta['input']).name if source==meta['input'] else meta['assets'][source]
        assert sha(run/target)==digest
    atoms=read(run/'initial.xyz')
    assert sha(run/'initial.xyz')==cfg['initial_structure_sha256']
    energy,gradient,engrad,digest=gradient_in_frame(meta,atoms)
    import re
    ranks=[int(x) for x in re.findall(r'running with\s+(\d+)\s+parallel MPI-processes',(run/'output.log').read_text(),re.I)]
    assert ranks and all(x==12 for x in ranks)
    return dict(job_id=cfg['job_id'],species=cfg['species'],grid_setting=cfg['grid_setting'],
        energy_Hartree=energy,gradient_Hartree_Bohr=gradient.tolist(),engrad_sha256=digest,
        same_geometry_verified=True,elapsed_seconds=meta['elapsed_seconds'],
        report=str(report.relative_to(ROOT)),report_sha256=sha(report),run_dir=meta['run_dir'])


def summarize(plan):
    rows=[validate(p) for p in plan['jobs']]
    lookup={(r['species'],r['grid_setting']):r for r in rows}
    binding={g:(lookup['dif2',g]['energy_Hartree']-2*lookup['dif',g]['energy_Hartree'])*HARTREE_KCAL
             for g in plan['grid_definitions']}
    differences={}
    for a,b in [('g2','g3'),('g3','dense')]:
        differences[f'{a}_to_{b}']=dict(
            monomer_total_delta_kcal_mol=(lookup['dif',b]['energy_Hartree']-lookup['dif',a]['energy_Hartree'])*HARTREE_KCAL,
            dimer_total_delta_kcal_mol=(lookup['dif2',b]['energy_Hartree']-lookup['dif2',a]['energy_Hartree'])*HARTREE_KCAL,
            self_binding_delta_kcal_mol=binding[b]-binding[a],
            max_dimer_gradient_component_delta_Hartree_Bohr=float(np.max(np.abs(np.array(lookup['dif2',b]['gradient_Hartree_Bohr'])-np.array(lookup['dif2',a]['gradient_Hartree_Bohr'])))))
    result=dict(status='diagnostic_completed_requires_protocol_decision',updated_at=now(),rows=rows,
        self_binding_kcal_mol=binding,differences=differences,
        g3_vs_dense_contrast_within_target=abs(binding['dense']-binding['g3'])<=plan['target_contrast_change_kcal_mol'],
        original_failed_gate_preserved=sha(ROOT/plan['original_failed_gate'])==plan['original_failed_gate_sha256'],
        scope=plan['reuse_scope'])
    save(OUT/'precision-diagnostic.json',result)
    lines=['# DIF 자기 결합 수치 정밀도 진단','',f"계산 완료: {result['updated_at']}",'',
        '동일 좌표와 VeryTightSCF에서 격자만 변경했다. 원래 실패 판정은 보존한다.','',
        '| 설정 | E(AA)−2E(A), kcal/mol |','|---|---:|']
    lines += [f'| {g} | {v:.6f} |' for g,v in binding.items()]
    lines += ['', '| 변경 | 단량체 총에너지 변화 | 이합체 총에너지 변화 | 결합에너지 변화 |','|---|---:|---:|---:|']
    for k,v in differences.items():lines.append(f"| {k} | {v['monomer_total_delta_kcal_mol']:.6f} | {v['dimer_total_delta_kcal_mol']:.6f} | {v['self_binding_delta_kcal_mol']:.6f} |")
    lines += ['', '단위는 kcal/mol. 전체 패널의 B/X 정밀도를 대신하지 않으며, 구조 자체의 재최적화 효과도 이번 고정 좌표 비교에 포함하지 않는다.']
    (OUT/'PRECISION_DIAGNOSTIC.ko.md').write_text('\n'.join(lines)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--execute',action='store_true');args=parser.parse_args()
    plan=prepare()
    if not args.execute:
        print(json.dumps(plan,indent=2));return
    frozen=sha(Path(__file__))
    state=dict(status='running',started_at=now(),pid=os.getpid(),plan_sha256=sha(PLAN),code_sha256=frozen,steps=[])
    save(OUT/'precision-run-state.json',state)
    try:
        for p in plan['jobs']:
            assert sha(Path(__file__))==frozen and sha(ROOT/p)==plan['input_hashes'][p]
            cfg=json.loads((ROOT/p).read_text());report=ROOT/f"reports/{cfg['job_id']}.json"
            if report.exists():
                validate(p);continue
            assert not list((ROOT/'runs').glob(f"*-{cfg['job_id']}")), 'Existing attempt without report requires inspection; do not duplicate'
            original=json.loads((ROOT/'reports/week2-summary.json').read_text())
            spent=0.
            for mp in (ROOT/'runs').glob('*/manifest.json'):
                m=json.loads(mp.read_text())
                if m.get('engine')=='orca' and m.get('week',0)>=3:spent+=m.get('elapsed_seconds',0.)
            assert original['remaining_budget_hours']['drug']*3600-spent>=cfg['timeout_seconds']
            state.update(current_job=cfg['job_id'],updated_at=now());save(OUT/'precision-run-state.json',state)
            with (OUT/'precision-controller.log').open('a') as log:
                result=subprocess.run([sys.executable,'scripts/research_run.py','--config',p],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
            assert result.returncode==0,f"Review {cfg['job_id']}"
            row=validate(p);state['steps'].append(row);save(OUT/'precision-run-state.json',state)
        result=summarize(plan)
        state.update(status='completed',ended_at=now(),current_job=None,
                     g3_vs_dense_contrast_within_target=result['g3_vs_dense_contrast_within_target'])
    except BaseException as exc:
        state.update(status='stopped_requires_review',ended_at=now(),reason=f'{type(exc).__name__}: {exc}')
        raise
    finally:save(OUT/'precision-run-state.json',state)


if __name__=='__main__':main()
