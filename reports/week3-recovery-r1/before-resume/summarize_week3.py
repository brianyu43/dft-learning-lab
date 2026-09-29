"""Live and final week-three report; incomplete/failed jobs remain visible."""
import json
from pathlib import Path

import numpy as np
from ase.io import read

from research_run import ROOT, save, sha, now

HARTREE_KCAL=627.509474


def geometry_groups(rows):
    groups=[]
    for row in sorted((r for r in rows if r['status']=='passed' and r['minimum_verified']),key=lambda r:r['energy_Hartree']):
        path=Path(row['optimized_xyz']);atoms=read(path if path.is_absolute() else ROOT/path)
        for group in groups:
            representative=group['_row'];ref=group['_atoms']
            if row['species']!=representative['species'] or atoms.get_chemical_symbols()!=ref.get_chemical_symbols():continue
            p=atoms.positions-atoms.positions.mean(axis=0);q=ref.positions-ref.positions.mean(axis=0)
            u,_,vt=np.linalg.svd(p.T@q);rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
            rmsd=float(np.sqrt(np.mean(np.sum((p@rotation-q)**2,axis=1))))
            delta=abs(row['energy_Hartree']-representative['energy_Hartree'])*HARTREE_KCAL
            if rmsd<=.05 and delta<=.01:
                group['members'].append(dict(job_id=row['job_id'],all_atom_RMSD_angstrom=rmsd,delta_energy_kcal_mol=delta));break
        else:
            groups.append(dict(species=row['species'],representative_job_id=row['job_id'],
                               members=[dict(job_id=row['job_id'],all_atom_RMSD_angstrom=0.,delta_energy_kcal_mol=0.)],
                               _row=row,_atoms=atoms))
    return [{k:v for k,v in g.items() if not k.startswith('_')} for g in groups]


def summarize():
    plan_path=ROOT/'experiments/week3/drug/panel-plan.json'
    plan=json.loads(plan_path.read_text())
    prior=json.loads((ROOT/'reports/week2-drug-verification.json').read_text())
    rows=[]
    for old in prior['checks']:
        rows.append(dict(species=old['species'],candidate_id=old['candidate_id'],role='primary',source_week=2,
                         status='passed',minimum_verified=old['minimum_verified'],energy_Hartree=old['energy_Hartree'],
                         job_id=old['job_id'],optimized_xyz=old['optimized_xyz'],
                         geometry_audit=old['geometry_audit'],elapsed_seconds=old['elapsed_seconds']))
    for cfgpath in [plan['first_job'],*plan['remaining_jobs']]:
        cfg=json.loads((ROOT/cfgpath).read_text());p=ROOT/f"reports/{cfg['job_id']}.json"
        row=dict(species=cfg['species'],candidate_id=cfg['candidate_id'],role=cfg['role'],source_week=3,
                 job_id=cfg['job_id'],status='not_started',minimum_verified=False,config=cfgpath)
        if p.exists():
            m=json.loads(p.read_text());r=m.get('result',{})
            assert m['config_sha256']==sha(ROOT/cfgpath)
            row.update(status=m['status'],minimum_verified=r.get('minimum_verified',False),
                       energy_Hartree=r.get('energy'),elapsed_seconds=m.get('elapsed_seconds'),
                       geometry_audit=r.get('geometry_audit'),run_dir=m['run_dir'],
                       optimized_xyz=str(Path(m['run_dir'])/(Path(m['input']).stem+'.xyz')),
                       lowest_frequency_cm1=r.get('lowest_vibrational_frequency_cm1'),
                       error=r.get('week3_validation_error') or m.get('stop_reason'))
        else:
            running=sorted((ROOT/'runs').glob(f"*-{cfg['job_id']}/manifest.json"))
            if running:
                m=json.loads(running[-1].read_text());row.update(status=m['status'],run_dir=str(running[-1].parent))
                progress=running[-1].parent/'progress.json'
                if progress.exists():row['progress']=json.loads(progress.read_text())
        rows.append(row)
    scores={}
    for scope in ['primary_only','all_available_poses']:
        pool=[r for r in rows if r['minimum_verified'] and r['status']=='passed' and (scope!='primary_only' or r['role']=='primary')]
        best={}
        for r in pool:
            if r['species'] not in best or r['energy_Hartree']<best[r['species']]['energy_Hartree']:best[r['species']]=r
        records=[]
        for b in ['adenine','caffeine','theobromine','xanthine']:
            names=['diflunisal',b,f'diflunisal--{b}','diflunisal--diflunisal',f'{b}--{b}']
            rec=dict(coformer=b,B_kcal_mol=None,X_kcal_mol=None,
                     available_verified_structure_counts={name:sum(r['species']==name for r in pool) for name in names},
                     missing=[name for name in names if name not in best])
            if all(name in best for name in names[:3]):
                rec['B_kcal_mol']=(best[names[2]]['energy_Hartree']-best[names[0]]['energy_Hartree']-best[names[1]]['energy_Hartree'])*HARTREE_KCAL
            if all(name in best for name in names[2:]):
                rec['X_kcal_mol']=(best[names[2]]['energy_Hartree']-.5*best[names[3]]['energy_Hartree']-.5*best[names[4]]['energy_Hartree'])*HARTREE_KCAL
            rec['selected_job_ids']={name:best[name]['job_id'] for name in names if name in best}
            records.append(rec)
        scores[scope]=records
    new=[r for r in rows if r['source_week']==3]
    complete=all(r['status']=='passed' and r['minimum_verified'] for r in new)
    drug=dict(updated_at=now(),status='bounded_panel_completed' if complete else 'partial',rows=rows,scores=scores,
              complete_new_jobs=sum(r['status']=='passed' and r['minimum_verified'] for r in new),expected_new_jobs=len(new),
              scope='Isolated vacuum electronic energies; counts are verified structures, not polymorphs or independent observations',
              comparison_warning='All-available-pose coverage is unequal. Report it separately from the primary-only panel; no exhaustive ranking or cocrystal free-energy claim.',
              conservative_geometry_groups=geometry_groups(rows),
              grouping_scope='Same species and fixed all-atom mapping, proper rotation, RMSD <=0.05 Angstrom and energy gap <=0.01 kcal/mol. No symmetry permutation; may split equivalent minima. Not a polymorph count.')
    save(ROOT/'reports/week3-drug-comparison.json',drug)
    state_path=ROOT/'reports/week3-run-state.json';state=json.loads(state_path.read_text()) if state_path.exists() else {'status':'first_job_running_preparing_pipeline'}
    material_path=ROOT/'reports/week3-materials.json';mat=json.loads(material_path.read_text()) if material_path.exists() else {'status':'not_started'}
    finished=complete and mat['status']=='passed' and state['status']=='completed'
    summary=dict(updated_at=now(),status='completed' if finished else 'in_progress_or_requires_review',
                 drug=drug['status'],drug_complete_new_jobs=drug['complete_new_jobs'],drug_expected_new_jobs=len(new),
                 materials=mat['status'],controller=state,week3_complete=finished,
                 cocrystal_polymorphs_determined=False,week4_started=False)
    save(ROOT/'reports/week3-summary.json',summary)
    lines=['# 3주차 실험 진행 및 결과','',f"갱신: {summary['updated_at']}. 전체 상태: **{summary['status']}**.",'',
           f"새 약물 계산 {drug['complete_new_jobs']}/{len(new)}건 최소점 검증 완료. 반도체: {mat['status']}. 제어기: {state['status']}.",'',
           '## 문헌 근거','',
           '정식 SI를 확보하고 실험법·스펙트럼을 확인했다. CAF의 COOH–N9 접촉과 ADE의 N1 양성자 이동은 분광학적 해석이다. DIF 공결정 CIF와 다형 수, 미형성 시도별 조건은 확인되지 않았다. [방법과 한계](../docs/WEEK3_PROTOCOL.ko.md), [정식 SI](../data/week3/literature/ANIE-64-e202507835-s001.pdf).','',
           '## 약물 계산 상태','',
           '| 대상 | 시작점 | 역할 | 상태 | 최소점 |','|---|---|---|---|---|']
    for r in rows:lines.append(f"| {r['species']} | {r['candidate_id']} | {r['role']} | {r['status']} | {r['minimum_verified']} |")
    for scope in scores:
        lines += ['',f'### {scope}','', '| 공형성제 | B kcal/mol | X kcal/mol | 미완료 대조군 |','|---|---:|---:|---|']
        for r in scores[scope]:
            b='미계산' if r['B_kcal_mol'] is None else f"{r['B_kcal_mol']:.3f}"
            x='미계산' if r['X_kcal_mol'] is None else f"{r['X_kcal_mol']:.3f}"
            lines.append(f"| {r['coformer']} | {b} | {x} | {', '.join(r['missing']) or '없음'} |")
    lines += ['', 'B=E(AB)−E(A)−E(B), X=E(AB)−½E(AA)−½E(BB). 전자 에너지이며 ZPE·열·용매·결정 격자 자유에너지가 아니다. 기본 후보와 추가 자세를 분리하고, 확보된 자세 수가 다른 한계와 중복 가능성을 유지한다.','',
              '## 반도체','',f"현재 상태: {mat['status']}."]
    if mat.get('nine_point_grid'):
        lines += ['','| 상대 위치 | Mo 평면 거리 Å | 층간 에너지 meV/Å² | 상층 z 힘 eV/Å |','|---|---:|---:|---:|']
        for r in mat['nine_point_grid']:lines.append(f"| r{r['registry_id']} | {r['Mo_plane_distance_angstrom']:.4f} | {r['gamma_meV_angstrom2']:.5f} | {r['upper_layer_force_z_eV_angstrom']:.6f} |")
    lines += ['', '[반도체 원시 비교·정밀도 기록](week3-materials.json). 24원자 계면/외측 치환 비교는 이번 범위에 포함하지 않는다.','',
              '## 실행 기록','', '[현재 실행 상태](week3-run-state.json), [약물 상세 결과](week3-drug-comparison.json), [입력과 선택 규칙](../experiments/week3/drug/panel-plan.json).',
              '실패·시간/메모리 한도·정밀도 미달이면 제어기는 중단하고 검토 사유를 기록한다. 자동 재시도하지 않는다.']
    (ROOT/'reports/WEEK3_RESULTS.ko.md').write_text('\n'.join(lines)+'\n')
    return summary


if __name__=='__main__':
    print(json.dumps(summarize(),indent=2))
