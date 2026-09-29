"""Close week two from actual engine outputs and independent verification records."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads((ROOT/p).read_text())
def save(p,x):(ROOT/p).write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

scan=read('reports/week2-monolayer-convergence.json')
combined=read('reports/week2-monolayer-combined.json')
bilayer=read('reports/week2-bilayer-analysis.json')
drug=read('reports/week2-drug-verification.json')
xtb=[read(f'reports/week2-xtb-{phase}.json') for phase in ('monomers','dimers')]
screen=read('reports/week2-xtb-analysis.json')
parallel_path=ROOT/'reports/week2-parallel-benchmark.json'
parallel=read('reports/week2-parallel-benchmark.json') if parallel_path.exists() else None
assert scan['completed']==7 and bilayer['completed_comparisons']==bilayer['expected_comparisons'] and drug['completed']==3
assert all(x['status']=='completed' for x in xtb)
jobs=[json.loads(p.read_text()) for p in sorted((ROOT/'runs').glob('*-w2-*/manifest.json'))]
assert all('ended_at' in j and 'elapsed_seconds' in j for j in jobs),'An engine is still running'
seconds={engine:sum(j['elapsed_seconds'] for j in jobs if j['engine']==engine) for engine in ('orca','qe')}
seconds['xtb']=sum(x['calculation_wall_seconds'] for x in xtb)
parallel_jobs=[json.loads(p.read_text()) for p in sorted((ROOT/'runs').glob('*-p2-*/manifest.json'))]
assert all('ended_at' in j for j in parallel_jobs),'A parallel benchmark engine is still running'
seconds['parallel_validation']=sum(j['elapsed_seconds'] for j in parallel_jobs)
seconds['total']=sum(seconds.values())
previous=read('reports/week1-summary.json')['computing_wall_seconds']
materials_passed=combined['status']=='passed' and bilayer['status']=='completed_checks_passed'
passed=drug['status']=='passed' and materials_passed
status='week2_core_checks_complete_with_open_research_preconditions' if passed else 'week2_executed_with_review_items'
pause=read('experiments/week2/pause-checkpoint.json')
user_paused=[j['job_id'] for j in jobs if str(Path(j['run_dir']).relative_to(ROOT))==pause['source_run_dir']]
engine_failures=[j['job_id'] for j in jobs if j['status']!='passed' and j['job_id'] not in user_paused]
xtb_failures=[r['candidate_id'] for x in xtb for r in x['results'] if r['result']['status']!='passed']
summary={'updated_at':datetime.now(timezone.utc).isoformat(),'status':status,'week':2,
         'engines':{'ORCA':'6.1.1','QE':'7.6','xTB':'6.7.1'},'job_counts':{'DFT_attempts':len(jobs)+len(parallel_jobs),'core_DFT_attempts':len(jobs),'parallel_validation_attempts':len(parallel_jobs),'xTB_attempts':sum(len(x['results']) for x in xtb)},
         'DFT_failed_or_review_job_ids':engine_failures,'xTB_failed_or_review_candidate_ids':xtb_failures,
         'user_interrupted_original_attempts':user_paused,
         'computing_wall_seconds':seconds,'timing_definition':'Sum of individually measured engine-wrapper wall durations, excluding preparation, analysis and downloads',
         'cumulative_week1_week2_hours':{'drug':(previous['drug']+seconds['orca']+seconds['xtb']+seconds['parallel_validation'])/3600,
                                       'semiconductor':(previous['semiconductor']+seconds['qe'])/3600},
         'monolayer_combined':combined,'bilayer_checks':{k:v for k,v in bilayer.items() if k!='rows'},
         'first_dimer_binding_energy_kcal_mol':drug['binding_energy_kcal_mol'],
         'drug_minima_verified':sum(r['minimum_verified'] for r in drug['checks']),
         'week3_engine_jobs_started':False,'future_work_autoscheduled':False,'automatic_retries':0,
         'parallel_validation_status':parallel['status'] if parallel else 'prepared_not_run',
         'parallel_failed_or_review_job_ids':[j['job_id'] for j in parallel_jobs if j['status']!='passed'],
         'open_gates':['final journal supplementary methods and negative-attempt conditions',
                       'final DFT state/pose coverage for the full control panel',
                       'self-association controls, method and solvent sensitivity',
                       'equilibrium interlayer distance beyond the initial rigid-layer samples',
                       'equilibrium in-plane lattice and heterobilayer/Se numerical convergence',
                       'exact reference pseudopotential revision for matched model benchmarking']}
summary['remaining_budget_hours']={name:max(0,budget-summary['cumulative_week1_week2_hours'][name])
                                   for name,budget in [('drug',80),('semiconductor',40)]}
if drug['status']=='passed':
    times={r['species']:r['elapsed_seconds'] for r in drug['checks']}
    mono=max(times['diflunisal'],times['caffeine']);dimer=times['diflunisal--caffeine']
    # Full panel: 5 monomer + 9 dimer species; already computed 2 monomers + 1 dimer.
    core=(3*mono+8*dimer)/3600
    capped=(13*mono+26*dimer)/3600
    summary['next_panel_cost_scenarios_hours']={'remaining_one_start_per_species':core,'remaining_core_plus_one_direct_ion_probe':core+dimer/3600,'remaining_up_to_three_starts_per_species':capped,
                                              'two_times_slowdown_remaining_three_starts':2*capped,
                                              'basis':'Observed monomer wall time and cumulative dimer effort including interrupted serial work and resumed parallel Freq; historical-cost planning scenario, not an uninterrupted same-core runtime prediction; SP/solvent/retries extra'}
if drug['status']=='passed':
    next_plan={'status':'plan_only_not_started','is_executable_scheduler':False,'execution_enabled':False,
               'remaining_primary_monomer_optimizations':3,'remaining_primary_dimer_optimizations':8,
               'additional_direct_DIF_ADE_proton_transfer_optimization':1,'remaining_optimization_start_count':12,
               'additional_frequency_calculation_upper_bound':12,
               'initial_MPI_processes':12,
               'early_cost_review':'Measure the first new uninterrupted 12-process Opt plus Freq before expanding the panel; the EnGrad speedup is not a full-job runtime prediction.',
               'frequency_policy':'Optimize first, then identify equivalent same-state/same-species minima before running redundant Hessians. Every distinct reported minimum still requires a successful frequency calculation.',
               'state_gate':'Select the exact direct proton-transfer start after reviewing acceptor-site chemistry; keep the original XYZ instead of forcing xTB relaxation first.',
               'further_pose_expansion':'Conditional, within the original maximum of 3 starts per species and the remaining 80 hour budget',
               'measured_cost_scenario_hours':summary['next_panel_cost_scenarios_hours']['remaining_core_plus_one_direct_ion_probe'],
               'remaining_drug_budget_hours':summary['remaining_budget_hours']['drug'],
               'three_starts_each_historical_scenario_within_remaining_budget':summary['next_panel_cost_scenarios_hours']['remaining_up_to_three_starts_per_species']<=summary['remaining_budget_hours']['drug'],
               'cost_scope':'Full Opt+Freq cost conservatively assigned even to the extra ion start; actual duplicate-Hessian savings not assumed'}
    save('experiments/week2/next-panel-plan.json',next_plan)
save('reports/week2-summary.json',summary)
state={'week':2,'status':status,'scope':'week2 only; week3 not launched','summary':'reports/week2-summary.json',
       'report':'reports/WEEK2_RESULTS.ko.md','open_gates':summary['open_gates']}
save('experiments/week2/status.json',state)
protocol=read('experiments/protocol-v1.json')
protocol.update(status=status,current_status_file='experiments/week2/status.json',manual_execution_authorized_scope='week1 and week2',
                execution_enabled=False,is_executable_scheduler=False)
task_status={'DR00':'week2: 93 monomer and 412 dimer starts screened with xTB; explicit state changes and provisional clusters recorded; full DFT state review remains',
             'DR01':'complete: first monomers/dimer and frequencies verified' if drug['status']=='passed' else 'executed: inspect week2 drug verification',
             'SR01':'week2 monolayer and pure-bilayer numerical checks passed at fixed a=3.18 A; equilibrium lattice and heterobilayer accuracy remain' if materials_passed else 'executed: inspect week2 numerical review items'}
for task in protocol['tasks']:
    if task['id'] in task_status:task['status']=task_status[task['id']]
protocol['drug'].update(week2_chemical_state_and_pose_inventory='experiments/week2/chemistry-candidates.json',
                        week2_first_DFT_selection='reports/week2-first-dft-selection.json')
protocol['semiconductor'].update(week2_convergence='reports/week2-monolayer-combined.json',week2_pure_bilayer='reports/week2-bilayer-analysis.json')
save('experiments/protocol-v1.json',protocol)

drug_rows=[]
for r in drug['checks']:
    freq=f"{r.get('vibrational_modes','—')}개, 최저 {r.get('lowest_vibrational_frequency_cm1',float('nan')):.2f} cm⁻¹"
    drug_rows.append(f"| {r['species']} | {r['elapsed_seconds']/60:.2f}분 | {freq} | {'최소점 확인' if r['minimum_verified'] else '검토 필요'} |")
screen_rows=[]
for name,s in screen['species'].items():
    screen_rows.append(f"| {name} | {s['started']} | {s['engine_passed']} | {s['geometry_usable']} | {s['cluster_count']} |")
bil_rows=[]
for r in bilayer['rows']:
    bil_rows.append(f"| {r['setting']} / r{r['registry_id']} | {r['Mo_plane_distance_angstrom']:.2f} | {r['gamma_meV_angstrom2']:.6f} |")
B=drug['binding_energy_kcal_mol'];Btext=f'{B:.3f} kcal/mol' if B is not None else '미확정'
r0_lowest=min((r for r in bilayer['rows'] if r['setting']=='selected' and r['registry_id']==0),key=lambda r:r['gamma_meV_angstrom2'])
c=combined['selected_setting'];estimate=summary.get('next_panel_cost_scenarios_hours',{})
cost_text=(f"종류별 1개로 남은 필수 대조군을 채우는 시나리오는 약 {estimate['remaining_one_start_per_species']:.1f}시간, 종류별 최대 3개까지 채우는 시나리오는 약 {estimate['remaining_up_to_three_starts_per_species']:.1f}시간이다. 후자가 2배 느려지면 {estimate['two_times_slowdown_remaining_three_starts']:.1f}시간이다. 남은 약물 예산은 {summary['remaining_budget_hours']['drug']:.1f}시간이며, 최대 3개 시나리오는 {'남은 예산 이내다' if estimate['remaining_up_to_three_starts_per_species']<=summary['remaining_budget_hours']['drug'] else '남은 예산을 초과한다'}. 첫 새 병렬 Opt+Freq의 실제 비용을 확인하기 전 전체 확장을 확정하지 않는다."
           if estimate else '최소점 검증 미완료로 후속 패널 비용 확정을 보류한다.')
outcome='2주차 핵심 계산과 검증을 완료했다' if passed else '2주차 실행을 마감했으며 추가 검토 항목이 남았다'
resumed_drug=next((r for r in drug['checks'] if 'resumed_frequency_wall_seconds' in r),None)
resume_note=(f"이합체는 87회 최적화 완료 후 사용자 요청으로 진동 계산 중 정지했다. 보존한 동일 좌표에서 {resumed_drug['frequency_MPI_processes']}프로세스로 Freq만 새로 실행했다. 위 이합체 시간에는 최초 시도 {resumed_drug['original_interrupted_wall_seconds']/60:.2f}분과 재개한 진동 계산 {resumed_drug['resumed_frequency_wall_seconds']/60:.2f}분이 모두 포함된다. 따라서 동일 코어 수로 연속 실행한 Opt+Freq 비용으로 해석하지 않는다."
             if resumed_drug else '이합체의 재개 진동 검증은 아직 완료되지 않았다.')
parallel_text='사용자의 코어 수 지적을 반영해 OpenMPI 4.1.1과 별도 실행기를 설치했다. 로컬 MPI 실행기 점검은 통과했으며 ORCA의 실제 병렬 에너지·힘 비교는 아직 미실행이다.'
if parallel and parallel['status']=='passed':
    times=parallel['median_wall_seconds_by_nprocs'];speed=parallel['speedup_vs_one']
    dimer_bench=[r for r in parallel['rows'] if r['nprocs'] in (1,8,12,16)]
    max_energy_delta=max(abs(r['energy_delta_reference_Hartree']) for r in dimer_bench)
    max_gradient_delta=max(r['max_gradient_delta_reference_Hartree_Bohr'] for r in dimer_bench)
    table='\n'.join(f"| {n} | {times[str(n)]:.3f} | {speed[str(n)]:.2f}배 |" for n in [1,8,12,16])
    parallel_text=f'''물의 2프로세스 Opt/Freq가 정상 완료됐고, 기존 직렬 결과와 에너지 1e-6 Hartree·진동수 0.1 cm⁻¹ 이내로 일치했다. 이어 동일한 50원자 이합체 좌표에서 EnGrad를 설정마다 두 번 실행해 에너지와 힘 차이 기준을 확인했다.

| MPI 프로세스 | 중앙 실행 시간(초) | 1프로세스 대비 속도 |
|---|---:|---:|
{table}

직렬 기준 최대 에너지 차이는 {max_energy_delta:.3g} Hartree, 최대 원자별 기울기 차이는 {max_gradient_delta:.3g} Hartree/Bohr로 사전 기준 1e-6 및 1e-5 이내였다.

시험한 병렬 설정 중 가장 빠르게 관측된 것은 **{parallel['fastest_observed_parallel_nprocs']}프로세스**다. 12프로세스가 가장 빠른 설정의 105% 이내에 들어가므로 재개에는 **12프로세스**를 선택했다. 두 반복의 결과이며 열·OS 변동과 다른 분자의 스케일링은 남는다. 전체 최적화 반복 수나 이합체 Hessian의 가속률을 측정한 것은 아니다. 기존 1코어 연구 결과와 실제 병렬 검증 비용을 분리해 보존했다.'''
report=f'''# 2주차 실제 실험 결과

작성 시각 {summary['updated_at']}. **{outcome}**. 초기 구조 생성, xTB 선별, 첫 약물 DFT 비용 측정, MoS₂ 정밀도 검사와 순수 이중층 대조를 실행했다. 전체 공결정 판정과 이종 이중층 본 실험은 다음 단계다.

## 약물: 후보 선별과 첫 최소점

31개 호변이성질체 후보 × 시드 3개로 단량체 93개, 9종 분자쌍으로 이합체 412개를 생성했다. 모든 초기 구조의 원자 수·순서·거리 기반 결합 연결을 검사했다. 이합체는 수소결합·방향족 접촉을 포함하며, 디플루니살–아데닌에는 양성자를 옮긴 초기 후보 12개를 포함했다. 모두 총 전하 0·단일항이다. 시작 형태는 상태 공간의 완전한 열거가 아니다.

| 종류 | xTB 시작 | 엔진 수렴 | 구조 사용 가능 | 보수적 구조 군집 |
|---|---:|---:|---:|---:|
{chr(10).join(screen_rows)}

엔진 수렴과 구조 사용 가능은 별도 판정이다. 무거운 원자 연결 변화·미배정 수소·분리는 제외 사유로 보존했고, 양성자 소유 원자 변화는 삭제하지 않고 따로 기록했다. 수소 위치는 원래 상태 그래프 대비 변화와 실제 입력 좌표 대비 순변화를 구분했다. 의도적으로 양성자를 옮긴 12개 시작점에서는 두 기준이 다르다. 거리로 추정한 수소 위치는 전자 전하나 결합차수의 검증이 아니다. 군집은 시작 상태·양성자 변화·기하학적 접촉을 나누고 무거운 원자 RMSD 0.75 Å로 묶었다. 원자 인덱스 접촉 서명은 대칭적으로 같은 접촉을 나눌 수 있어 군집 수를 실제 서로 다른 최소점 수로 해석하지 않는다. GFN2-xTB 선별에는 진동 계산이 없다. 양성자 이동 시작점 12개 중 수렴한 11개는 모두 원래 중성 상태 그래프와 같은 수소 소유 원자로 돌아왔다. 이는 진공 xTB의 기하 결과이며 결정의 양성자 상태를 판정하지 않는다. 다음 DFT 상태 검토에는 원래의 양성자 이동 XYZ도 남겨, xTB 이완만으로 이온쌍 가능성이 사라지지 않도록 한다.

![xTB 후보 수](figures/week2/xtb-screening-counts.png)

첫 DFT는 상태가 유지된 디플루니살·카페인 단량체와 사용 가능한 최저 xTB 에너지의 이합체를 골랐다. **r2SCAN-3c, TightSCF/TightOpt, 진공**이다. 최초 최적화는 직렬, 이합체의 재개 진동 계산은 검증된 MPI 환경을 사용했다. 각 분자쌍 안에서 자세를 생성하고 xTB 에너지로 정렬할 때 형성/미형성 라벨을 사용하지 않았다. 첫 비용 측정 분자쌍은 원 계획에 따라 알려진 양성 사례를 선택했으므로, 독립 블라인드 검증은 아니다.

| DFT 대상 | 실제 누적 계산 시간 | 병진·회전 제외 진동 | 판정 |
|---|---:|---|---|
{chr(10).join(drug_rows)}

{resume_note}

첫 이합체의 전자 결합 에너지 `B = E(AB) − E(A) − E(B)`는 **{Btext}**이다. 한 이합체 자세와 각 단량체의 최적화 에너지를 사용한 값이다. 전역 최소점·결정 형성 자유에너지·용해도·약효가 아니다. 자기 결합 `AA`, `BB` DFT 대조군을 아직 계산하지 않았으므로 경쟁 점수 `X = E(AB) − ½E(AA) − ½E(BB)`는 보고하지 않는다. ZPE·열·용매 보정도 이 수치에 포함하지 않는다.

## 반도체: 같은 구조에서 수치 정밀도 확인

1주차 내부 이완 좌표와 a=3.18 Å를 고정했다. vdW-DF-C09 + PseudoDojo nc-sr-04 PBE scalar NC, 전하밀도 cutoff=4×파동함수 cutoff, 전자 수렴 1e-10 Ry다. 원 논문과 정확한 의사퍼텐셜 파일 버전은 아직 맞추지 못했다.

첫 7개 SCF에서 cutoff 80/100/120 Ry, k 6/9/12, 실제 진공 간격 20/25/30 Å를 한 축씩 바꿨다. 120 Ry·k12·30 Å는 비교 기준이며 절대 정확한 값으로 가정하지 않는다. 단층 목표는 에너지 1 meV/atom, 최대 원자 힘 변화 0.005 eV/Å다.

선택한 조합은 **{c['cutoff']} Ry / k={c['k_grid']}×{c['k_grid']}×1 / 진공 {c['vacuum']} Å**이다. 세 축을 동시에 합친 재계산의 기준 대비 에너지 차이는 **{combined['energy_delta_meV_atom']:.6f} meV/atom**, 최대 힘 변화는 **{combined['max_force_delta_eV_angstrom']:.6g} eV/Å**였고, 조합 판정은 `{combined['status']}`다. 이 정밀도는 해당 고정 구조의 MoS₂에 한정한다. 허용오차를 만족하는 가장 작은 시험 설정을 후보로 골랐다. 실제 실행 시간은 SCF 반복 횟수의 영향도 받아 cutoff와 단조롭게 비례하지 않는다.

![단층 수렴](figures/week2/monolayer-convergence.png)

같은 방향의 R-type MoS₂/MoS₂ 이중층을 만들고 원시 셀 벡터로 r0=(0,0), r1=(1/3,2/3), r2=(2/3,1/3)을 비교했다. 각 이중층에서 동일 셀·동일 좌표·동일 계산 설정인 두 분리층을 빼고 면적으로 나눴다. 재사용은 동일한 입력을 다시 참조하는 경우만 허용했다.

| 설정 / 이동 | Mo 평면 거리 Å | 층간 에너지 meV/Å² |
|---|---:|---:|
{chr(10).join(bil_rows)}

대표 r0·D=6.3 Å의 두 설정 차이는 **{bilayer.get('representative_setting_difference_meV_angstrom2',float('nan')):.6g} meV/Å²**이다. 대칭상 동등한 r1/r2 차이는 **{bilayer.get('equivalent_registry_difference_meV_angstrom2',float('nan')):.6g} meV/Å²**다. 각각의 운영 기준은 0.1 meV/Å²다. 이 검사는 모든 거리·이동·Se 포함 조성에서 같은 오차 상한을 보장하지 않는다.

r0에서 가장 낮게 나온 표본 거리는 **{bilayer.get('lowest_sampled_r0_distance_angstrom')} Å**이며, 양쪽에 더 높은 에너지 표본이 있는지는 `{bilayer.get('r0_minimum_bracketed_by_samples')}`다. 그 표본에서 상층 전체의 z방향 힘은 **{r0_lowest['upper_layer_force_z_eV_angstrom']:+.6f} eV/Å**이다. 여기서 양의 힘은 층간 거리를 늘리는 방향이다. 음의 층간 에너지와 힘이 0인 평형은 별개이며, 이번 고정 층 거리 표본으로 평형 거리·면내 격자·포논 안정성을 확정하지 않는다.

![순수 이중층 대조](figures/week2/bilayer-controls.png)

## 계산 비용과 다음 단계

이번 엔진 점유 시간 합계는 **{seconds['total']/3600:.3f}시간**: xTB {seconds['xtb']/60:.2f}분, 핵심 ORCA {seconds['orca']/60:.2f}분, QE {seconds['qe']/60:.2f}분, 별도 병렬 검증 {seconds['parallel_validation']/60:.2f}분이다. 준비·다운로드·분석 시간은 제외하며 중단된 실제 계산 시간은 포함한다. 1주차를 합쳐 약물 **{summary['cumulative_week1_week2_hours']['drug']:.3f}/80시간**, 반도체 **{summary['cumulative_week1_week2_hours']['semiconductor']:.3f}/40시간**을 사용했다. 동시 엔진 1개, QE OpenMP 4/BLAS 1이며 ORCA의 직렬/병렬 비용은 위에서 구분했다. RSS는 프로세스 트리 샘플 값이며 전체 Mac 사용량·순간 최대치와 다르다.

{cost_text} 이는 단량체 두 개 중 느린 시간과 중단·병렬 재개를 포함한 이합체 누적 시간을 대입한 과거 비용 기반의 참고 시나리오다. 앞으로의 병렬 Opt+Freq 전체 실행 시간을 예측한 값이 아니며, 병렬 EnGrad의 가속률도 전체 Opt/Freq에 곱하지 않았다. 분자·배좌별 반복 횟수, 함수/용매 단일점, 재시도 비용은 별도다. 다음 약물 묶음은 남은 단량체 3종·이합체 8종과 아데닌의 직접 양성자 이동 시작점 1개, 총 12개 최적화 시작점으로 계획했다. 새 진동 계산은 최대 12개이며, 최적화 후 같은 상태·같은 종류의 동일 최소점이면 중복 Hessian을 줄인다. 구별되는 최소점마다 진동 검증은 유지한다. 그 뒤 상태·접촉 다양성과 비용을 보고 종류별 최대 3개까지 확장한다. 현재 군집별 최대 3개 목록은 상태 검토가 필요한 임시 목록이다.

## 병렬 실행 재검토

{parallel_text}

[코어 수 비교 계획과 설치·검증 범위](../docs/WEEK2_PARALLEL_REVIEW.ko.md)

MoSe₂와 이종 이중층에는 단층·층간 정밀도를 다시 확인해야 한다. 공통 면내 격자와 기준 층간 거리를 확정한 뒤 9점 기준 행렬로 진행한다. 24원자 계면 교환 27점 행렬은 그 이후의 조건부 단계다. **3주차 계산은 시작하지 않았다.**

최종 저널 보충자료와 미형성 시도별 조건은 미확정이다. 이를 확보하기 전 실험 실패 원인을 계산 에너지로 설명하지 않는다. 이번 계산은 D1/S1 가설의 지지·반박이나 학술적 신규성을 판정하는 데이터까지 도달하지 않았다.

## 증거와 재검증

핵심 DFT 검토 필요 작업: {engine_failures}. 사용자 요청으로 중단한 최초 시도: {user_paused}. xTB 검토 필요 후보: {len(xtb_failures)}개. 병렬 초기 검증 실패: {summary['parallel_failed_or_review_job_ids']}. 첫 물 MPI 시도는 라이브러리 검색 오류로 종료됐고, 작업 폴더에 필요한 라이브러리를 동봉한 별도 입력으로 수정 후 재검증했다. 모든 실패와 출력은 보존했고 자동 반복은 없었다. xTB는 91개 시점에서 3회 연속 미수렴으로 중단됐고, 누적 6개의 정상 SCC/구조 반복 한도 실패를 검토한 뒤 미실행 후보만 이어갔다. 이후 같은 실패 유형은 계속 제외하고 예상 밖 엔진/SCC 실패에 대한 중단 규칙을 유지했다. 최초 코드·중단 요약·검토 기록·후속 코드가 각각 보존돼 있다. ORCA 로그·Hessian·XYZ 에너지와 원자 순서를 교차 확인했다. QE는 XML의 전체 정밀도 좌표·셀·원자단위 에너지/힘으로 대조했다. 기존 파서/기하 검사 13개와 실제 Hessian을 활용한 재개 검증 회귀 검사 4개는 엔진 실행 횟수에 포함하지 않는다.

1주차 해시 115개 중 실제 계산 증거는 유지했다. 갱신한 README와 현재 프로토콜은 이전 바이트를 스냅샷으로 보존했다. 아래 명령은 계산 엔진을 시작하지 않고 저장된 결과를 다시 검사한다.

```bash
source scripts/activate.sh
python scripts/analyze_week2_monolayer.py
python scripts/check_week2_combined.py
python scripts/analyze_week2_xtb.py
python scripts/verify_week2_xtb.py
python scripts/verify_week2_drug.py
python scripts/verify_week2_bilayer_inputs.py
python scripts/analyze_week2_bilayers.py
python scripts/verify_week2_integrity.py
```

- [실행 계획과 사전 판정 기준](../docs/WEEK2_EXECUTION.ko.md)
- [수치·시간 요약](week2-summary.json)
- [후속 12개 시작점 계획, 미실행](../experiments/week2/next-panel-plan.json)
- [약물 독립 형식 검증](week2-drug-verification.json)
- [화학적 상태 검토와 이온쌍 후속 범위](../docs/WEEK2_STATE_REVIEW.ko.md)
- [xTB 전체 상태·군집 기록](week2-xtb-analysis.json)
- [첫 DFT 선택 근거](week2-first-dft-selection.json)
- [단층 수렴](week2-monolayer-convergence.json) · [조합 확인](week2-monolayer-combined.json)
- [이중층 좌표·기준층 대조](week2-bilayer-analysis.json)
- [1주차 보존 확인](week2-week1-preservation.json)
- [SHA-256 증거 목록](week2-evidence-sha256.json)

방법 근거: [ORCA 3c](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/3cmethods.html), [xTB 최적화](https://xtb-docs.readthedocs.io/en/latest/optimization.html), [층내/층간 분리 원 논문](https://arxiv.org/html/2503.15432v1). 수치 허용오차와 후보 개수는 이번 파일럿의 운영상 선택이다.
'''
(ROOT/'reports/WEEK2_RESULTS.ko.md').write_text(report)
(ROOT/'reports/WEEK2_PROGRESS.ko.md').write_text('# 2주차 진행 기록\n\n실행을 마감했다. [실제 결과와 남은 검토 항목](WEEK2_RESULTS.ko.md)을 참조한다.\n')
r=ROOT/'README.md';text=r.read_text()
text=text.replace('- [2주차 중단·재개 기록](reports/WEEK2_PAUSED.ko.md): 사용자 요청으로 정지, 최적화 좌표 보존',
                  '- [2주차 중단 당시 기록](reports/WEEK2_PAUSED.ko.md) · [재개 방법과 검증](docs/WEEK2_RESUME.ko.md): 87회 최적화 좌표 보존, 별도 12프로세스 진동 계산')
text=text.replace('- [2주차 실행 계획](docs/WEEK2_EXECUTION.ko.md) · [현재 진행 기록](reports/WEEK2_PROGRESS.ko.md)',
                  '- [2주차 실제 실험 결과](reports/WEEK2_RESULTS.ko.md): xTB 후보 선별, 첫 약물 DFT/진동, MoS₂ 수렴·순수 이중층 대조\n- [2주차 실행 계획](docs/WEEK2_EXECUTION.ko.md): 범위와 정밀도 판정 기준')
text=text.replace('Mo/S/Se 의사퍼텐셜은 확보했고 cutoff/k-grid/진공 수렴과 이중층 기준은 연구 프로토콜의 2주차 작업이다.',
                  'Mo/S/Se 의사퍼텐셜은 확보했다. MoS₂의 고정 구조 수렴과 순수 이중층 대조는 2주차 결과에 있고, MoSe₂/이종 이중층으로의 정밀도 일반화는 아직 검증하지 않았다.')
r.write_text(text)

paths=[]
for folder in ('experiments/week2','inputs/orca/week2','inputs/qe/week2','structures/week2','structures/optimized/week2','structures/checkpoints/week2','reports/figures/week2'):
    paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file())
paths.extend((ROOT/'reports').glob('week2-*.json'));paths.extend((ROOT/'reports').glob('week2-*.log'));paths.extend((ROOT/'reports').glob('w2-*.json'))
paths.extend((ROOT/'reports').glob('p2-*.json'))
paths.extend((ROOT/'scripts').glob('*week2*.py'));paths.extend((ROOT/'tests').glob('test_week2*.py'))
paths.extend(ROOT/p for p in ('docs/WEEK2_EXECUTION.ko.md','docs/WEEK2_STATE_REVIEW.ko.md','reports/WEEK2_RESULTS.ko.md','reports/week2-tests.log','README.md',
                              'docs/WEEK2_PARALLEL_REVIEW.ko.md','environment/orca-mpi-explicit.txt','scripts/orca_parallel.sh',
                              'docs/WEEK2_RESUME.ko.md',
                              'experiments/protocol-v1.json','experiments/archive/protocol-v1-before-week2.json','docs/archive/README-before-week2.md','runs/week2-finish-controller/runner_snapshot.py'))
paths.extend(p for p in (ROOT/'runs/week2-resume-20260911/pause_snapshot').rglob('*') if p.is_file())
paths.extend(p for p in (ROOT/'runs/week2-resume-20260911/prior-verifiers').rglob('*') if p.is_file())
paths.extend(p for p in (ROOT/'runs/week2-resume-20260911').glob('*snapshot.py'))
for j in [*jobs,*parallel_jobs]:
    run=Path(j['run_dir']);paths.extend(p for p in run.iterdir() if p.is_file() and p.suffix in ('.inp','.in','.xyz','.hess','.json','.log','.py'))
    paths.extend(run.glob('scratch/*.save/data-file-schema.xml'))
for phase in ('monomers','dimers'):
    root=ROOT/f'runs/week2-xtb-{phase}';paths.extend(p for p in root.iterdir() if p.is_file())
    for run in root.iterdir():
        if run.is_dir():paths.extend(p for p in run.iterdir() if p.is_file() and (p.suffix in ('.xyz','.json','.log') or p.name in ('charges','wbo','gradient','energy')))
index=ROOT/'reports/week2-evidence-sha256.json'
excluded={index,ROOT/'reports/week2-evidence-verification.json'}
paths=sorted(set(p for p in paths if p not in excluded))
save('reports/week2-evidence-sha256.json',{'algorithm':'SHA-256','scope':'compact inputs, outputs, Hessians, QE XML, structures, scripts and reports; large wavefunction binaries excluded',
                                        'files':{str(p.relative_to(ROOT)):sha(p) for p in paths}})
print(json.dumps({'status':status,'engine_hours':seconds['total']/3600,'report':'reports/WEEK2_RESULTS.ko.md','hashed_files':len(paths)}))
