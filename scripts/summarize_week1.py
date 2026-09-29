"""Build the week-one evidence index/report from completed, verified outputs."""
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
def read_json(path):
    return json.loads((ROOT / path).read_text())
def save_json(path, value):
    (ROOT / path).write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

verification = read_json('reports/week1-verification.json')
assert verification['status'] == 'passed'
ids = ('week1-water-opt-freq','week1-acetaminophen-opt-freq','week1-mos2-relax','week1-mos2-force-check')
jobs = [read_json(f'reports/{name}.json') for name in ids]
assert all(j['status']=='passed' for j in jobs)
checks = {item['species']: item for item in verification['checks']}
identities = read_json('reports/week1-molecule-identities.json')
literature = read_json('experiments/week1/literature-cases.json')
drug_seconds = sum(j['elapsed_seconds'] for j in jobs if j['engine']=='orca')
material_seconds = sum(j['elapsed_seconds'] for j in jobs if j['engine']=='qe')
total = drug_seconds + material_seconds
statuses = {
    'P00': 'partial: final pair labels verified; positive methods from author preprint; final SI and negative attempt conditions unresolved',
    'P01': 'complete for local serial ORCA/OpenMP QE; documented RSS visibility limits',
    'P02': 'complete: SCF/Opt/Freq and Hessian cross-check passed',
    'P03': 'complete: one conformer SCF/Opt/Freq and Hessian cross-check passed',
    'DR00': 'partial: five identities, 31 unreviewed tautomer candidates, 15 MMFF initial geometries; DFT state/pose selection pending week2',
    'SR00': 'week1 files/method frozen with explicit PBE-pseudo/C09 approximation; exact reference revision unverified',
    'SR01': 'week1 fixed-cell monolayer relaxation and targeted stricter SCF passed; full convergence/bilayer reference pending week2'
}
summary = {
    'updated_at': datetime.now(timezone.utc).isoformat(), 'protocol_id': 'dft-pilot-20260909-v1',
    'status': 'week1_core_calculations_complete_with_open_research_preconditions',
    'successful_engine_invocations': len(jobs), 'successful_systems': ['water','acetaminophen','MoS2'],
    'failed_before_engine_launch': ['runs/20260910T092726021436Z-week1-water-opt-freq'],
    'engine_failures_after_launch': 0, 'automatic_restarts': 0,
    'computing_wall_seconds': {'drug': drug_seconds,'semiconductor': material_seconds,'total':total},
    'budget_fraction_used': {'drug_80h':drug_seconds/(80*3600),'semiconductor_40h':material_seconds/(40*3600),'total_120h':total/(120*3600)},
    'timing_definition': 'Sum of engine-runner wall durations; excludes setup, downloads, analysis, and unlaunched failed wrapper',
    'tasks': statuses,
    'jobs': [{key:j[key] for key in ('job_id','run_dir','status','elapsed_seconds','output_sha256','result')} for j in jobs],
    'molecule_count': len(identities['molecules']),
    'tautomer_inventory_count': sum(m['tautomer_count'] for m in identities['molecules']),
    'MMFF_initial_geometry_count': sum(len(m['initial_geometries']) for m in identities['molecules']),
    'synthetic_parser_gate_tests': {'count':9,'status':'passed','log':'reports/week1-parser-tests.log','is_engine_validation':False},
    'independent_format_checks': 'reports/week1-verification.json',
    'experimental_evidence': 'experiments/week1/literature-cases.json',
    'week2_started': False, 'future_work_autoscheduled': False,
    'open_gates': ['final journal supplementary methods and negative attempt conditions',
                   'chemically reviewed tautomer/proton-transfer and dimer pose selection',
                   'cutoff/k/vacuum/lattice convergence and bilayer reference',
                   'exact reference pseudopotential revision for matched model benchmarking']}
save_json('reports/week1-summary.json', summary)
save_json('experiments/week1/status.json', {'week':1,'status':summary['status'],'tasks':statuses,
                                         'summary':'reports/week1-summary.json','report':'reports/WEEK1_RESULTS.ko.md'})

# Preserve the pre-execution plan before attaching live status. This JSON is not a scheduler.
archive = ROOT / 'experiments/archive/protocol-v1-before-week1.json'
archive.parent.mkdir(exist_ok=True)
protocol_path = ROOT / 'experiments/protocol-v1.json'
if not archive.exists():
    archive.write_bytes(protocol_path.read_bytes())
protocol = json.loads(protocol_path.read_text())
protocol.update(status=summary['status'], execution_enabled=False,
                execution_enabled_meaning='No automatic scheduling; week1 was run through explicit reviewed job configs',
                manual_execution_authorized_scope='week1 only', current_status_file='experiments/week1/status.json',
                pre_execution_snapshot=str(archive.relative_to(ROOT)))
for task in protocol['tasks']:
    if task['id'] in statuses:
        task['status'] = statuses[task['id']]
protocol['drug'].update(experimental_conditions_status=literature['overall_status'],
                        experimental_conditions_detail='experiments/week1/literature-cases.json',
                        molecular_identity_status='verified: five PubChem structures; DFT chemical states not yet frozen',
                        molecular_identity_detail='reports/week1-molecule-identities.json')
protocol['semiconductor'].update(week1_method_config='experiments/week1/mos2-relax.json',
                                 week1_pseudopotentials_frozen=True,
                                 matched_reference_pseudopotential_revision_verified=False)
save_json('experiments/protocol-v1.json',protocol)

rows=[]
names=['물 Opt+Freq','아세트아미노펜 Opt+Freq','MoS₂ 고정 셀 이완','MoS₂ 엄격 SCF 확인']
for name, job in zip(names,jobs):
    result=job['result']
    rms=job['max_single_child_rss_bytes']/1024**2
    result_text = f"양의 진동 {len(result['frequencies_cm1'])}개" if job['engine']=='orca' else '전자·좌표 수렴' if job['requested']['optimization'] else '동일 구조 힘 재확인'
    rows.append(f"| {name} | {job['elapsed_seconds']:.1f}초 | {result_text} | {rms:.1f} MiB |")

report = f'''# 1주차 실험 결과

실행일 2026-09-10. **핵심 계산 3종과 MoS₂ 추가 SCF 검사를 완료했다.** 물과 아세트아미노펜은 선택한 방법에서 국소 최소점임을 확인했고, MoS₂는 고정 셀의 원자 좌표 이완과 잔류 힘을 확인했다. 연구 결과를 해석하기 위한 문헌·화학적 상태·정밀도 검토는 일부 남아 있다.

| 실제 계산 | 소요 시간 | 확인 결과 | OS 최대 단일 자식 RSS |
|---|---:|---|---:|
{chr(10).join(rows)}

계산 엔진 점유 시간 합계는 **{total/60:.2f}분**이다. 약물 {drug_seconds/3600:.4f}시간/80시간, 반도체 {material_seconds/3600:.4f}시간/40시간을 사용했다. 준비·다운로드·분석 시간은 포함하지 않는다. 이 작은 계산의 시간을 큰 이합체·이중층에 그대로 외삽하지 않는다.

**계산에서 배운 것**

SCF는 주어진 원자 배치에서 전자 상태를 수렴시키는 과정이다. 구조 최적화는 원자를 움직여 힘을 줄인다. 진동 계산은 그 구조 주위의 에너지 곡률을 조사한다. 힘이 작아도 안장점일 수 있으므로 분자 최소점 판정에는 양의 진동 모드를 확인했다. 비선형 N원자 분자는 병진·회전 6개를 제외한 3N−6개 모드를 가지므로 물은 3개, 20원자 아세트아미노펜은 54개다.

- 물: r2SCAN-3c/진공, 에너지 {jobs[0]['result']['energy']:.12f} Ha. 진동수 1653.27, 3813.60, 3932.73 cm⁻¹. O–H {checks['water']['OH_bonds_angstrom'][0]:.4f} Å, H–O–H {checks['water']['HOH_angle_degrees']:.2f}°. 조화 진동 계산값이며 실험 스펙트럼과 직접 맞춘 값은 아니다.
- 아세트아미노펜: r2SCAN-3c/진공/TightOpt, 에너지 {jobs[1]['result']['energy']:.12f} Ha. 진동 54개 모두 양수, 최저 {checks['acetaminophen']['hessian_lowest_vibration_cm1']:.2f} cm⁻¹. **한 시작 배좌의 국소 최소점**이며 전역 최저 구조나 약효를 확인한 결과는 아니다.
- MoS₂: vdW-DF-C09 + PBE 생성 PseudoDojo NC 파일, 80/320 Ry, k=6×6×1. a=3.18 Å와 셀을 고정했다. 이완 후 S–S 두께 {checks['MoS2']['final_S_to_S_thickness_angstrom']:.4f} Å, 진공 간격 {checks['MoS2']['final_vacuum_gap_angstrom']:.4f} Å. 원시 전자 에너지는 {jobs[2]['result']['energy']:.8f} Ry다.

MoS₂ 출력의 SCF 힘 보정 경고를 따라 같은 좌표에서 수렴 기준을 1e-8 → 1e-10 Ry로 강화했다. 최대 힘은 {checks['MoS2-strict-SCF']['max_atomic_force_eV_angstrom']:.3g} eV/Å, 힘 변화는 {checks['MoS2-strict-SCF']['max_atomic_force_change_eV_angstrom']:.3g} eV/Å였다. 목표 0.01 eV/Å를 통과했다. **포논 안정성, 평형 격자, cutoff/k/진공 수렴은 아직 검증하지 않았다.** 다른 물질·함수의 절대 총에너지끼리 크기를 비교하지 않는다.

**연구 준비와 남은 판단**

| 항목 | 이번에 확보한 것 | 아직 필요한 것 |
|---|---|---|
| 문헌 | 최종 논문의 네 분자쌍 라벨; 저자 원고의 양성 제조 조건 | 최종 보충자료 대조, 미형성 시도의 용매·시간·몰비 등 |
| 약물 구조 | 5종 식별 검증, MMFF 초기 3D 구조 15개, 2D 호변이성질체 후보 31개 | 화학적 상태와 양성자 이동 검토, 이합체 자세·DFT 시작점 선정 |
| 반도체 | Mo/S/Se 공식 파일·해시, 함수 선택, MoS₂ 첫 이완 | 원 논문과 정확한 퍼텐셜 버전 대응, 수렴 연구, 이중층 기준 |

문헌에서 아데닌 조합의 양성자 이동 가능성을 확인해 이온쌍 후보를 다음 상태 검토에 포함했다. 성공 사례의 에탄올 조건을 네 사례 공통 조건으로 일반화하지 않는다. 본문은 최종 출판본에서 확인했고 세부 조건은 저자 원고를 사용했으며, 두 자료의 버전을 구분했다. [최종 논문](https://onlinelibrary.wiley.com/doi/10.1002/anie.202507835), [저자 원고 v3](https://doi.org/10.26434/chemrxiv-2024-vgvhk-v3)

2주차의 약물 작업은 상태·자세 검토 후 디플루니살 단량체와 첫 카페인 복합체의 비용 측정이다. 반도체 작업은 cutoff/k-grid/진공 수렴과 순수 이중층 대조군이다. 주차 2 이후의 계산은 이번 실행에 포함하지 않았다. 연구 가설 D1/S1의 지지·반박을 판단할 데이터는 아직 없다.

**검증 기록과 재현**

실제 ORCA 로그의 진동수를 `.hess` 파일과 대조했고 최종 XYZ의 원자 수·에너지·해시를 확인했다. QE는 ASE 판독 결과를 원시 Ry로 환산해 대조하고, 반올림된 화면 좌표 대신 전체 정밀도 XML로 셀 고정·좌표·힘·단위를 검증했다. cclib의 ORCA 6.1.1 판독 오류는 별도로 기록했다. 합성 출력 9개 검사는 실패를 성공으로 잘못 분류하지 않는지 확인하는 파서 검사이며 실제 계산 9회를 뜻하지 않는다.

처음 물 실행은 sandbox의 swap 조회 오류로 엔진을 시작하기 전에 중단됐다. 계측을 수정한 뒤 물을 실행했으며, 이후 엔진 실행 4건은 모두 통과했다. 초기 실패 디렉터리도 보존했다. 자동 재시작은 없었다.

ORCA 두 작업은 프로세스 트리 메모리 조회가 막혀 표의 최대 단일 자식 RSS만 확보했다. QE는 권한이 허용된 실행에서 트리도 샘플링했고 MoS₂ 이완의 피크는 {jobs[2]['peak_sampled_process_tree_rss_bytes']/1024**3:.3f} GiB였다. QE 두 작업 중 관측한 시스템 swap 증가량은 0바이트였다. 트리 샘플링은 0.5초 간격이라 순간 피크를 놓칠 수 있고, 이 값들은 전체 Mac 사용량이 아니다.

검증만 다시 수행하려면 프로젝트에서 다음을 실행한다. 엔진을 다시 돌리지 않는다.

```bash
source scripts/activate.sh
python -m unittest discover -s tests -v
python scripts/verify_week1.py
python scripts/summarize_week1.py
```

새 계산을 의도적으로 다시 실행하려면 `python scripts/research_run.py --config experiments/week1/<설정파일>.json`을 사용한다. 새 run 디렉터리가 만들어지고 이전 출력은 보존된다. QE의 프로세스 메모리 계측에는 macOS 실행 권한이 필요할 수 있다. 원자 파일 재다운로드 후에는 보관된 SHA-256과 비교한다.

- [수치·작업별 원시 기록](week1-summary.json)
- [독립 형식 검증](week1-verification.json)
- [방법 결정과 원문 접근 한계](../docs/WEEK1_METHOD_DECISIONS.ko.md)
- [실험 조건과 결측 필드](../experiments/week1/literature-cases.json)
- [5종 식별·초기 구조 기록](week1-molecule-identities.json)
- [최적화 구조 폴더](../structures/optimized/week1/)
- [실행 전 원 계획 스냅샷](../experiments/archive/protocol-v1-before-week1.json)
'''
(ROOT / 'reports/WEEK1_RESULTS.ko.md').write_text(report)

files = []
for directory in ('experiments/week1','structures/week1','structures/optimized/week1','data/pseudopotentials/pseudo-dojo-nc-sr-04-pbe-standard','data/literature/week1'):
    files.extend(p for p in (ROOT/directory).rglob('*') if p.is_file())
files += [ROOT/f'reports/{j}.json' for j in ids]
files += [ROOT/p for p in ('reports/week1-summary.json','reports/week1-verification.json','reports/week1-molecule-identities.json','reports/WEEK1_RESULTS.ko.md','docs/WEEK1_METHOD_DECISIONS.ko.md','scripts/research_run.py','scripts/verify_week1.py','scripts/summarize_week1.py','scripts/prepare_week1_molecules.py','scripts/prepare_week1_mos2.py','scripts/prepare_week1_force_check.py','scripts/fetch_week1_pseudos.py','tests/test_research_parser.py','reports/week1-parser-tests.log')]
files += [ROOT/p for p in ('README.md','docs/EXPERIMENT_PROTOCOL_V1.ko.md','experiments/protocol-v1.json','experiments/archive/protocol-v1-before-week1.json')]
for job in jobs:
    run=Path(job['run_dir'])
    for p in run.iterdir():
        if p.suffix in ('.inp','.in','.xyz','.hess','.py','.json','.log'):
            files.append(p)
    files.extend(run.glob('scratch/*.save/data-file-schema.xml'))
save_json('reports/week1-evidence-sha256.json', {'algorithm':'SHA-256','files':{str(p.relative_to(ROOT)):digest(p) for p in sorted(set(files))}})
print(json.dumps({'successful_jobs':len(jobs),'total_engine_wall_minutes':total/60,'report':'reports/WEEK1_RESULTS.ko.md','hashed_files':len(set(files))}))
