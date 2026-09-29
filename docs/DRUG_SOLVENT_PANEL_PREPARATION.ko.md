# 에탄올 고정 구조 단일점 패널 준비

2026-09-29. `scripts/prepare_completion_drug_solvent.py`로 공통 진공 최소점의 좌표를 유지한 에탄올 SMD 단일점 입력 생성 절차를 준비했다. **임시 자료의 입력 조립·검증 절차 8개 검사를 통과했으며 실제 생산 입력이나 용매 점수는 아직 없다.** [준비 기록](../reports/completion/drug-solvent-preparation-r1.json), [검사 로그](../reports/completion/drug-solvent-preparation-tests-r1.log).

## 이번 단계의 범위

원래 17개 출발점·14종을 유지하고 같은 r2SCAN-3c 진공 좌표에서 SMD 에탄올 EnGrad를 구한다. dense/strict 격자 각각 17건, 총 34개 단일점을 준비하도록 구성했다. 전자 적분 격자는 기존 AngularGrid 7/IntAcc 6.0과 AngularGrid 7/IntAcc 7.0·Unpruned·HGridReduced false다. VeryTightSCF, SMD18 false, DRACO false를 두 조건 모두에 명시한다.

에너지는 `orca_smd_total_energy_including_cds`로 표시한다. 실제 출력에서 최종 SCF 총에너지에 CDS가 이미 포함됐는지, 분산·gCP까지 합한 최종 값이 일치하는지 확인해야 한다. 진공의 `electronic_energy`와 다른 양으로 다루며 CDS를 한 번 더 더하지 않는다. [호환성 준비와 공식 출처](ORCA_METHOD_SOLVENT_COMPATIBILITY.ko.md).

이 단계는 좌표를 고정한 환경 비교다. 용매 안에서 최적화한 최소점이나 공결정 형성 자유에너지로 해석하지 않는다. 기존 주 방법·진공의 좌표 출처, 기본/추가 역할 및 자기 결합 대조군을 유지한다.

## 생산 입력을 만들기 전 필요한 근거

공통 진공 전자 패널의 실제 수치 검증과 물의 SMD 에탄올 Opt/Freq 상세 검토가 모두 있어야 한다. 두 문서의 통과 표기만 읽지 않고 원시 결과를 다시 확인한다.

- 공통 패널: 17개 최소점·Hessian·dense/strict 기울기 및 B/X 정밀도.
- 물 호환성: 실제 r2SCAN-3c·기저·에탄올·SMD 보정 분해, 동일 원자 대응·기울기와 양의 내부 진동 3개.
- 상세 검토: 실제 실행 부분을 인용한 용매 설정, CDS 회계, 해석적 Hessian과 기울기의 설명.

상세 검토 파일 `reports/completion/orca-smdethanol-compatibility-review-r1.json`은 아직 없다. 생성 시 상태 `smdethanol_runtime_reviewed`, 실제 호환성 계획·설정 해시 및 `audit_job`이 반환한 출처 해시를 기록한다. `composite_method`, `ethanol_smd_parameters`, `cds_energy_accounting`, `analytic_hessian_and_gradient` 각 항목에 판정·설명·실행 부분 인용을 남긴다. 인용의 존재 검사는 해석 자체를 자동 증명하지 않으므로 실제 출력 검토도 필요하다. 물의 통과가 큰 약물의 정확도나 비용을 보장하지 않는 범위도 명시한다.

17개 원래 이력과 좌표·후보·전하·스핀·역할을 모두 대조한 뒤, 좌표를 바이트 단위로 동일하게 복사하고 입력·의존성·부모 근거를 고정한다. 잘못된 용매, 다른 에너지 종류, CDS 중복 적용, 불완전하거나 비양수인 물 진동, 변경된 출처나 누락된 출발점은 진행을 막는다.

## 실행과 다음 단계

선행 검증 완료 후 다음 명령으로 준비하고 실행한다.

```sh
.venv/bin/python scripts/prepare_completion_drug_solvent.py
.venv/bin/python scripts/completion_batch.py --plan experiments/completion/smdethanol-at-uniform-r1.json --execute
```

각 작업은 ORCA 12 MPI·OpenMP 1, RSS 40 GiB, 30분 제한이며 한 엔진만 실행한다. 개별 시간 제한 합은 17시간이고 실측 예상 시간은 아니다. 기존 약물 파일럿 잔여 예산을 매 작업 전 확인하며 자동 재시도하지 않는다.

실제 출력 형식을 확인한 뒤 모든 SMD 에너지·기울기·좌표 대응과 dense/strict B/X 수치 변화를 검증해야 한다. 기본 자세만/추가 자세 포함, 진공 선택 자세 고정/용매 에너지 재선택, 동률 폭 0/0.5/1/2 kcal/mol의 구분을 유지한다. 단일점 집계·검증기가 실제 결과를 통과시킨 뒤 종류별 최저 시작점 1개씩 **최초 최대 14개 용매 Opt/Freq**를 준비한다. 가까운 대안이나 순위 역전이 있으면 기존 프로토콜에 따라 대안을 균형 있게 검토한다.

현재 코드에는 용매 최적화를 자동으로 활성화하는 경로가 없다. `solvent_geometry_optimization_enabled`는 false이며, 용매 최소점·보조 ZPE·대안 자세 및 최종 용매 순위 판정은 후속 단계다. 원래 D04의 범위를 고정 구조 비교만으로 축소하지 않는다.

## 검증한 범위

임시 폴더에서 원래 17개 좌표를 이용한 34개 입력 조립, 동일 좌표·원자 대응, 정확한 용매/격자 설정, 자원 제한, 잘못된 방법·용매·보정·Hessian 근거의 거부를 확인했다. 재실행은 동일 입력을 유지하고 변경된 의존성이나 비교 행을 거부한다. 시험의 선행 과학 판정은 명시적으로 대체했으므로 실제 공통 패널·SMD 엔진 검증을 의미하지 않는다. 실제 작업 폴더에서도 미완료 선행 조건 때문에 생산 입력이 생성되지 않음을 확인했다.

```sh
.venv/bin/python -m unittest discover -s tests -p 'test_completion_drug_solvent_preparation.py' -v
```
