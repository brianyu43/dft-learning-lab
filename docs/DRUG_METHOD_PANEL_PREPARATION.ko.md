# 같은 약물 좌표에서 wB97X-3c 방법 민감도 비교 준비

2026-09-29. 입력 생성 도구와 8개 검사를 준비했다. **실제 방법 패널 입력과 계산 결과는 아직 없다.** 공통 진공 패널의 정밀도 검증과 실제 wB97X-3c 물 계산의 상세 검토가 모두 끝나야 생산 입력을 생성한다. [준비 검증 기록](../reports/completion/drug-method-preparation-r1.json).

## 비교의 의미와 범위

이미 고정한 17개 시작점, 14종의 분자·분자쌍을 그대로 추적한다. 각 시작점에서 검증된 r2SCAN-3c 공통 dense 최소점의 좌표를 바꾸지 않고 wB97X-3c EnGrad로 에너지와 기울기를 구한다. 이는 주 방법의 구조에서 방법을 바꾸는 비교이며 wB97X-3c로 다시 최적화한 최소점 비교가 아니다. 새 구조 탐색이나 독립 표본 수에 추가하지 않는다. 기본 14개와 대안 3개의 역할도 바꾸지 않는다.

같은 좌표에 dense와 strict를 각각 적용해 17+17개 단일점을 계획한다. r2SCAN-3c에서 통과한 적분 격자가 다른 방법의 수치 정밀도도 보장한다고 가정하지 않는다. 두 조건은 기존 AngularGrid 7/IntAcc 6.0과 AngularGrid 7/IntAcc 7.0·Unpruned·HGridReduced false이며 모두 VeryTightSCF다. 방법 전용 `wB97X-3c` 키워드 외에 기저·분산·gCP를 임의로 추가하지 않는다. [방법 검토의 출처와 호환성 준비](ORCA_METHOD_SOLVENT_COMPATIBILITY.ko.md).

## 입력 생성 전 검증

`scripts/prepare_completion_drug_method.py`는 다음을 확인한다.

1. 공통 진공 점수 보고서와 상세 방법 검토 기록이 실제로 존재해야 한다.
2. 물 호환성 계산의 원시 파일을 다시 읽어 정상 종료·실제 방법·기저·진공·MPI·기울기와 해시를 재검증한다.
3. 공통 진공 검증기로 17개 최소점과 Hessian, 모든 strict 기울기 및 B/X 수치 차이를 다시 확인한다. `passed` 문자열만으로 진행하지 않는다.
4. 고정된 원래 이력과 종·기본/추가 역할·후보 그래프를 대조하고, 좌표 파일의 해시·원자 순서·유한성·화학적 연결을 확인한다.
5. 원래 좌표 파일을 바이트 단위로 동일하게 복사하고 모든 의존성·부모 결과를 해시로 고정한다. 기존 결과나 입력을 덮어쓰지 않는다.

상세 방법 검토는 `reports/completion/orca-wb97x3c-compatibility-review-r1.json`에 남길 예정이다. 이 파일은 현재 존재하지 않는다. 요구 형식은 상태 `wb97x3c_runtime_reviewed`, 실제 호환성 계획·설정 해시, `audit_job`이 다시 확인한 `source_sha256`, 그리고 `composite_method`, `basis_ecp`, `dispersion_parameters`, `energy_components` 각 항목의 판정·설명·실제 실행 부분 인용이다. 프로그램 저자 목록이나 입력 재출력은 인용 근거로 받지 않는다. 인용의 존재 검사는 설명의 과학적 타당성을 자동 증명하지 않으므로 원시 출력에 대한 검토도 필요하다.

물만으로 약물의 모든 원소·ECP 설정, 정확도나 비용을 확인했다고 간주하지 않는다. 검토 기록의 `water_only_not_drug_accuracy_or_cost`를 명시하고, 실제 약물 출력에서도 방법·기저·ECP·D4의 적용을 확인해야 한다. SMD 결과를 진공 방법 검증에 대체 사용하지 않는다.

## 실행과 후속 판정

선행 검증이 완료되면 다음 순서로 진행한다. 입력 생성만으로 엔진이 실행되지는 않는다.

```sh
.venv/bin/python scripts/prepare_completion_drug_method.py
.venv/bin/python scripts/completion_batch.py --plan experiments/completion/wb97x3c-at-uniform-r1.json --execute
```

각 작업은 ORCA 12 MPI·OpenMP 1, RSS 40 GiB, 30분 제한이다. 34개 개별 시간 제한의 합은 17시간이며 실측 예상 시간이나 전체 예산 확보를 뜻하지 않는다. 기존 제어기는 매 작업 전 약물 파일럿 잔여 예산과 한 엔진 제한을 적용한다. 실패 시 자동 재시도하지 않는다.

후속 결과 검증기는 실제 새 방법의 출력 형식을 확인한 뒤 연결한다. 각 점의 원시 방법·에너지·좌표 대응을 확인하고, dense/strict B와 X 변화 ≤0.1 kcal/mol 및 기울기 성분 차이 ≤1e−4 Eh/bohr를 점검한다. 기울기 차이는 수치 민감도이며 wB97X-3c에서 힘 자체가 0이라는 판정이 아니다.

방법 비교는 기본 자세만 사용한 범위와 모든 선택 자세 범위를 각각 유지한다. r2SCAN-3c의 선택 자세를 고정한 점수와 wB97X-3c 에너지로 다시 선택한 점수를 함께 보고하고, 자세 교체로 큰 차이가 가려지지 않게 한다. 동률 폭 0/0.5/1/2 kcal/mol도 모두 표시한다. 방법 간 절대 에너지를 섞지 않으며 ZPE·용매·형성 자유에너지와 별도 분석한다.

## 현재 검증의 한계

[8개 검사 로그](../reports/completion/drug-method-preparation-tests-r1.log)는 임시 폴더의 명시적 시험 자료로 17×2 입력 조립, 동일 좌표, 대조군·역할 보존, 자원 제한, 누락·변경된 근거 차단 및 재실행의 동일성을 확인한 것이다. 시험에 사용한 기존 좌표는 새 공통 최소점으로 인증하지 않았고 시험용 실행 파일은 실행 불가능한 표시 문자열이다. 실제 작업 폴더에서도 선행 자료가 없으면 아무 생산 입력도 쓰지 않는 것을 확인했다. 새 ORCA 실행은 0건이다.
