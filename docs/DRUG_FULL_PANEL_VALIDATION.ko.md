# 공통 약물 패널의 원시 검증과 수치 민감도 집계

2026-09-29. `scripts/completion_drug_panel.py`는 공통 dense 조건에서 완료된 17개 시작점의 결과를 검증하고, 각 최종 좌표의 strict EnGrad와 실제 B/X 점수를 연결한다. **공통 패널 및 strict 계산은 아직 미완료이며, 실제 네 사례의 B/X 결과 파일은 생성하지 않았다.**

## 원시 결과를 읽는 방법

기존 실행기의 `parse_output`은 과거 `result.json`을 다시 쓸 수 있으므로 사용하지 않는다. 새 검증기는 저장된 입력·자산·실행기·설정·로그 해시와 실제 정상 종료·SCF·MPI 12개·ORCA 6.1.1·첫 SCF 격자를 직접 읽는다. 공통 패널에는 사전에 정한 입력 문자열과 정확히 일치하는 결과만 허용한다.

최적화 결과에서는 최종 XYZ의 원자 대응·화학적 기하를 점검하고, 로그·`.engrad` 에너지와 기울기 좌표를 대조한다. 전체 3N−6 내부 진동을 원시 Hessian에서 읽어 모두 양수인지, Hessian 좌표가 최종 구조와 같은지 확인한다. 회전·병진만 허용하는 기존 좌표 대응을 재사용한다.

양성자 소유권 변화는 삭제하거나 전자적 전하로 단정하지 않고 기하 상태로 기록한다. 무거운 원자 연결이 바뀌거나 소유권이 불명확한 수소·분리된 분자 배치는 검토 대상으로 중단한다. 기존 아데닌 직접 양성자 이동 시작점의 이력도 유지한다.

## 전체 패널 연결

`uniform-vacuum-r4.json`의 파일 해시와 원래 시작점 이력을 먼저 확인한다. 시작점 17개와 종류·기본/추가 역할의 대응을 유지하고, 실제 원시 후보 그래프까지 고정된 이력과 대조한다. 과거 설정에서 검증된 최소점을 새 dense 결과로 대체 사용하지 않는다.

전체 dense 최소점이 검증되면 대표 두 구조의 기존 strict 작업을 재사용하고, 나머지 **15개 strict EnGrad 입력**을 생성한다. 모든 입력은 해당 dense 최종 XYZ와 정확히 같은 파일에서 출발한다. 코드·입력·부모 원시 자료의 해시를 포함한 별도 체크 계획을 고정하며, 준비 도구는 ORCA를 실행하지 않는다.

실제 전체 최소점이 아직 없으므로 이 입력 생성과 전체 원시 결과 집계 경로의 실제 자료 통합 검증은 후속 계산이 끝난 뒤 해야 한다.

## 점수와 통과 조건

전자 에너지 B/X를 다음 축으로 분리한다.

- 기본 시작점만 / 선택된 추가 자세까지 포함.
- dense 결과 / 같은 구조에서 구한 strict 결과. 수치 프로토콜 이름도 따로 기록한다.
- dense에서 선택한 자세를 고정한 변화 / strict 에너지로 최저 자세를 다시 선택한 변화.

두 선택 방식의 B/X 변화 모두 기존 **0.1 kcal/mol** 기준을 확인한다. 자세가 바뀌면서 큰 수치 차이가 가려지는 경우를 통과시키지 않는다. 각 시작점의 같은 좌표에서 기울기 성분 차이는 **1e−4 Eh/bohr** 이하여야 한다. 누락된 자기 결합 대조군이 있으면 점수는 null이고 전체 통과는 불가능하다. 동률 폭 0/0.5/1/2 kcal/mol의 기존 비교표도 각각 유지한다.

통과 상태 `uniform_vacuum_numerical_panel_passed`는 선택된 진공 자세의 수치 검증만 뜻한다. wB97X-3c, 에탄올, ZPE 보조 분석, 추가 사례 독립 검증과 결정 구조·실제 형성 확률은 이 보고서의 완료 범위가 아니다.

## 검증한 범위

8개 검사를 통과했다. 실제 과거 DIF–DIF p05 정밀 재최적화를 읽어 에너지 −1850.964363429242 Eh, 내부 진동 150개와 최저 +9.18442963 cm⁻¹를 다시 확인했고 원본 보고서·manifest·result·출력 파일 해시는 그대로였다. 동시에 이 결과를 새 dense 패널 입력으로 승인하려는 요청은 거부했다.

나머지 검사는 변조된 출력 해시 거부, 미완료 패널에서 파일 생성 차단, 알려진 0.08 kcal/mol 변화, 자세 변경에 가려질 수 있는 0.2 kcal/mol 실패, 누락된 자기 결합·혼합 수치 조건 및 시작점 대응 오류를 확인한다. 합성 에너지 검사는 실제 네 사례의 계산 결과가 아니다.

```sh
.venv/bin/python -m unittest discover -s tests -p 'test_completion_drug_panel.py' -v
# 전체 dense 최소점이 완료된 후:
.venv/bin/python scripts/completion_drug_panel.py prepare-checks
# 생성된 uniform-vacuum-r4-checks.json의 계산이 완료된 후:
.venv/bin/python scripts/completion_drug_panel.py evaluate
```

[검사 로그](../reports/completion/drug-panel-tests.log), [준비 증거](../reports/completion/drug-panel-preparation-r1.json), [산술 규칙](DRUG_SCORE_ASSEMBLY.ko.md), [공통 패널의 시작점](DRUG_UNIFORM_PANEL_PREPARATION.ko.md).
