# MACE CPU·Apple GPU 비교 실행 준비

2026-09-29. [실행 코드](../scripts/completion_mace_benchmark.py)와 [고정 입력](../experiments/completion/mace-device-benchmark-r1.json)을 준비했다. **실제 모델 추론과 가속률 측정은 아직 하지 않았다.** 현재 QE 기준 배치가 종료한 뒤 실행한다. 이번 준비는 MACE 비교에 관한 것이며 QE/ORCA 계산을 GPU로 전환한 결과가 아니다.

## 동일 조건 비교

공식 MACE-MP-0b3 medium의 동일 가중치를 사용한다. CPU float64·float32는 각각 PyTorch 4·8스레드, MPS float32는 CPU 준비 작업 4스레드로 비교한다. BLAS 관련 환경 변수는 1이고 PyTorch 실제 스레드 수도 기록한다. 첫 라운드는 CPU64/4→CPU32/4→CPU64/8→CPU32/8→MPS32/4, 두 번째는 역순이다. 모델은 각 작업에서 CPU로 읽은 뒤 필요한 자료형으로 변환하여 해당 장치에 올린다.

구조는 정상 종료한 과거 6원자 MoS₂ 이중층 SCF의 원시 입력·출력·힘·XML을 다시 확인한 뒤 만들었다. 위층 S를 Se로 치환한 6원자 구조와 이를 2×2로 늘리고 안쪽 S/Se 한 쌍을 교환한 24원자 구조다. 각각 Mo₂S₂Se₂와 Mo₈S₈Se₈이며, 두 분리층은 쌍과 **동일한 셀·좌표·주기 조건의 정확한 부분집합**이다. 이 구조들은 성능 비교용으로 명시한 미이완 구조이며, 현재 진행 중인 공통 격자 최적화나 24원자 본 실험의 결과로 세지 않는다.

각 장치·자료형·스레드 조합에서 6/24원자 쌍과 두 분리층, 총 6개 구조를 계산한다. 구조마다 준비 계산 2회와 측정 5회를 수행한다. 10개 작업 전체는 준비 계산 120회와 측정 300회다. `ASE Calculator.calculate`를 직접 호출하므로 이전 에너지 캐시를 읽는 시간을 성능으로 세지 않는다. 시간에는 이웃 그래프 생성과 계산기가 기본적으로 구하는 응력도 포함된다. 모델 import·읽기 시간과 각 작업 전체 시간은 따로 기록한다. GPU 작업은 전후 동기화하고, 실제 forward의 에너지·힘 텐서가 MPS/float32인지 검사한다. CPU 자동 대체와 MPS fast math는 끈다.

## 정밀도와 선택 기준

CPU float64/4스레드를 기준으로 모든 측정 반복의 다음 차이를 확인한다.

- 원자당 총에너지 차이 ≤1 meV.
- 원자별 힘 벡터 차이의 최댓값 ≤0.005 eV/Å.
- `(E쌍 − E아래층 − E위층) / 면적` 차이 ≤0.1 meV/Å².

CPU float32와 MPS float32의 차이도 별도로 남겨 자료형 변화와 장치 차이를 구별한다. 기준을 통과한 CPU 조합 중 중앙 실행 시간이 가장 짧은 조합과 GPU를 비교한다. GPU는 두 라운드 모두 1.1배 이상 빠르고 정밀도 기준도 통과할 때만 해당 크기의 후보로 추천한다. 쌍과 두 분리층의 중앙 실행 시간을 합산하여 비교하며, 준비·모델 읽기를 제외한 반복 사용 기준이라는 한계를 명시한다. 추천 결과는 기존 생산 설정을 자동 수정하지 않는다.

분산 보정은 이 장치 비교에 추가하지 않는다. 따라서 이 비교를 통과해도 C09 DFT 정확도나 PBE-D3 일치, 계면 모델 적합성, 장거리 분산의 검증이 되지 않는다. 실제 과학 비교는 [모델 준비도 검토](../reports/completion/MATERIALS_FOUNDATION_SCREEN.ko.md)의 별도 조건을 충족해야 한다. 여기서는 단일 구조 ASE 경로를 비교하며 묶음 추론의 가속률은 주장하지 않는다.

## 실행 제한과 검증

기존 제어기·연구 실행 잠금 두 개를 모두 확보한 뒤에만 모델 프로세스를 시작한다. 작업마다 상한은 180초·CPU 작업 프로세스 RSS 40 GiB, 전체 10개 작업의 실행 상한 합계는 30분이다. GPU는 PyTorch 권장 작업 메모리의 50%로 할당 한도를 설정한다. CPU RSS는 GPU 전체 메모리 사용량과 같은 수치가 아니다. 종료·실패·시간 초과 시 기록을 보존하고 자동 재시도하지 않는다. 단독 worker 명령도 제어기가 넘긴 잠금 없이 실행을 시작하지 못한다.

해석 가능한 인공값으로 단위 변환, 에너지 상쇄, 총에너지는 통과하지만 상호작용 에너지는 실패하는 경우, 힘 벡터 기준, 불완전 자료 거절, 동시 실행 차단을 확인한 7개 검사가 통과했다. 현재 실제 QE 잠금에 대한 실행 거절도 확인했다. GPU 비교 run 폴더는 생성되지 않았고 모델 추론은 시작하지 않았다. 13개 출처 해시와 입력 계획 해시를 [준비 검증 기록](../reports/completion/mace-device-benchmark-preparation-r1.json)에 보존한다.

```sh
.venv/bin/python scripts/completion_mace_benchmark.py prepare
# QE/ORCA 배치가 끝나고, 실제 MPS 장치 접근이 가능한 실행 환경에서 진행한다.
.venv/bin/python scripts/completion_mace_benchmark.py run
.venv/bin/python scripts/completion_mace_benchmark.py evaluate \
  --output reports/completion/mace-device-benchmark-r1.json
```

이 모델 비교의 실제 소요 시간은 별도 manifest로 기록한다. 현재 `completion_batch.spent()`가 합산하는 QE/ORCA 실행 시간에는 자동 포함되지 않으므로 최종 비용표에서 모델 비용을 별도 합산해야 한다.
