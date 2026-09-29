# QE 4·8·12 OpenMP 스레드 비교

2026-09-29. 사용자의 계산 시간 단축 요청에 따라 실제 실행 가능한 6건의 비교 입력을 준비했다. **입력만 준비됐고 아직 실행하지 않았다.** 현재 반도체 기준 배치가 종료한 뒤 한 엔진 제한과 잔여 예산을 확인하고 실행한다.

## 고정한 비교

기존 `w2-bilayer-selected-r0-d6p3-pair`는 정상 종료하고 원시 XML·에너지·힘 검증을 통과한 6원자 MoS₂ 이중층 SCF다. 이 실행의 입력 바이트를 그대로 복사했다. 함수·의사퍼텐셜·구조·셀·80/320 Ry cutoff·6×6×1 k 격자·SCF 기준 1e−10 Ry와 실행 파일을 유지하고 OpenMP 수만 바꾼다. BLAS 관련 스레드 수는 1이다.

첫 반복은 4→8→12, 두 번째 반복은 12→8→4다. 각 실행은 새 run 디렉터리와 빈 scratch에서 시작한다. 입력의 기존 prefix도 바이트 일치를 위해 유지하지만 scratch 경로는 각 run 내부이므로 이전 파동함수를 재사용하거나 덮어쓰지 않는다.

개별 상한은 15분·RSS 40 GiB, 전체 6건 상한 합계는 90분이다. 상한은 실제 예상 시간이 아니다. 기존 자료의 같은 입력 4스레드 시간은 133.097초였지만 현재 실행의 소요 시간을 보장하지 않는다. 비교 비용도 반도체 파일럿 예산에서 차감한다. 자동 재시도는 하지 않는다.

## 결과 판정

- 각 계산의 정상 종료와 원시 XML·로그, 실제 OpenMP 헤더, 입력·실행 파일·의사퍼텐셜·좌표·힘 단위 및 해시를 다시 검증한다.
- 모든 새 실행은 기존 기준 대비 총에너지 차이 ≤1e−5 eV, 최대 원자별 힘 벡터 차이 ≤1e−5 eV/Å를 만족해야 한다. 수치 일치가 깨지면 성능 채택을 보류한다.
- 8 또는 12스레드는 대응하는 두 반복 모두에서 4스레드 대비 1.1배 이상 빠를 때만 후보가 된다. 후보들의 중앙 실행 시간을 비교해 이 사례에서의 권장값을 기록한다. 이 조건을 만족하는 후보가 없으면 4스레드를 유지한다.
- 한 구조·두 반복으로 시스템 열 상태나 OS 변동을 완전히 통제했다고 보지 않는다. 이 결과가 24원자 구조·격자 최적화·모든 k 격자의 최적값을 입증하지 않는다. 실제 후속 구조의 비용도 확인한다.

## 실행 절차

```sh
.venv/bin/python scripts/completion_qe_thread_benchmark.py prepare
# 현재 materials 제어기의 종료와 예산 검토 후에만 다음 명령을 실행한다.
.venv/bin/python scripts/completion_batch.py --plan experiments/completion/qe-thread-benchmark-r1.json --execute
.venv/bin/python scripts/completion_qe_thread_benchmark.py evaluate
```

계획은 [qe-thread-benchmark-r1.json](../experiments/completion/qe-thread-benchmark-r1.json)에 26개 입력·의존 자료 해시와 함께 고정했다. 입력 생성 재실행의 동일성과 6개 입력의 원본 바이트 일치를 확인했으며, 미실행 상태에서 `evaluate`가 결과 파일을 만들지 않고 거절하는 것을 확인했다. [준비 검증 기록](../reports/completion/qe-thread-benchmark-preparation-r1.json).

현재 생산 설정은 바뀌지 않았다. 평가 함수도 권장값만 기록하고 기존 입력이나 실행 중 프로세스를 수정하지 않는다. 실제 채택에는 새로운 설정 버전과 해시를 사용하고 과거 결과를 보존한다.
