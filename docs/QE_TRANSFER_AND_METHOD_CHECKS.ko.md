# QE 결과 전달과 대체 분산 방법의 실행 전 점검

2026-09-29. 이 기록은 완료한 과거 SCF 2건의 재검증과 새 검증 입력 3건의 준비를 구분한다. 새 3건은 아직 실행하지 않았다. 진행 중인 반도체 기준 배치를 변경하거나 다른 연구 엔진을 동시에 실행하지 않았다.

## 실제 출력 2건에서 확인한 사항

`scripts/completion_qe_bridge.py`가 보관된 MoS₂ 단층과 6원자 이중층 SCF 결과를 ASE의 에너지(eV)·힘(eV/Å)으로 전달한다. 입력·출력·실행 코드·의사퍼텐셜의 해시, 정상 종료, XML 완료 표시, SCF 오차를 확인한 후 다음 항목을 교차 검사한다.

- 입력과 XML의 함수·분산 처리, cutoff, k-grid/offset, SCF 임계값, 스핀과 점유 설정.
- 입력과 최종 XML의 원자 순서·셀·좌표, 의사퍼텐셜별 실제 원자가 전자 수와 전하.
- XML 에너지와 원시 Ry 출력, XML 힘과 원시 Ry/au 힘.
- 요청한 원자 좌표가 실제 계산 좌표와 같은지, 좌표를 바꿨을 때 이전 힘을 재사용하지 않는지.

단위는 설치한 QE 7.6의 `Modules/constants.f90`와 대조했다. 소스 커밋은 `9f93ddec427d2b9a45bb72d828c6d324f62fcabd`다. 사용한 변환은 1 Hartree = 27.211386245988034 eV, 1 Bohr = 0.529177210903 Å이다. XML의 cutoff/SCF 임계값/에너지는 Hartree 기준, 입력과 화면 에너지는 Ry 기준으로 구분한다. XML 힘은 Hartree/Bohr, 화면 힘은 Ry/au다.

6원자 이중층의 XML과 화면 힘 차이 최댓값은 약 1.269e-7 eV/Å로, 화면의 소수점 8자리 반올림 범위 이내였다. 입력/출력의 정합성 확인이며 에너지 유한차분이나 실제 제약 이완 검증을 대신하지 않는다.

실제 출력 2건 확인과 별도로, 잘못된 함수·cutoff·SCF 상태·전자 수·스핀·k-grid·스레드 메타데이터·힘 부호·좌표·미완료 XML을 거부하는 **9개 검사**를 통과했다. 기존 원본은 변경하지 않고 임시 복사본을 사용했다. [결과 JSON](../reports/completion/qe-bridge-archived-validation.json), [검사 로그](../reports/completion/qe-bridge-tests.log).

## 스레드 표시의 차이

과거 SCF와 현재 실행의 표준 출력은 `Serial multi-threaded version, running on 4 processor cores`를 표시한다. 현재 `pw.x` 프로세스에서도 OMP_NUM_THREADS=4, BLAS=1 설정과 여러 실행 스레드를 확인했다. 한 번의 1초 CPU 표본은 160.9%였으며, 이것을 평균 사용률이나 병렬 효율로 일반화하지 않는다.

한편 XML `parallel_info/nthreads`는 과거 결과에서 1이었다. 설치 소스 `Modules/qexsd.f90`의 해당 기록은 `__OMP` 조건부 분기에 의존한다. 그러므로 이 XML 필드만으로 실제 계산이 한 스레드였다고 판단하지 않고, 원시 로그·빌드 옵션·실행 프로세스 근거를 따로 남긴다. 실행 중인 바이너리나 소스는 고치지 않았다.

## 준비된 새 검증 3건

| 계산 | 변경 | 완료 후 판정 |
|---|---|---|
| C09 이중층 −변위 | 기존 6원자 r0/D=6.3 Å 구조에서 원자 1(S)의 z를 −0.002 Å 이동 | 중앙 에너지 미분의 한쪽 값 |
| C09 이중층 +변위 | 같은 원자의 z를 +0.002 Å 이동 | `−[E(+h)−E(−h)]/(2h)`와 기존 XML 힘의 차이 ≤1e-4 eV/Å |
| PBE-D3(BJ) 단층 | 기존 MoS₂ 좌표·PBE 의사퍼텐셜 유지, 함수/분산 처리 변경 | 실제 XML에 PBE, D3 version 4, three-body=true가 적용되고 힘/에너지 판독 검증 통과 |

세 계산 모두 cutoff 80/320 Ry, Γ 중심 6×6×1 k-grid, SCF 임계값 1e-12 Ry다. 변위 비교는 같은 셀·원자 순서를 사용한다. 기존 중심 힘은 SCF 1e-10 Ry 결과이므로 잔여 SCF 및 유한 변위 오차를 포함한 점검이다. 기준을 넘으면 변위 크기·전자 수렴·힘 전달을 검토하며 자동으로 기준을 완화하지 않는다.

QE 7.6의 실제 `PW/Doc/INPUT_PW.def`에서 `vdw_corr='grimme-d3'`, `dftd3_version=4`가 BJ damping임을 확인했다. three-body 항은 `.true.`로 명시한다. 현재 문서 확인과 입력 저장·재읽기까지만 검증되었으며, 이 세 옵션이 실제 엔진 출력에 나타나야 런타임 호환성 확인이 끝난다. C09와 D3를 동시에 적용하지 않는다.

각 계산은 30분/40 GiB 상한, QE OpenMP 4/BLAS 1, 동시 연구 엔진 1개 조건이다. 원래 대표 SCF 시간은 이중층 133.097초, 단층 76.467초였고, 새 계산의 더 엄격한 SCF 및 다른 함수 실행 시간은 아직 모른다. 3건 합계 최대 90분도 기존 반도체 예산에 포함한다.

20개 입력·코드·출처 파일 해시와 판정 기준을 [고정 계획](../experiments/completion/qe-force-and-method-checks-r1.json)에 저장했다. 현재 컨트롤러가 종료하고 예산이 확인된 뒤 실행한다.

```sh
.venv/bin/python scripts/completion_batch.py --plan experiments/completion/qe-force-and-method-checks-r1.json --execute
.venv/bin/python scripts/evaluate_completion_qe_checks.py
```

판독 코드는 9개 검사를 통과했지만, 후속 평가 스크립트는 문법 확인만 완료했고 실제 새 결과로는 아직 실행하지 않았다. 이 점검이 끝나도 CI/CO의 방법 민감도와 실제 QE 제약 이완 9개가 끝나는 것은 아니다. 각각 본 구조·전 궤적에 대한 검증이 별도로 필요하다.
