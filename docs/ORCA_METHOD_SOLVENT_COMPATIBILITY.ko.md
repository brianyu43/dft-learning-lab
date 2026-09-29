# 약물 방법·용매 계산의 실행 검증 준비

2026-09-29. D03/D04의 실제 약물 비교 전에 설치된 ORCA 6.1.1이 필요한 설정과 출력을 처리하는지 확인한다. 현재는 입력 준비와 검증기 점검만 완료했으며 두 계산은 실행 전이다. [준비 기록](../reports/completion/orca-compatibility-preparation-r1.json)에 28개 파일 해시, 원본 물 구조의 연결, 반복 준비의 동일성 및 5개 검사 결과를 남겼다.

| 계산 | 구조 | 확인할 출력 | 최대 시간 |
|---|---|---|---|
| wB97X-3c / 진공 EnGrad | 1주차 검증된 물 최소점의 같은 좌표 | 실제 복합 방법·기저·ECP·D4, 에너지와 기울기 좌표 대응 | 30분 |
| r2SCAN-3c / SMD 에탄올 Opt+Freq | 같은 물 좌표에서 시작 | 실제 용매, 보정 항, 최적화, 해석적 Hessian 및 양의 내부 진동 3개 | 30분 |

두 작업 모두 조밀 적분 격자 AngularGrid 7/IntAcc 6.0, VeryTightSCF, 12 MPI·OpenMP 1, 프로세스 트리 RSS 40 GiB 제한을 사용한다. 두 번째는 VeryTightOpt이며 SMD18과 DRACO를 끈다. 물의 실행 시간이나 메모리를 큰 약물 이합체 비용 추정에 그대로 사용하지 않는다. 이 두 작업은 파일럿 B/X나 독립 표본 수에 포함하지 않는다.

## 방법과 에너지의 구별

wB97X-3c의 전용 키워드는 vDZP/ECP와 조정된 D4를 포함하는 복합 방법을 지정한다. 다른 기저나 gCP를 따로 붙이지 않는다. 실제 적용 여부는 입력 문자열 외에 실행 출력에서도 확인한다. [ORCA 복합 방법 설명](https://www.faccts.de/docs/orca/6.1/manual/contents/modelchemistries/3cmethods.html)

SMD의 CDS 항은 출력의 SCF 총에너지에 이미 포함된다. r2SCAN-3c의 최종 에너지는 여기에 출력된 분산 및 gCP 항이 반영된 값인지 확인하며 CDS를 다시 더하지 않는다. 이 양은 `orca_smd_total_energy_including_cds`로 구분해 저장한다. 진공 전자 에너지와 같은 비교 묶음에 섞거나 공결정 형성 자유에너지로 부르지 않는다. [ORCA 용매 설명](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/solvationmodels.html)

## 실행과 판정

현재 반도체 제어기가 종료되고 실행 잠금이 풀린 후, 잔여 약물 예산을 확인하여 실행한다. 동시 실행이나 자동 재시도는 하지 않는다.

```sh
.venv/bin/python scripts/completion_batch.py --plan experiments/completion/orca-method-solvent-compatibility-r1.json --execute
.venv/bin/python scripts/evaluate_completion_orca_compatibility.py
```

검증기는 원래 보고서를 수정하지 않는다. 해시·설정·정상 종료·전자 수렴·MPI·실제 적분 격자·원자 순서·기울기 좌표·Hessian과 진동을 재확인한다. 방법명을 포함하는 프로그램 저자 목록이나 입력 재출력만으로 방법 적용을 인정하지 않는다. SMD 에너지의 마지막 단계 분해가 맞아야 하며 앞선 단계가 정상이라는 이유로 마지막 오류를 넘기지 않는다.

현재 5개 검사는 기존 r2SCAN 물 출력과 명시적으로 만든 가상 SMD 숫자를 사용했다. 실제 wB97X-3c 및 SMD 출력 형식과 엔진 호환성을 확인한 검사는 아니다. 새 출력의 형식이 예상과 다르면 원시 자료를 확인하고 검증기를 검토한다. 기준을 조용히 생략하지 않는다.

자동 검증 결과만으로 전체 약물 패널을 활성화하지 않는다. 실제 ECP/D4 매개변수, 에탄올 설정 및 해석적 Hessian을 원시 출력에서 추가 검토한 뒤 호환성 판정을 별도로 남긴다. 실제 D03/D04는 공통 진공 패널의 정밀도 통과 후 모든 비교 성분에 대한 방법 단일점, 에탄올 단일점, 종류별 최적화·진동 및 필요한 대안 자세까지 완료해야 한다.
