# Apple GPU 활용 확인과 후속 적용

2026-09-29. 사용자는 가능한 단계에 Mac GPU를 적극적으로 사용해 계산 시간을 줄이도록 요청했다. 시스템에서 **Apple M5 Pro, CPU 18코어, GPU 20코어, 메모리 64GB, Metal 지원**을 확인했다.

## 지금 확인된 결과

| 대상 | 확인 결과 | 현재 적용 상태 |
|---|---|---|
| 실행 중인 QE 7.6 | 설치 빌드는 `QE_GPU`가 비어 있고 OpenMP ON·MPI OFF다. 소스의 GPU 선택지는 CUDA·ROCm·oneAPI와 OpenACC/OpenMP이며 Metal 경로가 없다. | 현재 C09 계산은 CPU에서 계속 실행 |
| ORCA 6.1.1 macOS arm64 | 이 프로젝트에서 실제로 검증한 경로는 CPU MPI 병렬 실행이다. 확인한 배포·공식 병렬 문서에서 Apple GPU 실행 경로를 확보하지 못했다. | 검증된 12 MPI 프로세스 구성 유지 |
| PyTorch MPS | 실제 `mps:0`에서 float32 행렬 연산·미분 성공. CPU 자동 대체를 끈 상태에서 확인했다. | GPU 설치 점검 통과 |
| MACE | 별도 환경에 설치하고 계산기 import를 확인한다. 공식 프로젝트는 MPS 가속을 안내한다. | 실제 모델의 에너지·힘·속도 비교는 아직 미실행 |

QE 빌드 근거는 `.tools/build/qe-7.6/CMakeCache.txt`, GPU 선택 구현은 `.tools/src/qe-7.6/cmake/DetermineGPUFeatures.cmake`다. 현재 DFT 엔진이나 실행 중인 입력을 GPU 코드로 바꾸지 않았다.

## 분리한 실행 환경

`.tools/mps-env`에 PyTorch 2.14.0, mace-torch 0.3.16과 의존성을 설치했다. 기존 `.venv`와 ORCA/QE 설치에는 섞지 않았다. 의존성은 [요구 목록](../environment/mps-requirements.in), [버전·해시 고정 목록](../environment/mps-requirements.lock), [설치 결과](../environment/mps-installed-r1.json)에 기록한다.

기본 샌드박스에서는 MPS가 설치되어 있어도 장치 사용 가능 상태가 false였다. [첫 시도](../reports/completion/apple-gpu-preflight-sandbox-r1.json)를 보존하고, 허용된 호스트 장치 접근으로 같은 점검을 다시 실행해 [성공 기록](../reports/completion/apple-gpu-preflight-host-r1.json)을 남겼다. 이는 비과학적 128×128 행렬 한 번과 미분을 확인한 작은 설치 검사다.

GPU 계산값의 CPU float64 기준 대비 차이는 무차원 목적함수 3.50e−6, 최대 미분 성분 3.71e−8이었다. 메모리 할당은 약 128KiB였다. 이 수치를 약물·반도체 에너지 오차나 연구 계산 가속률로 해석하지 않는다. MPS float64 할당은 실제로 거절됐다.

## 실제 계산에 적용하는 기준

1. MACE 후보는 기존 최대 2개 모델 비교 범위와 적합성 검토를 유지한다. 먼저 사용 범위가 확인된 후보의 같은 가중치·원자·좌표·주기 셀에서 CPU float64, CPU float32, MPS float32 에너지와 힘을 비교한다. 가중치는 별도로 출처·버전·해시를 고정한다.
2. MPS가 실제로 쓰였는지 확인하고, 자동 CPU 대체를 끈다. float64 가중치를 MPS로 바로 읽으면 실패할 수 있으므로 CPU에서 읽어 float32로 바꾼 뒤 옮기는 경로를 점검한다. 이 형 변환을 DFT의 수치 조건 변경과 혼동하지 않는다.
3. 총에너지 오차뿐 아니라 쌍−분리층의 상호작용 에너지 차이와 원자별 힘 차이를 확인한다. 비교하려는 물리량의 기존 허용 오차를 만족해야 채택한다. MACE의 학습 방법·분산 처리와 C09 DFT의 차이도 [모델 준비도 검토](../reports/completion/MATERIALS_FOUNDATION_SCREEN.ko.md)에 따라 분리한다.
4. 실제 사용할 6원자·24원자 크기 및 묶음 처리에서 준비 시간을 제외한 반복 시간과 전체 시간을 함께 측정한다. 작은 구조에서 GPU가 더 느리면 빠른 CPU 경로를 사용한다. GPU를 썼다는 사실만으로 가속됐다고 보고하지 않는다.
5. 계산 엔진을 한 번에 하나 실행하는 기존 규칙은 유지한다. 이번 작은 GPU 설치 검사 외의 모델 계산·성능 비교는 현재 QE 배치와 겹치지 않게 실행한다. GPU가 빠르고 필요한 정밀도를 만족하는 단계는 후속 입력·설정에 명시해 사용한다.

실제 모델 가중치·추론 결과나 연구 계산의 GPU 가속률은 아직 없다. MACE를 사용할 수 있다는 사실은 기존 r2SCAN-3c/C09 DFT를 같은 계산으로 대체할 수 있다는 뜻이 아니다. 최종 DFT 검증은 원래 방법과 정밀도를 유지한다.

## CPU 시간 단축 점검

현재 MoSe₂ 참조 SCF를 10.01초 동안 관측했을 때 `pw.x`의 평균 CPU 사용량은 1.704코어 상당이었다. [실측 기록](../reports/completion/qe-live-cpu-use-r1.json). OpenMP 4개 설정이 매 순간 4코어를 가득 사용한다는 뜻은 아니며, 이 짧은 구간을 전체 실행의 평균이나 최적 코어 수로 일반화하지 않는다.

현재 배치 종료 후 대표 구조에서 4·8·12 OpenMP 스레드의 실제 에너지·힘 일치와 실행 시간을 비교한다. 기존 빌드는 MPI가 꺼져 있으므로 k-point pool MPI를 명령 옵션만으로 활성화할 수 없다. OpenMP 증설 효과가 작으면 별도 MPI 빌드와 k-point 분할을 후속 검토한다. 검증된 실행 파일과 이전 결과는 보존하며, 실행 중 의존 파일을 바꾸지 않는다.

## 근거

- [QE 공식 설치 안내](https://www.quantum-espresso.org/Doc/user_guide/node11.html). 이 문서보다 현재 설치 소스·빌드 설정을 버전별 판단에 우선했다.
- [ORCA 6.1 병렬 실행 문서](https://www.faccts.de/docs/orca/6.1/manual/contents/essentialelements/parallel.html).
- [MACE 공식 저장소의 Apple Silicon 안내](https://github.com/ACEsuit/mace).
- [PyTorch MPS와 배정밀도 제한](https://docs.pytorch.org/docs/main/notes/mps.html).
