# DFT 학습·탐색 실험실

Apple Silicon Mac에서 반도체와 약물 분자 계산을 함께 배우는 10주 프로젝트다. 목표는 기준 계산 재현과 후속 검증 후보 선정이다.

**2026-09-30 00:08 KST, 사용자 요청으로 실험을 중간 종료했다. 전체 연구 계획은 미완료다.** 진행 중이던 마지막 SCF가 정상 종료된 뒤 제어기를 종료했고, 계산 프로세스와 실행 잠금이 남지 않았음을 확인했다. 자동 재개하지 않는다. 완성된 예측 모델이나 검증된 공결정·반도체 발견을 주장하지 않으며, 과거 실패·중단 기록도 보존한다.

- [중단 지점·확보한 결과·미완료 항목](reports/completion/USER_STOP_CHECKPOINT.ko.md)
- [전체 완료 계획과 현재 단계](docs/COMPLETION_PLAN.ko.md)
- [공개 범위와 검증 방법](docs/PUBLIC_RELEASE.ko.md)
- [약물 정밀도 진단](reports/completion/PRECISION_DIAGNOSTIC.ko.md) · [전체 약물 패널 검증 절차](docs/DRUG_FULL_PANEL_VALIDATION.ko.md)
- [두 단층 격자·수치 검증](reports/completion/MONOLAYER_EQUILIBRIA.ko.md) · [공통 격자 내부 이완](reports/completion/COMMON_LATTICE_LAYERS.ko.md): 첫 이종 이중층 단일점까지 종료했으나 분리층 기준과 결합에너지는 미완료다.

- [2주차 중단 당시 기록](reports/WEEK2_PAUSED.ko.md) · [재개 방법과 검증](docs/WEEK2_RESUME.ko.md): 87회 최적화 좌표 보존, 별도 12프로세스 진동 계산
- [2주차 실제 실험 결과](reports/WEEK2_RESULTS.ko.md): xTB 후보 선별, 첫 약물 DFT/진동, MoS₂ 수렴·순수 이중층 대조
- [2주차 실행 계획](docs/WEEK2_EXECUTION.ko.md): 범위와 정밀도 판정 기준
- [1주차 실제 실험 결과](reports/WEEK1_RESULTS.ko.md): 물·아세트아미노펜 최적화/진동, MoS₂ 단층 이완과 추가 SCF 확인 완료. 문헌·화학적 상태·정밀도 검토의 남은 항목 포함
- [1주차 방법 결정](docs/WEEK1_METHOD_DECISIONS.ko.md): 원문 조건, 양성자 이동 후보, 퍼텐셜·단위·파서 검증
- [10주 작업 계획](docs/PLAN.ko.md): 주차별 산출물, 계산 예산, 비교 방법, 통과·보류 기준
- [계산 실험 계획 v1](docs/EXPERIMENT_PROTOCOL_V1.ko.md): 첫 4주 파일럿, 약물 대조군·반도체 구조 행렬·판정·예산·실행 순서. 연구 실행 시 이 계획을 우선 적용
- [학술 기여 가능성 조사와 계획 수정 제안](docs/RESEARCH_OPPORTUNITIES_2026-09-09.ko.md): 2026-09-09 문헌 조사, 양 분야 파일럿과 진행·전환 기준. 새 계산은 미실행
- [설치 결과](reports/INSTALLATION.ko.md): 설치한 버전과 실제 실행 검증
- [ORCA 설치·실행 안내](docs/ORCA_SETUP.ko.md): 설치 완료, 물 분자 직렬 실행 검증

원본 대화와 로컬 설치물은 공개 저장소에 포함하지 않는다. 내려받은 논문·보충자료·외부 모델·의사퍼텐셜은 해당 출처에서 별도로 확보해야 한다. 배포 권한과 실행 환경의 범위는 [제3자 자료 안내](NOTICE.md)에 명시했다.

## 바로 실행

```bash
git clone https://github.com/brianyu43/dft-learning-lab.git
cd dft-learning-lab
# 아래 활성화·계산 명령은 Python 환경 및 해당 엔진 설치 후 사용한다.
source scripts/activate.sh
python scripts/doctor.py
python scripts/week2_status.py  # 계산을 시작하지 않고 현재 상태만 읽기
```

물 분자와 실리콘의 작은 실제 계산을 다시 실행하려면:

```bash
python scripts/smoke.py pyscf
python scripts/smoke.py xtb
python scripts/smoke.py qe
python scripts/smoke.py orca
```

각 명령은 최대 10분으로 제한한다. OpenMP 스레드 환경 변수의 기본값은 4이며 ORCA 예제는 `%pal` 없는 직렬 계산이다. `runs/<시각>-<엔진>/`에 입력·출력·메타데이터를 새로 저장한다. `reports/smoke-*.json`은 가장 최근 결과를 가리킨다. 원시 실행 디렉터리는 덮어쓰지 않는다. RSS는 하위 프로세스의 최대 상주 메모리이며 전체 Mac의 통합 메모리 사용량이 아니다.

JupyterLab을 열려면 같은 환경에서:

```bash
jupyter lab --ServerApp.ip=127.0.0.1
```

`structures/molecules`, `structures/materials`, `structures/week1`은 초기 구조·후보 목록이다. 실제 최적화 결과는 `structures/optimized/week1`에 구분해 보관했다. 물/Si 설치 검증의 단일점 결과와 1주차 연구 실행을 혼동하지 않는다. Mo/S/Se 의사퍼텐셜은 로컬에 확보했다. MoS₂의 고정 구조 수렴과 순수 이중층 대조는 2주차 결과에 있다. 이후 두 조성의 평형 단층에서 공통 수치 설정을 검증했으며, 이종 이중층 결합에너지의 정밀도 일반화는 미완료다.

1주차 완료 출력의 판독·해시·구조 확인만 다시 하려면, 원래 로컬 실행 환경과 의사퍼텐셜 등 필요한 자산이 있어야 한다. 과거 보고서에는 원본 절대 경로와 해시를 그대로 남겼으므로 새 복제본에서 전체 원시 검증 명령이 바로 실행된다고 가정하지 않는다. 공개 복제본의 파일 무결성 및 엔진 없는 검사는 [공개 범위](docs/PUBLIC_RELEASE.ko.md)를 따른다.

원래 실행 환경에서의 1주차 검증 명령:

```bash
python scripts/verify_week1.py
```

연구용 새 계산은 `python scripts/research_run.py --config experiments/week1/설정파일.json`으로 명시적으로 실행한다. 이 명령은 실제 계산을 다시 수행하므로 결과 검증 명령과 구분한다. 1주차 기록은 [1주차 요약](reports/week1-summary.json), 현재 중단 상태는 [최종 체크포인트](reports/completion/user-stop-checkpoint-r1.json)에 있다. 과거 상태 파일의 `running`은 당시 스냅샷일 수 있다. `summarize_week1.py`와 `verify_week1.py`는 과거 보고서를 다시 쓸 수 있으므로 이번 중단 검증에서 실행하지 않았다. 자동 예약 작업은 없다.

## 환경 구성과 재현

Python은 `.venv`, QE는 `.tools/qe-7.6`, xTB는 `.tools/xtb`, ORCA는 `.tools/orca-6.1.1`에 있다. QE는 로컬 Homebrew의 GCC/FFTW/OpenBLAS 라이브러리에 연결된다. 가상환경이나 엔진 폴더를 다른 Mac으로 복사해 실행하는 방식은 지원하지 않는다.

Python 환경은 다음과 같이 복원한다. Homebrew의 `uv`와 Python 3.11이 필요하다.

```bash
UV_CACHE_DIR="$PWD/.cache/uv" uv sync --frozen
```

QE 빌드 의존성은 `environment/Brewfile` 및 설치 버전 기록에 남아 있다. 해당 GCC 16 경로가 있는 Mac에서 `bash scripts/build_qe.sh`로 QE를 다시 빌드할 수 있다. xTB는 프로젝트 내부 micromamba와 `environment/xtb-explicit.txt`로 복원한다.

```bash
MAMBA_ROOT_PREFIX="$PWD/.tools/mamba-root" .tools/micromamba/bin/micromamba create \
  -y -p "$PWD/.tools/xtb" --file environment/xtb-explicit.txt
```

기본 스레드 수를 바꾸려면 활성화 전에 `export DFT_THREADS=4`를 설정한다. 전역 `.zshrc`는 수정하지 않는다. 재생성 스크립트 `scripts/prepare_examples.py`는 결정적인 시드로 예제 구조와 ORCA 입력을 다시 쓰므로, 사용자가 수정한 예제는 먼저 별도 이름으로 보관한다.
