# DFT 학습·탐색 실험실

Apple Silicon Mac에서 반도체와 약물 분자 계산을 함께 배우는 10주 프로젝트다. 목표는 기준 계산 재현과 후속 검증 후보 선정이다.

- [1주차 실제 실험 결과](reports/WEEK1_RESULTS.ko.md): 물·아세트아미노펜 최적화/진동, MoS₂ 단층 이완과 추가 SCF 확인 완료. 문헌·화학적 상태·정밀도 검토의 남은 항목 포함
- [1주차 방법 결정](docs/WEEK1_METHOD_DECISIONS.ko.md): 원문 조건, 양성자 이동 후보, 퍼텐셜·단위·파서 검증
- [10주 작업 계획](docs/PLAN.ko.md): 주차별 산출물, 계산 예산, 비교 방법, 통과·보류 기준
- [계산 실험 계획 v1](docs/EXPERIMENT_PROTOCOL_V1.ko.md): 첫 4주 파일럿, 약물 대조군·반도체 구조 행렬·판정·예산·실행 순서. 연구 실행 시 이 계획을 우선 적용
- [학술 기여 가능성 조사와 계획 수정 제안](docs/RESEARCH_OPPORTUNITIES_2026-09-09.ko.md): 2026-09-09 문헌 조사, 양 분야 파일럿과 진행·전환 기준. 새 계산은 미실행
- [설치 결과](reports/INSTALLATION.ko.md): 설치한 버전과 실제 실행 검증
- [ORCA 설치·실행 안내](docs/ORCA_SETUP.ko.md): 설치 완료, 물 분자 직렬 실행 검증
- [첨부 대화 원문](docs/source/original-conversation.ko.txt)

## 바로 실행

```bash
cd /Users/xavier/Documents/dev/generaldrugs/dft
source scripts/activate.sh
python scripts/doctor.py
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

`structures/molecules`, `structures/materials`, `structures/week1`은 초기 구조·후보 목록이다. 실제 최적화 결과는 `structures/optimized/week1`에 구분해 보관했다. 물/Si 설치 검증의 단일점 결과와 1주차 연구 실행을 혼동하지 않는다. Mo/S/Se 의사퍼텐셜은 확보했고 cutoff/k-grid/진공 수렴과 이중층 기준은 연구 프로토콜의 2주차 작업이다.

1주차 완료 출력의 판독·해시·구조 확인만 다시 하려면:

```bash
python scripts/verify_week1.py
python scripts/summarize_week1.py
```

연구용 새 계산은 `python scripts/research_run.py --config experiments/week1/설정파일.json`으로 명시적으로 실행한다. 이 명령은 실제 계산을 다시 수행하므로 결과 검증 명령과 구분한다. 작업별 설정·실행 시간·수렴 상태는 [실행 요약](reports/week1-summary.json), 앞으로 필요한 조건은 [현재 상태](experiments/week1/status.json)에 있다. 자동 예약 작업은 없다.

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
