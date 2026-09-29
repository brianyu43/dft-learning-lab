# 설치 및 검증 결과

확인일: 2026-09-09. 장비: Apple M5 Pro / CPU 18코어 / 메모리 64 GiB / macOS 26.6.2 / arm64.

ORCA를 포함한 계산·분석 환경을 설치했고 작은 실제 계산으로 검증했다. 2026-09-09 다운로드된 ORCA 6.1.1 arm64 패키지를 추가 설치하고 직렬 물 분자 계산을 확인했다. 10주 탐색 계산 전체를 실행한 상태는 아니다.

## 설치 완료

| 구성 | 버전 | 설치 위치·역할 |
|---|---|---|
| Quantum ESPRESSO | 7.6 | `.tools/qe-7.6/bin`, arm64 CPU/OpenMP, MPI 없음 |
| ORCA | 6.1.1 | `.tools/orca-6.1.1`, macOS arm64, 직렬 실행 검증 |
| xTB | 6.7.1 | `.tools/xtb`, conda-forge arm64 환경 |
| Python | 3.11.15 | 기존 uv Python을 이용한 프로젝트 `.venv` |
| ASE / RDKit | 3.29.0 / 2026.3.6 | 구조 생성·입출력·배좌 전처리 |
| PySCF / geomeTRIC | 2.14.0 / 1.1.1 | 공개 분자 DFT 학습·최적화 라이브러리 |
| pymatgen / spglib | 2026.5.4 / 2.7.0 | 결정·대칭 처리 |
| JupyterLab | 4.6.3 | localhost 웹 환경 |
| cclib / py3Dmol | 1.8.1 / 2.5.5 | 출력 파싱·분자 시각화 라이브러리 |

NumPy·SciPy·pandas·matplotlib·ipykernel도 설치했다. 전체 버전과 import 결과는 `doctor.json`, Python 잠금은 `../uv.lock`, xTB 명시적 패키지 목록은 `../environment/xtb-explicit.txt`에 있다. Homebrew 빌드 도구는 GCC 16.2.0, CMake 4.4.3, FFTW 3.3.11, OpenBLAS 0.3.34이다.

QE는 공식 `qe-7.6` 릴리스의 commit `9f93ddec427d2b9a45bb72d828c6d324f62fcabd`에서 빌드했다. 설치 실행 파일은 `pw.x`, `ph.x`, `pp.x`, `bands.x`, `dos.x`, `projwfc.x`, `q2r.x`, `matdyn.x`, `dynmat.x`이다. 이 중 실제 계산 검증은 `pw.x`에서 수행했다. 나머지 실행 파일의 실제 물성 계산과 전체 upstream 회귀 테스트는 수행하지 않았다.

## 실제 설치 검증

| 엔진 | 계산 | 경과 시간 | 최대 RSS | 결과 |
|---|---|---:|---:|---|
| qe 7.6 | Si2 diamond / PBE/ONCV | 0.349초 | 24.8 MiB | 통과 |
| pyscf 2.14.0 | H2O / PBE/def2-SVP | 0.451초 | 132.9 MiB | 통과 |
| xtb 6.7.1 | H2O / GFN2-xTB (semiempirical, not DFT) | 0.021초 | 13.3 MiB | 통과 |
| orca 6.1.1 | H2O / PBE/def2-SVP | 4.287초 | 122.2 MiB | 통과 |

모두 종료 코드 0, 유한한 에너지, 전자 계산 수렴을 확인했다. 시간은 이 Mac에서 최종 재실행한 소규모 작업의 측정값이며 초기 로딩과 작업 규모에 따라 달라진다. 약물 최적화나 12원자 MoSSe 탐색의 시간 예측으로 사용하지 않는다. 메모리는 `getrusage(RUSAGE_CHILDREN).ru_maxrss`의 macOS 바이트 값을 MiB로 환산한 것으로, 프로세스 트리 메모리 합이나 Mac 전체 메모리 최고값이 아니다.

- PySCF 물: -76.272014472667 Hartree, PBE/def2-SVP, 진공, charge 0, multiplicity 1.
- QE 실리콘: -15.57645486 Ry, Si2, PBE/ONCV, cutoff 30/120 Ry, 2×2×2 k-grid. 이 설정의 정밀도 수렴은 검증하지 않았다.
- ORCA 물: -76.272118859779 Hartree, PBE/def2-SVP, 진공, charge 0, multiplicity 1. `ORCA TERMINATED NORMALLY` 확인. PySCF와 수치 적분 등 세부 설정이 달라 에너지가 동일하다고 가정하지 않는다.
- xTB 물: -5.070364532464 Hartree, GFN2-xTB. DFT 에너지와 서로 직접 비교하지 않는다.
- JupyterLab: localhost에서 실제 서버를 잠시 시작해 `/api/status`와 `/lab`의 HTTP 200 응답을 확인하고 종료했다. 상시 서버는 남겨두지 않았다.
- Python 130개 설치 패키지에 대해 `uv pip check` 통과. `uv sync --frozen --offline` 성공. 주요 라이브러리 import 및 작성 스크립트 문법 점검 통과.
- 6개 분자의 원자 수·XYZ 해시·SDF 파싱, Si2·MoS2/MoSe2의 3/12원자 CIF·20 Å 층간 진공 확인. 원문을 바이트 단위로 동일하게 보존했다.

초기 계측에서는 `/usr/bin/time -l`이 샌드박스의 `sysctl` 제한으로 종료 코드 1을 반환했다. 계산 자체는 수렴했으나 초기 manifest는 실패로 보존했다. 계측을 Python의 OS resource API로 바꾸고 세 엔진을 다시 실행해 통과 결과를 기록했다. QE stderr의 Fortran 부동소수점 상태 플래그 알림은 원시 로그에 남아 있다. 이 확인은 소규모 설치 검증이며 연구 정밀도의 검증을 대신하지 않는다.

## 만들어 둔 작업물

- `docs/PLAN.ko.md`: 주차별 연구 질문·실행·산출물, 자원 예산, 진행/중단 기준, 비교 규칙.
- `structures/molecules/`: 아세트아미노펜과 5종 공형성제의 MMFF 초기 XYZ·배좌 SDF, 물.
- `structures/materials/`: Si와 MoS2/MoSe2 원시·2×2 초기 CIF.
- `inputs/orca/`: 물 단일점·최적화/진동, 아세트아미노펜 진공·물 최적화 입력. 물 단일점은 실행 검증했고 최적화·진동·아세트아미노펜은 미실행.
- `inputs/qe/`: 실리콘 입력과 공식 QE 예제에서 가져온 Si 의사퍼텐셜, 출처·SHA-256.
- `scripts/activate.sh`, `doctor.py`, `smoke.py`: 환경 활성화, 설치 확인, 제한시간 있는 실제 실행·기록.
- `environment/`: 설치 메타데이터·잠금 목록·QE 소스 commit/하위 모듈·일반 계산 기록 양식.

## 설치 과정의 조정

xTB Homebrew 설치에서 추가 외부 formula를 영구 신뢰하는 작업이 자동 승인 검토에 거부됐다. 시스템 신뢰 설정 변경 없이 개발팀이 안내하는 conda-forge 경로를 사용해 프로젝트 내부에 설치하고 검증했다. 미사용 Homebrew tap은 정리했다.

Jupyter의 Node 스크립트가 상위 폴더의 `type: module` 설정을 물려받는 오류가 있었다. 이 DFT 폴더에 `type: commonjs`를 지정하는 최소 `package.json`을 추가해 해결했다. 상위 프로젝트 설정은 수정하지 않았다.

## 남은 단계

ORCA 설치와 직렬 물 계산 검증은 완료했다. 설치 파일 내 MD5 무결성 검사를 통과했으며 별도로 SHA-256을 기록했다. 패키지는 Open MPI 4.1.1용 병렬 실행 파일도 포함하지만, 해당 MPI 런타임 설치와 병렬 계산 검증은 아직 수행하지 않았다. 현재 사용법은 `../docs/ORCA_SETUP.ko.md`에 있다.

SCF 설치 검증에서 기하 최적화·진동·안정 최소점은 확인하지 않았다. 실제 아세트아미노펜 최적화, Mo/S/Se 의사퍼텐셜 선택, 문헌 재현, 수렴 시험, 후보 탐색은 10주 계획에 따라 별도로 진행한다.
