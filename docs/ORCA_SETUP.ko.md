# ORCA 설치·실행 안내

2026-09-09 ORCA 6.1.1 macOS arm64 설치 완료. 물 분자의 PBE/def2-SVP 직렬 계산에서 SCF 수렴, 종료 코드 0, `ORCA TERMINATED NORMALLY`를 확인했다.

공식 문서: [설치 안내](https://www.faccts.de/docs/orca/6.1/tutorials/first_steps/install.html), [ORCA 이용 안내](https://www.faccts.de/orca/). 공식 포럼 또는 FACCTs 다운로드 계정이 필요하며, 배포처에서 이용 조건에 동의해야 한다. 로그인·비밀번호는 대화에 보내지 않는다.

## 이제 실행하는 방법

```bash
cd /Users/xavier/Documents/dev/generaldrugs/dft
source scripts/activate.sh
python scripts/smoke.py orca
```

입력·출력·측정값은 실행별 `runs/<시각>-orca/`에 저장한다. 최근 결과는 `reports/smoke-orca.json`이다. 설치 `.run` 파일을 다시 실행할 필요는 없다.

자신의 입력 파일을 직접 계산할 때는 활성화 후 별도의 계산 디렉터리에서 다음 형식으로 실행한다.

```bash
"$ORCA_DIR/orca" molecule.inp > molecule.out 2> molecule.err
```

이 명령에는 시간 제한이 없다. `molecule.inp`는 해당 디렉터리에 실제로 준비한 입력 파일명으로 바꾼다. ORCA는 GUI 앱을 여는 방식이 아니라 입력 파일을 읽고 출력 파일을 만드는 명령줄 프로그램이다.

## 설치 기록

- 원본: `/Users/xavier/Downloads/orca_6_1_1_macosx_arm64_openmpi411.run`
- 설치 위치: `/Users/xavier/Documents/dev/generaldrugs/dft/.tools/orca-6.1.1`
- 실행 파일의 Mach-O `arm64` 확인, 내장 MD5 무결성 검사 통과. 내장 SHA-256·서명은 없었으며 별도의 파일 SHA-256을 `environment/orca-installer.sha256`에 기록했다.
- `.run`은 shell 기반 Makeself 설치기다. 원본을 `sh <파일>.run --noexec --target <프로젝트 내 경로>`로 풀고 설치 디렉터리 이름과 실행 경로를 설정했다.
- `environment/orca.local.sh`를 `scripts/activate.sh`에서 읽도록 연결했다. 전역 `.zshrc`는 수정하지 않았다. 원본 `setup`은 파일 이동과 전역 shell 설정 추가를 수행한다.
- 패키지에 포함된 EULA·제3자 라이선스 파일은 설치 폴더에 보존했다.

## 검증 범위와 다음 단계

이번 확인은 물의 직렬 단일점 계산이다. 패키지 이름의 `openmpi411`은 Open MPI 4.1.1용 병렬 구성이라는 뜻이며, 해당 MPI 런타임 설치·연결·병렬 실행은 아직 확인하지 않았다. 다음 단계에서 `inputs/orca/water_opt_freq.inp`의 r2SCAN-3c 최적화·진동 계산을 진행한다.

준비된 `acetaminophen_vacuum_opt.inp`와 `acetaminophen_water_opt.inp`는 동일 MMFF 시작 기하에서 진공/물의 최적화를 시작하는 예제다. 아직 실행하지 않았다. 최종 주파수 계산은 각각 최적화된 기하를 사용해야 한다.

`python scripts/smoke.py pyscf`로 공개 분자 DFT 실습도 실행할 수 있다. 이 PySCF PBE 계산은 ORCA r2SCAN-3c와 동일한 방법이 아니다.
