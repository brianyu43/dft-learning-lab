# GitHub 공개 범위와 검증

2026-09-29. 공개 저장소: `brianyu43/dft-learning-lab`.

## 공개하는 내용

- 프로젝트 코드·검사·환경 버전·입력·초기 및 최적화 구조.
- 실험 계획, 문헌에 근거한 사실 표, 결과 요약, 그림, 실패·중단 및 정밀도 검토 기록.
- 종료된 계산에서 선택한 입력·로그·좌표·기울기·Hessian·QE XML과 실행 메타데이터. 파동함수·실행 프로그램·동적 라이브러리·의사퍼텐셜은 제외한다.
- `publication-manifest.json`에 공개 파일의 경로·크기·SHA-256을 고정한다. manifest 자체는 자기 해시 순환을 피하기 위해 목록에서 제외한다.

원본 파일은 해시와 과거 기록을 보존하기 위해 임의로 고치지 않는다. 따라서 실행 보고서의 `run_dir`에는 원래 머신의 절대 경로가 남아 있다. 이 저장소는 원본 환경 전체를 복제한 실행 이미지가 아니다. 해시가 가리키는 외부 설치물·논문·제3자 자산이 저장소에 모두 포함돼 있다는 뜻도 아니다.

## 현재 과학적 상태

MoS₂의 C09 면내 격자 최적화는 종료돼 원시 검증을 통과했다. 공통 격자와 이종 이중층 기준 계산은 진행 중이다. 약물은 기존 최소점과 수치 정밀도 진단을 확보했지만, 공통 조건의 전체 패널·방법·용매 민감도·독립 검증은 미완료다.

실행 중인 작업의 원시 출력은 첫 공개본에서 제외한다. 상태 JSON은 커밋 시점의 스냅샷이며 GitHub가 로컬 계산의 실시간 모니터라는 뜻이 아니다. 자세한 완료 기준은 [전체 계획](COMPLETION_PLAN.ko.md)을 따른다.

## 공개 복제본에서 가능한 검사

파일 해시는 Python 표준 라이브러리만으로 검증한다.

```sh
python3 scripts/verify_publication.py
```

Python 3.11과 uv로 의존성을 복원한 뒤, 엔진 실행이나 라이선스 자료 없이 다음 검사를 실행할 수 있다.

```sh
uv sync --frozen
uv run python -m unittest discover -s tests -p 'test_research_parser.py' -v
uv run python -m unittest discover -s tests -p 'test_completion_drug_contrasts.py' -v
uv run python -m unittest discover -s tests -p 'test_completion_alloy_geometry.py' -v
uv run python -m unittest discover -s tests -p 'test_completion_alloy_probe.py' -v
```

이 검사는 파서의 실패 처리, 점수 산술, 구조·제약 및 입력 생성을 확인한다. 실제 DFT 재실행이나 전체 연구 완료를 검증하는 검사는 아니다. 다른 검사에는 원래 실행 환경과 제외된 자산이 필요할 수 있다.

첫 공개 준비에서는 Git에 올릴 파일만 추출한 깨끗한 복제본에서 위 **36개 검사(9+9+12+6)**와 공개 파일의 SHA-256 검증을 통과했다. 검사는 기존 Python 3.11 환경으로 실행했으며 새 머신에서 엔진 설치와 DFT 재실행까지 확인했다는 뜻은 아니다. `.gitattributes`는 자동 줄바꿈 변환을 끄므로 원시 파일의 바이트와 해시를 보존한다.

계산 재실행에는 ORCA/QE/xTB의 적절한 설치와 출처별 자료 확보, 경로 및 원시 해시의 명시적인 이관이 필요하다. 원본 메타데이터를 조용히 덮어쓰거나 예전 결과를 새 계산으로 표시하지 않는다.
