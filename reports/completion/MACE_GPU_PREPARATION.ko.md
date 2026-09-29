# MACE-MP-0b3 가중치 확보와 Apple GPU 비교 준비

2026-09-29 14:23 UTC 기준. 사용자 GPU 활용 요청에 따라 공식 ACEsuit 배포에서 MIT 표기 **MACE-MP-0b3 medium** 가중치를 확보했다. 이전 `apple-gpu-readiness-r1.json`과 `MATERIALS_FOUNDATION_SCREEN.ko.md`는 다운로드 이전 조사 기록으로 보존한다.

| 확인 항목 | 실제 확인 결과 |
|---|---|
| 파일 크기 | 79,472,952 bytes |
| SHA-256 | `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012` |
| 모델 클래스 | ScaleShiftMACE |
| 매개변수 수 | 9,063,204 |
| 원소 목록 | 89종, Mo·S·Se 포함 |
| 모델 cutoff | 6.0 Å |
| 저장된 부동소수점 가중치·버퍼 | float64 |
| 읽기 검증 환경 | CPU, PyTorch 2.14.0 / mace-torch 0.3.16 |
| 실제 에너지·힘 추론 | 아직 미실행 |

[공식 배포](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_mp_0b3)에서 받은 파일의 크기를 릴리스 메타데이터와 대조하고 로컬 SHA-256을 계산했다. GitHub 메타데이터의 upstream digest는 없었으므로 이 해시는 로컬 바이트의 재현성을 고정하며 별도의 서명 검증을 뜻하지 않는다. [다운로드 기록](mace-mp0b3-checkpoint-r1.json)과 [실제 가중치 읽기 결과](mace-mp0b3-checkpoint-metadata-r1.json)를 보존했다.

원소 지원 목록에 들어 있다는 사실은 MoS₂/MoSe₂ 계면 정확도나 C09 조건과의 일치를 입증하지 않는다. 6 Å cutoff 모델의 장거리 분산, 학습 데이터 중복, PBE 계열 학습 조건과 기준 C09의 차이는 [기존 모델 검토](MATERIALS_FOUNDATION_SCREEN.ko.md)의 별도 검증 대상으로 남는다. 이번 읽기에서는 분산 보정을 추가하거나 에너지·힘을 계산하지 않았다.

Apple GPU 설치 검사는 이미 통과했지만 MPS는 float64를 지원하지 않는다. 실제 GPU 비교에서는 이 가중치를 CPU에서 읽고 float32로 변환한 뒤 옮겨야 한다. CPU float64·CPU float32·MPS float32의 에너지·힘·상호작용 에너지 차이와 실제 소요 시간을 비교한 후 사용 범위를 결정한다. 단순히 메모리에 읽힌 것을 GPU 추론 성공이나 시간 단축으로 세지 않는다.

현재 QE와 실제 모델 추론을 동시에 실행하지 않는다. 현재 기준 배치 종료 후 CPU/MPS 비교를 진행할 준비가 된 상태이며, 기존 GPU 사용 기준은 [Apple GPU 적용 문서](../../docs/APPLE_GPU_ACCELERATION.ko.md)를 따른다. 실제 모델 비교 최대 2개의 기존 범위와 독립 검증 기준은 유지한다.
