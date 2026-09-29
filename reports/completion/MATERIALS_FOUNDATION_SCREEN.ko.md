# 반도체 공개 기반 모델 후보와 비교 조건

2026-09-29. 공개 배포 메타데이터와 계산 소스를 읽어 후속 모델 비교의 가능성을 추가로 조사했다. **가중치 다운로드·설치·추론은 하지 않았으며, CI/CO 계면에 대한 정확도나 학습 중복은 아직 검증하지 않았다.** 기존 C09 DFT 실행과 판정 기준은 그대로 유지한다.

## 실제 등록된 두 후보

| 후보 | 저자가 명시한 학습 기준 | 등록된 가중치 크기 | 배포 조건과 현재 판정 |
|---|---|---:|---|
| MACE-MP-0b3 medium | MPTrj, PBE+U 계열 | 79,472,952바이트 | 공식 모델 표는 MIT로 표시한다. 후속 비교 후보이며, Mo/S/Se 학습 계산의 상세 조건·층간 적용 범위·학습 중복은 추가 확인해야 한다. |
| MACE-MATPES-PBE-0 | MATPES PBE, +U 없음 | 79,471,284바이트 | 공식 표와 릴리즈는 ASL로 표시한다. 사용 자격이 확인되지 않아 실행 대상으로 확정하지 않았다. |

근거는 [고정한 공식 모델 표](https://github.com/ACEsuit/mace-foundations/blob/16a9f178706ce053f3ca8531efbab00a305d0854/README.md), [MP-0b3 릴리즈](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_mp_0b3), [MATPES 릴리즈](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_matpes_0)다. 두 모델의 GitHub 자산 `digest`는 조회 시 `null`이었다. 표의 크기는 API에 등록된 값이며 실제 가중치 파일을 내려받아 무결성을 확인한 결과가 아니다. 조사한 후보 수와 실제 비교할 모델 수를 구별하며, 실행 비교는 원래 계획대로 적합한 최대 2개다.

ASL 13(a)는 지정한 기관에서 수행하는 학술·비영리 연구로 사용 범위를 정의한다. 현재 대화에서 그 요건을 확인하지 않았으므로 개인 학습이라는 이유만으로 적용 가능하다고 가정하지 않는다. 이는 MATPES 후보의 실행 조건이며, 다른 DFT 실험을 중단할 이유는 아니다. [확인한 ASL 원문](https://github.com/gabor1/ASL/blob/52733adb881efea63e0a63c5b8534faf95db475b/ASL.md)

## 분산 옵션의 이름만으로 조건을 맞출 수 없다

MACE의 `mace_mp` 함수는 기본적으로 분산 보정이 꺼져 있다. 이를 켠 경우의 실제 소스와, 프로젝트가 준비한 PBE-D3(BJ) 입력을 비교했다.

| 항목 | 조사한 MACE 보조 함수에서 D3를 켰을 때 | 프로젝트의 QE 7.6 대체 방법 준비 |
|---|---|---|
| 함수·감쇠 | PBE, BJ | PBE, D3 version 4 |
| ATM 세 원자 보정 | 연결된 TorchDFTD3의 기본 `abc=False` | 입력에서 `dftd3_threebody=.true.` |
| 분산 거리 절단 | 40 bohr | API 기본 √9000 ≈94.86833 bohr |
| 배위수 거리 절단 | 40 bohr | API 기본 40 bohr |
| 모델 계산 자료형 | 기본 float32 | 모델을 사용한다면 float64 검증을 별도로 수행해야 함 |

TorchDFTD3 계산기를 직접 만들 때의 기본 감쇠는 `zero`, 거리 절단은 95 bohr다. MACE 보조 함수는 이 값을 BJ와 40 bohr로 덮어쓴다. 따라서 라이브러리 이름이나 어느 한 생성자의 기본값만 읽어 전체 계산 조건이 같다고 판단하면 안 된다. 위 내용은 [MACE 소스](https://github.com/ACEsuit/mace/blob/c5c65fe77aa6fdc0e2c587a95abc52734d9ddb8e/mace/calculators/foundations_models.py)와 [TorchDFTD3 소스](https://github.com/pfnet-research/torch-dftd/blob/28b43557e6d950502941990fa5a1c9a3aa9e8339/torch_dftd/torch_dftd3_calculator.py)의 정적 확인이다. 실제 두 구현의 에너지·힘 일치는 아직 시험하지 않았다.

QE 쪽은 설치에 사용한 `dft-d3/api.f90`의 기본값과 `PW/src/input.f90`의 초기화 호출, 이미 준비한 `w3-qe-pbe-d3bj-compatibility.in`을 읽었다. 현재 실행 중인 C09 계산에는 이 D3 설정을 추가하지 않았다.

## 후속 실행에 필요한 순서

1. 원래 기준·비용 검증과 PBE-D3(BJ) 호환성 계산을 완료한다. 모델 후보의 공개 등록은 이 단계의 통과를 대신하지 않는다.
2. 사용 가능한 후보를 고정하고, 가중치·환경·학습 자료 출처와 Mo/S/Se 지원을 확인한다. 학습 중복이 불명확하면 독립 검증이라고 부르지 않는다. MATPES 후보의 사용 요건이 불명확한 동안에는 실행하지 않는다.
3. 동일한 구조·셀에서 분산 보정만 비교한다. PBE/BJ·세 원자 항·거리 절단·주기 경계·단위를 명시적으로 맞추고, 에너지·힘 차이를 측정한다. 소스 기본값을 맞추는 것만으로 구현 일치를 선언하지 않는다.
4. 후속 PBE-D3 DFT와 동일한 pair/bottom/top 좌표에서 모델의 상호작용 에너지와 상대 위치 차이를 계산한다. 원자 에너지 기준을 섞지 않고, 같은 모델의 세 성분으로 차이를 구한다. 실제 Mac에서의 시간·메모리·유한값 검증도 이 단계에 필요하다.
5. C09와의 차이는 모델 근사와 DFT 기준 방법 차이로 나누어 보고한다. 이 구분에 필요한 PBE-D3 자료를 확보하지 못하면 C09 대비 차이를 순수 모델 오차로 해석하지 않는다.

동일 구조에서 이 구분은 다음 항등식으로 표현된다.

\[
\gamma_{\mathrm{MLP+D3}}-\gamma_{\mathrm{C09}}
=\left(\gamma_{\mathrm{MLP+D3}}-\gamma_{\mathrm{PBE-D3}}\right)
+\left(\gamma_{\mathrm{PBE-D3}}-\gamma_{\mathrm{C09}}\right).
\]

첫 괄호에도 학습 DFT와 이번 QE 계산의 의사퍼텐셜·세부 조건 차이가 남을 수 있다. 정확한 기준 일치를 확인하기 전에는 그 한계도 표시한다. 현재 이 식의 모델 예측값이나 새로운 PBE-D3 결과는 산출하지 않았다.

## 보존한 근거

공식 저장소 커밋·릴리즈 API·소스·라이선스는 `data/completion/materials-models/foundation-screen-r1/`에 읽기 자료로 보존했다. [최초 취득 목록](../../data/completion/materials-models/foundation-screen-r1/manifest.json), [MP-0b3 추가 목록](../../data/completion/materials-models/foundation-screen-r1/addendum-manifest.json), [정적 비교와 파일 해시](materials-foundation-screen-r1.json)에 출처와 범위를 기록했다. 다운로드한 Python 소스를 import하거나 실행하지 않고 AST로 함수의 기본값을 읽었다. 이 자료 묶음에는 모델 가중치가 없다.
