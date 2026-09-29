# 3주차 실험 진행 및 결과

갱신: 2026-09-11T10:41:11.075987+00:00. 전체 상태: **in_progress_or_requires_review**.

새 약물 계산 5/14건 최소점 검증 완료. 반도체: not_started. 제어기: stopped_requires_review.

## 문헌 근거

정식 SI를 확보하고 실험법·스펙트럼을 확인했다. CAF의 COOH–N9 접촉과 ADE의 N1 양성자 이동은 분광학적 해석이다. DIF 공결정 CIF와 다형 수, 미형성 시도별 조건은 확인되지 않았다. [방법과 한계](../docs/WEEK3_PROTOCOL.ko.md), [정식 SI](../data/week3/literature/ANIE-64-e202507835-s001.pdf).

## 약물 계산 상태

| 대상 | 시작점 | 역할 | 상태 | 최소점 |
|---|---|---|---|---|
| diflunisal | diflunisal-t02-s43 | primary | passed | True |
| caffeine | caffeine-t00-s17 | primary | passed | True |
| diflunisal--caffeine | diflunisal--caffeine-p36 | primary | passed | True |
| diflunisal--caffeine | diflunisal--caffeine-p14 | first_alternative_pose_and_full_job_cost | passed | True |
| diflunisal--caffeine | diflunisal--caffeine-p06 | literature_contact_COOH_to_N9 | passed | True |
| adenine | adenine-t06-s89 | primary | passed | True |
| theobromine | theobromine-t00-s89 | primary | passed | True |
| xanthine | xanthine-t00-s17 | primary | passed | True |
| diflunisal--diflunisal | diflunisal--diflunisal-p05 | primary | needs_review | False |
| caffeine--caffeine | caffeine--caffeine-p19 | primary | not_started | False |
| diflunisal--adenine | diflunisal--adenine-p39 | primary | not_started | False |
| adenine--adenine | adenine--adenine-p06 | primary | not_started | False |
| diflunisal--theobromine | diflunisal--theobromine-p28 | primary | not_started | False |
| theobromine--theobromine | theobromine--theobromine-p43 | primary | not_started | False |
| diflunisal--xanthine | diflunisal--xanthine-p10 | primary | not_started | False |
| xanthine--xanthine | xanthine--xanthine-p27 | primary | not_started | False |
| diflunisal--adenine | diflunisal--adenine-p36 | literature_N1_proton_transfer_probe | not_started | False |

### primary_only

| 공형성제 | B kcal/mol | X kcal/mol | 미완료 대조군 |
|---|---:|---:|---|
| adenine | 미계산 | 미계산 | diflunisal--adenine, diflunisal--diflunisal, adenine--adenine |
| caffeine | -5.420 | 미계산 | diflunisal--diflunisal, caffeine--caffeine |
| theobromine | 미계산 | 미계산 | diflunisal--theobromine, diflunisal--diflunisal, theobromine--theobromine |
| xanthine | 미계산 | 미계산 | diflunisal--xanthine, diflunisal--diflunisal, xanthine--xanthine |

### all_available_poses

| 공형성제 | B kcal/mol | X kcal/mol | 미완료 대조군 |
|---|---:|---:|---|
| adenine | 미계산 | 미계산 | diflunisal--adenine, diflunisal--diflunisal, adenine--adenine |
| caffeine | -5.420 | 미계산 | diflunisal--diflunisal, caffeine--caffeine |
| theobromine | 미계산 | 미계산 | diflunisal--theobromine, diflunisal--diflunisal, theobromine--theobromine |
| xanthine | 미계산 | 미계산 | diflunisal--xanthine, diflunisal--diflunisal, xanthine--xanthine |

B=E(AB)−E(A)−E(B), X=E(AB)−½E(AA)−½E(BB). 전자 에너지이며 ZPE·열·용매·결정 격자 자유에너지가 아니다. 기본 후보와 추가 자세를 분리하고, 확보된 자세 수가 다른 한계와 중복 가능성을 유지한다.

## 반도체

현재 상태: not_started.

[반도체 원시 비교·정밀도 기록](week3-materials.json). 24원자 계면/외측 치환 비교는 이번 범위에 포함하지 않는다.

## 실행 기록

[현재 실행 상태](week3-run-state.json), [약물 상세 결과](week3-drug-comparison.json), [입력과 선택 규칙](../experiments/week3/drug/panel-plan.json).
실패·시간/메모리 한도·정밀도 미달이면 제어기는 중단하고 검토 사유를 기록한다. 자동 재시도하지 않는다.
