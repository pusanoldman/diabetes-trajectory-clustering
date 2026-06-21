# Project Summary

## 프로젝트명

당뇨병 trajectory clustering

## 핵심 연구 주제

당뇨병이 발생하기까지 당뇨병 관련 주요 임상변수의 trajectory를 이용해 환자 아형을 나눌 수 있는지 확인한다. 나아가 어떤 아형에서 혈당 상승이 빠르게 일어나고 조기 중재가 필요한지, 여러 임상변수 조합으로 아형을 예측할 수 있는지 검토한다.

## 데이터 개요

- 전체 10020명 중 기존 당뇨병 환자를 제외한 7850명
- AS1부터 AS10까지 약 2년 간격 반복 측정
- 최대 약 20년 추적
- 추적 중 당뇨병 발생자: 2066명
- `status_당뇨병 = 1`: 추적 중 당뇨병 발생
- `tt_당뇨병`: 당뇨병 발생 시점 또는 마지막 추적 시점까지의 기간
- `마지막추적기수`: 마지막으로 관찰된 AS 기수

## 1차 분석 아이디어

당뇨병 발생자 2066명을 대상으로 먼저 clustering을 수행한다.

첫 접근은 진단 시점 또는 진단 직전 기수의 임상변수를 사용한다. 예를 들어 어떤 대상자가 `마지막추적기수 = 8`에서 당뇨병이 발생했다면 `AS8_BMI`, `AS8_허리둘레` 같은 AS8 변수를 진단 시점 feature로 사용한다.

## 후보 clustering 변수

연구 메모에서 우선 고려한 변수는 다음과 같다.

- HOMA2-B 또는 HOMA-B
- HOMA2-IR 또는 HOMA-IR
- 공복인슐린
- BMI
- 허리둘레
- HDL
- TG
- HbA1c
- 당뇨병 진단 시 연령

현재 코드 예시는 실제 데이터에 이미 존재하는 변수명을 고려해 `HBA1C`, `BMI`, `HOMA_IR`, `HOMA_B`를 기본값으로 둔다. 데이터에 `HOMA2_B`, `HOMA2_IR`가 계산되어 있다면 config의 `cluster_base_vars`에 추가하면 된다.

## Trajectory 분석 아이디어

진단 전 일정 기간, 예를 들어 -10년부터 진단 시점까지의 반복 측정값으로 feature를 만든다.

현재 코드가 지원하는 summary feature:

- `last`: 진단 시점 또는 진단 직전 마지막 관찰값
- `delta`: window 안 첫 관찰값과 마지막 관찰값의 차이
- `slope`: 전체 window의 선형 기울기
- `slope_last5`: 진단 전 5년 이내 기울기
- `auc_mean`: 관찰 기간 평균 AUC

grid feature도 지원한다. 예를 들어 -10, -8, -6, -4, -2, 0년의 보간값을 feature로 만들 수 있다.

## 기대 subtype

문헌과 연구 메모를 바탕으로 다음과 같은 profile을 탐색한다.

- Beta-cell deficient type: HOMA-B 또는 HOMA2-B가 상대적으로 낮은 유형
- Obesity/insulin resistance type: BMI, 허리둘레, 공복인슐린, HOMA-IR 또는 HOMA2-IR이 높은 유형
- Lipodystrophy-like type: BMI/허리둘레는 높지 않지만 인슐린저항성 및 TG가 높은 유형
- Liver/lipid type: 비비만 또는 상대적 비비만 상태에서 지질 및 간/대사 축이 두드러지는 유형

아시아 비비만 당뇨 특성을 고려하면 obesity type과 liver/lipid type이 주요하게 관찰될 가능성이 있다.

## 고려사항

- BMI와 허리둘레는 성별 차이가 커서 남녀 분리 clustering도 검토한다.
- HOMA2 계산값이 없을 경우, 기존 HOMA-IR/HOMA-B를 우선 사용하거나 별도 계산 후 상관 및 분포를 검증한다.
- HOMA, insulin, TG 등 오른쪽 꼬리가 긴 변수는 log 또는 log1p 변환을 검토한다.
- 결측과 이상값이 많은 변수는 clustering 안정성을 크게 흔들 수 있으므로 결측률, winsorization, 민감도 분석을 남겨야 한다.
- k 선택은 silhouette 하나만 보지 말고 cluster 크기, 임상 해석 가능성, trajectory plot을 함께 본다.

## 참고 문헌 방향

- Ahlqvist et al., The Lancet Diabetes & Endocrinology, 2018: adult-onset diabetes data-driven subgrouping
- Udler et al., PLOS Medicine, 2018: genetic loci and metabolic trait 기반 soft clustering
- Diabetes Care, 2026: T2D subtype 발생 전 pathophysiological risk factor trajectory
- 한국/아시아 cohort 기반 T2D subtype 연구
