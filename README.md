# Diabetes Trajectory Clustering

당뇨병 발생 전후 임상변수 trajectory와 진단 시점 변수를 이용해 환자 아형을 clustering하는 분석 저장소입니다.

이 저장소에는 실행 코드와 문서만 포함합니다. 개인/연구 데이터(`data.csv`), 원본 `.docx`, 논문 PDF, 결과물은 공개 GitHub에 올리지 않도록 `.gitignore`에 제외했습니다.

## 연구 질문

- 당뇨병이 발생하기까지 주요 임상변수의 trajectory로 DM progressor subtype을 나눌 수 있는가?
- 어떤 subtype에서 혈당 상승이 빠르게 일어나고, 조기 중재가 필요한가?
- baseline 또는 진단 시점의 임상변수 조합으로 subtype을 예측할 수 있는가?

## 데이터 구조

예상 입력 파일은 wide-format CSV입니다.

- 개인 ID: `RID`
- 당뇨병 발생 여부: `status_당뇨병`
- 당뇨병 발생 또는 추적 기간: `tt_당뇨병`
- 반복 측정 변수: `AS1_BMI`, `AS2_BMI`, ..., `AS10_BMI`처럼 `AS{기수}_{변수명}` 구조
- 추적 시간 변수: 기수별 추적 기간을 나타내는 컬럼

현재 연구 메모 기준으로 전체 10020명 중 기존 당뇨병 환자를 제외한 7850명을 최대 20년 추적했고, 추적 중 당뇨병 발생자는 2066명입니다.

## 설치

```bash
git clone https://github.com/pusanoldman/diabetes-trajectory-clustering.git
cd diabetes-trajectory-clustering
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

macOS/Linux에서는 가상환경 활성화 명령만 다릅니다.

```bash
source .venv/bin/activate
```

Jupyter Notebook/Lab은 분석 필수 의존성에서 제외했습니다. 이미 쓰는 Jupyter 환경에서 이 저장소를 열거나, 필요할 때만 별도로 설치하세요.

```bash
pip install notebook ipykernel
```

## 데이터 넣기

로컬에서만 다음 위치에 CSV를 넣습니다.

```text
data/data.csv
```

`data/` 폴더는 `.gitignore` 처리되어 있으므로 실제 데이터는 GitHub에 올라가지 않습니다.

## 바로 실행

Trajectory feature 기반 clustering:

```bash
python scripts/run_clustering.py --config configs/example_trajectory_config.json
```

진단 시점 또는 진단 직전 값 기반 clustering:

```bash
python scripts/run_clustering.py --config configs/example_diagnosis_config.json
```

데이터 경로와 결과 폴더를 직접 지정할 수도 있습니다.

```bash
python scripts/run_clustering.py --data "C:\Users\me\Desktop\data.csv" --out outputs/my_run --mode trajectory --k 4
```

## Jupyter에서 실행

```python
import sys
from pathlib import Path

repo = Path.cwd()
sys.path.insert(0, str(repo / "src"))

from diabetes_trajectory_clustering import run_from_config

result = run_from_config(
    "configs/example_trajectory_config.json",
    data_path="data/data.csv",
    out_dir="outputs/trajectory"
)

result["cluster_sizes"]
result["k_evaluation"]
```

## 주요 설정

설정 파일은 `configs/` 아래 JSON 파일을 수정하면 됩니다.

- `analysis_mode`: `trajectory` 또는 `diagnosis`
- `cluster_base_vars`: clustering에 넣을 AS 반복측정 변수의 base name
- `k_primary`: 최종 cluster 개수
- `k_range`: k 비교 범위
- `cluster_method`: `kmeans`, `gmm`, `hierarchical`
- `window_years`: 진단 전 몇 년까지 볼지
- `min_visits`: trajectory 분석에 필요한 최소 관찰 기수 수
- `transforms`: log/log1p 변환할 변수
- `nonpositive_as_na`: 0 이하 값을 결측 처리할 변수

## `cluster_base_vars`에 넣을 수 있는 base 변수

`cluster_base_vars`에는 `AS1_BMI`, `AS2_BMI`처럼 반복 측정 컬럼에서 `AS{기수}_`를 뺀 base name을 넣습니다. 예를 들어 데이터에 `AS1_TG`, `AS2_TG`가 있으면 config에는 `TG`만 적습니다.

기본 추천 조합:

```python
CLUSTER_VARS = ["HBA1C", "BMI", "HOMA_IR", "HOMA_B"]
```

임상적으로 해석하기 좋은 추가 후보:

```python
["TG", "HDL", "공복혈당", "공복인슐린", "허리둘레", "수축기혈압", "이완기혈압"]
```

OGTT/인슐린 반응 관련 후보:

```python
["60분혈당", "120분혈당", "60분인슐린", "120분인슐린", "AUCglucose", "AUCinsulin"]
```

간/신장/지질 관련 후보:

```python
["ALT", "AST", "BUN", "CREATININE", "TCHL", "TG", "HDL"]
```

파생 지표 후보:

```python
["HOMA_IR", "HOMA_B", "IGI60", "ISI", "ISSI2", "DI", "meanPG", "meanPI"]
```

실제로 clustering에 우선 고려할 base 변수 목록:

```text
BMI
HBA1C
HOMA_B
HOMA_IR
TG
HDL
공복혈당
공복인슐린
허리둘레
수축기혈압
이완기혈압
60분혈당
120분혈당
60분인슐린
120분인슐린
AUCglucose
AUCinsulin
ALT
AST
BUN
CREATININE
TCHL
IGI60
ISI
ISSI2
DI
meanPG
meanPI
```

오른쪽 꼬리가 긴 변수는 config의 `transforms` 설정에 따라 log 변환됩니다. 현재 trajectory config에서는 다음 변수가 `log1p` 변환됩니다.

```json
"transforms": {
  "HOMA_IR": "log1p",
  "TG": "log1p",
  "공복인슐린": "log1p",
  "60분인슐린": "log1p",
  "120분인슐린": "log1p",
  "AUCinsulin": "log1p"
}
```

예를 들어 `CLUSTER_VARS`에 `HOMA_IR`를 넣으면 실제 feature matrix에는 `log1p_HOMA_IR`가 만들어지고, 그 값이 clustering에 들어갑니다. `TG`, `공복인슐린`, `60분인슐린`, `120분인슐린`, `AUCinsulin`도 같은 방식입니다.

0 이하 값 처리는 `nonpositive_as_na`에서 지정합니다.

```json
"nonpositive_as_na": ["HOMA_IR", "HOMA_B"]
```

즉 `HOMA_IR`, `HOMA_B`의 0 이하 값은 결측으로 처리됩니다. `HOMA_B`는 현재 기본 config에서 log 변환하지 않고 0 이하 결측 처리만 합니다. `HOMA_B`도 오른쪽 꼬리가 심해 log 변환하고 싶으면 `transforms`에 `"HOMA_B": "log1p"`를 추가하면 됩니다.

이상치는 IQR 기준으로 제거하지 않고, 변수별 1-99% 분위수 기준으로 winsorization합니다.

```json
"winsorize": true,
"winsorize_q": [0.01, 0.99]
```

즉 하위 1%보다 작은 값은 1% 분위수 값으로 올리고, 상위 99%보다 큰 값은 99% 분위수 값으로 내립니다. 대상자를 삭제하는 방식이 아니라 극단값을 경계값으로 잘라내는 방식입니다.

전처리 순서는 다음과 같습니다.

1. `nonpositive_as_na`에 지정된 변수의 0 이하 값을 결측 처리
2. `transforms`에 지정된 변수에 log/log1p 변환 적용
3. 변수별 1-99% 분위수 winsorization
4. 결측 대체
5. 표준화
6. clustering

주의: 당뇨병 진단 여부, 진단 기준 충족 여부, 방문 시점, 검진일처럼 outcome 또는 시간 정보를 직접 나타내는 변수는 clustering 변수로 넣으면 결과 해석이 왜곡될 수 있습니다. 처음에는 4-8개 정도의 연속형 임상변수로 시작하는 것이 좋습니다.

## 결과물

실행 후 `outputs/` 아래에 다음 파일들이 생성됩니다.

- `available_as_base_variables.csv`: 데이터에서 감지된 AS base 변수 목록
- `feature_matrix_before_imputation.csv`: 결측 대체 전 feature matrix
- `feature_matrix_imputed.csv`: 결측 대체 후 feature matrix
- `k_evaluation.csv`: k별 silhouette, Calinski-Harabasz, Davies-Bouldin 지표
- `cluster_assignments.csv`: 개인별 cluster 배정
- `cluster_sizes.csv`: cluster별 대상자 수
- `cluster_feature_means.csv`: cluster별 feature 평균
- `cluster_feature_zmeans.csv`: 표준화 scale의 cluster profile
- `long_data_with_clusters.csv`: long-format 반복측정 자료와 cluster label
- `fig_*.png`: cluster size, PCA, heatmap, trajectory plot

## 분석 방향

1. 당뇨병 발생자만 선택합니다.
2. 진단 시점 기준 변수 또는 진단 전 trajectory summary feature를 만듭니다.
3. 결측률이 높은 feature와 대상자를 제외합니다.
4. median/KNN/iterative imputation 중 하나로 결측을 대체합니다.
5. feature를 표준화하고 k-means/GMM/hierarchical clustering을 수행합니다.
6. cluster profile과 trajectory plot으로 subtype 해석 가능성을 확인합니다.

자세한 연구 정리는 [docs/project_summary.md](docs/project_summary.md)를 참고하세요.
