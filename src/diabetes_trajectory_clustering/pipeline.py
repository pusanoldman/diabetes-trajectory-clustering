from __future__ import annotations

import json
import os
import re
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", str(Path.cwd() / ".matplotlib_cache"))
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer, SimpleImputer
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

try:
    from sklearn.experimental import enable_iterative_imputer  # noqa: F401
    from sklearn.impute import IterativeImputer
except Exception:  # pragma: no cover - depends on sklearn build
    IterativeImputer = None

warnings.filterwarnings("ignore", category=FutureWarning)


DEFAULT_CONFIG: dict[str, Any] = {
    "data_path": "data/data.csv",
    "out_dir": "outputs",
    "id_col": "RID",
    "status_col": "status_당뇨병",
    "dx_time_col": "tt_당뇨병",
    "incident_dm_value": 1,
    "analysis_mode": "trajectory",
    "cluster_base_vars": ["HBA1C", "BMI", "HOMA_IR", "HOMA_B"],
    "static_vars": [],
    "include_age_at_dx": False,
    "age_base_var": "연령",
    "k_primary": 4,
    "k_range": [2, 3, 4, 5, 6, 7],
    "cluster_method": "kmeans",
    "kmeans_n_init": 10,
    "kmeans_max_iter": 300,
    "kmeans_tol": 1e-4,
    "gmm_n_init": 5,
    "silhouette_sample_size": 600,
    "window_years": 10,
    "min_visits": 3,
    "include_diagnosis_visit": True,
    "max_post_dx_years": 0.25,
    "feature_mode": "summary",
    "trajectory_features": ["last", "delta", "slope", "slope_last5", "auc_mean"],
    "grid_times": [-10, -8, -6, -4, -2, 0],
    "allow_grid_extrapolation": False,
    "transforms": {
        "HOMA_IR": "log1p",
        "TG": "log1p",
        "공복인슐린": "log1p",
        "60분인슐린": "log1p",
        "120분인슐린": "log1p",
        "AUCinsulin": "log1p",
    },
    "nonpositive_as_na": ["HOMA_IR", "HOMA_B"],
    "winsorize": True,
    "winsorize_q": [0.01, 0.99],
    "max_row_missing": 0.50,
    "max_col_missing": 0.80,
    "imputer_method": "median",
    "time_bin_width": 2.0,
    "random_state": 2026,
    "save_outputs": True,
}


def load_config(config_path: str | Path | None = None, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    config = DEFAULT_CONFIG.copy()
    if config_path:
        with Path(config_path).open("r", encoding="utf-8") as f:
            config.update(json.load(f))
    if overrides:
        config.update({k: v for k, v in overrides.items() if v is not None})
    return config


def make_json_ready(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {k: make_json_ready(v) for k, v in value.items()}
    if isinstance(value, list):
        return [make_json_ready(item) for item in value]
    return value


def read_csv_safely(path: str | Path) -> pd.DataFrame:
    path = Path(path).expanduser()
    last_error: Exception | None = None
    for encoding in ["utf-8-sig", "cp949", "euc-kr", "latin1"]:
        try:
            return pd.read_csv(path, encoding=encoding, low_memory=False)
        except Exception as exc:
            last_error = exc
    raise RuntimeError(f"Could not read CSV file: {path}") from last_error


def clean_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out = out.replace(r"^\s*$", np.nan, regex=True)
    for col in out.columns:
        if out[col].dtype == "object":
            series = out[col].astype(str).str.strip().replace({"nan": np.nan, "None": np.nan, "": np.nan})
            converted = pd.to_numeric(series, errors="coerce")
            denom = series.notna().sum()
            if denom > 0 and converted.notna().sum() / denom >= 0.70:
                out[col] = converted
            else:
                out[col] = series
    return out


def infer_as_structure(df: pd.DataFrame) -> tuple[dict[int, dict[str, str]], list[str]]:
    pattern = re.compile(r"^AS(\d+)_(.+)$")
    visit_map: defaultdict[int, dict[str, str]] = defaultdict(dict)
    bases: set[str] = set()
    for col in df.columns:
        match = pattern.match(str(col))
        if match:
            visit_num = int(match.group(1))
            base = match.group(2)
            visit_map[visit_num][base] = col
            bases.add(base)
    return dict(visit_map), sorted(bases)


def resolve_base_vars(requested: list[str], available_bases: list[str]) -> list[str]:
    lower_map = {x.lower(): x for x in available_bases}
    resolved: list[str] = []
    missing: list[str] = []
    for raw_var in requested:
        var = str(raw_var)
        match = re.match(r"^AS\d+_(.+)$", var)
        if match:
            var = match.group(1)
        if var in available_bases:
            resolved.append(var)
        elif var.lower() in lower_map:
            resolved.append(lower_map[var.lower()])
        else:
            missing.append(raw_var)
    if missing:
        preview = ", ".join(available_bases[:80])
        raise ValueError(f"Missing AS base variables: {missing}\nAvailable examples: {preview}")
    return list(dict.fromkeys(resolved))


def get_fu_year(df: pd.DataFrame, visit_num: int) -> pd.Series:
    if visit_num == 1:
        return pd.Series(0.0, index=df.index)
    candidates = [f"AS{visit_num}_FU_YEAR", f"AS{visit_num}_검진개월수", f"AS{visit_num}_검진개월"]
    for col in candidates:
        if col in df.columns:
            values = pd.to_numeric(df[col], errors="coerce")
            if "개월" in col:
                values = values / 12.0
            return values
    return pd.Series(np.nan, index=df.index)


def make_long(df: pd.DataFrame, base_vars: list[str], visit_map: dict[int, dict[str, str]], config: dict[str, Any]) -> pd.DataFrame:
    id_col = config["id_col"]
    static_cols = [id_col]
    for col in [config["status_col"], config["dx_time_col"]]:
        if col in df.columns and col not in static_cols:
            static_cols.append(col)

    pieces = []
    for visit_num in sorted(visit_map):
        temp = df[static_cols].copy()
        temp["visit"] = f"AS{visit_num}"
        temp["visit_num"] = visit_num
        temp["FU_YEAR"] = get_fu_year(df, visit_num)
        for base in base_vars:
            col = visit_map.get(visit_num, {}).get(base)
            temp[base] = pd.to_numeric(df[col], errors="coerce") if col else np.nan
        pieces.append(temp)

    long = pd.concat(pieces, axis=0, ignore_index=True)
    dx_time_col = config["dx_time_col"]
    if dx_time_col in long.columns:
        long["time_to_dx"] = pd.to_numeric(long["FU_YEAR"], errors="coerce") - pd.to_numeric(long[dx_time_col], errors="coerce")
    else:
        long["time_to_dx"] = np.nan
    return long


def transformed_name(var: str, transforms: dict[str, str]) -> str:
    transform = transforms.get(var)
    if transform == "log":
        return f"log_{var}"
    if transform == "log1p":
        return f"log1p_{var}"
    return var


def apply_var_qc_and_transform(long: pd.DataFrame, base_vars: list[str], config: dict[str, Any]) -> tuple[pd.DataFrame, list[str]]:
    out = long.copy()
    transforms = config.get("transforms", {}) or {}
    nonpositive_as_na = set(config.get("nonpositive_as_na", []) or [])
    q_low, q_high = config.get("winsorize_q", [0.01, 0.99])
    feature_vars: list[str] = []

    for var in base_vars:
        if var not in out.columns:
            continue
        x = pd.to_numeric(out[var], errors="coerce").astype(float)
        if var in nonpositive_as_na:
            x = x.mask(x <= 0)

        transform = transforms.get(var)
        out_name = transformed_name(var, transforms)
        if transform == "log":
            x = np.log(x.mask(x <= 0))
        elif transform == "log1p":
            x = np.log1p(x.mask(x <= -1))
        elif transform not in [None, "none", ""]:
            raise ValueError(f"Unsupported transform for {var}: {transform}")

        if config.get("winsorize", True):
            lo, hi = x.quantile(q_low), x.quantile(q_high)
            if np.isfinite(lo) and np.isfinite(hi) and lo < hi:
                x = x.clip(lo, hi)
        out[out_name] = x
        feature_vars.append(out_name)

    return out, feature_vars


def safe_slope(t: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(t) & np.isfinite(y)
    t, y = t[mask], y[mask]
    if len(t) < 2 or np.nanstd(t) == 0:
        return np.nan
    return float(np.polyfit(t, y, 1)[0])


def safe_auc_mean(t: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(t) & np.isfinite(y)
    t, y = t[mask], y[mask]
    if len(t) < 2:
        return np.nan
    order = np.argsort(t)
    t, y = t[order], y[order]
    span = t.max() - t.min()
    if span <= 0:
        return np.nan
    if hasattr(np, "trapezoid"):
        auc = np.trapezoid(y, t)
    else:
        auc = np.trapz(y, t)
    return float(auc / span)


def summarize_one_variable(group: pd.DataFrame, var: str, features: list[str]) -> dict[str, float]:
    x = group[["time_to_dx", var]].dropna().sort_values("time_to_dx")
    output: dict[str, float] = {}
    if x.empty:
        return {f"{var}__{feature}": np.nan for feature in features}

    t = x["time_to_dx"].to_numpy(dtype=float)
    y = x[var].to_numpy(dtype=float)
    for feature in features:
        key = f"{var}__{feature}"
        if feature == "first":
            output[key] = y[0]
        elif feature == "last":
            output[key] = y[-1]
        elif feature == "mean":
            output[key] = float(np.nanmean(y))
        elif feature == "sd":
            output[key] = float(np.nanstd(y, ddof=1)) if len(y) >= 2 else np.nan
        elif feature == "min":
            output[key] = float(np.nanmin(y))
        elif feature == "max":
            output[key] = float(np.nanmax(y))
        elif feature == "delta":
            output[key] = y[-1] - y[0] if len(y) >= 2 else np.nan
        elif feature == "slope":
            output[key] = safe_slope(t, y)
        elif feature == "auc_mean":
            output[key] = safe_auc_mean(t, y)
        elif feature == "delta_last5":
            mask = t >= -5
            output[key] = y[mask][-1] - y[mask][0] if mask.sum() >= 2 else np.nan
        elif feature == "slope_last5":
            mask = t >= -5
            output[key] = safe_slope(t[mask], y[mask]) if mask.sum() >= 2 else np.nan
        else:
            raise ValueError(f"Unsupported trajectory feature: {feature}")
    return output


def grid_one_variable(group: pd.DataFrame, var: str, grid_times: list[float], allow_extrapolation: bool) -> dict[str, float]:
    x = group[["time_to_dx", var]].dropna().sort_values("time_to_dx")
    output = {f"{var}__grid_{gt:g}y": np.nan for gt in grid_times}
    if x.empty:
        return output
    tmp = x.groupby("time_to_dx", as_index=False)[var].mean().sort_values("time_to_dx")
    t = tmp["time_to_dx"].to_numpy(dtype=float)
    y = tmp[var].to_numpy(dtype=float)
    if len(t) == 1:
        for gt in grid_times:
            if abs(gt - t[0]) < 1e-8:
                output[f"{var}__grid_{gt:g}y"] = y[0]
        return output
    for gt in grid_times:
        if allow_extrapolation or (gt >= t.min() and gt <= t.max()):
            output[f"{var}__grid_{gt:g}y"] = float(np.interp(gt, t, y))
    return output


def build_trajectory_features(long: pd.DataFrame, feature_vars: list[str], config: dict[str, Any]) -> pd.DataFrame:
    rows = []
    for rid, group in long.groupby(config["id_col"]):
        row: dict[str, Any] = {config["id_col"]: rid}
        if config["feature_mode"] in ["summary", "both"]:
            for var in feature_vars:
                row.update(summarize_one_variable(group, var, config["trajectory_features"]))
        if config["feature_mode"] in ["grid", "both"]:
            for var in feature_vars:
                row.update(grid_one_variable(group, var, config["grid_times"], config["allow_grid_extrapolation"]))
        rows.append(row)
    return pd.DataFrame(rows).set_index(config["id_col"])


def build_diagnosis_features(long: pd.DataFrame, feature_vars: list[str], config: dict[str, Any]) -> pd.DataFrame:
    rows = []
    for rid, group in long.groupby(config["id_col"]):
        group = group.dropna(subset=["time_to_dx"]).copy()
        if group.empty:
            continue
        group["abs_time"] = group["time_to_dx"].abs()
        row: dict[str, Any] = {config["id_col"]: rid}
        for var in feature_vars:
            nearest = group.dropna(subset=[var]).sort_values("abs_time")
            row[f"{var}__diagnosis"] = nearest[var].iloc[0] if not nearest.empty else np.nan
        rows.append(row)
    return pd.DataFrame(rows).set_index(config["id_col"])


def extract_static_features(df: pd.DataFrame, ids: pd.Index, config: dict[str, Any]) -> pd.DataFrame:
    static_vars = config.get("static_vars", []) or []
    id_col = config["id_col"]
    if not static_vars:
        return pd.DataFrame(index=pd.Index(ids, name=id_col))
    base = df.set_index(id_col)
    frames = []
    for var in static_vars:
        col = var if var in df.columns else f"AS1_{var}" if f"AS1_{var}" in df.columns else None
        if col is None:
            print(f"Warning: static variable not found: {var}")
            continue
        series = base.loc[ids, col]
        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().sum() / max(1, series.notna().sum()) >= 0.70:
            frames.append(numeric.rename(var).to_frame())
        else:
            frames.append(pd.get_dummies(series.astype("category"), prefix=str(var), dummy_na=True))
    return pd.concat(frames, axis=1) if frames else pd.DataFrame(index=pd.Index(ids, name=id_col))


def extract_age_at_dx(long: pd.DataFrame, config: dict[str, Any]) -> pd.Series:
    age_col = transformed_name(config.get("age_base_var", "연령"), config.get("transforms", {}) or {})
    if age_col not in long.columns:
        return pd.Series(dtype=float, name="age_dx")
    group = long.dropna(subset=["time_to_dx", age_col]).copy()
    if group.empty:
        return pd.Series(dtype=float, name="age_dx")
    group["abs_time"] = group["time_to_dx"].abs()
    idx = group.sort_values([config["id_col"], "abs_time"]).groupby(config["id_col"]).head(1).index
    series = group.loc[idx, [config["id_col"], age_col]].set_index(config["id_col"])[age_col]
    series.name = "age_dx"
    return series


def preprocess_feature_matrix(features: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    x = features.copy().apply(pd.to_numeric, errors="coerce")
    col_missing = x.isna().mean().sort_values(ascending=False)
    keep_cols = col_missing[col_missing <= config["max_col_missing"]].index.tolist()
    dropped_cols = [col for col in x.columns if col not in keep_cols]
    x = x[keep_cols]

    row_missing = x.isna().mean(axis=1)
    keep_rows = row_missing[row_missing <= config["max_row_missing"]].index
    dropped_rows = [idx for idx in x.index if idx not in keep_rows]
    x = x.loc[keep_rows]
    if x.shape[1] == 0:
        raise ValueError("No features remain after column-missing filtering.")
    if x.shape[0] < 5:
        raise ValueError("Too few participants remain after row-missing filtering.")

    method = config["imputer_method"]
    if method == "median":
        imputer = SimpleImputer(strategy="median")
    elif method == "knn":
        imputer = KNNImputer(n_neighbors=5)
    elif method == "iterative":
        if IterativeImputer is None:
            raise ImportError("IterativeImputer is unavailable. Use imputer_method='median'.")
        imputer = IterativeImputer(random_state=config["random_state"], max_iter=20)
    else:
        raise ValueError("imputer_method must be one of: median, knn, iterative.")

    x_imputed = pd.DataFrame(imputer.fit_transform(x), index=x.index, columns=x.columns)
    scaler = StandardScaler()
    x_scaled = pd.DataFrame(scaler.fit_transform(x_imputed), index=x.index, columns=x.columns)
    report = {
        "n_initial": int(features.shape[0]),
        "p_initial": int(features.shape[1]),
        "n_final": int(x_scaled.shape[0]),
        "p_final": int(x_scaled.shape[1]),
        "dropped_rows": [str(x) for x in dropped_rows],
        "dropped_cols": dropped_cols,
        "col_missing": col_missing,
    }
    return x_scaled, x_imputed, report


def fit_clustering(x_scaled: pd.DataFrame, k: int, config: dict[str, Any], random_state: int | None = None) -> tuple[np.ndarray, Any]:
    random_state = config["random_state"] if random_state is None else random_state
    method = config["cluster_method"]
    x_arr = x_scaled.to_numpy()
    if method == "kmeans":
        model = KMeans(
            n_clusters=k,
            n_init=config["kmeans_n_init"],
            max_iter=config["kmeans_max_iter"],
            tol=config["kmeans_tol"],
            random_state=random_state,
        )
        labels = model.fit_predict(x_arr)
    elif method == "gmm":
        model = GaussianMixture(n_components=k, covariance_type="full", n_init=config["gmm_n_init"], random_state=random_state)
        labels = model.fit_predict(x_arr)
    elif method == "hierarchical":
        model = AgglomerativeClustering(n_clusters=k, linkage="ward")
        labels = model.fit_predict(x_arr)
    else:
        raise ValueError("cluster_method must be one of: kmeans, gmm, hierarchical.")
    return labels.astype(int) + 1, model


def evaluate_k_range(x_scaled: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    rows = []
    n = x_scaled.shape[0]
    for k in config["k_range"]:
        if k < 2 or k >= n:
            continue
        labels, model = fit_clustering(x_scaled, int(k), config, random_state=config["random_state"] + int(k))
        counts = pd.Series(labels).value_counts().sort_index()
        row = {
            "k": int(k),
            "method": config["cluster_method"],
            "silhouette": silhouette_score(
                x_scaled,
                labels,
                sample_size=min(config["silhouette_sample_size"], n),
                random_state=config["random_state"],
            ),
            "calinski_harabasz": calinski_harabasz_score(x_scaled, labels),
            "davies_bouldin": davies_bouldin_score(x_scaled, labels),
            "min_cluster_n": int(counts.min()),
            "max_cluster_n": int(counts.max()),
        }
        if config["cluster_method"] == "kmeans":
            row["inertia"] = float(model.inertia_)
        if config["cluster_method"] == "gmm":
            row["bic"] = float(model.bic(x_scaled))
            row["aic"] = float(model.aic(x_scaled))
        rows.append(row)
    return pd.DataFrame(rows)


def save_basic_figures(
    x_scaled: pd.DataFrame,
    x_imputed: pd.DataFrame,
    assignments: pd.DataFrame,
    cluster_sizes: pd.DataFrame,
    cluster_zmeans: pd.DataFrame,
    long_qc: pd.DataFrame,
    feature_vars: list[str],
    out_dir: Path,
    config: dict[str, Any],
) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(cluster_sizes.index.astype(str), cluster_sizes["n"])
    ax.set_xlabel("Cluster")
    ax.set_ylabel("N")
    ax.set_title("Cluster size")
    for idx, row in cluster_sizes.iterrows():
        ax.text(int(idx) - 1, row["n"], f"{int(row['n'])}\n({row['percent']:.1f}%)", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_cluster_sizes.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    pca = PCA(n_components=2, random_state=config["random_state"])
    coords = pca.fit_transform(x_scaled)
    pc = pd.DataFrame(coords, index=x_scaled.index, columns=["PC1", "PC2"]).join(assignments)
    fig, ax = plt.subplots(figsize=(7, 5))
    for cluster, group in pc.groupby("cluster"):
        ax.scatter(group["PC1"], group["PC2"], s=20, alpha=0.75, label=f"Cluster {cluster}")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.set_title("PCA visualization")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_pca_clusters.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    pc.to_csv(out_dir / "pca_coordinates.csv", encoding="utf-8-sig")

    heatmap = cluster_zmeans.copy()
    if heatmap.shape[1] > 40:
        heatmap = heatmap[heatmap.var(axis=0).sort_values(ascending=False).head(40).index]
    fig, ax = plt.subplots(figsize=(max(8, 0.35 * heatmap.shape[1]), max(4, 0.25 * heatmap.shape[1])))
    image = ax.imshow(heatmap.values, aspect="auto")
    ax.set_yticks(np.arange(heatmap.shape[0]))
    ax.set_yticklabels([f"Cluster {c}" for c in heatmap.index])
    ax.set_xticks(np.arange(heatmap.shape[1]))
    ax.set_xticklabels(heatmap.columns, rotation=90, fontsize=8)
    ax.set_title("Cluster profile heatmap: standardized feature means")
    fig.colorbar(image, ax=ax, fraction=0.02, pad=0.02)
    fig.tight_layout()
    fig.savefig(out_dir / "fig_feature_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    trajectory_data = long_qc.merge(assignments.reset_index(), on=config["id_col"], how="inner")
    trajectory_data = trajectory_data.dropna(subset=["time_to_dx"])
    trajectory_data["time_bin"] = (np.round(trajectory_data["time_to_dx"] / config["time_bin_width"]) * config["time_bin_width"]).astype(float)
    trajectory_data = trajectory_data[
        (trajectory_data["time_bin"] >= -config["window_years"]) & (trajectory_data["time_bin"] <= config["max_post_dx_years"])
    ]
    for var in feature_vars:
        if var not in trajectory_data.columns:
            continue
        summary = trajectory_data.dropna(subset=[var]).groupby(["cluster", "time_bin"])[var].agg(["mean", "count", "std"]).reset_index()
        if summary.empty:
            continue
        summary["se"] = summary["std"] / np.sqrt(summary["count"].clip(lower=1))
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for cluster, group in summary.groupby("cluster"):
            group = group.sort_values("time_bin")
            ax.plot(group["time_bin"], group["mean"], marker="o", label=f"Cluster {cluster}")
            ax.fill_between(group["time_bin"], group["mean"] - 1.96 * group["se"], group["mean"] + 1.96 * group["se"], alpha=0.15)
        ax.axvline(0, linestyle="--", linewidth=1)
        ax.set_xlabel("Years before diagnosis")
        ax.set_ylabel(var)
        ax.set_title(f"Trajectory by cluster: {var}")
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.25)
        fig.tight_layout()
        safe = re.sub(r"[^A-Za-z0-9가-힣_\-]+", "_", var)[:80]
        fig.savefig(out_dir / f"fig_trajectory_{safe}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)


def run_pipeline(config: dict[str, Any]) -> dict[str, Any]:
    data_path = Path(config["data_path"]).expanduser()
    out_dir = Path(config["out_dir"]).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = read_csv_safely(data_path)
    df = clean_dataframe(raw)
    visit_map, available_bases = infer_as_structure(df)
    cluster_base_vars = resolve_base_vars(config["cluster_base_vars"], available_bases)

    long_base_vars = list(cluster_base_vars)
    if config.get("include_age_at_dx") and config.get("age_base_var") in available_bases:
        long_base_vars.append(config["age_base_var"])
    long_base_vars = list(dict.fromkeys(long_base_vars))

    long0 = make_long(df, long_base_vars, visit_map, config)
    if config["status_col"] in long0.columns:
        long0 = long0[long0[config["status_col"]] == config["incident_dm_value"]].copy()

    if config["analysis_mode"] == "trajectory":
        upper = config["max_post_dx_years"] if config["include_diagnosis_visit"] else 0
        long_window = long0[(long0["time_to_dx"] >= -config["window_years"]) & (long0["time_to_dx"] <= upper)].copy()
    elif config["analysis_mode"] == "diagnosis":
        long_window = long0[(long0["time_to_dx"] >= -config["window_years"]) & (long0["time_to_dx"] <= config["max_post_dx_years"])].copy()
    else:
        raise ValueError("analysis_mode must be 'trajectory' or 'diagnosis'.")

    long_qc, transformed_vars = apply_var_qc_and_transform(long_window, long_base_vars, config)
    feature_vars = [transformed_name(var, config.get("transforms", {}) or {}) for var in cluster_base_vars]
    long_qc["has_any_feature"] = long_qc[feature_vars].notna().any(axis=1)
    visit_counts = long_qc[long_qc["has_any_feature"]].groupby(config["id_col"])["visit_num"].nunique()
    eligible_ids = visit_counts[visit_counts >= config["min_visits"]].index
    long_qc = long_qc[long_qc[config["id_col"]].isin(eligible_ids)].copy()

    if config["analysis_mode"] == "trajectory":
        feature_matrix = build_trajectory_features(long_qc, feature_vars, config)
    else:
        feature_matrix = build_diagnosis_features(long_qc, feature_vars, config)

    static_features = extract_static_features(df, feature_matrix.index, config)
    if static_features.shape[1] > 0:
        feature_matrix = feature_matrix.join(static_features, how="left")
    if config.get("include_age_at_dx"):
        feature_matrix = feature_matrix.join(extract_age_at_dx(long_qc, config), how="left")

    x_scaled, x_imputed, preprocess_report = preprocess_feature_matrix(feature_matrix, config)
    k_eval = evaluate_k_range(x_scaled, config)
    labels, _ = fit_clustering(x_scaled, int(config["k_primary"]), config)
    assignments = pd.DataFrame({config["id_col"]: x_scaled.index, "cluster": labels}).set_index(config["id_col"]).sort_values("cluster")
    cluster_sizes = assignments["cluster"].value_counts().sort_index().rename("n").to_frame()
    cluster_sizes["percent"] = cluster_sizes["n"] / cluster_sizes["n"].sum() * 100
    features_labeled = x_imputed.join(assignments, how="inner")
    scaled_labeled = x_scaled.join(assignments, how="inner")
    cluster_feature_means = features_labeled.groupby("cluster").mean(numeric_only=True)
    cluster_feature_zmeans = scaled_labeled.groupby("cluster").mean(numeric_only=True)

    if config.get("save_outputs", True):
        with (out_dir / "run_config_resolved.json").open("w", encoding="utf-8") as f:
            json.dump(make_json_ready(config), f, ensure_ascii=False, indent=2)
        pd.Series(available_bases, name="available_base_variable").to_csv(out_dir / "available_as_base_variables.csv", index=False, encoding="utf-8-sig")
        feature_matrix.to_csv(out_dir / "feature_matrix_before_imputation.csv", encoding="utf-8-sig")
        x_imputed.to_csv(out_dir / "feature_matrix_imputed.csv", encoding="utf-8-sig")
        assignments.reset_index().to_csv(out_dir / "cluster_assignments.csv", index=False, encoding="utf-8-sig")
        cluster_sizes.to_csv(out_dir / "cluster_sizes.csv", encoding="utf-8-sig")
        cluster_feature_means.to_csv(out_dir / "cluster_feature_means.csv", encoding="utf-8-sig")
        cluster_feature_zmeans.to_csv(out_dir / "cluster_feature_zmeans.csv", encoding="utf-8-sig")
        k_eval.to_csv(out_dir / "k_evaluation.csv", index=False, encoding="utf-8-sig")
        long_qc.merge(assignments.reset_index(), on=config["id_col"], how="inner").to_csv(
            out_dir / "long_data_with_clusters.csv", index=False, encoding="utf-8-sig"
        )
        save_basic_figures(
            x_scaled,
            x_imputed,
            assignments,
            cluster_sizes,
            cluster_feature_zmeans,
            long_qc,
            feature_vars,
            out_dir,
            config,
        )

    return {
        "raw_shape": raw.shape,
        "clean_shape": df.shape,
        "available_bases": available_bases,
        "cluster_base_vars": cluster_base_vars,
        "feature_vars": feature_vars,
        "feature_matrix": feature_matrix,
        "x_scaled": x_scaled,
        "x_imputed": x_imputed,
        "assignments": assignments,
        "cluster_sizes": cluster_sizes,
        "cluster_feature_means": cluster_feature_means,
        "cluster_feature_zmeans": cluster_feature_zmeans,
        "k_evaluation": k_eval,
        "preprocess_report": preprocess_report,
        "out_dir": out_dir,
    }


def run_from_config(config_path: str | Path | None = None, **overrides: Any) -> dict[str, Any]:
    config = load_config(config_path, overrides)
    return run_pipeline(config)
