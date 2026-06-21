from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from diabetes_trajectory_clustering import run_from_config


def main() -> int:
    parser = argparse.ArgumentParser(description="Run diabetes trajectory or diagnosis-time clustering.")
    parser.add_argument("--config", default="configs/example_trajectory_config.json", help="Path to JSON config file.")
    parser.add_argument("--data", default=None, help="Optional CSV path override.")
    parser.add_argument("--out", default=None, help="Optional output directory override.")
    parser.add_argument("--mode", choices=["trajectory", "diagnosis"], default=None, help="Optional analysis mode override.")
    parser.add_argument("--k", type=int, default=None, help="Optional primary k override.")
    args = parser.parse_args()

    overrides = {
        "data_path": args.data,
        "out_dir": args.out,
        "analysis_mode": args.mode,
        "k_primary": args.k,
    }
    result = run_from_config(args.config, **overrides)
    summary = {
        "raw_shape": result["raw_shape"],
        "clean_shape": result["clean_shape"],
        "cluster_base_vars": result["cluster_base_vars"],
        "feature_vars": result["feature_vars"],
        "cluster_sizes": result["cluster_sizes"].reset_index().to_dict(orient="records"),
        "out_dir": str(result["out_dir"]),
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
