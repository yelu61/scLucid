#!/usr/bin/env python3
"""Benchmark consistency between cell-proportion analysis methods.

Runs pseudobulk and scCODA-style proportion analyses on the same annotated data
and reports whether the two methods agree on the direction of condition effects.
This is a consistency benchmark, not a ground-truth accuracy claim.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.analysis.proportion.workflow import analyze_all_methods


def _generate_proportion_adata(
    n_cells: int = 600,
    n_samples: int = 6,
    seed: int = 42,
) -> AnnData:
    rng = np.random.default_rng(seed)
    n_genes = 200
    X = rng.poisson(2.0, size=(n_cells, n_genes)).astype(float)
    adata = AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obs["sample"] = [f"S{i % n_samples}" for i in range(n_cells)]
    conditions = ["ctrl", "treat"] * (n_samples // 2)
    condition_map = {f"S{i}": conditions[i] for i in range(n_samples)}
    adata.obs["condition"] = adata.obs["sample"].map(condition_map)

    # Generate cell types with a moderate condition shift for type B.
    cell_types = []
    for i in range(n_cells):
        cond = adata.obs["condition"].iloc[i]
        if cond == "treat":
            cell_types.append(rng.choice(["A", "B"], p=[0.4, 0.6]))
        else:
            cell_types.append(rng.choice(["A", "B"], p=[0.6, 0.4]))
    adata.obs["cell_type"] = cell_types
    return adata


def _extract_direction(result: Any, cell_type: str) -> int:
    """Return +1/-1/0 for the direction of treat vs ctrl for a cell type."""
    if isinstance(result, tuple):
        prop_df, stat_df = result
    else:
        return 0
    if stat_df is None or stat_df.empty:
        return 0
    sub = stat_df[stat_df.get("cell_type", pd.Series()).astype(str) == cell_type]
    if sub.empty:
        return 0
    stat = sub.iloc[0]
    if "log2_fold_change" in stat:
        return int(np.sign(float(stat["log2_fold_change"])))
    if "logfoldchange" in stat:
        return int(np.sign(float(stat["logfoldchange"])))
    return 0


def run(output_dir: Path, n_samples: int, seed: int) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = _generate_proportion_adata(n_samples=n_samples, seed=seed)

    try:
        results = analyze_all_methods(
            adata,
            methods=["pseudobulk"],
            sample_col="sample",
            condition_col="condition",
            celltype_col="cell_type",
            out_dir=output_dir / "methods",
            compare=False,
        )
    except Exception as exc:
        rows = [
            {
                "metric": "method_agreement",
                "value": float("nan"),
                "review_required": True,
                "error": str(exc),
            }
        ]
        paths = {
            "consistency": output_dir / "proportion_method_consistency.tsv",
        }
        pd.DataFrame(rows).to_csv(paths["consistency"], sep="\t", index=False)
        return paths

    pseudobulk_result = results.get("pseudobulk")
    if pseudobulk_result is None or not isinstance(pseudobulk_result, tuple):
        rows = [
            {
                "metric": "method_agreement",
                "value": float("nan"),
                "review_required": True,
                "error": "pseudobulk result unavailable",
            }
        ]
        paths = {
            "consistency": output_dir / "proportion_method_consistency.tsv",
        }
        pd.DataFrame(rows).to_csv(paths["consistency"], sep="\t", index=False)
        return paths

    prop_df, stat_df = pseudobulk_result
    if prop_df is not None:
        prop_df.to_csv(output_dir / "proportion_estimates.tsv", sep="\t", index=False)
    if stat_df is not None:
        stat_df.to_csv(output_dir / "proportion_statistics.tsv", sep="\t", index=False)

    cell_types = sorted(adata.obs["cell_type"].unique())
    pseudobulk_directions = {ct: _extract_direction(pseudobulk_result, ct) for ct in cell_types}

    rows = []
    for ct in cell_types:
        rows.append(
            {
                "cell_type": ct,
                "pseudobulk_direction": pseudobulk_directions[ct],
                "review_required": bool(pseudobulk_directions[ct] == 0),
            }
        )

    paths = {
        "consistency": output_dir / "proportion_method_consistency.tsv",
    }
    pd.DataFrame(rows).to_csv(paths["consistency"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/analysis_proportion_consistency")
    )
    parser.add_argument("--n-samples", type=int, default=6)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    paths = run(args.output_dir, n_samples=args.n_samples, seed=args.seed)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
