#!/usr/bin/env python3
"""Benchmark pseudobulk DE type-I error control on synthetic null data.

Generates data with no true condition effect (same cell-type composition and
expression distribution across conditions) but with multiple biological samples.
Then runs scLucid's pseudobulk DE and reports the fraction of genes with
adjusted p-value < 0.05. Under the null this fraction should be close to 0.05.
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

from scLucid.analysis.config import PseudobulkDEConfig
from scLucid.analysis.differential_expression.de_core import run_pseudobulk_de


def _generate_null_adata(
    n_cells: int = 600,
    n_genes: int = 500,
    n_samples: int = 6,
    seed: int = 42,
) -> AnnData:
    rng = np.random.default_rng(seed)
    X = rng.poisson(2.0, size=(n_cells, n_genes)).astype(float)
    adata = AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obs["sample"] = [f"S{i % n_samples}" for i in range(n_cells)]
    conditions = ["ctrl", "treat"] * (n_samples // 2)
    condition_map = {f"S{i}": conditions[i] for i in range(n_samples)}
    adata.obs["condition"] = adata.obs["sample"].map(condition_map)
    adata.obs["cell_type"] = rng.choice(["A", "B"], size=n_cells)
    adata.layers["counts"] = X.copy()
    return adata


def _fdr_at_alpha(pvals: pd.Series, alpha: float = 0.05) -> float:
    valid = pd.to_numeric(pvals, errors="coerce").dropna()
    if valid.empty:
        return float("nan")
    return float((valid < alpha).mean())


def run(output_dir: Path, n_genes: int, n_samples: int, seed: int) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = _generate_null_adata(n_genes=n_genes, n_samples=n_samples, seed=seed)

    config = PseudobulkDEConfig(
        sample_col="sample",
        condition_key="condition",
        contrasts=[["ctrl", "treat"]],
        groupby=None,
        method="welch_logcpm",
        min_samples_per_condition=2,
    )
    result = run_pseudobulk_de(adata, config=config)

    observed_fdr = _fdr_at_alpha(result.get("padj"), alpha=0.05)
    rows = [
        {
            "metric": "observed_fdr_at_alpha_0.05",
            "value": observed_fdr,
            "expected_under_null": 0.05,
            "review_required": not (0.03 <= observed_fdr <= 0.10) if pd.notna(observed_fdr) else True,
            "n_genes_tested": int(result.shape[0]) if not result.empty else 0,
            "n_samples": n_samples,
        }
    ]

    paths = {
        "de_results": output_dir / "pseudobulk_de_null_results.tsv",
        "fdr_summary": output_dir / "pseudobulk_de_fdr_summary.tsv",
    }
    result.to_csv(paths["de_results"], sep="\t", index=False)
    pd.DataFrame(rows).to_csv(paths["fdr_summary"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/analysis_pseudobulk_de")
    )
    parser.add_argument("--n-genes", type=int, default=500)
    parser.add_argument("--n-samples", type=int, default=6)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    paths = run(args.output_dir, n_genes=args.n_genes, n_samples=args.n_samples, seed=args.seed)
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
