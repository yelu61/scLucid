#!/usr/bin/env python3
"""Benchmark MT% threshold strategies against tumor program preservation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scLucid.qc.policy.intelligent_qc import (
    IntelligentQCConfig,
    IntelligentQCRecommender,
    StrategyType,
)
from validation.dataset_registry import DATASETS
from validation.gene_panels import MARKER_PANELS, TUMOR_PROGRAM_PANELS, present_genes

STRATEGIES = (
    "fixed_20",
    "fixed_5",
    "bimodal_gmm",
    "sample_aware",
    "multicomponent",
)

MT_STRATEGY_CONFIG = {
    "fixed_20": ("percentile", {}),
    "fixed_5": ("percentile", {}),
    "bimodal_gmm": ("bimodal_gmm", {}),
    "sample_aware": ("sample_aware", {}),
    "multicomponent": ("multicomponent", {"mt_max_components": 4}),
}


def _matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _subset(adata: ad.AnnData, max_cells: int | None, seed: int) -> ad.AnnData:
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    return adata[keep].copy()


def _ensure_qc_metrics(adata: ad.AnnData) -> None:
    X = _matrix(adata)
    if "total_counts" not in adata.obs:
        adata.obs["total_counts"] = np.asarray(X.sum(axis=1)).ravel()
    if "n_genes_by_counts" not in adata.obs:
        if sp.issparse(X):
            adata.obs["n_genes_by_counts"] = np.asarray((X > 0).sum(axis=1)).ravel()
        else:
            adata.obs["n_genes_by_counts"] = (np.asarray(X) > 0).sum(axis=1)
    if "pct_counts_mt" not in adata.obs:
        names = pd.Index(adata.var_names.astype(str))
        mt_mask = names.str.startswith("MT-") | names.str.startswith("mt-")
        if mt_mask.any():
            mt_counts = np.asarray(X[:, mt_mask].sum(axis=1)).ravel()
            total = np.asarray(adata.obs["total_counts"], dtype=float)
            adata.obs["pct_counts_mt"] = np.divide(
                mt_counts * 100.0,
                total,
                out=np.zeros_like(total, dtype=float),
                where=total > 0,
            )
        else:
            adata.obs["pct_counts_mt"] = np.nan


def _program_retention(
    adata: ad.AnnData, keep: pd.Series, program_genes: tuple[str, ...]
) -> float:
    present = present_genes(adata.var_names, program_genes)
    if len(present) < 2:
        return float("nan")
    X = _matrix(adata)
    idx = pd.Index(adata.var_names.astype(str)).get_indexer(present)
    before = float(np.asarray(X[:, idx].mean()).mean())
    after = float(np.asarray(X[keep.to_numpy(), :][:, idx].mean()).mean()) if keep.any() else 0.0
    if before <= 0:
        return float("nan")
    return after / before


def _run_strategy(
    adata: ad.AnnData, dataset: str, strategy: str, sample_key: str | None
) -> dict[str, Any]:
    mt_pct = np.asarray(adata.obs["pct_counts_mt"].values, dtype=float)

    if strategy == "fixed_20":
        threshold = 20.0
    elif strategy == "fixed_5":
        threshold = 5.0
    else:
        model, kwargs = MT_STRATEGY_CONFIG[strategy]
        cfg_kwargs = {"mt_model": model, **kwargs}
        if model == "sample_aware" and sample_key is not None:
            cfg_kwargs["sample_key"] = sample_key
        cfg = IntelligentQCConfig(**cfg_kwargs)
        recommender = IntelligentQCRecommender(
            strategy=StrategyType.TUMOR_AWARE, config=cfg
        )
        rec = recommender._recommend_max_mt(
            adata, tissue_type="tumor", strategy=StrategyType.TUMOR_AWARE, plot=False
        )
        threshold = float(rec.threshold)

    keep = pd.Series(mt_pct <= threshold, index=adata.obs_names)
    removed = (~keep).sum()

    program_scores = {
        name: _program_retention(adata, keep, genes)
        for name, genes in TUMOR_PROGRAM_PANELS.items()
    }
    marker_scores = {
        name: _program_retention(adata, keep, genes)
        for name, genes in MARKER_PANELS.items()
    }

    return {
        "dataset": dataset,
        "strategy": strategy,
        "threshold": round(threshold, 3),
        "retention_rate": float(keep.mean()),
        "removed_cells": int(removed),
        "mean_program_retention": float(
            np.nanmean([v for v in program_scores.values() if np.isfinite(v)])
        ),
        "min_program_retention": float(
            np.nanmin([v for v in program_scores.values() if np.isfinite(v)] or [np.nan])
        ),
        "mean_marker_retention": float(
            np.nanmean([v for v in marker_scores.values() if np.isfinite(v)])
        ),
        "program_scores": program_scores,
        "marker_scores": marker_scores,
    }


def run(
    output_dir: Path, datasets: set[str] | None, max_cells: int | None, seed: int
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []

    for spec in DATASETS:
        if datasets and spec.key not in datasets:
            continue
        if not spec.path.exists():
            continue
        is_tumor = "tumor" in spec.modality_role or any("tumor" in role for role in spec.qc_roles)
        if not is_tumor:
            continue
        adata = ad.read_h5ad(spec.path)
        adata = _subset(adata, max_cells=max_cells, seed=seed)
        _ensure_qc_metrics(adata)
        sample_key = (
            "sample" if "sample" in adata.obs.columns
            else "sampleID" if "sampleID" in adata.obs.columns
            else None
        )
        for strategy in STRATEGIES:
            rows.append(_run_strategy(adata, spec.key, strategy, sample_key))

    scorecard_rows = []
    for row in rows:
        scorecard_rows.append(
            {
                "dataset": row["dataset"],
                "strategy": row["strategy"],
                "threshold": row["threshold"],
                "retention_rate": row["retention_rate"],
                "removed_cells": row["removed_cells"],
                "mean_program_retention": row["mean_program_retention"],
                "min_program_retention": row["min_program_retention"],
                "mean_marker_retention": row["mean_marker_retention"],
                "review_required": row["mean_program_retention"] < 0.8
                or row["retention_rate"] < 0.35,
            }
        )

    paths = {
        "scorecard": output_dir / "mt_threshold_scorecard.tsv",
        "details": output_dir / "mt_threshold_details.json",
    }
    pd.DataFrame(scorecard_rows).to_csv(paths["scorecard"], sep="\t", index=False)
    import json

    paths["details"].write_text(json.dumps(rows, indent=2, default=str), encoding="utf-8")
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("validation_outputs/qc_mt_threshold_benchmark")
    )
    parser.add_argument("--datasets", nargs="*", default=None)
    parser.add_argument(
        "--max-cells",
        type=int,
        default=5000,
        help="Deterministic subset size per dataset for pilot runs. Use 0 for full data.",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    paths = run(
        args.output_dir,
        datasets=set(args.datasets) if args.datasets else None,
        max_cells=None if args.max_cells == 0 else args.max_cells,
        seed=args.seed,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
