"""Tests for tumor biological fidelity benchmark runner."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from validation.qc.run_tumor_biological_fidelity_benchmark import (
    _adaptive_thresholds,
    _malignant_like_retention,
    _retention_bias_rows,
)


class TestTumorBiologicalFidelityBenchmark:
    def test_adaptive_thresholds_uses_production_recommender(self):
        rng = np.random.default_rng(0)
        X = rng.poisson(2, (300, 50)).astype(float)
        adata = AnnData(X)
        adata.obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
        adata.obs["total_counts"] = X.sum(axis=1)
        adata.obs["pct_counts_mt"] = rng.uniform(0, 30, size=300)
        adata.obs["sampleID"] = "s1"
        adata.layers["counts"] = X.copy()

        thresholds = _adaptive_thresholds(adata, tumor_aware=True)
        assert "min_genes" in thresholds
        assert "max_mt_percent" in thresholds
        assert "recommendation_method" in thresholds
        assert thresholds["overall_strategy"] == "tumor_aware"

    def test_retention_bias_rows_decompose_review_reason(self):
        rng = np.random.default_rng(0)
        X = rng.poisson(2, (100, 20)).astype(float)
        adata = AnnData(X)
        adata.obs["sample"] = ["A"] * 60 + ["B"] * 40
        adata.obs["cell_type"] = ["T"] * 100
        keep = pd.Series([True] * 20 + [False] * 40 + [True] * 40, index=adata.obs_names)

        rows = _retention_bias_rows(adata, "demo", "sclucid_tumor_aware", keep)
        by_group = {r["group"]: r for r in rows}
        assert by_group["A"]["review_required"] is True
        assert by_group["A"]["review_reason"] == "strategy_bias_low_retention"
        assert by_group["B"]["review_required"] is False

    def test_retention_bias_rows_skips_small_groups(self):
        rng = np.random.default_rng(0)
        X = rng.poisson(2, (30, 20)).astype(float)
        adata = AnnData(X)
        adata.obs["sample"] = ["A"] * 28 + ["B"] * 2
        keep = pd.Series([True] * 30, index=adata.obs_names)

        rows = _retention_bias_rows(adata, "demo", "sclucid_tumor_aware", keep)
        by_group = {r["group"]: r for r in rows}
        assert by_group["B"]["review_required"] is False
        assert by_group["B"]["review_reason"] == "small_group_skipped"

    def test_malignant_like_retention_computes_signature_cells(self):
        rng = np.random.default_rng(0)
        X = rng.poisson(2, (100, 20)).astype(float)
        adata = AnnData(X)
        adata.var_names = [f"gene_{i}" for i in range(20)]
        # Make first 6 genes the malignant-like panel, and upregulate in last 30 cells.
        X[:70, :6] = rng.poisson(1, (70, 6)).astype(float)
        X[70:, :6] = rng.poisson(10, (30, 6)).astype(float)
        adata.X = X
        adata.layers["counts"] = X.copy()
        genes = tuple(adata.var_names[:6])
        keep = pd.Series([True] * 80 + [False] * 20, index=adata.obs_names)

        result = _malignant_like_retention(adata, keep, genes)
        assert result["malignant_like_cells"] > 0
        assert 0.0 <= result["malignant_like_retention_rate"] <= 1.0
        assert result["malignant_like_genes_present"] == 6
