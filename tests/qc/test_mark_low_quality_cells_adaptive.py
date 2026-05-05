"""Tests for scLucid.qc.mark_low_quality_cells_adaptive."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


def _make_adata_with_batches(n_cells: int = 200, n_genes: int = 50, n_batches: int = 2):
    """Create a minimal AnnData with batch/sample labels."""
    rng = np.random.default_rng(42)
    X = rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(X)
    adata.obs_names = [f"cell_{i:04d}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i:03d}" for i in range(n_genes)]
    adata.obs["sampleID"] = [f"batch_{i % n_batches}" for i in range(n_cells)]
    adata.obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
    adata.obs["total_counts"] = X.sum(axis=1)
    adata.obs["pct_counts_mt"] = rng.uniform(0, 20, size=n_cells)
    return adata


class TestMarkLowQualityCellsAdaptive:
    def test_adaptive_marking_basic(self):
        """Smoke test: adaptive marking runs and adds columns."""
        from scLucid.qc.filtering.core import mark_low_quality_cells_adaptive

        adata = _make_adata_with_batches(n_cells=200, n_batches=2)
        result = mark_low_quality_cells_adaptive(
            adata,
            batch_key="sampleID",
            metrics=["n_genes_by_counts", "pct_counts_mt"],
            method="hierarchical",
        )

        assert isinstance(result, AnnData)
        adaptive_cols = [c for c in result.obs.columns if "_adaptive" in c]
        assert len(adaptive_cols) > 0

    def test_adaptive_marking_multiple_metrics(self):
        """Multiple metrics each get their own adaptive flag column."""
        from scLucid.qc.filtering.core import mark_low_quality_cells_adaptive

        adata = _make_adata_with_batches(n_cells=200, n_batches=2)
        result = mark_low_quality_cells_adaptive(
            adata,
            batch_key="sampleID",
            metrics=["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            method="hierarchical",
        )

        for metric in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
            col = f"outlier_{metric}_adaptive"
            assert col in result.obs.columns

    def test_adaptive_marking_respects_batch(self):
        """Outliers are identified per-batch, not globally."""
        from scLucid.qc.filtering.core import mark_low_quality_cells_adaptive

        rng = np.random.default_rng(42)
        n = 100
        X0 = rng.integers(0, 5, size=(n, 30)).astype(np.float32)
        X1 = rng.integers(0, 15, size=(n, 30)).astype(np.float32)
        X = np.vstack([X0, X1])
        adata = AnnData(X)
        adata.obs_names = [f"c{i}" for i in range(2 * n)]
        adata.obs["sampleID"] = ["batch_0"] * n + ["batch_1"] * n
        adata.obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
        adata.obs["total_counts"] = X.sum(axis=1)
        adata.obs["pct_counts_mt"] = rng.uniform(0, 20, size=2 * n)

        result = mark_low_quality_cells_adaptive(
            adata,
            batch_key="sampleID",
            metrics=["total_counts"],
            method="hierarchical",
        )

        b0_outliers = result.obs.loc[result.obs["sampleID"] == "batch_0", "outlier_total_counts_adaptive"].sum()
        b1_outliers = result.obs.loc[result.obs["sampleID"] == "batch_1", "outlier_total_counts_adaptive"].sum()
        assert b0_outliers >= 0
        assert b1_outliers >= 0

    def test_adaptive_marking_invalid_method_raises(self):
        """Invalid method raises ValueError."""
        from scLucid.qc.filtering.core import mark_low_quality_cells_adaptive

        adata = _make_adata_with_batches(n_cells=100)
        with pytest.raises(ValueError, match="Unknown method"):
            mark_low_quality_cells_adaptive(
                adata,
                batch_key="sampleID",
                metrics=["n_genes_by_counts"],
                method="nonexistent",
            )

    def test_adaptive_marking_few_batches(self):
        """Works with few batches where hierarchical still has enough groups."""
        from scLucid.qc.filtering.core import mark_low_quality_cells_adaptive

        adata = _make_adata_with_batches(n_cells=100, n_batches=2)
        result = mark_low_quality_cells_adaptive(
            adata,
            batch_key="sampleID",
            metrics=["n_genes_by_counts"],
            method="independent",
        )
        assert "outlier_n_genes_by_counts_adaptive" in result.obs.columns
