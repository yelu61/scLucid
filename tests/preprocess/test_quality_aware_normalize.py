"""Tests for scLucid.preprocess.quality_aware_normalize."""

import numpy as np
import pytest
from anndata import AnnData


def _make_adata_with_qc(n_cells: int = 100, n_genes: int = 50):
    """Create an AnnData with QC metrics in .obs."""
    rng = np.random.default_rng(42)
    X = rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(X)
    adata.obs_names = [f"cell_{i:04d}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i:03d}" for i in range(n_genes)]
    adata.layers["counts"] = X.copy()
    adata.obs["pct_counts_mt"] = rng.uniform(0, 30, size=n_cells)
    adata.obs["doublet_score"] = rng.uniform(0, 1, size=n_cells)
    adata.obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
    return adata


class TestQualityAwareNormalize:
    def test_quality_aware_normalize_basic(self):
        """Smoke test: creates a normalized layer."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=100)
        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt", "doublet_score"],
            input_layer="counts",
            output_layer="quality_normalized",
            target_sum=1e4,
            log_transform=True,
        )
        assert "quality_normalized" in result.layers
        assert result.layers["quality_normalized"].shape == (100, 50)

    def test_quality_aware_normalize_output_is_non_negative(self):
        """Normalized output should be non-negative (log1p)."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=100)
        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt"],
            input_layer="counts",
            output_layer="quality_normalized",
            target_sum=1e4,
            log_transform=True,
        )
        assert (result.layers["quality_normalized"] >= 0).all()

    def test_quality_aware_normalize_no_log(self):
        """log_transform=False skips log1p."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=100)
        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt"],
            input_layer="counts",
            output_layer="quality_normalized",
            log_transform=False,
        )
        assert "quality_normalized" in result.layers

    def test_quality_aware_normalize_missing_metric(self):
        """Missing quality metric raises ValueError."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=100)
        with pytest.raises(ValueError, match="Missing quality metrics"):
            quality_aware_normalize(
                adata,
                quality_metrics=["nonexistent_metric"],
                input_layer="counts",
            )

    def test_quality_aware_normalize_different_bins(self):
        """Different n_bins values should all work."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=100)
        for n_bins in [2, 5, 10]:
            result = quality_aware_normalize(
                adata,
                quality_metrics=["pct_counts_mt"],
                input_layer="counts",
                n_bins=n_bins,
            )
            assert "quality_normalized" in result.layers
