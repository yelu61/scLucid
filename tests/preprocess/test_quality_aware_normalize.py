"""Tests for scLucid.preprocess.quality_aware_normalize."""

import numpy as np
import pytest
import scipy.sparse
from anndata import AnnData


def _make_adata_with_qc(n_cells=100, n_genes=50):
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

    # --- Edge case tests ---

    def test_stores_side_effect_columns(self):
        """Verify quality_score, quality_bin, and quality_weight are created."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=100)
        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt", "doublet_score"],
            input_layer="counts",
            n_bins=5,
            target_sum=1e4,
        )
        assert "quality_score" in result.obs
        assert "quality_bin" in result.obs
        assert "quality_weight" in result.obs
        assert result.obs["quality_score"].between(0, 1).all()
        assert result.obs["quality_bin"].nunique() >= 2  # at least 2 bins were created

    def test_quality_score_direction_heuristic(self):
        """pct_ / mt_ metrics treated as lower-is-better; n_genes as higher-is-better."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=100)
        # Create cells with clearly distinct QC profiles
        adata.obs["pct_counts_mt"] = 0.0
        adata.obs["n_genes_by_counts"] = 500.0
        # Cell 0: high MT (bad) + low n_genes (bad) → should get low quality score
        adata.obs.loc[adata.obs_names[0], "pct_counts_mt"] = 95.0
        adata.obs.loc[adata.obs_names[0], "n_genes_by_counts"] = 5.0
        # Cell 1: low MT (good) + high n_genes (good) → should get high quality score
        adata.obs.loc[adata.obs_names[1], "pct_counts_mt"] = 1.0
        adata.obs.loc[adata.obs_names[1], "n_genes_by_counts"] = 900.0
        # Remaining cells: moderate values
        adata.obs.loc[adata.obs_names[2:], "pct_counts_mt"] = 20.0
        adata.obs.loc[adata.obs_names[2:], "n_genes_by_counts"] = 400.0

        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt", "n_genes_by_counts"],
            input_layer="counts",
            n_bins=3,
            target_sum=1e4,
        )
        scores = result.obs["quality_score"]
        # The bad cell (0) should have lower score than the good cell (1)
        assert scores.iloc[0] < scores.iloc[1]
        # And should also be below the average
        assert scores.iloc[0] < scores.iloc[2:].mean()

    def test_zero_count_cell_is_handled(self):
        """A cell with total count = 0 should not produce inf/nan."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=20)
        adata.layers["counts"][0, :] = 0
        adata.obs.loc[adata.obs_names[0], "n_genes_by_counts"] = 0

        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt"],
            input_layer="counts",
            n_bins=2,
            target_sum=1e4,
        )
        normalized = result.layers["quality_normalized"]
        assert not np.any(np.isinf(normalized))
        assert not np.any(np.isnan(normalized))

    def test_all_identical_quality_metric(self):
        """A metric where all values are identical should not crash pd.qcut."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=50)
        adata.obs["uniform_metric"] = 5.0  # all identical

        result = quality_aware_normalize(
            adata,
            quality_metrics=["uniform_metric"],
            input_layer="counts",
            n_bins=5,
        )
        assert "quality_normalized" in result.layers
        # All quality scores should be identical (1e-10 guard span → single bin)
        scores = result.obs["quality_score"]
        np.testing.assert_allclose(scores, scores.iloc[0], rtol=1e-5)
        assert result.obs["quality_bin"].nunique() == 1

    def test_nan_in_quality_metric_not_poisoning_all_cells(self):
        """A single NaN metric value should not produce NaN quality scores for all cells."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=30)
        adata.obs.loc[adata.obs_names[5], "pct_counts_mt"] = np.nan

        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt"],
            input_layer="counts",
        )
        scores = result.obs["quality_score"]
        # The NaN cell should get a mid-score, not poison others
        assert not np.isnan(scores.iloc[[i for i in range(30) if i != 5]]).any()
        assert not np.isnan(scores.iloc[5])  # NaN cell gets 0.5
        assert 0.4 <= scores.iloc[5] <= 0.6

    def test_sparse_matrix_input(self):
        """Sparse counts layer should work correctly."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=50)
        sparse_X = scipy.sparse.csr_matrix(adata.layers["counts"])
        sparse_adata = AnnData(
            X=sparse_X.copy(),
            obs=adata.obs.copy(),
            var=adata.var.copy(),
        )
        sparse_adata.layers["counts"] = sparse_X.copy()

        result = quality_aware_normalize(
            sparse_adata,
            quality_metrics=["pct_counts_mt"],
            input_layer="counts",
        )
        assert "quality_normalized" in result.layers
        # Output should be sparse since input was sparse
        assert scipy.sparse.issparse(result.layers["quality_normalized"])

    def test_input_layer_fallback_to_X(self):
        """When input_layer not in layers, falls back to adata.X."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=30)
        # No "counts" layer — should use adata.X
        del adata.layers["counts"]

        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt"],
        )
        assert "quality_normalized" in result.layers

    def test_fewer_cells_than_bins(self):
        """When cells < bins, pd.qcut with duplicates='drop' handles it."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=5)
        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt"],
            input_layer="counts",
            n_bins=10,  # more bins than cells
        )
        assert "quality_normalized" in result.layers
        # Verify fewer bins than requested were created
        assert result.obs["quality_bin"].nunique() <= 5  # can't have more bins than cells

    def test_target_sum_auto_uses_per_bin_median(self):
        """When target_sum is None, each bin normalizes to its own median."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=60)
        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt"],
            input_layer="counts",
            target_sum=None,
        )
        assert "quality_normalized" in result.layers
        # Output exists and is well-formed
        assert result.layers["quality_normalized"].shape == (60, 50)

    def test_single_cell_does_not_crash(self):
        """A single-cell dataset should not crash."""
        from scLucid.preprocess import quality_aware_normalize

        adata = _make_adata_with_qc(n_cells=1)
        result = quality_aware_normalize(
            adata,
            quality_metrics=["pct_counts_mt"],
            input_layer="counts",
            n_bins=1,
        )
        assert "quality_normalized" in result.layers
