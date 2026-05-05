"""Tests for scLucid.preprocess.plot_normalization_effect."""

import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest
from anndata import AnnData


def _make_adata_with_layers(n_cells: int = 50, n_genes: int = 20):
    """Create an AnnData with 'counts' and 'normalized' layers."""
    rng = np.random.default_rng(42)
    X = rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    adata = AnnData(X)
    adata.obs_names = [f"cell_{i:04d}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i:03d}" for i in range(n_genes)]
    adata.layers["counts"] = X.copy()
    # Simple normalization for testing
    adata.layers["normalized"] = np.log1p(X / X.sum(axis=1, keepdims=True) * 1e4)
    return adata


class TestPlotNormalizationEffect:
    def test_plot_normalization_effect_basic(self):
        """Smoke test: returns a matplotlib Figure."""
        from scLucid.preprocess import plot_normalization_effect

        adata = _make_adata_with_layers(n_cells=50)
        fig = plot_normalization_effect(
            adata,
            original_layer="counts",
            normalized_layer="normalized",
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_normalization_effect_with_save_dir(self, tmp_path):
        """Saving to directory produces a file."""
        from scLucid.preprocess import plot_normalization_effect

        adata = _make_adata_with_layers(n_cells=50)
        save_dir = tmp_path / "norm_plots"
        fig = plot_normalization_effect(
            adata,
            original_layer="counts",
            normalized_layer="normalized",
            save_dir=str(save_dir),
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
        # Directory should be created
        assert save_dir.exists()

    def test_plot_normalization_effect_missing_layer(self):
        """Missing layer raises ValueError."""
        from scLucid.preprocess import plot_normalization_effect

        adata = _make_adata_with_layers()
        with pytest.raises(ValueError, match="not found"):
            plot_normalization_effect(
                adata,
                original_layer="nonexistent",
                normalized_layer="normalized",
            )

    def test_plot_normalization_effect_gene_subset(self):
        """gene_subset limits the plotted genes."""
        from scLucid.preprocess import plot_normalization_effect

        adata = _make_adata_with_layers(n_cells=50)
        fig = plot_normalization_effect(
            adata,
            original_layer="counts",
            normalized_layer="normalized",
            gene_subset=["gene_000", "gene_001"],
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_normalization_effect_no_log(self):
        """log_transformed=False adjusts plot titles/scales."""
        from scLucid.preprocess import plot_normalization_effect

        adata = _make_adata_with_layers(n_cells=50)
        fig = plot_normalization_effect(
            adata,
            original_layer="counts",
            normalized_layer="normalized",
            log_transformed=False,
        )
        assert isinstance(fig, plt.Figure)
        plt.close(fig)
