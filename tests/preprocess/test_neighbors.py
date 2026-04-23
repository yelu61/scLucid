"""Unit tests for scLucid.preprocess.neighbors."""

import numpy as np
import pytest

from scLucid.preprocess.config import NeighborsConfig
from scLucid.preprocess.neighbors import optimize_neighbors_pcs


@pytest.mark.unit
class TestOptimizeNeighborsPcs:
    """Tests for optimize_neighbors_pcs function."""

    def test_basic_run(self, minimal_adata):
        adata = minimal_adata.copy()
        # Create a simple PCA embedding
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 20))

        config = NeighborsConfig(
            n_neighbors_list=[5, 10],
            n_pcs_list=[5, 10],
            resolution=0.5,
            plot=False,
            report=False,
            verbose=False,
        )

        result_df = optimize_neighbors_pcs(adata, config=config)

        assert result_df is not None
        assert len(result_df) == 4  # 2 n_neighbors x 2 n_pcs
        assert "n_neighbors" in result_df.columns
        assert "n_pcs" in result_df.columns
        assert "silhouette_score" in result_df.columns

    def test_stores_best_params_in_uns(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 20))

        config = NeighborsConfig(
            n_neighbors_list=[5, 10],
            n_pcs_list=[5, 10],
            resolution=0.5,
            plot=False,
            report=False,
            verbose=False,
        )

        optimize_neighbors_pcs(adata, config=config)

        assert "sclucid" in adata.uns
        assert "preprocess" in adata.uns["sclucid"]
        meta = adata.uns["sclucid"]["preprocess"].get("neighbors_optimization", {})
        assert "best_params" in meta
        assert meta["best_params"] is not None
        assert "n_neighbors" in meta["best_params"]
        assert "n_pcs" in meta["best_params"]

    def test_missing_pca_raises(self, minimal_adata):
        adata = minimal_adata.copy()
        # No PCA embedding

        config = NeighborsConfig(
            n_neighbors_list=[5, 10],
            n_pcs_list=[5, 10],
            plot=False,
            report=False,
            verbose=False,
        )

        with pytest.raises(ValueError, match="not found"):
            optimize_neighbors_pcs(adata, config=config)

    def test_single_grid_point(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))

        config = NeighborsConfig(
            n_neighbors_list=[5],
            n_pcs_list=[5],
            resolution=0.5,
            plot=False,
            report=False,
            verbose=False,
        )

        result_df = optimize_neighbors_pcs(adata, config=config)

        assert len(result_df) == 1
        assert result_df.iloc[0]["n_neighbors"] == 5
        assert result_df.iloc[0]["n_pcs"] == 5

    def test_config_not_mutated(self, minimal_adata):
        adata = minimal_adata.copy()
        adata.obsm["X_pca"] = np.random.default_rng(0).normal(size=(adata.n_obs, 10))

        config = NeighborsConfig(
            n_neighbors_list=[5],
            n_pcs_list=[5],
            resolution=0.5,
            plot=False,
            report=False,
            verbose=False,
        )
        original_dict = config.to_dict()
        optimize_neighbors_pcs(adata, config=config)
        assert config.to_dict() == original_dict
