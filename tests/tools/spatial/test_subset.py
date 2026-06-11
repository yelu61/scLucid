"""Tests for spatial subsetting."""

import numpy as np
from anndata import AnnData

from scLucid.tools.spatial import subset_spatial_window, SpatialWindowConfig


def test_subset_spatial_window():
    n_spots = 20
    coords = np.column_stack([np.arange(n_spots), np.zeros(n_spots)])
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 5)))
    adata.obsm["spatial"] = coords
    config = SpatialWindowConfig(xlim=(5, 10), ylim=(-1, 1))
    subset = subset_spatial_window(adata, config=config)
    assert subset.n_obs == 6  # spots 5,6,7,8,9,10 -> 6 spots
    assert "spatial" in subset.obsm


def test_subset_spatial_window_empty():
    n_spots = 5
    coords = np.array([[0, 0], [1, 1], [2, 2], [3, 3], [4, 4]])
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 3)))
    adata.obsm["spatial"] = coords
    config = SpatialWindowConfig(xlim=(10, 20))
    subset = subset_spatial_window(adata, config=config)
    assert subset.n_obs == 0
