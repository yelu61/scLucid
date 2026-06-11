"""Tests for spatial neighbor graph construction."""

import numpy as np
from anndata import AnnData

from scLucid.tools.spatial import build_spatial_neighbors, SpatialNeighborsConfig


def test_build_spatial_neighbors_knn():
    n_spots = 10
    coords = np.random.rand(n_spots, 2)
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 5)))
    adata.obsm["spatial"] = coords
    out = build_spatial_neighbors(adata, config=SpatialNeighborsConfig(n_neigh=3))
    assert "spatial_connectivities" in out.obsp
    assert "spatial_distances" in out.obsp
    assert out.obsp["spatial_connectivities"].shape == (n_spots, n_spots)


def test_build_spatial_neighbors_radius():
    n_spots = 10
    coords = np.random.rand(n_spots, 2)
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 5)))
    adata.obsm["spatial"] = coords
    out = build_spatial_neighbors(
        adata, config=SpatialNeighborsConfig(method="radius", radius=0.5)
    )
    assert "spatial_connectivities" in out.obsp
    assert out.obsp["spatial_connectivities"].shape == (n_spots, n_spots)
