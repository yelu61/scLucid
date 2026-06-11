"""Tests for tissue zone detection."""

import numpy as np
from anndata import AnnData

from scLucid.tools.spatial import find_tissue_zones, TissueZonesConfig


def test_find_tissue_zones():
    n_spots = 30
    coords = np.column_stack([np.arange(n_spots), np.zeros(n_spots)])
    # Two spatially separated expression patterns
    expr = np.zeros((n_spots, 10))
    expr[:15, :5] = np.random.poisson(10, size=(15, 5))
    expr[15:, 5:] = np.random.poisson(10, size=(15, 5))
    adata = AnnData(X=expr)
    adata.obsm["spatial"] = coords
    find_tissue_zones(adata, config=TissueZonesConfig(n_components=2, input="expression"))
    assert "tissue_zones" in adata.obs.columns
    assert len(adata.obs["tissue_zones"].unique()) <= 2
    assert "X_tissue_zones" in adata.obsm


def test_find_tissue_zones_with_n_components():
    n_spots = 20
    coords = np.random.rand(n_spots, 2)
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 5)))
    adata.obsm["spatial"] = coords
    find_tissue_zones(adata, config=TissueZonesConfig(n_components=3, input="expression"))
    assert "tissue_zones" in adata.obs.columns
    assert len(adata.obs["tissue_zones"].unique()) <= 3
