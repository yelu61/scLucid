"""Tests for spatial autocorrelation and SVG detection."""

import numpy as np
from anndata import AnnData

from scLucid.tools.spatial import (
    SpatialNeighborsConfig,
    build_spatial_neighbors,
    compute_moran_i,
    find_spatially_variable_genes,
)


def test_compute_moran_i_positive_autocorr():
    n_spots = 20
    coords = np.column_stack([np.arange(n_spots), np.zeros(n_spots)])
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 5)))
    adata.obsm["spatial"] = coords
    build_spatial_neighbors(adata, config=SpatialNeighborsConfig(n_neigh=2))
    # Values increasing along the line -> positive spatial autocorrelation
    values = np.arange(n_spots).astype(float)
    result = compute_moran_i(adata, values)
    assert result["moran_i"] > 0
    assert result["mode"] == "moran"


def test_find_spatially_variable_genes():
    n_spots = 20
    coords = np.column_stack([np.arange(n_spots), np.zeros(n_spots)])
    # Gene 0 has strong spatial pattern; gene 1 is random
    expr = np.zeros((n_spots, 2))
    expr[:, 0] = np.arange(n_spots)
    expr[:, 1] = np.random.poisson(5, size=n_spots)
    adata = AnnData(X=expr, var={"gene": ["gene0", "gene1"]})
    adata.var_names = ["gene0", "gene1"]
    adata.obsm["spatial"] = coords
    result = find_spatially_variable_genes(adata)
    assert not result.empty
    assert "moran_i" in result.columns
    assert "pvals_adj" in result.columns
    assert "spatially_variable" in result.columns
    # Gene 0 should have higher Moran's I than gene 1
    moran_gene0 = result.loc[result["gene"] == "gene0", "moran_i"].iloc[0]
    moran_gene1 = result.loc[result["gene"] == "gene1", "moran_i"].iloc[0]
    assert moran_gene0 > moran_gene1


def test_find_spatially_variable_genes_handles_constant_genes():
    n_spots = 20
    coords = np.column_stack([np.arange(n_spots), np.zeros(n_spots)])
    expr = np.zeros((n_spots, 3))
    expr[:, 0] = np.arange(n_spots)
    expr[:, 1] = 5.0
    expr[:, 2] = np.random.poisson(5, size=n_spots)
    adata = AnnData(X=expr)
    adata.var_names = ["patterned", "constant", "random"]
    adata.obsm["spatial"] = coords

    result = find_spatially_variable_genes(adata)

    assert "constant" not in set(result["gene"])
    assert "spatially_variable" in adata.var
    assert bool(adata.var.loc["constant", "spatially_variable"]) is False
