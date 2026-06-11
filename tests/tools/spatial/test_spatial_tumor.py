"""Tests for tumor-focused spatial utilities."""

import numpy as np
from anndata import AnnData

from scLucid.tools.spatial import (
    analyze_spatial_niches,
    compute_immune_infiltration_score,
    find_tumor_stroma_boundary,
    spatial_ici_response_signature,
)


def test_find_tumor_stroma_boundary():
    n_spots = 40
    coords = np.column_stack([np.tile(np.arange(8), 5), np.repeat(np.arange(5), 8)])
    labels = ["tumor"] * 16 + ["stroma"] * 24
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 10)))
    adata.obsm["spatial"] = coords
    adata.obs["tumor_label"] = labels
    out = find_tumor_stroma_boundary(
        adata, tumor_label="tumor", label_key="tumor_label", key_added="boundary"
    )
    assert "boundary" in out.obs.columns
    assert set(out.obs["boundary"].unique()).issubset({"core", "boundary", "margin", "stroma"})
    assert "tumor_boundary" in out.uns["sclucid"]["tools"]["spatial"]


def test_compute_immune_infiltration_score():
    n_spots = 20
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 10)))
    adata.obsm["spatial"] = np.random.rand(n_spots, 2)
    out = compute_immune_infiltration_score(adata, immune_markers=None)
    assert "immune_infiltration_score" in out.obs.columns
    assert "immune_infiltration_score" in out.uns["sclucid"]["tools"]["spatial"]


def test_analyze_spatial_niches():
    n_spots = 30
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 10)))
    adata.obsm["spatial"] = np.random.rand(n_spots, 2)
    out = analyze_spatial_niches(adata, n_components=2)
    assert "spatial_niches" in out.obs.columns
    assert out.obs["spatial_niches"].nunique() <= 2
    assert "spatial_niches" in out.uns["sclucid"]["tools"]["spatial"]


def test_spatial_ici_response_signature():
    n_spots = 20
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 10)))
    result = spatial_ici_response_signature(adata)
    assert not result.empty
    assert any("ICI_" in c for c in result.columns)
    assert not result["valid_for_publication_inference"].any()
