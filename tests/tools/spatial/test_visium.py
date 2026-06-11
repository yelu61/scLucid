"""Tests for Visium I/O helpers."""

import numpy as np
from anndata import AnnData

from scLucid.tools.spatial import crop_visium, rotate_visium, VisiumIOConfig


def test_crop_visium():
    n_spots = 20
    coords = np.column_stack([np.arange(n_spots), np.zeros(n_spots)])
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 5)))
    adata.obsm["spatial"] = coords
    cropped = crop_visium(adata, xlim=(5, 10))
    assert cropped.n_obs == 6
    assert "spatial" in cropped.obsm


def test_rotate_visium():
    n_spots = 10
    coords = np.column_stack([np.arange(n_spots), np.zeros(n_spots)])
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 5)))
    adata.obsm["spatial"] = coords.copy()
    rotated = rotate_visium(adata, angle_degrees=90)
    assert rotated.n_obs == n_spots
    assert "spatial" in rotated.obsm
    # 90-degree rotation around centroid: verify distances preserved
    orig = coords.copy()
    rot = rotated.obsm["spatial"]
    orig_dists = np.linalg.norm(orig - orig.mean(axis=0), axis=1)
    rot_dists = np.linalg.norm(rot - rot.mean(axis=0), axis=1)
    np.testing.assert_allclose(rot_dists, orig_dists, atol=1e-6)
