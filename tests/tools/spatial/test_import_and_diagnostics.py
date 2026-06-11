"""Tests for spatial diagnostics and import without squidpy."""

import numpy as np
import pytest
from anndata import AnnData


def test_import_without_squidpy():
    """scLucid.tools.spatial should import even when squidpy is absent."""
    import importlib

    import scLucid.tools.spatial as spmod

    # Re-import to ensure it is fresh
    importlib.reload(spmod)
    assert hasattr(spmod, "run_spatial_analysis")
    assert hasattr(spmod, "build_spatial_neighbors")
    assert hasattr(spmod, "compute_moran_i")


def test_diagnose_spatial_data_quality_passes():
    from scLucid.tools.spatial import diagnose_spatial_data_quality

    n_spots = 20
    coords = np.column_stack([np.arange(n_spots), np.arange(n_spots)])
    adata = AnnData(X=np.random.poisson(5, size=(n_spots, 10)))
    adata.obsm["spatial"] = coords
    result = diagnose_spatial_data_quality(adata)
    assert result["passed"]
    assert result["n_spots"] == n_spots
    assert result["n_duplicate_coords"] == 0


def test_diagnose_spatial_data_quality_missing_key():
    from scLucid.tools.spatial import diagnose_spatial_data_quality

    adata = AnnData(X=np.random.poisson(5, size=(5, 10)))
    result = diagnose_spatial_data_quality(adata)
    assert not result["passed"]
    assert any("not found" in w for w in result["warnings"])
