"""Spatial subsetting by coordinate window."""

from __future__ import annotations

from typing import Optional

import numpy as np
from anndata import AnnData

from .config import SpatialWindowConfig


def subset_spatial_window(
    adata: AnnData,
    config: Optional[SpatialWindowConfig] = None,
) -> AnnData:
    """Subset an AnnData to spots within a spatial coordinate window.

    Parameters
    ----------
    adata
        AnnData with spatial coordinates.
    config
        Window configuration.

    Returns
    -------
    AnnData
        Subsetted AnnData.
    """
    if config is None:
        config = SpatialWindowConfig()

    if config.spatial_key not in adata.obsm:
        raise KeyError(f"Spatial key '{config.spatial_key}' not found in adata.obsm")

    coords = np.asarray(adata.obsm[config.spatial_key])[:, :2]
    mask = np.ones(adata.n_obs, dtype=bool)

    if config.xlim is not None:
        x_min, x_max = config.xlim
        mask &= (coords[:, 0] >= x_min) & (coords[:, 0] <= x_max)
    if config.ylim is not None:
        y_min, y_max = config.ylim
        mask &= (coords[:, 1] >= y_min) & (coords[:, 1] <= y_max)

    return adata[mask].copy()
