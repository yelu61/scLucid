"""Spatial tool utility helpers."""

from __future__ import annotations

from typing import Optional

import numpy as np
from anndata import AnnData


def infer_spatial_platform(adata: AnnData, spatial_key: str = "spatial") -> Optional[str]:
    """Infer spatial platform from AnnData structure."""
    if "spatial" in adata.uns:
        return "visium"
    return None


def validate_spatial_coords(adata: AnnData, spatial_key: str = "spatial") -> None:
    """Validate that spatial coordinates exist and are well-formed."""
    if spatial_key not in adata.obsm:
        raise KeyError(f"Spatial key '{spatial_key}' not found in adata.obsm")
    coords = np.asarray(adata.obsm[spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError(
            f"Spatial coordinates have shape {coords.shape}; expected (n_spots, >=2)"
        )
