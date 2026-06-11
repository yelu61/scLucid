"""Visium IO and coordinate manipulation utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np
from anndata import AnnData

from .config import VisiumIOConfig


def read_visium_10x(
    path: Union[str, Path],
    config: Optional[VisiumIOConfig] = None,
) -> AnnData:
    """Read 10x Visium data, delegating to scanpy or squidpy if available.

    Parameters
    ----------
    path
        Path to the spaceranger output directory.
    config
        IO configuration.

    Returns
    -------
    AnnData
        Loaded Visium AnnData.
    """
    if config is None:
        config = VisiumIOConfig()

    try:
        import squidpy as sq

        return sq.read.visium(path, library_id=config.library_id, load_images=config.load_images)
    except ImportError:
        pass

    try:
        import scanpy as sc

        return sc.read_visium(path, library_id=config.library_id, load_images=config.load_images)
    except ImportError:
        raise ImportError(
            "Reading Visium data requires squidpy or scanpy. Install with: pip install scLucid[spatial]"
        )


def crop_visium(
    adata: AnnData,
    xlim: tuple[float, float],
    ylim: Optional[tuple[float, float]] = None,
    spatial_key: str = "spatial",
) -> AnnData:
    """Crop a Visium AnnData to a coordinate rectangle."""
    from .subset import subset_spatial_window
    from .config import SpatialWindowConfig

    config = SpatialWindowConfig(spatial_key=spatial_key, xlim=xlim, ylim=ylim)
    return subset_spatial_window(adata, config=config)


def rotate_visium(
    adata: AnnData,
    angle_degrees: float,
    spatial_key: str = "spatial",
    center: Optional[tuple[float, float]] = None,
) -> AnnData:
    """Rotate spatial coordinates around a center point.

    Parameters
    ----------
    adata
        Spatial AnnData.
    angle_degrees
        Rotation angle in degrees (counter-clockwise).
    spatial_key
        Key for spatial coordinates in ``adata.obsm``.
    center
        Rotation center (x, y). If None, use the centroid of coordinates.

    Returns
    -------
    AnnData
        AnnData with rotated coordinates.
    """
    if spatial_key not in adata.obsm:
        raise KeyError(f"Spatial key '{spatial_key}' not found in adata.obsm")

    coords = np.asarray(adata.obsm[spatial_key])[:, :2].copy().astype(float)
    if center is None:
        center = coords.mean(axis=0)

    angle = np.deg2rad(angle_degrees)
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    R = np.array([[cos_a, -sin_a], [sin_a, cos_a]])

    rotated = (coords - center) @ R.T + center
    adata.obsm[spatial_key] = rotated
    return adata
