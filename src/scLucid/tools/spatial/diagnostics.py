"""Spatial data quality diagnostics."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

from .config import SpatialDiagnosticsConfig


def diagnose_spatial_data_quality(
    adata: AnnData,
    config: Optional[SpatialDiagnosticsConfig] = None,
) -> Dict[str, Any]:
    """Diagnose whether spatial transcriptomics data is ready for analysis.

    Parameters
    ----------
    adata
        AnnData with spatial coordinates in ``obsm[config.spatial_key]``.
    config
        Diagnostic configuration.

    Returns
    -------
    dict
        Diagnostic report with ``passed``, ``warnings``, and spatial metadata.
    """
    if config is None:
        config = SpatialDiagnosticsConfig()

    warnings: List[str] = []

    if config.spatial_key not in adata.obsm:
        warnings.append(f"Spatial key '{config.spatial_key}' not found in adata.obsm.")
        return {
            "passed": False,
            "warnings": warnings,
            "n_spots": adata.n_obs,
            "n_duplicate_coords": 0,
            "spatial_extent": (np.nan, np.nan, np.nan, np.nan),
            "platform_hint": None,
            "image_key_present": False,
        }

    coords = np.asarray(adata.obsm[config.spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        warnings.append(
            f"Spatial coordinates have shape {coords.shape}; expected (n_spots, >=2)."
        )

    n_spots = int(coords.shape[0])
    nan_coords = int(np.isnan(coords[:, :2]).any(axis=1).sum())
    if nan_coords > 0:
        warnings.append(f"{nan_coords} spots have NaN spatial coordinates.")

    if n_spots < config.min_spots:
        warnings.append(
            f"Only {n_spots} spots; minimum requested is {config.min_spots}."
        )

    n_duplicate_coords = 0
    if config.check_duplicate_coords and coords.shape[0] > 0:
        unique = np.unique(coords[:, :2], axis=0)
        n_duplicate_coords = int(coords.shape[0] - unique.shape[0])
        if n_duplicate_coords > 0:
            warnings.append(f"{n_duplicate_coords} duplicate (x, y) coordinates found.")

    x_min, x_max = float(np.nanmin(coords[:, 0])), float(np.nanmax(coords[:, 0]))
    y_min, y_max = float(np.nanmin(coords[:, 1])), float(np.nanmax(coords[:, 1]))

    image_key_present = False
    if config.require_image:
        image_key_present = "spatial" in adata.uns and any(
            isinstance(v, dict) and "images" in v for v in adata.uns.get("spatial", {}).values()
        )
        if not image_key_present:
            warnings.append("Required spatial image key not found in adata.uns['spatial'].")

    platform_hint = None
    if "platform" in adata.uns:
        platform_hint = str(adata.uns["platform"])
    elif "spatial" in adata.uns:
        # Try to infer from first library entry
        first_lib = next(iter(adata.uns["spatial"].values()), {})
        if "scalefactors" in first_lib:
            platform_hint = "visium"

    passed = not warnings
    return {
        "passed": passed,
        "warnings": warnings,
        "n_spots": n_spots,
        "n_duplicate_coords": n_duplicate_coords,
        "spatial_extent": (x_min, x_max, y_min, y_max),
        "platform_hint": platform_hint,
        "image_key_present": image_key_present,
        "spatial_key": config.spatial_key,
    }
