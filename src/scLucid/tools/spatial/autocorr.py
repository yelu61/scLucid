"""Spatial autocorrelation statistics (Moran's I, Geary's c)."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse as sp
from scipy import stats

from .config import SpatialAutocorrConfig
from .neighbors import build_spatial_neighbors

log = logging.getLogger(__name__)


def _moran_i(values: np.ndarray, weights: sp.spmatrix) -> Tuple[float, float]:
    """Compute Moran's I and a permutation-based p-value.

    Returns (I, pval_two_sided) where pval is approximate and based on a small
    number of permutations when requested by the caller.
    """
    n = len(values)
    z = values - np.mean(values)
    z = np.asarray(z).ravel()

    W = weights.astype(float)
    if sp.issparse(W):
        num = float(z @ (W @ z))
        S0 = float(W.sum())
    else:
        num = float(z @ W @ z)
        S0 = float(np.asarray(W).sum())

    denom = float(np.sum(z**2))
    if denom == 0 or S0 == 0:
        return np.nan, np.nan

    I = (n / S0) * (num / denom)
    return I, np.nan


def _geary_c(values: np.ndarray, weights: sp.spmatrix) -> float:
    """Compute Geary's c spatial autocorrelation statistic."""
    n = len(values)
    z = np.asarray(values).ravel()
    diff = z[:, None] - z[None, :]

    W = weights.astype(float)
    if sp.issparse(W):
        W = W.toarray()
    W = np.asarray(W)
    S0 = W.sum()
    if S0 == 0:
        return np.nan

    num = float(np.sum(W * (diff**2)))
    denom = float(2 * np.sum((z - z.mean()) ** 2))
    if denom == 0:
        return np.nan

    return ((n - 1) / (2 * S0)) * (num / denom)


def compute_spatial_autocorr(
    adata: AnnData,
    values: np.ndarray,
    config: Optional[SpatialAutocorrConfig] = None,
) -> Dict[str, float]:
    """Compute Moran's I or Geary's c for a single numeric vector.

    Parameters
    ----------
    adata
        AnnData with a spatial neighbor graph already built, or with spatial
        coordinates from which neighbors will be built.
    values
        Numeric values aligned with ``adata.obs_names``.
    config
        Autocorrelation configuration.

    Returns
    -------
    dict
        Statistic and optional permutation p-value.
    """
    if config is None:
        config = SpatialAutocorrConfig()

    if "spatial_connectivities" not in adata.obsp:
        adata = build_spatial_neighbors(adata)

    W = adata.obsp["spatial_connectivities"]

    if config.mode == "moran":
        I, _ = _moran_i(values, W)
        result = {"moran_i": I, "pval": np.nan}
        if config.n_permutations > 0:
            observed = I
            permuted = []
            rng = np.random.default_rng(42)
            for _ in range(config.n_permutations):
                perm = rng.permutation(values)
                pi, _ = _moran_i(perm, W)
                permuted.append(pi)
            permuted = np.asarray(permuted)
            result["pval"] = float(np.mean(np.abs(permuted) >= np.abs(observed)))
            result["pval"] = max(result["pval"], 1.0 / (config.n_permutations + 1))
    elif config.mode == "geary":
        result = {"geary_c": _geary_c(values, W), "pval": np.nan}
    else:
        raise ValueError(f"Unknown autocorrelation mode: {config.mode}")

    result["mode"] = config.mode
    result["n_permutations"] = config.n_permutations
    result["inference_level"] = "exploratory_spatial"
    result["valid_for_publication_inference"] = False
    return result


def compute_moran_i(
    adata: AnnData,
    values: np.ndarray,
    n_permutations: int = 0,
) -> Dict[str, float]:
    """Convenience wrapper for Moran's I."""
    config = SpatialAutocorrConfig(mode="moran", n_permutations=n_permutations)
    return compute_spatial_autocorr(adata, values, config=config)
