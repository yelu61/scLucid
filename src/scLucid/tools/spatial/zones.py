"""Tissue zone detection via NMF on spatial features."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.decomposition import NMF

from .config import TissueZonesConfig

log = logging.getLogger(__name__)


def find_tissue_zones(
    adata: AnnData,
    config: Optional[TissueZonesConfig] = None,
) -> AnnData:
    """Identify tissue zones by NMF on spatial expression or deconvolution features.

    Parameters
    ----------
    adata
        Spatial AnnData.
    config
        Tissue zone configuration.

    Returns
    -------
    AnnData
        ``adata`` with zone loadings in ``obsm[config.key_added]`` and zone
        assignments in ``obs[config.key_added]``.
    """
    if config is None:
        config = TissueZonesConfig()

    if config.input == "expression":
        X = adata.X
        X = np.asarray(X.toarray() if hasattr(X, "toarray") else X)
        features = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    elif config.input == "deconvolution":
        deconv_key = "deconvolution"
        spatial_uns = adata.uns.get("sclucid", {}).get("tools", {}).get("spatial", {})
        if deconv_key not in spatial_uns:
            raise ValueError("input='deconvolution' requires prior spatial deconvolution results")
        proportions = spatial_uns[deconv_key].get("proportions")
        if proportions is None:
            raise ValueError("Deconvolution proportions not found")
        features = proportions.loc[adata.obs_names]
    else:
        raise ValueError(f"Unknown input type: {config.input}")

    features = features.fillna(0.0)
    features = features.loc[:, features.std() > 0]
    if features.shape[1] == 0:
        raise ValueError("No variable features available for NMF")

    n_components = min(config.n_components, features.shape[0], features.shape[1])
    model = NMF(n_components=n_components, init="nndsvda", random_state=42, max_iter=500)
    W = model.fit_transform(np.maximum(features.values, 0.0))
    H = model.components_

    obs_key = config.key_added
    obsm_key = f"X_{config.key_added}"

    adata.obsm[obsm_key] = W
    adata.obs[obs_key] = pd.Categorical(np.argmax(W, axis=1).astype(str))

    adata.uns.setdefault("sclucid", {}).setdefault("tools", {}).setdefault("spatial", {})
    adata.uns["sclucid"]["tools"]["spatial"]["tissue_zones"] = {
        "params": config.to_dict(),
        "n_components": n_components,
        "feature_names": features.columns.tolist(),
        "reconstruction_err": float(model.reconstruction_err_),
        "inference_level": "exploratory_spatial",
        "valid_for_publication_inference": False,
    }

    return adata
