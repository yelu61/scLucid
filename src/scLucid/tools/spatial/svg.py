"""Spatially variable gene detection."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from .autocorr import compute_moran_i
from .config import SVGConfig
from .neighbors import build_spatial_neighbors

log = logging.getLogger(__name__)


def find_spatially_variable_genes(
    adata: AnnData,
    config: Optional[SVGConfig] = None,
) -> pd.DataFrame:
    """Detect spatially variable genes using Moran's I.

    Parameters
    ----------
    adata
        AnnData with spatial coordinates and gene expression.
    config
        SVG configuration.

    Returns
    -------
    pd.DataFrame
        Per-gene Moran's I, p-value, and adjusted p-value.
    """
    if config is None:
        config = SVGConfig()

    if config.spatial_key not in adata.obsm:
        raise KeyError(f"Spatial key '{config.spatial_key}' not found in adata.obsm")

    if "spatial_connectivities" not in adata.obsp:
        build_spatial_neighbors(adata)

    X = adata.X
    if config.layer is not None:
        if config.layer not in adata.layers:
            raise KeyError(f"Layer '{config.layer}' not found in adata.layers")
        X = adata.layers[config.layer]
    X = np.asarray(X.toarray() if hasattr(X, "toarray") else X)

    results = []
    for i, gene in enumerate(adata.var_names):
        vals = X[:, i]
        if vals.std() == 0:
            continue
        moran = compute_moran_i(adata, vals, n_permutations=config.n_permutations)
        results.append(
            {
                "gene": gene,
                "moran_i": moran["moran_i"],
                "pval": moran["pval"],
                "n_permutations": moran["n_permutations"],
            }
        )

    df = pd.DataFrame(results)
    if df.empty:
        return df

    from statsmodels.stats.multitest import multipletests

    valid = df["pval"].notna() & np.isfinite(df["pval"])
    df["pvals_adj"] = np.nan
    if valid.any():
        _, qvals, _, _ = multipletests(df.loc[valid, "pval"].values, method="fdr_bh")
        df.loc[valid, "pvals_adj"] = qvals

    df["spatially_variable"] = df["pvals_adj"] < config.alpha
    df["inference_level"] = "exploratory_spatial"
    df["valid_for_publication_inference"] = False
    df["method"] = config.method

    adata.var[config.key_added] = df.set_index("gene").loc[adata.var_names, "spatially_variable"].fillna(False).values
    adata.uns.setdefault("sclucid", {}).setdefault("tools", {}).setdefault("spatial", {})
    adata.uns["sclucid"]["tools"]["spatial"]["svg"] = df.to_dict(orient="records")

    return df.sort_values("moran_i", ascending=False)
