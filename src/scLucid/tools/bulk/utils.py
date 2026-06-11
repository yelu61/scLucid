"""Bulk RNA-seq utility helpers."""

from __future__ import annotations

from typing import List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData


def deduplicate_var_names(adata: AnnData, aggr: str = "sum") -> AnnData:
    """Collapse duplicate gene symbols by summing counts across columns.

    Parameters
    ----------
    adata
        Bulk AnnData.
    aggr
        Aggregation method; currently only ``sum`` is supported.

    Returns
    -------
    AnnData
        New AnnData with unique var_names.
    """
    if aggr != "sum":
        raise ValueError("Only aggr='sum' is supported")

    var_counts = pd.Series(adata.var_names).value_counts()
    duplicates = var_counts[var_counts > 1].index.tolist()
    if not duplicates:
        return adata.copy()

    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    deduped = df.groupby(level=0, axis=1).sum()

    new_adata = adata[:, deduped.columns].copy()
    new_adata.X = deduped.values
    return new_adata


def filter_bulk_genes(
    adata: AnnData,
    min_counts: float = 10.0,
    min_samples: int = 2,
    top_n: Optional[int] = None,
) -> List[str]:
    """Return gene names passing bulk expression filters."""
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    gene_totals = X.sum(axis=0)
    n_samples_expressing = (X > 0).sum(axis=0)
    mask = (gene_totals >= min_counts) & (n_samples_expressing >= min_samples)
    genes = adata.var_names[mask].tolist()
    if top_n is not None:
        genes = genes[:top_n]
    return genes
