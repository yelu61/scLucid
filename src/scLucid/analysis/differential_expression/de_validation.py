"""Cross-validation between bulk DE and single-cell pseudobulk DE."""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats


def _safe_log2fc(col: pd.Series) -> pd.Series:
    return pd.to_numeric(col, errors="coerce")


def compare_bulk_vs_pseudobulk_de(
    bulk_de: pd.DataFrame,
    pseudobulk_de: pd.DataFrame,
    *,
    gene_col: str = "gene",
    log2fc_col: str = "log2fc",
    pval_col: str = "pvals_adj",
    alpha: float = 0.05,
) -> Dict[str, object]:
    """Compare bulk DE and pseudobulk DE results for concordance.

    Parameters
    ----------
    bulk_de
        Bulk DE result table.
    pseudobulk_de
        Single-cell pseudobulk DE result table.
    gene_col
        Column containing gene identifiers.
    log2fc_col
        Column containing log2 fold changes.
    pval_col
        Column containing adjusted p-values.
    alpha
        Significance threshold.

    Returns
    -------
    dict
        Concordance metrics including log2FC correlation, significant-gene overlap,
        directional concordance, Jaccard index, and hypergeometric test p-value.
    """
    bulk = bulk_de.copy()
    pseudo = pseudobulk_de.copy()

    for df in (bulk, pseudo):
        if gene_col not in df.columns:
            raise KeyError(f"Column '{gene_col}' not found in DE result")
        if log2fc_col not in df.columns:
            raise KeyError(f"Column '{log2fc_col}' not found in DE result")

    bulk["__log2fc"] = _safe_log2fc(bulk[log2fc_col])
    pseudo["__log2fc"] = _safe_log2fc(pseudo[log2fc_col])

    merged = pd.merge(
        bulk[[gene_col, "__log2fc"]],
        pseudo[[gene_col, "__log2fc"]],
        on=gene_col,
        suffixes=("_bulk", "_pseudo"),
    )
    merged = merged.dropna(subset=["__log2fc_bulk", "__log2fc_pseudo"])

    n_shared = len(merged)
    log2fc_corr = float(merged["__log2fc_bulk"].corr(merged["__log2fc_pseudo"], method="spearman"))

    # Directional concordance: same sign and both magnitude > 0
    same_direction = float((merged["__log2fc_bulk"] * merged["__log2fc_pseudo"] > 0).mean())

    # Significant gene overlap
    bulk_sig = set(bulk.loc[bulk[pval_col] < alpha, gene_col]) if pval_col in bulk.columns else set()
    pseudo_sig = set(pseudo.loc[pseudo[pval_col] < alpha, gene_col]) if pval_col in pseudo.columns else set()

    intersection = bulk_sig & pseudo_sig
    union = bulk_sig | pseudo_sig
    jaccard = float(len(intersection) / len(union)) if union else 0.0

    # Hypergeometric overlap p-value
    hypergeom_pval = np.nan
    if bulk_sig and pseudo_sig and n_shared > 0:
        M = n_shared
        N = len(bulk_sig & set(merged[gene_col]))
        n = len(pseudo_sig & set(merged[gene_col]))
        k = len(intersection & set(merged[gene_col]))
        if N > 0 and n > 0 and k >= 0:
            hypergeom_pval = float(stats.hypergeom.sf(k - 1, M, N, n))

    return {
        "schema_version": "bulk_pseudobulk_concordance_v1",
        "n_shared_genes": n_shared,
        "log2fc_spearman": log2fc_corr,
        "directional_concordance": same_direction,
        "bulk_significant_genes": len(bulk_sig),
        "pseudobulk_significant_genes": len(pseudo_sig),
        "shared_significant_genes": len(intersection),
        "jaccard_index": jaccard,
        "hypergeometric_overlap_pval": hypergeom_pval,
        "alpha": alpha,
        "shared_significant_gene_list": sorted(intersection),
    }


def store_bulk_pseudobulk_concordance(
    adata,
    concordance: Dict[str, object],
    key_added: str = "bulk_pseudobulk_concordance",
) -> None:
    """Store concordance result in ``adata.uns['sclucid']['analysis']['de']``."""
    from ...utils import sanitize_for_hdf5

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("de", {})
    adata.uns["sclucid"]["analysis"]["de"][key_added] = sanitize_for_hdf5(concordance)
