"""Bulk RNA-seq normalization utilities."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from .config import BulkNormalizationConfig

log = logging.getLogger(__name__)


def _require_var_column(adata: AnnData, col: str) -> pd.Series:
    if col not in adata.var.columns:
        raise KeyError(f"Gene length column '{col}' not found in adata.var")
    return pd.to_numeric(adata.var[col], errors="coerce")


def normalize_bulk_counts(
    adata: AnnData,
    config: Optional[BulkNormalizationConfig] = None,
) -> pd.DataFrame:
    """Normalize bulk RNA-seq counts to CPM, TPM, FPKM, or RPKM.

    Parameters
    ----------
    adata
        AnnData with samples as observations and genes as variables. ``X`` should
        contain non-negative count-like values.
    config
        Normalization configuration.

    Returns
    -------
    pd.DataFrame
        Normalized expression matrix (samples x genes) with the same index and
        columns as ``adata``.
    """
    if config is None:
        config = BulkNormalizationConfig()

    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    if np.min(X) < 0:
        raise ValueError("normalize_bulk_counts requires non-negative input.")

    method = config.method.upper()
    lib_sizes = X.sum(axis=1, keepdims=True) + config.pseudocount

    if method == "CPM":
        normalized = X / lib_sizes * config.target_sum
    elif method in {"RPKM", "FPKM"}:
        if config.gene_length_col is None:
            raise ValueError(f"{method} normalization requires gene_length_col")
        lengths = _require_var_column(adata, config.gene_length_col).values
        lengths = np.where(lengths > 0, lengths, np.nan)
        with np.errstate(invalid="ignore"):
            normalized = (X / lib_sizes) / (lengths / 1e3) * config.target_sum
        normalized = np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)
    elif method == "TPM":
        if config.gene_length_col is None:
            raise ValueError("TPM normalization requires gene_length_col")
        lengths = _require_var_column(adata, config.gene_length_col).values
        lengths = np.where(lengths > 0, lengths, np.nan)
        with np.errstate(invalid="ignore"):
            rpk = X / (lengths / 1e3)
        rpk = np.nan_to_num(rpk, nan=0.0, posinf=0.0, neginf=0.0)
        rpk_sums = rpk.sum(axis=1, keepdims=True) + config.pseudocount
        normalized = rpk / rpk_sums * config.target_sum
    elif method == "DESEQ2_SIZE_FACTORS":
        size_factors = estimate_size_factors_median_ratio(adata)
        normalized = X / size_factors[:, np.newaxis]
    else:
        raise ValueError(f"Unknown normalization method: {config.method}")

    return pd.DataFrame(
        normalized,
        index=adata.obs_names,
        columns=adata.var_names,
    )


def estimate_size_factors_median_ratio(
    adata: AnnData,
    pseudocount: float = 1.0,
) -> np.ndarray:
    """Estimate DESeq2-style size factors using the median-of-ratios method.

    For each gene, compute the geometric mean across samples. Each sample's size
    factor is the median of its ratios to the gene geometric means.

    Parameters
    ----------
    adata
        AnnData with non-negative count-like values.
    pseudocount
        Pseudocount added before taking logs to avoid log(0).

    Returns
    -------
    np.ndarray
        Size factor per sample (length ``n_obs``).
    """
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    if np.min(X) < 0:
        raise ValueError("Size factor estimation requires non-negative input.")

    # Add pseudocount and compute log geometric mean per gene
    log_counts = np.log(X + pseudocount)
    log_geo_means = np.mean(log_counts, axis=0)

    # Ratio of each sample to geometric mean, per gene
    log_ratios = log_counts - log_geo_means[np.newaxis, :]

    # Median ratio per sample, ignoring genes with -inf geometric mean
    size_factors = np.exp(np.median(log_ratios[:, np.isfinite(log_geo_means)], axis=1))
    return size_factors
