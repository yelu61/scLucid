"""Bulk RNA-seq data quality diagnostics."""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from .config import BulkDiagnosticsConfig


def diagnose_bulk_data_quality(
    adata: AnnData,
    config: Optional[BulkDiagnosticsConfig] = None,
    condition_col: Optional[str] = None,
) -> Dict[str, Any]:
    """Diagnose whether bulk RNA-seq data is suitable for downstream inference.

    Parameters
    ----------
    adata
        AnnData with samples as observations and genes as variables. ``X`` should
        contain count-like values (non-negative).
    config
        Diagnostic configuration.
    condition_col
        Optional column in ``adata.obs`` defining biological conditions. If provided,
        replicate balance is checked per condition.

    Returns
    -------
    dict
        Diagnostic report with ``passed``, ``warnings``, ``replicate_requirement_met``,
        and descriptive statistics.
    """
    if config is None:
        config = BulkDiagnosticsConfig()

    warnings: List[str] = []
    n_samples = int(adata.n_obs)
    n_genes = int(adata.n_vars)

    # Basic size checks
    if n_samples < config.min_samples_total:
        warnings.append(
            f"Only {n_samples} samples available; minimum requested is {config.min_samples_total}."
        )

    # Library size statistics
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)

    if np.min(X) < 0:
        warnings.append("Expression matrix contains negative values; expected count-like non-negative input.")

    lib_sizes = np.asarray(X.sum(axis=1)).ravel()
    zero_gene_fraction = float(np.mean((X > 0).sum(axis=1) == 0))

    if zero_gene_fraction > 0:
        warnings.append(f"{zero_gene_fraction:.1%} of samples have zero expressed genes.")

    if config.max_zero_gene_fraction is not None and zero_gene_fraction > config.max_zero_gene_fraction:
        warnings.append(
            f"Zero-gene fraction {zero_gene_fraction:.1%} exceeds threshold {config.max_zero_gene_fraction:.1%}."
        )

    library_size_cv = float(np.std(lib_sizes) / (np.mean(lib_sizes) + 1e-12))
    if config.max_library_size_cv is not None and library_size_cv > config.max_library_size_cv:
        warnings.append(
            f"Library size CV ({library_size_cv:.2f}) exceeds threshold ({config.max_library_size_cv:.2f})."
        )

    # Condition/replicate checks
    replicate_requirement_met = False
    n_conditions = 1
    min_replicates = n_samples
    max_replicates = n_samples

    if condition_col is not None and condition_col in adata.obs.columns:
        cond_counts = adata.obs[condition_col].value_counts(dropna=False)
        n_conditions = int(cond_counts.shape[0])
        min_replicates = int(cond_counts.min())
        max_replicates = int(cond_counts.max())

        if n_conditions < 2:
            warnings.append(f"Only {n_conditions} condition level found; need at least 2 for DE.")

        if config.require_replicates and min_replicates < config.min_samples_per_condition:
            warnings.append(
                f"Minimum replicate count is {min_replicates}; "
                f"requested at least {config.min_samples_per_condition} per condition."
            )
        else:
            replicate_requirement_met = True
    else:
        if condition_col is not None:
            warnings.append(f"Condition column '{condition_col}' not found in adata.obs.")
        if n_samples >= 2:
            replicate_requirement_met = not config.require_replicates or n_samples >= config.min_samples_per_condition

    # Normalization state heuristics
    fraction_integer = float(np.mean(np.abs(X - np.rint(X)) < 1e-6)) if X.size else 0.0
    if fraction_integer < 0.95:
        warnings.append(
            f"Only {fraction_integer:.1%} of values are integer-like; input may already be normalized."
        )

    passed = not warnings
    recommended_method = "welch" if replicate_requirement_met else "descriptive"

    return {
        "passed": passed,
        "warnings": warnings,
        "replicate_requirement_met": replicate_requirement_met,
        "n_samples": n_samples,
        "n_genes": n_genes,
        "n_conditions": n_conditions,
        "min_replicates": min_replicates,
        "max_replicates": max_replicates,
        "library_size_mean": float(np.mean(lib_sizes)),
        "library_size_std": float(np.std(lib_sizes)),
        "library_size_cv": library_size_cv,
        "zero_gene_fraction": zero_gene_fraction,
        "fraction_integer": fraction_integer,
        "recommended_method": recommended_method,
        "condition_col": condition_col,
    }
