"""Bulk differential expression analysis."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from .config import BulkDEConfig
from .diagnostics import diagnose_bulk_data_quality

log = logging.getLogger(__name__)


def _benjamini_hochberg(pvals: pd.Series, method: str = "fdr_bh") -> pd.Series:
    """Adjust p-values using Benjamini-Hochberg or Bonferroni."""
    from statsmodels.stats.multitest import multipletests

    valid = pvals.notna() & np.isfinite(pvals)
    adjusted = pd.Series(np.nan, index=pvals.index)
    if valid.any():
        method_map = {"fdr_bh": "fdr_bh", "bonferroni": "bonferroni"}
        sm_method = method_map.get(method, "fdr_bh")
        _, qvals, _, _ = multipletests(pvals[valid].values, method=sm_method)
        adjusted.loc[valid] = qvals
    return adjusted


def _filter_low_expression(
    X: np.ndarray,
    gene_names: pd.Index,
    min_counts_per_gene: int,
    min_samples_expressing: int,
) -> np.ndarray:
    """Return boolean mask of genes passing expression filters."""
    gene_totals = X.sum(axis=0)
    n_samples_expressing = (X > 0).sum(axis=0)
    mask = (gene_totals >= min_counts_per_gene) & (n_samples_expressing >= min_samples_expressing)
    return np.asarray(mask).ravel()


def _run_welch_de(
    counts: np.ndarray,
    gene_names: pd.Index,
    condition_mask: pd.Series,
    condition1: str,
    condition2: str,
) -> pd.DataFrame:
    """Run Welch t-test DE on log2(CPM+1) values."""
    # Library-size normalize + log transform
    lib_sizes = counts.sum(axis=1, keepdims=True) + 1e-6
    cpm = counts / lib_sizes * 1e6
    logcpm = np.log2(cpm + 1.0)

    g1 = logcpm[condition_mask.values]
    g2 = logcpm[~condition_mask.values]

    records = []
    for idx, gene in enumerate(gene_names):
        a = g1[:, idx]
        b = g2[:, idx]
        if a.std() == 0 and b.std() == 0:
            continue
        stat, pval = stats.ttest_ind(b, a, equal_var=False)
        mean1 = float(a.mean())
        mean2 = float(b.mean())
        log2fc = mean2 - mean1
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logfoldchanges": log2fc,
                "log2fc": log2fc,
                "mean_logcpm_condition1": mean1,
                "mean_logcpm_condition2": mean2,
                "scores": stat,
                "statistic": stat,
                "pvals": pval,
                "pval": pval,
                "condition1": condition1,
                "condition2": condition2,
                "direction": f"{condition2} - {condition1}",
            }
        )

    return pd.DataFrame(records)


def _run_ttest_de(
    counts: np.ndarray,
    gene_names: pd.Index,
    condition_mask: pd.Series,
    condition1: str,
    condition2: str,
) -> pd.DataFrame:
    """Run Student's t-test DE on log2(CPM+1) values."""
    lib_sizes = counts.sum(axis=1, keepdims=True) + 1e-6
    cpm = counts / lib_sizes * 1e6
    logcpm = np.log2(cpm + 1.0)

    g1 = logcpm[condition_mask.values]
    g2 = logcpm[~condition_mask.values]

    records = []
    for idx, gene in enumerate(gene_names):
        a = g1[:, idx]
        b = g2[:, idx]
        if a.std() == 0 and b.std() == 0:
            continue
        stat, pval = stats.ttest_ind(b, a, equal_var=True)
        mean1 = float(a.mean())
        mean2 = float(b.mean())
        log2fc = mean2 - mean1
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logfoldchanges": log2fc,
                "log2fc": log2fc,
                "mean_logcpm_condition1": mean1,
                "mean_logcpm_condition2": mean2,
                "scores": stat,
                "statistic": stat,
                "pvals": pval,
                "pval": pval,
                "condition1": condition1,
                "condition2": condition2,
                "direction": f"{condition2} - {condition1}",
            }
        )

    return pd.DataFrame(records)


def _run_descriptive_de(
    counts: np.ndarray,
    gene_names: pd.Index,
    condition_mask: pd.Series,
    condition1: str,
    condition2: str,
) -> pd.DataFrame:
    """Return effect sizes only, without p-values, for single-sample contrasts."""
    lib_sizes = counts.sum(axis=1, keepdims=True) + 1e-6
    cpm = counts / lib_sizes * 1e6
    logcpm = np.log2(cpm + 1.0)

    g1 = logcpm[condition_mask.values]
    g2 = logcpm[~condition_mask.values]

    records = []
    for idx, gene in enumerate(gene_names):
        a = g1[:, idx]
        b = g2[:, idx]
        mean1 = float(a.mean())
        mean2 = float(b.mean())
        log2fc = mean2 - mean1
        records.append(
            {
                "names": gene,
                "gene": gene,
                "logfoldchanges": log2fc,
                "log2fc": log2fc,
                "mean_logcpm_condition1": mean1,
                "mean_logcpm_condition2": mean2,
                "scores": np.nan,
                "statistic": np.nan,
                "pvals": np.nan,
                "pval": np.nan,
                "condition1": condition1,
                "condition2": condition2,
                "direction": f"{condition2} - {condition1}",
            }
        )

    return pd.DataFrame(records)


def run_bulk_de(
    adata: AnnData,
    config: Optional[BulkDEConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """Run bulk differential expression analysis.

    Supports Welch t-test, Student's t-test, and optional pydeseq2/limma backends.
    Results are annotated with inference-level semantics.

    Parameters
    ----------
    adata
        Bulk RNA-seq AnnData with samples as observations.
    config
        DE configuration.
    **kwargs
        Overrides for config fields.

    Returns:
    -------
    pd.DataFrame
        DE result table with ``log2fc``, ``pvals``, ``pvals_adj``, and inference tags.
    """
    if config is None:
        config = BulkDEConfig(**kwargs)
    else:
        config = config.model_copy(update=kwargs)

    condition_col = config.condition_col
    if condition_col not in adata.obs.columns:
        raise KeyError(f"Condition column '{condition_col}' not found in adata.obs")

    cond_series = adata.obs[condition_col].astype(str)
    if config.condition1 not in cond_series.values:
        raise ValueError(f"Condition '{config.condition1}' not found in column '{condition_col}'")
    if config.condition2 not in cond_series.values:
        raise ValueError(f"Condition '{config.condition2}' not found in column '{condition_col}'")

    compare_mask = cond_series.isin([config.condition1, config.condition2])
    adata_compare = adata[compare_mask].copy()
    cond_series = adata_compare.obs[condition_col].astype(str)

    # Diagnostics are scoped to the requested two-condition contrast, so unrelated
    # condition levels cannot affect replicate counts or DE inputs.
    diag = diagnose_bulk_data_quality(adata_compare, condition_col=condition_col)

    X = adata_compare.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    counts = np.asarray(X, dtype=float)

    # Filter low-expression genes
    keep_genes = _filter_low_expression(
        counts,
        adata.var_names,
        config.min_counts_per_gene,
        config.min_samples_expressing,
    )
    counts = counts[:, keep_genes]
    gene_names = adata_compare.var_names[keep_genes]

    condition_mask = cond_series == config.condition1
    n1 = int((cond_series == config.condition1).sum())
    n2 = int((cond_series == config.condition2).sum())

    has_replicates = n1 >= config.min_samples_per_condition and n2 >= config.min_samples_per_condition
    if not has_replicates and not config.fallback_to_descriptive:
        raise ValueError(
            f"Insufficient replicates for DE ({n1} vs {n2}). "
            "Set fallback_to_descriptive=True to return effect sizes only."
        )

    # Dispatch
    if not has_replicates and config.fallback_to_descriptive:
        result = _run_descriptive_de(counts, gene_names, condition_mask, config.condition1, config.condition2)
        inference_level = "descriptive_sample_level"
        valid_for_publication = False
        result_warning = "Only one sample in at least one condition; returned descriptive effect sizes without formal p-values."
    elif config.method in {"welch", "ttest"}:
        if config.method == "welch":
            result = _run_welch_de(counts, gene_names, condition_mask, config.condition1, config.condition2)
        else:
            result = _run_ttest_de(counts, gene_names, condition_mask, config.condition1, config.condition2)
        inference_level = "sample_level"
        valid_for_publication = n1 >= 2 and n2 >= 2
        result_warning = None
    elif config.method == "pydeseq2":
        result = _run_pydeseq2_de(adata_compare[:, keep_genes].copy(), config)
        inference_level = "sample_level"
        valid_for_publication = n1 >= 2 and n2 >= 2
        result_warning = None
    elif config.method == "limma":
        raise NotImplementedError("limma backend not yet implemented; use 'welch' or 'pydeseq2'")
    else:
        raise ValueError(f"Unknown DE method: {config.method}")

    if result.empty:
        return result

    result["pvals_adj"] = _benjamini_hochberg(result["pvals"], method=config.p_adjust_method)
    result["padj"] = result["pvals_adj"]
    result["inference_level"] = inference_level
    result["valid_for_publication_inference"] = bool(valid_for_publication)
    result["replicate_requirement_met"] = has_replicates
    result["diagnostic_status"] = "passed" if diag["passed"] else "warning"
    result["n_samples_condition1"] = n1
    result["n_samples_condition2"] = n2
    result["method"] = config.method
    result["result_warning"] = result_warning

    bulk_ns = adata.uns.setdefault("sclucid", {}).setdefault("tools", {}).setdefault("bulk", {})
    bulk_ns["diagnostics"] = diag
    bulk_ns["de"] = {
        "params": config.to_dict(),
        "n_genes_tested": int(result.shape[0]),
        "n_samples_condition1": n1,
        "n_samples_condition2": n2,
        "condition1": config.condition1,
        "condition2": config.condition2,
        "method": config.method,
        "diagnostic_status": "passed" if diag["passed"] else "warning",
        "inference_level": inference_level,
        "valid_for_publication_inference": bool(valid_for_publication),
        "replicate_requirement_met": bool(has_replicates),
        "result_warning": result_warning,
    }

    return result


def _run_pydeseq2_de(adata: AnnData, config: BulkDEConfig) -> pd.DataFrame:
    """Run pydeseq2 on bulk counts."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise ImportError("pydeseq2 is required for method='pydeseq2'. Install with: pip install scLucid[bulk]") from exc

    meta = adata.obs[[config.condition_col]].copy()
    design_col = "__condition"
    meta[design_col] = pd.Categorical(
        meta[config.condition_col].astype(str),
        categories=[config.condition1, config.condition2],
    )

    counts = pd.DataFrame(
        np.asarray(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X).astype(int),
        index=adata.obs_names,
        columns=adata.var_names,
    )

    dds = DeseqDataSet(counts=counts, metadata=meta, design=f"~{design_col}", quiet=True, n_cpus=1)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=[design_col, config.condition2, config.condition1], quiet=True, n_cpus=1)
    ds.summary()
    df = ds.results_df.copy()

    df = df.rename(
        columns={
            "log2FoldChange": "logfoldchanges",
            "pvalue": "pvals",
            "padj": "pvals_adj",
            "stat": "scores",
            "baseMean": "base_mean",
        }
    )
    df["names"] = df.index.astype(str)
    df["gene"] = df["names"]
    df["log2fc"] = df["logfoldchanges"]
    df["statistic"] = df["scores"]
    df["pval"] = df["pvals"]
    df["padj"] = df["pvals_adj"]
    df["condition1"] = config.condition1
    df["condition2"] = config.condition2
    df["direction"] = f"{config.condition2} - {config.condition1}"
    return df
