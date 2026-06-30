"""Cluster-level marker annotation helpers."""

from __future__ import annotations

import logging
from importlib.metadata import version
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.stats import wilcoxon

from ...utils import Manager, sanitize_for_hdf5

log = logging.getLogger(__name__)


def _benjamini_hochberg(pvalues: pd.Series) -> pd.Series:
    """Apply Benjamini-Hochberg FDR correction to a Series of p-values."""
    p = pvalues.copy()
    n = len(p)
    if n == 0:
        return p
    sorted_idx = p.argsort().values
    sorted_p = p.iloc[sorted_idx].values
    adjusted = sorted_p * n / np.arange(1, n + 1)
    # enforce monotonicity from the bottom up
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    fdr = pd.Series(np.nan, index=p.index)
    fdr.iloc[sorted_idx] = adjusted
    return fdr


def annotate_clusters(
    adata: AnnData,
    cluster_key: str,
    marker_config: Union[str, Manager],
    method: Literal["max_score", "enrichment", "combined"] = "max_score",
    use_raw: bool = False,
    key_added: Optional[str] = None,
    min_confidence: float = 0.3,
    confidence_key: Optional[str] = None,
    min_score: float = 0.1,
    significance_threshold: float = 0.05,
    min_score_margin: float = 0.0,
    n_genes: int = 100,
    score_weight: float = 0.6,
    enrichment_weight: float = 0.4,
    plot: bool = False,
    copy: bool = False,
) -> AnnData:
    """
    Assign cell type labels to clusters using various evidence.

    Enhancements:
    - Robust to missing score columns and non-categorical cluster keys.
    - Stable Unknown handling and category order preservation.
    - Parameter trace including scanpy version and marker stats.
    """
    if copy:
        adata = adata.copy()
    if key_added is None:
        key_added = f"{cluster_key}_annotated"

    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError("marker_config must be a file path or Manager instance.")
    mgr.intersect_with(adata.raw if use_raw and adata.raw is not None else adata)

    # Ensure categorical
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")
    if not isinstance(adata.obs[cluster_key].dtype, pd.CategoricalDtype):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    # 1. Score-based
    def annotate_by_max_score():
        score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
        if not score_cols:
            raise RuntimeError("No *_score columns found. Please run score_cell_types first.")
        means = adata.obs.groupby(cluster_key, observed=False)[score_cols].mean()
        global_medians = adata.obs[score_cols].median()
        result = {}
        evidence = {}
        for cluster in means.index:
            cluster_means = means.loc[cluster]
            # Guard against all-NaN for this cluster
            if cluster_means.isna().all():
                result[str(cluster)] = "Unknown"
                evidence[str(cluster)] = {"assigned_label": "Unknown", "reason": "all_nan"}
                continue

            cluster_scores = adata.obs.loc[adata.obs[cluster_key] == cluster, score_cols]
            pvals = {}
            for col in score_cols:
                scores = cluster_scores[col].dropna()
                median = global_medians[col]
                if len(scores) < 2:
                    pvals[col] = 1.0
                    continue
                try:
                    diff = scores - median
                    if diff.nunique() <= 1:
                        pvals[col] = 1.0
                    else:
                        _, p = wilcoxon(diff, zero_method="zsplit", alternative="greater")
                        pvals[col] = float(p)
                except Exception:
                    pvals[col] = 1.0

            pval_series = pd.Series(pvals)
            fdr = _benjamini_hochberg(pval_series)

            sorted_means = cluster_means.sort_values(ascending=False)
            best = sorted_means.index[0]
            best_score = float(sorted_means.iloc[0])
            runner_up_score = float(sorted_means.iloc[1]) if len(sorted_means) > 1 else 0.0
            margin = best_score - runner_up_score
            cell_type = best[:-6] if best.endswith("_score") else best

            passed = (
                best_score >= min_score
                and fdr[best] <= significance_threshold
                and margin >= min_score_margin
            )
            result[str(cluster)] = cell_type if passed else "Unknown"
            evidence[str(cluster)] = {
                "pvalues": {k: float(v) for k, v in pvals.items()},
                "fdr": {k: float(v) for k, v in fdr.items()},
                "margin": float(margin),
                "best_score": float(best_score),
                "runner_up_score": float(runner_up_score),
                "assigned_label": result[str(cluster)],
            }
        return result, evidence

    # 2. Enrichment-based
    def annotate_by_enrichment():
        rgg_key = f"rank_genes_{cluster_key}"
        # Avoid overwriting existing rank_genes_groups results
        base_rgg_key = rgg_key
        counter = 1
        while rgg_key in adata.uns:
            rgg_key = f"{base_rgg_key}_{counter}"
            counter += 1
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw,
            key_added=rgg_key,
        )
        markers_df = sc.get.rank_genes_groups_df(adata, group=None, key=rgg_key)
        result = {}
        categories = list(adata.obs[cluster_key].cat.categories)
        for cluster in categories:
            genes = (
                markers_df.loc[markers_df["group"] == cluster, "names"]
                .head(n_genes)
                .astype(str)
                .tolist()
            )
            best_score, best_type = -1.0, "Unknown"
            for cell_type, cell in mgr.CELLS.items():
                if not cell.markers:
                    continue
                denom = max(1, len(cell.markers))
                overlap = len(set(genes) & set(cell.markers)) / denom
                if overlap > best_score:
                    best_score, best_type = overlap, cell_type
            result[str(cluster)] = best_type if best_score >= min_score else "Unknown"
        return result

    # 3. Combined
    def annotate_by_combined():
        score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
        if not score_cols:
            raise RuntimeError("No *_score columns found. Please run score_cell_types first.")
        means = adata.obs.groupby(cluster_key, observed=False)[score_cols].mean()
        rgg_key = f"rank_genes_{cluster_key}"
        base_rgg_key = rgg_key
        counter = 1
        while rgg_key in adata.uns:
            rgg_key = f"{base_rgg_key}_{counter}"
            counter += 1
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw,
            key_added=rgg_key,
        )
        markers_df = sc.get.rank_genes_groups_df(adata, group=None, key=rgg_key)
        categories = list(adata.obs[cluster_key].cat.categories)
        result = {}
        for cluster in categories:
            genes = (
                markers_df.loc[markers_df["group"] == cluster, "names"]
                .head(n_genes)
                .astype(str)
                .tolist()
            )
            combined_scores = {}
            for cell_type, cell in mgr.CELLS.items():
                score_col = f"{cell_type}_score"
                score_val = (
                    float(means.loc[cluster, score_col]) if score_col in means.columns else 0.0
                )
                denom = max(1, len(cell.markers))
                overlap_val = (len(set(genes) & set(cell.markers)) / denom) if cell.markers else 0.0
                combined_scores[cell_type] = (
                    score_weight * score_val + enrichment_weight * overlap_val
                )
            best_type = max(combined_scores, key=combined_scores.get)
            best_score = combined_scores[best_type]
            result[str(cluster)] = best_type if best_score >= min_score else "Unknown"
        return result

    # Select method
    if method == "max_score":
        mapping, max_score_evidence = annotate_by_max_score()
    elif method == "enrichment":
        mapping = annotate_by_enrichment()
        max_score_evidence = {}
    elif method == "combined":
        mapping = annotate_by_combined()
        max_score_evidence = {}
    else:
        raise ValueError(f"Unknown annotation method: {method}")

    # Assign labels
    cluster_codes = adata.obs[cluster_key].astype(str)
    assigned = cluster_codes.map(mapping)
    assigned = assigned.fillna("Unknown")
    adata.obs[key_added] = pd.Categorical(assigned)

    # Save and optional plot
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    params_dict = sanitize_for_hdf5(
        {
            "method": method,
            "min_score": min_score,
            "significance_threshold": significance_threshold,
            "min_score_margin": min_score_margin,
            "score_weight": score_weight,
            "enrichment_weight": enrichment_weight,
            "mapping": mapping,
            "scanpy_version": version("scanpy"),
            "n_markers": {k: len(v.markers) for k, v in getattr(mgr, "CELLS", {}).items()},
            "max_score_evidence": max_score_evidence,
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_params"] = params_dict
    if plot:
        if "X_umap" not in adata.obsm:
            sc.tl.umap(adata)
        sc.pl.umap(adata, color=[cluster_key, key_added], wspace=0.4)
    return adata
