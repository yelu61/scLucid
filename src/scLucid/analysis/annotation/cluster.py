"""Cluster-level marker annotation helpers."""

from __future__ import annotations

from importlib.metadata import version
from typing import Literal, Optional, Union

import logging
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ...utils import Manager, sanitize_for_hdf5

log = logging.getLogger(__name__)


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
    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    # 1. Score-based
    def annotate_by_max_score():
        score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
        if not score_cols:
            raise RuntimeError("No *_score columns found. Please run score_cell_types first.")
        means = adata.obs.groupby(cluster_key)[score_cols].mean()
        result = {}
        for cluster in means.index:
            best = means.loc[cluster].idxmax()
            best_score = float(means.loc[cluster, best])
            cell_type = best[:-6] if best.endswith("_score") else best
            result[str(cluster)] = cell_type if best_score >= min_score else "Unknown"
        return result

    # 2. Enrichment-based
    def annotate_by_enrichment():
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw,
            key_added=f"rank_genes_{cluster_key}",
        )
        markers_df = sc.get.rank_genes_groups_df(adata, group=None, key=f"rank_genes_{cluster_key}")
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
        means = adata.obs.groupby(cluster_key)[score_cols].mean()
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw,
            key_added=f"rank_genes_{cluster_key}",
        )
        markers_df = sc.get.rank_genes_groups_df(adata, group=None, key=f"rank_genes_{cluster_key}")
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
        mapping = annotate_by_max_score()
    elif method == "enrichment":
        mapping = annotate_by_enrichment()
    elif method == "combined":
        mapping = annotate_by_combined()
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
            "score_weight": score_weight,
            "enrichment_weight": enrichment_weight,
            "mapping": mapping,
            "scanpy_version": version("scanpy"),
            "n_markers": {k: len(v.markers) for k, v in getattr(mgr, "CELLS", {}).items()},
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_params"] = params_dict
    if plot:
        if "X_umap" not in adata.obsm:
            sc.tl.umap(adata)
        sc.pl.umap(adata, color=[cluster_key, key_added], wspace=0.4)
    return adata
