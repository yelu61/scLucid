"""
Cell type annotation functions for single-cell RNA-seq data.

This module provides functions to annotate cell clusters with cell type labels
based on marker gene expression.
"""

from typing import Literal, Optional

import scanpy as sc

from .manager import Manager

# --- Helper Annotation Functions ---


def _annotate_by_max_score(adata, cluster_key, mgr):
    """Helper to annotate clusters based on the highest average marker score."""
    score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
    if not score_cols:
        raise RuntimeError("No score columns found. Run `score_cell_types` first.")

    cluster_means = adata.obs.groupby(cluster_key)[score_cols].mean()
    best_cell_type = cluster_means.idxmax(axis=1)
    return best_cell_type.to_dict()


def _annotate_by_enrichment(adata, cluster_key, mgr):
    """Helper to annotate clusters based on marker gene enrichment."""
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        key_added=f"rank_genes_{cluster_key}",
    )
    markers_df = sc.get.rank_genes_groups_df(
        adata, key=f"rank_genes_{cluster_key}", group=None
    )

    annotations = {}
    for cluster in adata.obs[cluster_key].cat.categories:
        cluster_markers = (
            markers_df[markers_df["group"] == cluster]["names"].head(100).tolist()
        )
        best_score, best_cell_type = -1, "Unknown"

        for cell_type, cell in mgr.CELLS.items():
            known_markers = set(cell.markers)
            overlap = len(set(cluster_markers).intersection(known_markers))

            # Simple overlap score, can be replaced with more complex stats
            score = overlap / len(known_markers) if known_markers else 0
            if score > best_score:
                best_score, best_cell_type = score, cell_type
        annotations[cluster] = best_cell_type
    return annotations


# --- Main Annotation Functions ---


def score_cell_types(
    adata: sc.AnnData,
    marker_config: str,
    layer: Optional[str] = "log1p_norm",
    min_genes: int = 3,
) -> sc.AnnData:
    """
    Score cells for multiple cell types using `sc.tl.score_genes`.

    Args:
        adata: AnnData object.
        marker_config: Path to the marker TOML configuration file.
        layer: Layer to use for scoring.
        min_genes: Minimum number of marker genes required to compute a score.

    Returns:
        AnnData object with score columns added to `adata.obs`.
    """
    mgr = Manager(marker_config)
    mgr.intersect_with(adata)

    with use_layer_as_X(adata, layer):
        for cell_type, cell in mgr.CELLS.items():
            markers = cell.markers
            if len(markers) >= min_genes:
                print(f"Scoring for cell type: {cell_type} ({len(markers)} markers)")
                sc.tl.score_genes(adata, markers, score_name=f"{cell_type}_score")
            else:
                print(
                    f"Skipping '{cell_type}': only {len(markers)} markers found (min: {min_genes})."
                )
    return adata


def annotate_clusters(
    adata: sc.AnnData,
    cluster_key: str,
    marker_config: str,
    method: Literal["max_score", "enrichment"] = "max_score",
    key_added: Optional[str] = None,
) -> sc.AnnData:
    """
    Annotate clusters with cell type labels using pre-computed scores or enrichment.

    Args:
        adata: AnnData object with clustering results.
        cluster_key: Key in `adata.obs` for the clustering to annotate.
        marker_config: Path to the marker TOML configuration file.
        method: Annotation method ('max_score' or 'enrichment').
        key_added: Key in `adata.obs` to store the new annotations. Defaults to f"{cluster_key}_annotated".

    Returns:
        AnnData object with cell type annotations.
    """
    if key_added is None:
        key_added = f"{cluster_key}_annotated"

    mgr = Manager(marker_config)
    mgr.intersect_with(adata)

    print(f"Annotating clusters in '{cluster_key}' using '{method}' method...")
    if method == "max_score":
        mapping = _annotate_by_max_score(adata, cluster_key, mgr)
    elif method == "enrichment":
        mapping = _annotate_by_enrichment(adata, cluster_key, mgr)
    else:
        raise ValueError(f"Unknown annotation method: {method}")

    adata.obs[key_added] = adata.obs[cluster_key].map(mapping).astype("category")

    # Add colors for plotting
    adata.uns[f"{key_added}_colors"] = [
        mgr[cat].color
        for cat in adata.obs[key_added].cat.categories
        if cat in mgr.CELLS and mgr[cat].color
    ]
    print("Annotation complete.")
    return adata
