"""
Cell type annotation functions for single-cell RNA-seq data.

This module provides functions to annotate cell clusters with cell type labels
based on marker gene expression.
"""

import sys

sys.path.append("..")

from typing import Literal, Optional

import scanpy as sc

from ..utils import use_layer_as_X
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


def _annotate_by_enrichment(adata, cluster_key, mgr, use_raw):
    """Helper to annotate clusters based on marker gene enrichment."""
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        key_added=f"rank_genes_{cluster_key}",
        use_raw=use_raw
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

# In scRNA/analysis/annotation.py

def score_cell_types(
    adata: sc.AnnData,
    marker_config: str | Manager,
    layer: Optional[str] = "log1p_norm",
    use_raw: bool = False, 
    min_genes: int = 3,
) -> sc.AnnData:
    """
    Score cells for multiple cell types using `sc.tl.score_genes`.

    Args:
        adata: AnnData object.
        marker_config: A Manager instance or path to a marker TOML file.
        layer: Layer to use for scoring if `use_raw` is False.
        use_raw: If True, use `adata.raw` for scoring. This is recommended after
                 HVG selection to score against all genes. Defaults to False.
        min_genes: Minimum number of marker genes required to compute a score.

    Returns:
        AnnData object with score columns added to `adata.obs`.
    """
    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError(
            "marker_config must be a file path (str) or a Manager instance."
        )

    if use_raw:
        if adata.raw is None:
            raise ValueError(
                "adata.raw is not set. Please set it before using use_raw=True, "
                "e.g., after normalization: `adata.raw = adata`."
            )
        print(f"Scoring using adata.raw. The 'layer' parameter ('{layer}') will be ignored.")
        # Intersect markers with genes present in the raw data
        mgr.intersect_with(adata.raw)
        
        # Loop and score using the use_raw flag
        for cell_type, cell in mgr.CELLS.items():
            markers = cell.markers
            if len(markers) >= min_genes:
                print(f"Scoring for cell type: {cell_type} ({len(markers)} markers)")
                sc.tl.score_genes(
                    adata, 
                    markers, 
                    score_name=f"{cell_type}_score", 
                    use_raw=True # Pass use_raw=True to scanpy
                )
            else:
                print(
                    f"Skipping '{cell_type}': only {len(markers)} markers found in adata.raw (min: {min_genes})."
                )
    else:
        # This is the original logic for using a specified layer
        print(f"Scoring using data from layer: '{layer}'")
        mgr.intersect_with(adata)
        
        with use_layer_as_X(adata, layer):
            for cell_type, cell in mgr.CELLS.items():
                markers = cell.markers
                if len(markers) >= min_genes:
                    print(f"Scoring for cell type: {cell_type} ({len(markers)} markers)")
                    sc.tl.score_genes(adata, markers, score_name=f"{cell_type}_score")
                else:
                    print(
                        f"Skipping '{cell_type}': only {len(markers)} markers found in adata.var (min: {min_genes})."
                    )
    return adata

def annotate_clusters(
    adata: sc.AnnData,
    cluster_key: str,
    marker_config: str | Manager,
    method: Literal["max_score", "enrichment"] = "max_score",
    use_raw: bool = False,
    key_added: Optional[str] = None,
) -> sc.AnnData:
    """
    Annotate clusters with cell type labels using pre-computed scores or enrichment.

    Args:
        adata: AnnData object with clustering results.
        cluster_key: Key in `adata.obs` for the clustering to annotate.
        marker_config: Path to the marker TOML configuration file.
        method: Annotation method ('max_score' or 'enrichment').
        use_raw: If True, use `adata.raw` for the 'enrichment' method.
                 This should match the `use_raw` setting in `score_cell_types`.
        key_added: Key in `adata.obs` to store the new annotations. Defaults to f"{cluster_key}_annotated".

    Returns:
        AnnData object with cell type annotations.
    """
    if key_added is None:
        key_added = f"{cluster_key}_annotated"

    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError(
            "marker_config must be a file path (str) or a Manager instance."
        )
    mgr.intersect_with(adata.raw.copy())

    source_adata = adata.raw if use_raw else adata
    if source_adata is None:
        raise ValueError("adata.raw must be set to use `use_raw=True`.")
    mgr.intersect_with(source_adata)

    print(f"Annotating clusters in '{cluster_key}' using '{method}' method...")
    if method == "max_score":
        # This method relies on scores already calculated, so it implicitly
        # respects the use_raw flag from the score_cell_types step.
        mapping = _annotate_by_max_score(adata, cluster_key, mgr)
    elif method == "enrichment":
        # Pass use_raw to the helper
        mapping = _annotate_by_enrichment(adata, cluster_key, mgr, use_raw=use_raw)
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
