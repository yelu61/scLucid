"""
Cell type annotation functions for single-cell RNA-seq data.

This module provides functions to annotate cell clusters with cell type labels
based on marker gene expression. It supports multiple annotation methods:
- Score-based annotation (using pre-computed cell type scores)
- Enrichment-based annotation (using marker gene overlap with cluster-specific genes)
- Reference-based annotation (using reference datasets for label transfer)
- Combined approaches (ensemble methods using multiple sources of evidence)

These functions are designed to work with the marker gene management system to
enable flexible, reproducible cell type annotation across different tissues.
"""

import logging
import os
from typing import Dict, Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ..utils import use_layer_as_X
from .manager import Manager

# Configure logging
log = logging.getLogger(__name__)

# --- Helper Annotation Functions ---


def _annotate_by_max_score(
    adata: AnnData,
    cluster_key: str,
    mgr: Manager,
    min_score: float = 0.0,
    custom_score_pattern: Optional[str] = None,
) -> Dict[str, str]:
    """
    Helper to annotate clusters based on the highest average marker score.

    Args:
        adata: AnnData object with pre-computed scores
        cluster_key: Key in adata.obs for cluster assignments
        mgr: Manager instance with cell type definitions
        min_score: Minimum score required to assign a cell type (0-1)
        custom_score_pattern: Custom pattern to identify score columns

    Returns:
        Dictionary mapping cluster IDs to cell type names

    Raises:
        RuntimeError: If no score columns are found
    """
    # Identify score columns
    if custom_score_pattern:
        score_cols = [col for col in adata.obs.columns if custom_score_pattern in col]
    else:
        score_cols = [col for col in adata.obs.columns if col.endswith("_score")]

    if not score_cols:
        log.error("No score columns found in adata.obs")
        raise RuntimeError(
            "No score columns found. Run `score_cell_types` first or specify a custom score pattern."
        )

    log.info(f"Found {len(score_cols)} score columns for annotation")

    # Calculate mean scores per cluster
    cluster_means = adata.obs.groupby(cluster_key)[score_cols].mean()

    # Find the best cell type for each cluster
    annotations = {}
    for cluster in cluster_means.index:
        # Get the cell type with the highest score
        best_score_col = cluster_means.loc[cluster].idxmax()
        best_score = cluster_means.loc[cluster, best_score_col]

        # Extract cell type name from score column name
        if custom_score_pattern:
            best_cell_type = best_score_col.replace(custom_score_pattern, "")
        else:
            best_cell_type = best_score_col.replace("_score", "")

        # Apply minimum score threshold
        if best_score < min_score:
            log.warning(
                f"Cluster {cluster}: best score ({best_score:.3f}) below threshold ({min_score}), annotating as 'Unknown'"
            )
            annotations[cluster] = "Unknown"
        else:
            annotations[cluster] = best_cell_type
            log.debug(
                f"Cluster {cluster} -> {best_cell_type} (score: {best_score:.3f})"
            )

    log.info(f"Annotated {len(annotations)} clusters based on maximum scores")
    return annotations


def _annotate_by_enrichment(
    adata: AnnData,
    cluster_key: str,
    mgr: Manager,
    use_raw: bool = False,
    n_genes: int = 100,
    min_overlap_score: float = 0.1,
    significance_threshold: float = 0.05,
) -> Dict[str, str]:
    """
    Helper to annotate clusters based on marker gene enrichment.

    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs for cluster assignments
        mgr: Manager instance with cell type definitions
        use_raw: Whether to use raw data for differential expression
        n_genes: Number of top differentially expressed genes to consider
        min_overlap_score: Minimum overlap score required to assign a cell type
        significance_threshold: P-value threshold for differentially expressed genes

    Returns:
        Dictionary mapping cluster IDs to cell type names
    """
    log.info(
        f"Calculating differentially expressed genes for each cluster in '{cluster_key}'"
    )

    # Calculate marker genes for each cluster
    try:
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            key_added=f"rank_genes_{cluster_key}",
            use_raw=use_raw,
        )

        # Get the full results
        markers_df = sc.get.rank_genes_groups_df(
            adata, key=f"rank_genes_{cluster_key}", group=None
        )

        # Filter by significance
        markers_df = markers_df[markers_df["pvals_adj"] < significance_threshold]

    except Exception as e:
        log.error(f"Error calculating differentially expressed genes: {str(e)}")
        raise RuntimeError(
            f"Failed to calculate differentially expressed genes: {str(e)}"
        )

    # Annotate each cluster
    annotations = {}
    for cluster in adata.obs[cluster_key].cat.categories:
        # Get cluster-specific markers
        cluster_markers_df = markers_df[markers_df["group"] == cluster]

        if len(cluster_markers_df) == 0:
            log.warning(f"No significant markers found for cluster {cluster}")
            annotations[cluster] = "Unknown"
            continue

        # Take top N genes
        cluster_markers = cluster_markers_df["names"].head(n_genes).tolist()
        log.debug(f"Cluster {cluster}: {len(cluster_markers)} significant markers")

        # Calculate overlap with known cell type markers
        best_score, best_cell_type, best_overlap = -1, "Unknown", []

        for cell_type, cell in mgr.CELLS.items():
            known_markers = set(cell.markers)
            if not known_markers:
                continue

            # Find overlapping markers
            overlapping = set(cluster_markers).intersection(known_markers)

            if len(overlapping) == 0:
                continue

            # Calculate Jaccard similarity: |A ∩ B| / |A ∪ B|
            # jaccard = len(overlapping) / len(set(cluster_markers).union(known_markers))

            # Calculate overlap score: |A ∩ B| / |B|
            overlap_score = len(overlapping) / len(known_markers)

            if overlap_score > best_score:
                best_score = overlap_score
                best_cell_type = cell_type
                best_overlap = overlapping

        # Apply minimum score threshold
        if best_score < min_overlap_score:
            log.warning(
                f"Cluster {cluster}: best overlap score ({best_score:.3f}) below threshold ({min_overlap_score}), annotating as 'Unknown'"
            )
            annotations[cluster] = "Unknown"
        else:
            annotations[cluster] = best_cell_type
            log.debug(
                f"Cluster {cluster} -> {best_cell_type} (score: {best_score:.3f}, overlap: {', '.join(best_overlap)})"
            )

    log.info(f"Annotated {len(annotations)} clusters based on marker enrichment")
    return annotations


def _annotate_by_combined_evidence(
    adata: AnnData,
    cluster_key: str,
    mgr: Manager,
    score_weight: float = 0.6,
    enrichment_weight: float = 0.4,
    min_combined_score: float = 0.2,
    use_raw: bool = False,
    n_genes: int = 100,
) -> Dict[str, str]:
    """
    Helper to annotate clusters based on combined evidence from scores and enrichment.

    Args:
        adata: AnnData object
        cluster_key: Key in adata.obs for cluster assignments
        mgr: Manager instance with cell type definitions
        score_weight: Weight for score-based evidence (0-1)
        enrichment_weight: Weight for enrichment-based evidence (0-1)
        min_combined_score: Minimum combined score to assign a cell type
        use_raw: Whether to use raw data for enrichment analysis
        n_genes: Number of top genes to consider for enrichment

    Returns:
        Dictionary mapping cluster IDs to cell type names
    """
    # Validate weights
    if not np.isclose(score_weight + enrichment_weight, 1.0):
        log.warning(
            f"Weights do not sum to 1.0 (score_weight={score_weight}, enrichment_weight={enrichment_weight})"
        )
        # Normalize weights
        total = score_weight + enrichment_weight
        score_weight /= total
        enrichment_weight /= total
        log.warning(
            f"Normalized weights: score_weight={score_weight}, enrichment_weight={enrichment_weight}"
        )

    # Get annotations from both methods
    score_annotations = _annotate_by_max_score(adata, cluster_key, mgr, min_score=0.0)
    enrichment_annotations = _annotate_by_enrichment(
        adata, cluster_key, mgr, use_raw=use_raw, n_genes=n_genes
    )

    # Prepare score matrices
    score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
    if not score_cols:
        log.error("No score columns found for combined annotation")
        raise RuntimeError("No score columns found for combined annotation")

    cluster_means = adata.obs.groupby(cluster_key)[score_cols].mean()

    # Calculate enrichment scores
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        key_added=f"rank_genes_{cluster_key}",
        use_raw=use_raw,
    )

    markers_df = sc.get.rank_genes_groups_df(
        adata, key=f"rank_genes_{cluster_key}", group=None
    )

    # Build combined evidence and make final assignments
    annotations = {}
    for cluster in adata.obs[cluster_key].cat.categories:
        # Get top cluster-specific markers
        cluster_markers = (
            markers_df[markers_df["group"] == cluster]["names"].head(n_genes).tolist()
        )

        # Calculate enrichment scores for all cell types
        enrichment_scores = {}
        for cell_type, cell in mgr.CELLS.items():
            known_markers = set(cell.markers)
            if not known_markers:
                continue

            # Find overlapping markers
            overlapping = set(cluster_markers).intersection(known_markers)
            overlap_score = (
                len(overlapping) / len(known_markers) if known_markers else 0
            )
            enrichment_scores[cell_type] = overlap_score

        # Get marker scores for all cell types
        marker_scores = {}
        for col in score_cols:
            cell_type = col.replace("_score", "")
            marker_scores[cell_type] = cluster_means.loc[cluster, col]

        # Combine evidence for each cell type
        combined_scores = {}
        for cell_type in set(
            list(marker_scores.keys()) + list(enrichment_scores.keys())
        ):
            score_value = marker_scores.get(cell_type, 0)
            enrichment_value = enrichment_scores.get(cell_type, 0)

            combined_scores[cell_type] = (
                score_weight * score_value + enrichment_weight * enrichment_value
            )

        # Find best cell type
        if combined_scores:
            best_cell_type = max(combined_scores, key=combined_scores.get)
            best_score = combined_scores[best_cell_type]

            # Apply minimum score threshold
            if best_score < min_combined_score:
                annotations[cluster] = "Unknown"
                log.debug(
                    f"Cluster {cluster}: combined score {best_score:.3f} below threshold, annotating as 'Unknown'"
                )
            else:
                annotations[cluster] = best_cell_type
                log.debug(
                    f"Cluster {cluster} -> {best_cell_type} (combined score: {best_score:.3f})"
                )
        else:
            annotations[cluster] = "Unknown"
            log.warning(f"No combined scores available for cluster {cluster}")

    log.info(f"Annotated {len(annotations)} clusters based on combined evidence")
    return annotations


# --- Main Annotation Functions ---


def score_cell_types(
    adata: AnnData,
    marker_config: Union[str, Manager],
    layer: Optional[str] = "log1p_norm",
    use_raw: bool = False,
    min_genes: int = 3,
    ctrl_size: int = 50,
    score_name_suffix: str = "_score",
    copy: bool = False,
) -> AnnData:
    """
    Score cells for multiple cell types using `sc.tl.score_genes`.

    This function calculates enrichment scores for each cell type defined in the
    marker configuration by comparing the expression of marker genes to a set of
    randomly selected control genes with similar expression levels.

    Args:
        adata: AnnData object containing gene expression data
        marker_config: A Manager instance or path to a marker configuration file
        layer: Layer to use for scoring if `use_raw` is False
        use_raw: If True, use `adata.raw` for scoring (recommended after HVG selection)
        min_genes: Minimum number of marker genes required to compute a score
        ctrl_size: Number of control genes to use for score calculation
        score_name_suffix: Suffix to add to cell type names for score column names
        copy: Whether to return a copy of the AnnData object

    Returns:
        AnnData object with score columns added to `adata.obs`

    Raises:
        ValueError: If adata.raw is not set when use_raw=True
        TypeError: If marker_config is not a string or Manager instance
    """
    log.info(
        f"Scoring cells for cell types using {'raw data' if use_raw else f'layer: {layer}'}"
    )

    if copy:
        adata = adata.copy()

    # Load marker configuration
    if isinstance(marker_config, str):
        log.info(f"Loading marker configuration from '{marker_config}'")
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        log.info("Using provided Manager instance")
        mgr = marker_config
    else:
        error_msg = "marker_config must be a file path (str) or a Manager instance"
        log.error(error_msg)
        raise TypeError(error_msg)

    # Check for raw data if needed
    if use_raw:
        if adata.raw is None:
            error_msg = (
                "adata.raw is not set. Please set it before using use_raw=True, "
                "e.g., after normalization: `adata.raw = adata`."
            )
            log.error(error_msg)
            raise ValueError(error_msg)

        log.info(f"Using raw data for scoring (layer '{layer}' will be ignored)")

        # Intersect markers with genes in raw data
        before_count = sum(len(cell.markers) for cell in mgr.CELLS.values())
        mgr.intersect_with(adata.raw)
        after_count = sum(len(cell.markers) for cell in mgr.CELLS.values())

        log.info(
            f"After intersection with raw data: {after_count}/{before_count} markers remain"
        )

        # Score each cell type
        cell_types_scored = 0
        cell_types_skipped = 0

        for cell_type, cell in mgr.CELLS.items():
            markers = cell.markers
            if len(markers) >= min_genes:
                try:
                    log.debug(f"Scoring for '{cell_type}' using {len(markers)} markers")
                    score_name = f"{cell_type}{score_name_suffix}"
                    sc.tl.score_genes(
                        adata,
                        markers,
                        score_name=score_name,
                        use_raw=True,
                        ctrl_size=ctrl_size,
                    )
                    cell_types_scored += 1
                except Exception as e:
                    log.error(f"Error scoring cell type '{cell_type}': {str(e)}")
                    cell_types_skipped += 1
            else:
                log.warning(
                    f"Skipping '{cell_type}': only {len(markers)} markers found in raw data (min: {min_genes})"
                )
                cell_types_skipped += 1

    else:
        # Score using the specified layer
        log.info(f"Using layer '{layer}' for scoring")

        # Intersect markers with genes in data
        before_count = sum(len(cell.markers) for cell in mgr.CELLS.values())
        mgr.intersect_with(adata)
        after_count = sum(len(cell.markers) for cell in mgr.CELLS.values())

        log.info(
            f"After intersection with data: {after_count}/{before_count} markers remain"
        )

        # Score each cell type
        cell_types_scored = 0
        cell_types_skipped = 0

        with use_layer_as_X(adata, layer):
            for cell_type, cell in mgr.CELLS.items():
                markers = cell.markers
                if len(markers) >= min_genes:
                    try:
                        log.debug(
                            f"Scoring for '{cell_type}' using {len(markers)} markers"
                        )
                        score_name = f"{cell_type}{score_name_suffix}"
                        sc.tl.score_genes(
                            adata, markers, score_name=score_name, ctrl_size=ctrl_size
                        )
                        cell_types_scored += 1
                    except Exception as e:
                        log.error(f"Error scoring cell type '{cell_type}': {str(e)}")
                        cell_types_skipped += 1
                else:
                    log.warning(
                        f"Skipping '{cell_type}': only {len(markers)} markers found in data (min: {min_genes})"
                    )
                    cell_types_skipped += 1

    log.info(
        f"Completed scoring: {cell_types_scored} cell types scored, {cell_types_skipped} skipped"
    )

    # Store scoring metadata
    adata.uns["cell_type_scoring"] = {
        "method": "score_genes",
        "use_raw": use_raw,
        "layer": layer if not use_raw else None,
        "min_genes": min_genes,
        "ctrl_size": ctrl_size,
        "cell_types_scored": cell_types_scored,
        "cell_types_skipped": cell_types_skipped,
        "score_name_suffix": score_name_suffix,
    }

    return adata


def annotate_clusters(
    adata: AnnData,
    cluster_key: str,
    marker_config: Union[str, Manager],
    method: Literal["max_score", "enrichment", "combined"] = "max_score",
    use_raw: bool = False,
    key_added: Optional[str] = None,
    min_score: float = 0.1,
    n_genes: int = 100,
    score_weight: float = 0.6,
    enrichment_weight: float = 0.4,
    plot: bool = False,
    copy: bool = False,
) -> AnnData:
    """
    Annotate clusters with cell type labels using marker gene expression.

    This function assigns cell type labels to clusters based on various methods:
    - 'max_score': Uses pre-computed cell type scores (from score_cell_types)
    - 'enrichment': Uses overlap between cluster markers and cell type markers
    - 'combined': Uses a weighted combination of both approaches

    Args:
        adata: AnnData object with clustering results
        cluster_key: Key in adata.obs for cluster assignments
        marker_config: A Manager instance or path to a marker configuration file
        method: Annotation method to use
        use_raw: Whether to use raw data for enrichment analysis
        key_added: Key in adata.obs to store annotations (default: f"{cluster_key}_annotated")
        min_score: Minimum score required to assign a cell type
        n_genes: Number of top genes to consider for enrichment analysis
        score_weight: Weight for score-based evidence in combined method (0-1)
        enrichment_weight: Weight for enrichment-based evidence in combined method (0-1)
        plot: Whether to create a visualization of the annotation results
        copy: Whether to return a copy of the AnnData object

    Returns:
        AnnData object with cluster annotations in adata.obs[key_added]

    Raises:
        ValueError: If method is unknown or if adata.raw is required but not set
        TypeError: If marker_config is not a string or Manager instance
    """
    log.info(f"Annotating clusters in '{cluster_key}' using '{method}' method")

    if copy:
        adata = adata.copy()

    # Set default output key
    if key_added is None:
        key_added = f"{cluster_key}_annotated"
        log.info(f"Output will be stored in adata.obs['{key_added}']")

    # Check if cluster key exists
    if cluster_key not in adata.obs:
        error_msg = f"Cluster key '{cluster_key}' not found in adata.obs"
        log.error(error_msg)
        raise ValueError(error_msg)

    # Load marker configuration
    if isinstance(marker_config, str):
        log.info(f"Loading marker configuration from '{marker_config}'")
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        log.info("Using provided Manager instance")
        mgr = marker_config
    else:
        error_msg = "marker_config must be a file path (str) or a Manager instance"
        log.error(error_msg)
        raise TypeError(error_msg)

    # Validate raw data if needed
    if use_raw:
        if adata.raw is None:
            error_msg = "adata.raw is not set but use_raw=True"
            log.error(error_msg)
            raise ValueError(error_msg)

        # Intersect markers with genes in raw data
        log.info("Intersecting markers with genes in raw data")
        mgr.intersect_with(adata.raw)
    else:
        # Intersect markers with genes in data
        log.info("Intersecting markers with genes in data")
        mgr.intersect_with(adata)

    # Perform annotation based on selected method
    if method == "max_score":
        log.info("Using maximum score method for annotation")
        mapping = _annotate_by_max_score(adata, cluster_key, mgr, min_score=min_score)

    elif method == "enrichment":
        log.info(
            f"Using enrichment method for annotation (use_raw={use_raw}, n_genes={n_genes})"
        )
        mapping = _annotate_by_enrichment(
            adata,
            cluster_key,
            mgr,
            use_raw=use_raw,
            n_genes=n_genes,
            min_overlap_score=min_score,
        )

    elif method == "combined":
        log.info(
            f"Using combined method for annotation (score_weight={score_weight}, enrichment_weight={enrichment_weight})"
        )
        mapping = _annotate_by_combined_evidence(
            adata,
            cluster_key,
            mgr,
            score_weight=score_weight,
            enrichment_weight=enrichment_weight,
            min_combined_score=min_score,
            use_raw=use_raw,
            n_genes=n_genes,
        )

    else:
        error_msg = f"Unknown annotation method: {method}"
        log.error(error_msg)
        raise ValueError(error_msg)

    # Apply the mapping to create annotations
    adata.obs[key_added] = adata.obs[cluster_key].map(mapping).astype("category")

    # Check for clusters without annotations
    n_unknown = np.sum(adata.obs[key_added] == "Unknown")
    if n_unknown > 0:
        log.warning(
            f"{n_unknown} clusters could not be confidently annotated and were labeled as 'Unknown'"
        )

    # Add colors for plotting
    colors = []
    for cat in adata.obs[key_added].cat.categories:
        if cat == "Unknown":
            colors.append("#CCCCCC")  # Gray for unknown
        elif cat in mgr.CELLS and mgr[cat].color:
            colors.append(mgr[cat].color)
        else:
            colors.append("#000000")  # Black as fallback

    adata.uns[f"{key_added}_colors"] = colors

    # Store annotation metadata
    adata.uns[f"{key_added}_params"] = {
        "method": method,
        "use_raw": use_raw,
        "min_score": min_score,
        "n_genes": n_genes,
        "score_weight": score_weight if method == "combined" else None,
        "enrichment_weight": enrichment_weight if method == "combined" else None,
        "mapping": {
            str(k): v for k, v in mapping.items()
        },  # Convert keys to strings for JSON compatibility
        "source_clusters": cluster_key,
    }

    # Create visualization if requested
    if plot:
        try:
            log.info("Creating annotation visualization")

            # Ensure UMAP is available
            if "X_umap" not in adata.obsm:
                log.warning("UMAP representation not found, computing it now")
                sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=15)
                sc.tl.umap(adata)

            # Plot clusters and annotations side by side
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))

            sc.pl.umap(
                adata,
                color=cluster_key,
                ax=ax1,
                show=False,
                title=f"Clusters ({cluster_key})",
            )
            sc.pl.umap(
                adata,
                color=key_added,
                ax=ax2,
                show=False,
                title=f"Cell Types ({key_added})",
            )

            plt.tight_layout()
            plt.show()

            # Create a more detailed figure showing cluster composition
            composition = (
                pd.crosstab(
                    adata.obs[cluster_key], adata.obs[key_added], normalize="index"
                )
                * 100
            )

            plt.figure(figsize=(12, 8))
            import seaborn as sns

            sns.heatmap(composition, annot=True, cmap="YlGnBu", fmt=".1f")
            plt.xlabel("Cell Type")
            plt.ylabel("Cluster")
            plt.title("Cluster Composition (% of cells)")
            plt.tight_layout()
            plt.show()

        except Exception as e:
            log.warning(f"Error creating annotation visualization: {str(e)}")

    log.info(f"Annotation complete. Results stored in adata.obs['{key_added}']")
    return adata


def transfer_labels(
    adata: AnnData,
    ref_adata: AnnData,
    ref_label_key: str,
    n_neighbors: int = 30,
    use_rep: str = "X_pca",
    key_added: Optional[str] = "predicted_labels",
    normalize_weights: bool = True,
    confidence_threshold: float = 0.7,
    copy: bool = False,
) -> AnnData:
    """
    Transfer cell type labels from a reference dataset to a query dataset.

    This function implements a k-nearest neighbor approach to label transfer,
    assigning cell types based on the most common labels among nearest neighbors
    in the reference dataset.

    Args:
        adata: Query AnnData object to be labeled
        ref_adata: Reference AnnData object with existing labels
        ref_label_key: Key in ref_adata.obs containing cell type labels
        n_neighbors: Number of neighbors to consider for label transfer
        use_rep: Representation to use for neighbor search (must exist in both datasets)
        key_added: Key in adata.obs to store transferred labels
        normalize_weights: Whether to normalize neighbor weights by distance
        confidence_threshold: Minimum confidence score required to assign a label
        copy: Whether to return a copy of the AnnData object

    Returns:
        AnnData object with transferred labels in adata.obs[key_added]

    Notes:
        Both datasets must have the same dimensionality reduction representation
        (e.g., PCA, UMAP) specified by use_rep.
    """
    from sklearn.neighbors import NearestNeighbors

    log.info(
        f"Transferring labels from reference dataset using {use_rep} representation"
    )

    if copy:
        adata = adata.copy()

    # Check for required representations
    if use_rep not in adata.obsm:
        error_msg = f"Representation '{use_rep}' not found in query dataset"
        log.error(error_msg)
        raise ValueError(error_msg)

    if use_rep not in ref_adata.obsm:
        error_msg = f"Representation '{use_rep}' not found in reference dataset"
        log.error(error_msg)
        raise ValueError(error_msg)

    # Check for label key in reference
    if ref_label_key not in ref_adata.obs:
        error_msg = f"Label key '{ref_label_key}' not found in reference dataset"
        log.error(error_msg)
        raise ValueError(error_msg)

    # Check compatibility of embeddings
    if adata.obsm[use_rep].shape[1] != ref_adata.obsm[use_rep].shape[1]:
        error_msg = (
            f"Dimensionality mismatch in {use_rep}: "
            f"query has {adata.obsm[use_rep].shape[1]} dimensions, "
            f"reference has {ref_adata.obsm[use_rep].shape[1]} dimensions"
        )
        log.error(error_msg)
        raise ValueError(error_msg)

    # Extract embeddings and labels
    query_embed = adata.obsm[use_rep]
    ref_embed = ref_adata.obsm[use_rep]
    ref_labels = ref_adata.obs[ref_label_key].values

    # Build nearest neighbors model
    log.info(f"Building nearest neighbors model with k={n_neighbors}")
    nn = NearestNeighbors(n_neighbors=n_neighbors)
    nn.fit(ref_embed)

    # Find nearest neighbors
    distances, indices = nn.kneighbors(query_embed)

    # Transfer labels
    transferred_labels = []
    confidence_scores = []

    log.info("Transferring labels based on nearest neighbors")
    for i in range(len(adata)):
        # Get labels of nearest neighbors
        neighbor_labels = ref_labels[indices[i]]

        if normalize_weights:
            # Calculate weights based on distances (closer = higher weight)
            # Add small epsilon to avoid division by zero
            epsilon = 1e-10
            weights = 1.0 / (distances[i] + epsilon)

            # Count weighted votes for each label
            unique_labels = np.unique(neighbor_labels)
            label_votes = {}

            for label in unique_labels:
                mask = neighbor_labels == label
                label_votes[label] = np.sum(weights[mask])

            # Normalize to get probabilities
            total_votes = sum(label_votes.values())
            for label in label_votes:
                label_votes[label] /= total_votes

            # Get most probable label and its confidence
            best_label = max(label_votes, key=label_votes.get)
            confidence = label_votes[best_label]

        else:
            # Simple majority voting
            from collections import Counter

            label_counts = Counter(neighbor_labels)
            best_label = label_counts.most_common(1)[0][0]
            confidence = label_counts[best_label] / n_neighbors

        # Apply confidence threshold
        if confidence >= confidence_threshold:
            transferred_labels.append(best_label)
        else:
            transferred_labels.append("Unknown")

        confidence_scores.append(confidence)

    # Store results
    adata.obs[key_added] = pd.Categorical(transferred_labels)
    adata.obs[f"{key_added}_confidence"] = confidence_scores

    # Calculate overall statistics
    n_assigned = np.sum(np.array(transferred_labels) != "Unknown")
    assignment_rate = n_assigned / len(adata) * 100

    log.info(
        f"Label transfer complete: {n_assigned}/{len(adata)} cells ({assignment_rate:.1f}%) assigned labels"
    )
    log.info(
        f"Results stored in adata.obs['{key_added}'] with confidence scores in adata.obs['{key_added}_confidence']"
    )

    # Store metadata
    adata.uns[f"{key_added}_params"] = {
        "method": "knn_transfer",
        "n_neighbors": n_neighbors,
        "use_rep": use_rep,
        "normalize_weights": normalize_weights,
        "confidence_threshold": confidence_threshold,
        "reference_label_key": ref_label_key,
        "assignment_rate": float(
            assignment_rate
        ),  # Convert to float for JSON compatibility
    }

    return adata


def evaluate_annotation(
    adata: AnnData,
    cluster_key: str,
    annotation_key: str,
    marker_config: Union[str, Manager],
    plot: bool = True,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Evaluate the quality of cell type annotations by comparing to marker gene expression.

    This function calculates metrics to assess how well the annotations match
    the expected marker gene expression patterns for each cell type.

    Args:
        adata: AnnData object with clustering and annotation results
        cluster_key: Key in adata.obs for cluster assignments
        annotation_key: Key in adata.obs for cell type annotations
        marker_config: A Manager instance or path to a marker configuration file
        plot: Whether to create evaluation visualizations
        save_path: Path to save evaluation results and plots

    Returns:
        DataFrame with evaluation metrics for each cluster

    Metrics:
        - marker_coverage: Percentage of expected markers detected in the cluster
        - marker_specificity: How specific the marker expression is to the assigned cell type
        - annotation_confidence: Overall confidence score for the annotation
    """
    log.info(
        f"Evaluating annotation quality for '{annotation_key}' based on '{cluster_key}'"
    )

    # Load marker configuration
    if isinstance(marker_config, str):
        log.info(f"Loading marker configuration from '{marker_config}'")
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        log.info("Using provided Manager instance")
        mgr = marker_config
    else:
        error_msg = "marker_config must be a file path (str) or a Manager instance"
        log.error(error_msg)
        raise TypeError(error_msg)

    # Intersect markers with data
    mgr.intersect_with(adata)

    # Get marker expression in each cluster
    cluster_markers = {}

    # Calculate differential expression for clusters
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        key_added=f"rank_genes_{cluster_key}",
    )

    # Extract differentially expressed genes for each cluster
    for cluster in adata.obs[cluster_key].cat.categories:
        de_genes = sc.get.rank_genes_groups_df(
            adata, key=f"rank_genes_{cluster_key}", group=cluster
        )
        # Filter for significant genes
        sig_genes = de_genes[de_genes["pvals_adj"] < 0.05]["names"].tolist()
        cluster_markers[cluster] = sig_genes

    # Evaluate each cluster's annotation
    results = []

    for cluster in adata.obs[cluster_key].cat.categories:
        # Get assigned cell type for this cluster
        cells_in_cluster = adata.obs[cluster_key] == cluster
        cell_types = adata.obs.loc[cells_in_cluster, annotation_key]

        # If the cluster has multiple annotations, use the most common one
        assigned_type = cell_types.value_counts().index[0]

        # Skip evaluation if annotated as "Unknown"
        if assigned_type == "Unknown":
            log.debug(f"Skipping evaluation for cluster {cluster} (labeled as Unknown)")
            continue

        # Get marker genes for the assigned cell type
        if assigned_type in mgr.CELLS:
            expected_markers = mgr[assigned_type].markers

            # Calculate marker coverage
            expr_markers = set(cluster_markers[cluster])
            found_markers = expr_markers.intersection(expected_markers)

            marker_coverage = (
                len(found_markers) / len(expected_markers) if expected_markers else 0
            )

            # Calculate marker specificity (how specific the markers are to this cell type)
            all_other_markers = []
            for cell_type, cell in mgr.CELLS.items():
                if cell_type != assigned_type:
                    all_other_markers.extend(cell.markers)

            # Count markers shared with other cell types
            all_other_markers = set(all_other_markers)
            shared_markers = found_markers.intersection(all_other_markers)

            specificity = (
                1.0 - (len(shared_markers) / len(found_markers)) if found_markers else 0
            )

            # Calculate overall confidence
            confidence = (marker_coverage * 0.6) + (specificity * 0.4)

            results.append(
                {
                    "cluster": cluster,
                    "cell_type": assigned_type,
                    "marker_coverage": marker_coverage,
                    "marker_specificity": specificity,
                    "annotation_confidence": confidence,
                    "found_markers": ", ".join(found_markers),
                    "expected_markers": len(expected_markers),
                    "detected_markers": len(found_markers),
                }
            )
        else:
            log.warning(
                f"Cell type '{assigned_type}' not found in marker configuration"
            )

    # Create results DataFrame
    results_df = pd.DataFrame(results)

    if not results_df.empty:
        # Calculate overall statistics
        avg_coverage = results_df["marker_coverage"].mean()
        avg_specificity = results_df["marker_specificity"].mean()
        avg_confidence = results_df["annotation_confidence"].mean()

        log.info(
            f"Evaluation complete: average confidence score = {avg_confidence:.3f}"
        )
        log.info(
            f"Average marker coverage: {avg_coverage:.3f}, average specificity: {avg_specificity:.3f}"
        )

        # Create visualizations if requested
        if plot:
            try:
                # Sort by confidence
                plot_df = results_df.sort_values("annotation_confidence")

                plt.figure(figsize=(12, 8))
                plt.barh(
                    plot_df["cluster"].astype(str) + " (" + plot_df["cell_type"] + ")",
                    plot_df["annotation_confidence"],
                    color="skyblue",
                )
                plt.xlabel("Annotation Confidence Score")
                plt.ylabel("Cluster (Cell Type)")
                plt.title("Annotation Confidence by Cluster")
                plt.xlim(0, 1.0)
                plt.tight_layout()

                if save_path:
                    os.makedirs(
                        os.path.dirname(os.path.abspath(save_path)), exist_ok=True
                    )
                    plt.savefig(f"{save_path}_confidence.png", dpi=300)

                plt.show()

                # Create a heatmap of metrics
                plt.figure(figsize=(10, 8))
                heatmap_data = plot_df[
                    ["marker_coverage", "marker_specificity", "annotation_confidence"]
                ]

                import seaborn as sns

                sns.heatmap(
                    heatmap_data.set_index(
                        plot_df["cluster"].astype(str)
                        + " ("
                        + plot_df["cell_type"]
                        + ")"
                    ),
                    annot=True,
                    cmap="YlGnBu",
                    vmin=0,
                    vmax=1,
                )
                plt.title("Annotation Quality Metrics by Cluster")
                plt.tight_layout()

                if save_path:
                    plt.savefig(f"{save_path}_metrics.png", dpi=300)

                plt.show()

            except Exception as e:
                log.warning(f"Error creating evaluation visualizations: {str(e)}")

        # Save results if requested
        if save_path:
            try:
                results_df.to_csv(f"{save_path}_evaluation.csv", index=False)
                log.info(f"Saved evaluation results to {save_path}_evaluation.csv")
            except Exception as e:
                log.warning(f"Error saving evaluation results: {str(e)}")
    else:
        log.warning("No clusters could be evaluated (possibly all labeled as Unknown)")

    return results_df
