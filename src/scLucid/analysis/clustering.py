"""
Clustering functions for single-cell RNA-seq data.

This module provides comprehensive clustering solutions for single-cell RNA-seq data,
including marker-guided clustering, optimal parameter selection, and cluster merging.
It supports various clustering algorithms (Leiden, Louvain, K-means, HDBSCAN) with
automatic parameter optimization.
"""

import gc
import logging
import os
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy.stats import entropy
from sklearn import metrics
from dataclasses import dataclass
from ..utils.marker_manager import Manager

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = [
    "ClusteringConfig",
    "ResolutionSearchConfig",
    "find_resolution",
    "cluster_cells",
    "merge_clusters",
]

@dataclass
class ClusteringConfig:
    """Configuration for a single clustering run."""
    method: Literal["leiden", "louvain"] = "leiden"
    resolution: float = 1.0
    use_rep: str = "X_pca"
    key_added: Optional[str] = None
    random_state: int = 42
    
    def __post_init__(self):
        if self.key_added is None:
            self.key_added = f"{self.method}_res{self.resolution}"

@dataclass
class ResolutionSearchConfig:
    """Configuration for optimizing clustering resolution."""
    resolution_range: Tuple[float, float, int] = (0.1, 2.0, 10)
    metric: Literal["marker_separation", "silhouette"] = "silhouette"
    marker_config: Optional[Union[str, Manager]] = None
    use_raw_for_markers: bool = False


# --- Helper Functions ---
def _evaluate_marker_separation(
    adata: sc.AnnData,
    cluster_key: str,
    marker_genes: Dict[str, List[str]],
    use_raw: bool = False,
) -> float:
    """
    Calculates how well clusters separate known marker gene sets.

    This metric measures how specifically marker genes are expressed in clusters.
    Higher scores indicate better separation of cell types based on known markers.

    Args:
        adata: AnnData object with clustering results
        cluster_key: Key in adata.obs containing cluster assignments
        marker_genes: Dictionary mapping cell types to marker gene lists
        use_raw: Whether to use adata.raw for expression data

    Returns:
        Separation score (0-1, higher is better)
    """
    # Determine which expression data to use
    if use_raw:
        if adata.raw is None:
            log.error("adata.raw is not set but use_raw=True")
            raise ValueError("adata.raw must be set to use `use_raw=True`.")
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        # Fallback to layer or .X
        if "log1p_norm" in adata.layers:
            X = adata.layers["log1p_norm"]
            log.info("Using 'log1p_norm' layer for marker separation score")
        else:
            log.warning(
                "'log1p_norm' layer not found. Using adata.X for marker separation score."
            )
            X = adata.X
        var_names = adata.var_names

    clusters = adata.obs[cluster_key].cat.categories
    n_clusters = len(clusters)
    if n_clusters <= 1:
        log.warning("Only one cluster found. Marker separation score will be 0.")
        return 0.0

    scores = []
    for cell_type, markers in marker_genes.items():
        marker_indices = [i for i, gene in enumerate(var_names) if gene in markers]
        if not marker_indices:
            log.debug(f"No markers for '{cell_type}' found in data. Skipping.")
            continue

        # Mean expression of marker set in each cluster
        try:
            cluster_means = np.array(
                [
                    X[(adata.obs[cluster_key] == c).values][:, marker_indices].mean()
                    for c in clusters
                ]
            )

            if np.sum(cluster_means) == 0:
                log.debug(f"Zero expression for '{cell_type}' markers. Skipping.")
                continue

            # Normalize to a probability distribution
            cluster_probs = cluster_means / np.sum(cluster_means)

            # Calculate score as 1 - normalized_entropy
            marker_entropy = entropy(cluster_probs)
            max_entropy = np.log(n_clusters)
            score = 1 - (marker_entropy / max_entropy) if max_entropy > 0 else 0.0
            scores.append(score)
            log.debug(f"Cell type '{cell_type}': separation score = {score:.4f}")
        except Exception as e:
            log.warning(
                f"Error calculating marker separation for '{cell_type}': {str(e)}"
            )

    if not scores:
        log.warning("No valid marker sets found for separation score calculation")
        return 0.0

    final_score = np.mean(scores)
    log.info(
        f"Overall marker separation score: {final_score:.4f} (from {len(scores)} cell types)"
    )
    return final_score


def _evaluate_silhouette(
    adata: sc.AnnData, cluster_key: str, use_rep: str = "X_pca"
) -> float:
    """
    Calculates silhouette score for clustering quality assessment.

    Args:
        adata: AnnData object with clustering results
        cluster_key: Key in adata.obs containing cluster assignments
        use_rep: Representation to use for distance calculation

    Returns:
        Silhouette score (-1 to 1, higher is better)
    """
    if use_rep not in adata.obsm:
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")

    try:
        labels = adata.obs[cluster_key].cat.codes
        unique_labels = np.unique(labels)

        if len(unique_labels) <= 1:
            log.warning("Only one cluster found. Silhouette score will be 0.")
            return 0.0

        # Use a sample for very large datasets
        if adata.n_obs > 20000:
            sample_size = min(20000, adata.n_obs)
            log.info(
                f"Using {sample_size} cells (subsample) for silhouette calculation"
            )
            score = metrics.silhouette_score(
                adata.obsm[use_rep], labels, sample_size=sample_size, random_state=42
            )
        else:
            score = metrics.silhouette_score(adata.obsm[use_rep], labels)

        log.info(f"Silhouette score: {score:.4f}")
        return score
    except Exception as e:
        log.error(f"Error calculating silhouette score: {str(e)}")
        return 0.0


# --- Main Functions ---
def find_resolution(
    adata: sc.AnnData,
    resolution_range: Tuple[float, float, int] = (0.1, 2.0, 10),
    clustering_method: Literal["leiden", "louvain"] = "leiden",
    metric: Literal["marker_separation", "silhouette"] = "marker_separation",
    marker_config: Optional[Union[str, Manager]] = None,
    use_rep: str = "X_pca",
    neighbors_key: Optional[str] = None,
    plot: bool = True,
    use_raw: bool = False,
    random_state: int = 42,
    copy: bool = False,
) -> sc.AnnData:
    """
    Performs clustering across a range of resolutions and finds the optimal one.

    This function evaluates clustering quality using either marker gene separation
    or silhouette score. If a marker_config is provided, it defaults to using marker
    gene separation. Otherwise, it automatically falls back to silhouette score.

    Args:
        adata: AnnData object
        resolution_range: Tuple of (start, end, steps) for resolution search
        clustering_method: 'leiden' or 'louvain' algorithm
        metric: Metric for evaluation ('marker_separation' or 'silhouette')
        marker_config: A Manager instance or path to a marker TOML file
        use_rep: Representation to use for silhouette score (e.g., 'X_pca')
        neighbors_key: Key in `adata.uns` for the neighbors graph
        plot: Whether to plot the evaluation metrics vs. resolution
        use_raw: Whether to use `adata.raw` for marker_separation metric
        random_state: Seed for clustering reproducibility
        copy: If True, return a copy of the AnnData object

    Returns:
        AnnData object with the optimal clustering stored in `adata.obs`

    Examples:
        >>> # Using silhouette score for optimization
        >>> adata = find_resolution(
        ...     adata,
        ...     resolution_range=(0.1, 2.0, 10),
        ...     metric="silhouette"
        ... )
        >>>
        >>> # Using marker separation with a cell type marker manager
        >>> from scRNA.analysis import get_marker_manager
        >>> marker_mgr = get_marker_manager(species="human", tissue="pbmc")
        >>> adata = find_resolution(
        ...     adata,
        ...     metric="marker_separation",
        ...     marker_config=marker_mgr
        ... )
    """
    log.info(f"Starting resolution optimization for {clustering_method} clustering")

    # Handle copy if requested
    if copy:
        adata = adata.copy()

    # Validate parameters
    if use_rep not in adata.obsm and metric == "silhouette":
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")

    if "neighbors" not in adata.uns and neighbors_key is None:
        log.info(f"Computing neighbors using '{use_rep}'")
        sc.pp.neighbors(adata, use_rep=use_rep, random_state=random_state)

    # Determine which metric to use
    effective_metric = metric
    if marker_config is None and metric == "marker_separation":
        effective_metric = "silhouette"
        log.warning(
            f"No marker_config provided. Falling back to '{effective_metric}' metric."
        )

    if effective_metric == "marker_separation" and marker_config is None:
        log.error("marker_config must be provided for marker_separation metric")
        raise ValueError(
            "`marker_config` must be provided when explicitly using 'marker_separation' metric."
        )

    # Process marker configuration if provided
    if isinstance(marker_config, Manager):
        log.info("Using provided Manager instance for marker genes")
        mgr = marker_config
    elif isinstance(marker_config, str):
        log.info(f"Loading marker configuration from '{marker_config}'")
        mgr = Manager(marker_config)
    else:
        mgr = None

    if mgr:
        log.info("Intersecting marker genes with dataset")
        if use_raw and adata.raw is not None:
            # Use raw data for intersection if specified
            mgr.intersect_with(adata.raw)
        else:
            mgr.intersect_with(adata)

    # Extract marker genes from Manager
    marker_genes = (
        {cell.name: cell.markers for cell in mgr.CELLS.values()} if mgr else {}
    )
    if mgr:
        log.info(
            f"Using {len(marker_genes)} cell types with marker genes for evaluation"
        )

    # Generate resolution range
    start, end, steps = resolution_range
    resolutions = np.linspace(start, end, steps)
    log.info(f"Testing {steps} resolutions from {start} to {end}")

    eval_results = []

    # Evaluate each resolution
    log.info(f"Searching for optimal resolution using '{effective_metric}' metric")

    for res in resolutions:
        key = f"{clustering_method}_res_{res:.2f}"

        # Perform clustering
        try:
            if clustering_method == "leiden":
                sc.tl.leiden(
                    adata,
                    resolution=res,
                    key_added=key,
                    neighbors_key=neighbors_key,
                    random_state=random_state,
                )
            else:
                sc.tl.louvain(
                    adata,
                    resolution=res,
                    key_added=key,
                    neighbors_key=neighbors_key,
                    random_state=random_state,
                )
        except Exception as e:
            log.error(
                f"Error in {clustering_method} clustering with resolution {res}: {str(e)}"
            )
            continue

        # Evaluate clustering quality
        if effective_metric == "marker_separation":
            score = _evaluate_marker_separation(
                adata, key, marker_genes, use_raw=use_raw
            )
        else:  # Silhouette
            score = _evaluate_silhouette(adata, key, use_rep)

        n_clusters = adata.obs[key].nunique()
        eval_results.append(
            {"resolution": res, "n_clusters": n_clusters, "score": score}
        )
        log.info(
            f"Resolution {res:.2f}: {n_clusters} clusters, {effective_metric} score = {score:.4f}"
        )

    # Process evaluation results
    eval_df = pd.DataFrame(eval_results).dropna()
    if eval_df.empty:
        log.error("Clustering evaluation failed for all resolutions")
        raise RuntimeError("Clustering evaluation failed for all resolutions.")

    # Find optimal resolution
    optimal_idx = eval_df["score"].idxmax()
    optimal_res = eval_df.loc[optimal_idx, "resolution"]
    optimal_clusters = eval_df.loc[optimal_idx, "n_clusters"]
    optimal_score = eval_df.loc[optimal_idx, "score"]

    log.info(
        f"Optimal resolution: {optimal_res:.2f} with {optimal_clusters} clusters "
        f"and {effective_metric} score of {optimal_score:.4f}"
    )

    # Create final clustering with optimal resolution
    key_added = f"{clustering_method}_optimal"
    try:
        if clustering_method == "leiden":
            sc.tl.leiden(
                adata,
                resolution=optimal_res,
                key_added=key_added,
                neighbors_key=neighbors_key,
                random_state=random_state,
            )
        else:
            sc.tl.louvain(
                adata,
                resolution=optimal_res,
                key_added=key_added,
                neighbors_key=neighbors_key,
                random_state=random_state,
            )
    except Exception as e:
        log.error(f"Error in final {clustering_method} clustering: {str(e)}")
        raise RuntimeError(f"Final clustering failed: {str(e)}")

    # --- Store results in unified namespace ---
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('clustering', {})
    adata.uns['scrnatk']['analysis'][f"{key_added}_evaluation"] = {
        "metric": effective_metric,
        "optimal_resolution": float(optimal_res),
        "optimal_score": float(optimal_score),
        "resolutions": eval_df["resolution"].tolist(),
        "scores": eval_df["score"].tolist(),
        "n_clusters": eval_df["n_clusters"].tolist(),
        "parameters": {
            "method": clustering_method,
            "use_rep": use_rep,
            "use_raw": use_raw,
            "random_state": random_state,
        },
    }

    log.info(f"Final clustering stored in `adata.obs['{key_added}']`")

    # Plot evaluation results
    if plot:
        try:
            log.info("Creating resolution optimization plot")
            fig, ax1 = plt.subplots(figsize=(10, 6))
            ax2 = ax1.twinx()

            # Plot score
            ax1.plot(
                eval_df["resolution"],
                eval_df["score"],
                "o-",
                color="b",
                label=f"{effective_metric.replace('_', ' ').title()} Score",
            )

            # Plot cluster counts
            ax2.plot(
                eval_df["resolution"],
                eval_df["n_clusters"],
                "o-",
                color="g",
                label="# Clusters",
            )

            # Mark optimal resolution
            ax1.axvline(
                x=optimal_res,
                color="r",
                linestyle="--",
                label=f"Optimal Res: {optimal_res:.2f}",
            )

            # Add labels and legend
            ax1.set_xlabel("Resolution Parameter", fontsize=12)
            ax1.set_ylabel(
                f"{effective_metric.replace('_', ' ').title()} Score",
                color="b",
                fontsize=12,
            )
            ax2.set_ylabel("Number of Clusters", color="g", fontsize=12)

            # Create combined legend
            lines1, labels1 = ax1.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax1.legend(lines1 + lines2, labels1 + labels2, loc="best", frameon=True)

            plt.title(
                f"Clustering Resolution Optimization ({clustering_method.capitalize()})",
                fontsize=14,
            )
            plt.tight_layout()
            plt.show()

        except Exception as e:
            log.warning(f"Error creating optimization plot: {str(e)}")

    return adata



def cluster_cells(
    adata: sc.AnnData,
    config: Optional[ClusteringConfig] = None,
    # Backward compatibility params
    method: Literal["leiden", "louvain", "kmeans", "hdbscan"] = "leiden",
    resolution: float = 1.0,
    n_clusters: Optional[int] = None,
    use_rep: str = "X_pca",
    key_added: str = "cluster",
    random_state: int = 42,
    plot: bool = True,
    copy: bool = False,
    **kwargs,
) -> sc.AnnData:
    """
    Cluster cells using various clustering algorithms.

    This function provides a unified interface to multiple clustering algorithms,
    with appropriate parameter handling for each method.

    Args:
        adata: AnnData object
        config: ClusteringConfig object containing clustering parameters
        method: Clustering method to use:
            - 'leiden': Leiden community detection (resolution-based)
            - 'louvain': Louvain community detection (resolution-based)
            - 'kmeans': K-means clustering (requires n_clusters)
            - 'hdbscan': Density-based clustering (optional min_cluster_size)
        resolution: Resolution parameter for Leiden/Louvain (higher=more clusters)
        n_clusters: Number of clusters for K-means
        use_rep: Representation to use for clustering (e.g., 'X_pca', 'X_umap')
        key_added: Key in adata.obs to store the resulting cluster assignments
        random_state: Random seed for reproducibility
        plot: Whether to plot the clustering results
        copy: If True, return a copy of the AnnData object
        **kwargs: Additional parameters for the clustering algorithms
            - For K-means: n_init, max_iter, etc.
            - For HDBSCAN: min_cluster_size, min_samples, etc.

    Returns:
        AnnData object with clustering results in adata.obs[key_added]

    Examples:
        >>> # Leiden clustering with default parameters
        >>> adata = cluster_cells(adata, method="leiden", resolution=0.8)
        >>>
        >>> # K-means clustering with 10 clusters
        >>> adata = cluster_cells(adata, method="kmeans", n_clusters=10)
        >>>
        >>> # HDBSCAN clustering with custom parameters
        >>> adata = cluster_cells(
        ...     adata,
        ...     method="hdbscan",
        ...     min_cluster_size=50,
        ...     min_samples=10
        ... )
    """
    """Cluster cells using Leiden or Louvain community detection."""
    if config is None:
        config = ClusteringConfig(
            method=method, resolution=resolution, use_rep=use_rep, key_added=key_added, **kwargs
        )
    
    log.info(f"Clustering cells using '{config.method}' with resolution {config.resolution}")
    if "neighbors" not in adata.uns:
        log.info(f"Computing neighbors graph using '{config.use_rep}'")
        sc.pp.neighbors(adata, use_rep=config.use_rep, random_state=config.random_state)
    
    if config.method == "leiden":
        sc.tl.leiden(adata, resolution=config.resolution, key_added=config.key_added, random_state=config.random_state)
    elif config.method == "louvain":
        sc.tl.louvain(adata, resolution=config.resolution, key_added=config.key_added, random_state=config.random_state)
    
    log.info(f"Clustering cells using '{method}' algorithm")

    # Handle copy if requested
    if copy:
        adata = adata.copy()

    # Validate parameters
    if use_rep not in adata.obsm:
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")

    if method in ["leiden", "louvain"]:
        # Check resolution parameter
        if resolution <= 0:
            log.error(f"Invalid resolution: {resolution}. Must be positive.")
            raise ValueError(f"Resolution must be positive, got {resolution}")

        # Compute neighbors if needed
        if "neighbors" not in adata.uns:
            log.info(f"Computing neighbors using '{use_rep}'")
            sc.pp.neighbors(adata, use_rep=use_rep, random_state=random_state)

    elif method == "kmeans":
        if n_clusters is None:
            log.error("n_clusters parameter is required for K-means clustering")
            raise ValueError("n_clusters must be specified for k-means clustering")

        if n_clusters <= 0:
            log.error(f"Invalid n_clusters: {n_clusters}. Must be positive.")
            raise ValueError(f"n_clusters must be positive, got {n_clusters}")

    # Execute clustering based on method
    if method == "leiden":
        log.info(f"Running Leiden clustering with resolution={resolution}")
        sc.tl.leiden(
            adata, resolution=resolution, key_added=key_added, random_state=random_state
        )

    elif method == "louvain":
        log.info(f"Running Louvain clustering with resolution={resolution}")
        sc.tl.louvain(
            adata, resolution=resolution, key_added=key_added, random_state=random_state
        )

    elif method == "kmeans":
        log.info(f"Running K-means clustering with k={n_clusters}")
        try:
            from sklearn.cluster import KMeans

            # Get the embedding data
            X = adata.obsm[use_rep]

            # Set default parameters with better defaults than sklearn
            kmeans_params = {"n_init": 10, "max_iter": 300, "algorithm": "auto"}
            # Update with user-provided parameters
            kmeans_params.update(kwargs)

            # Run kmeans
            kmeans = KMeans(
                n_clusters=n_clusters, random_state=random_state, **kmeans_params
            )

            # Fit and assign clusters
            labels = kmeans.fit_predict(X)
            adata.obs[key_added] = labels.astype(str)
            adata.obs[key_added] = adata.obs[key_added].astype("category")

            # Store additional information
            adata.uns[f"{key_added}_params"] = {
                "method": "kmeans",
                "n_clusters": n_clusters,
                "random_state": random_state,
                "use_rep": use_rep,
                **kmeans_params,
            }

        except Exception as e:
            log.error(f"Error in K-means clustering: {str(e)}")
            raise RuntimeError(f"K-means clustering failed: {str(e)}")

    elif method == "hdbscan":
        log.info("Running HDBSCAN clustering")
        try:
            import hdbscan
        except ImportError:
            log.error("HDBSCAN is not installed")
            raise ImportError("Please install hdbscan: pip install hdbscan")

        try:
            # Get the embedding data
            X = adata.obsm[use_rep]

            # Set default parameters
            hdbscan_params = {
                "min_cluster_size": 50,
                "min_samples": 10,
                "metric": "euclidean",
                "cluster_selection_method": "eom",  # excess of mass
            }
            # Update with user-provided parameters
            hdbscan_params.update(kwargs)

            log.info(
                f"HDBSCAN parameters: min_cluster_size={hdbscan_params['min_cluster_size']}, "
                f"min_samples={hdbscan_params['min_samples']}"
            )

            # Run HDBSCAN
            clusterer = hdbscan.HDBSCAN(**hdbscan_params)
            labels = clusterer.fit_predict(X)

            # Store results
            adata.obs[key_added] = labels.astype(str)
            adata.obs[key_added] = adata.obs[key_added].astype("category")

            # Store additional information
            adata.uns[f"{key_added}_params"] = {
                "method": "hdbscan",
                "use_rep": use_rep,
                **hdbscan_params,
            }

            # Check for noise points
            n_noise = np.sum(labels == -1)
            if n_noise > 0:
                noise_percent = n_noise / adata.n_obs * 100
                log.warning(
                    f"HDBSCAN identified {n_noise} cells ({noise_percent:.1f}%) as noise points"
                )

        except Exception as e:
            log.error(f"Error in HDBSCAN clustering: {str(e)}")
            raise RuntimeError(f"HDBSCAN clustering failed: {str(e)}")

    else:
        available_methods = ["leiden", "louvain", "kmeans", "hdbscan"]
        log.error(f"Unknown clustering method: '{method}'")
        raise ValueError(f"Unknown method: '{method}'. Choose from {available_methods}")

    # --- Store results in unified namespace ---
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('clustering', {})
    adata.uns['scrnatk']['analysis']['clustering'][config.key_added] = {
        'config': config.__dict__,
        'n_clusters': adata.obs[config.key_added].nunique()
    }
    
    log.info(f"Found {adata.obs[config.key_added].nunique()} clusters. Result in .obs['{config.key_added}']")
    
    # Count clusters and log result
    n_clusters_found = adata.obs[key_added].nunique()
    log.info(f"Found {n_clusters_found} clusters")

    # Create visualization
    if plot:
        try:
            # Ensure UMAP exists for visualization
            if "X_umap" not in adata.obsm:
                log.info("Computing UMAP for visualization")
                sc.tl.umap(adata)

            plt.figure(figsize=(10, 8))
            sc.pl.umap(
                adata,
                color=key_added,
                legend_loc="on data",
                title=f"{method.capitalize()} Clustering",
            )
        except Exception as e:
            log.warning(f"Error creating clustering plot: {str(e)}")

    return adata


def merge_clusters(
    adata: sc.AnnData,
    cluster_key: str,
    similarity_threshold: float = 0.8,
    method: Literal["marker_overlap", "expression_correlation"] = "marker_overlap",
    key_added: Optional[str] = None,
    copy: bool = False,
) -> sc.AnnData:
    """
    Merges similar clusters based on marker overlap or expression correlation.

    This function identifies and combines clusters that are highly similar according
    to the specified similarity metric, reducing over-clustering.

    Args:
        adata: AnnData object with clustering results
        cluster_key: Key in adata.obs containing cluster assignments
        similarity_threshold: Threshold for similarity to merge clusters (0-1)
        method: Method to calculate similarity:
            - 'marker_overlap': Jaccard similarity of top marker genes
            - 'expression_correlation': Correlation of mean expression profiles
        key_added: Key in adata.obs to store merged cluster assignments
        copy: If True, return a copy of the AnnData object

    Returns:
        AnnData object with merged clusters in adata.obs[key_added]

    Examples:
        >>> # Merge clusters based on marker gene overlap
        >>> adata = merge_clusters(
        ...     adata,
        ...     cluster_key="leiden_optimal",
        ...     similarity_threshold=0.7,
        ...     method="marker_overlap"
        ... )
        >>>
        >>> # Merge clusters based on expression correlation
        >>> adata = merge_clusters(
        ...     adata,
        ...     cluster_key="leiden",
        ...     similarity_threshold=0.9,
        ...     method="expression_correlation",
        ...     key_added="leiden_merged_corr"
        ... )
    """
    log.info(f"Starting cluster merging using '{method}' method")

    # Handle copy if requested
    if copy:
        adata = adata.copy()

    # Validate parameters
    if cluster_key not in adata.obs:
        log.error(f"Cluster key '{cluster_key}' not found in adata.obs")
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")

    if not 0 <= similarity_threshold <= 1:
        log.error(
            f"Invalid similarity threshold: {similarity_threshold}. Must be between 0 and 1."
        )
        raise ValueError(
            f"Similarity threshold must be between 0 and 1, got {similarity_threshold}"
        )

    # Set output key
    if key_added is None:
        key_added = f"{cluster_key}_merged"
        log.info(f"Results will be stored in adata.obs['{key_added}']")

    # Get cluster categories
    if not hasattr(adata.obs[cluster_key], "cat"):
        log.warning(f"Converting {cluster_key} to categorical")
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    clusters = adata.obs[cluster_key].cat.categories
    n_clusters = len(clusters)

    log.info(
        f"Analyzing {n_clusters} clusters with similarity threshold {similarity_threshold}"
    )

    # Calculate similarity matrix based on method
    if method == "marker_overlap":
        log.info("Computing marker genes for each cluster")

        # Identify marker genes for each cluster
        try:
            sc.tl.rank_genes_groups(
                adata,
                groupby=cluster_key,
                method="wilcoxon",
                key_added=f"rank_genes_{cluster_key}",
            )
            markers_df = sc.get.rank_genes_groups_df(
                adata, key=f"rank_genes_{cluster_key}", group=None
            )

            # Extract top markers for each cluster
            top_markers = {}
            for group in clusters:
                group_markers = markers_df[markers_df["group"] == group]
                # Filter by significance and fold change
                sig_markers = group_markers[
                    (group_markers["pvals_adj"] < 0.05)
                    & (group_markers["logfoldchanges"] > 0.5)
                ]
                # Take top 50 markers, or all significant ones if fewer
                n_markers = min(50, len(sig_markers))
                top_markers[group] = sig_markers["names"].head(n_markers).tolist()
                log.debug(f"Cluster {group}: {len(top_markers[group])} marker genes")

            # Compute Jaccard similarity between marker sets
            log.info("Computing Jaccard similarity between marker gene sets")
            sim_matrix = np.zeros((n_clusters, n_clusters))
            for i, c1 in enumerate(clusters):
                for j, c2 in enumerate(clusters):
                    if i >= j:  # Skip redundant calculations
                        continue
                    set1, set2 = set(top_markers[c1]), set(top_markers[c2])
                    if not set1 or not set2:  # Skip if either set is empty
                        continue
                    # Jaccard similarity: intersection size / union size
                    sim = len(set1.intersection(set2)) / len(set1.union(set2))
                    sim_matrix[i, j] = sim_matrix[j, i] = sim

        except Exception as e:
            log.error(f"Error computing marker-based similarity: {str(e)}")
            raise RuntimeError(f"Marker overlap calculation failed: {str(e)}")

    elif method == "expression_correlation":
        log.info("Computing mean expression profiles for each cluster")

        try:
            # Get appropriate expression data
            if scipy.sparse.issparse(adata.X):
                # For sparse matrices, we need to handle the computation differently
                mean_profiles = {}
                for c in clusters:
                    # Get cells in this cluster
                    mask = adata.obs[cluster_key] == c
                    # Compute mean profile
                    if mask.sum() > 0:
                        cluster_mean = adata[mask].X.mean(axis=0)
                        if hasattr(cluster_mean, "A"):  # Convert from matrix to array
                            cluster_mean = cluster_mean.A1
                        else:
                            cluster_mean = cluster_mean.flatten()
                        mean_profiles[c] = cluster_mean

                # Create DataFrame from mean profiles
                mean_df = pd.DataFrame(mean_profiles, index=adata.var_names)

            else:
                # For dense matrices, we can use a more direct approach
                mean_profiles = {
                    c: adata[adata.obs[cluster_key] == c].X.mean(axis=0)
                    for c in clusters
                }
                mean_df = pd.DataFrame(mean_profiles, index=adata.var_names)

            # Compute correlation matrix
            log.info("Computing correlation between mean expression profiles")
            sim_matrix = np.corrcoef(mean_df.values.T)

        except Exception as e:
            log.error(f"Error computing expression correlation: {str(e)}")
            raise RuntimeError(f"Expression correlation calculation failed: {str(e)}")

    else:
        available_methods = ["marker_overlap", "expression_correlation"]
        log.error(f"Unknown method: {method}")
        raise ValueError(f"Unknown method: {method}. Choose from {available_methods}")

    # Build graph where edges connect clusters with similarity > threshold
    log.info("Building similarity graph for merging")
    G = nx.from_numpy_array(sim_matrix > similarity_threshold)

    # Identify connected components (these will be our merged clusters)
    components = list(nx.connected_components(G))
    log.info(f"Found {len(components)} connected components")

    # Create mapping from original cluster to merged cluster
    mapping = {}
    for i, comp in enumerate(components):
        for node in comp:
            mapping[clusters[node]] = f"M{i + 1}"

    # Apply mapping to create new cluster assignments
    adata.obs[key_added] = adata.obs[cluster_key].map(mapping).astype("category")

    # Store merging metadata
    adata.uns[f"{key_added}_params"] = {
        "source_clusters": cluster_key,
        "method": method,
        "similarity_threshold": similarity_threshold,
        "original_clusters": n_clusters,
        "merged_clusters": len(components),
        "mapping": {
            str(k): v for k, v in mapping.items()
        },  # Convert keys to strings for JSON compatibility
    }

    log.info(
        f"Merged {n_clusters} clusters into {len(components)} super-clusters "
        f"using '{method}' with threshold {similarity_threshold}"
    )
    log.info(f"Results stored in adata.obs['{key_added}']")

    return adata
