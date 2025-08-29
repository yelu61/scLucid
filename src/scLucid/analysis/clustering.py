"""
Clustering functions for single-cell RNA-seq analysis.

This module provides:
- Unsupervised clustering (Leiden, Louvain, K-means, HDBSCAN)
- Resolution optimization providing comprehensive metrics to guide user choice.
- Cluster merging based on marker overlap or expression correlation.
- Standardized config dataclasses, logging, and results traceability.
"""

import dataclasses
import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from sklearn import metrics

from .config import ClusteringConfig, MergeClustersConfig, ResolutionSearchConfig

log = logging.getLogger(__name__)

__all__ = [
    "find_resolution",
    "cluster_cells",
    "merge_clusters",
]

# ====================== Clustering Evaluation Helpers ======================


def _get_marker_abundance(
    adata: AnnData, cluster_key: str, de_method: str, min_log2fc: float, min_pct: float
) -> float:
    """
    Calculate the average number of significant markers per cluster.
    A higher number suggests better biological separability.
    """
    try:
        # Use a temporary key to avoid overwriting important results
        rank_key = f"rank_genes_{cluster_key}_temp"
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=de_method,
            key_added=rank_key,
            n_genes=100,  # Only need top genes to check for existence
            use_raw=False,  # Use normalized data for marker finding
        )

        markers_df = sc.get.rank_genes_groups_df(adata, key=rank_key, group=None)

        # Filter for significant markers
        sig_markers = markers_df[
            (markers_df["logfoldchanges"] > min_log2fc)
            & (markers_df["pvals_adj"] < 0.05)
            & (markers_df["pct_nz_group"] > min_pct)
        ]

        # Calculate average number of markers per cluster
        marker_counts = sig_markers.groupby("group").size()
        n_clusters = adata.obs[cluster_key].nunique()

        # For clusters with 0 markers, their count will be missing. Fill with 0.
        avg_markers = marker_counts.reindex(
            adata.obs[cluster_key].cat.categories, fill_value=0
        ).mean()

        # Clean up temporary results
        del adata.uns[rank_key]

        return float(avg_markers)

    except Exception as e:
        log.warning(f"Could not compute marker abundance for {cluster_key}: {e}")
        return 0.0


def _get_clustering_stability(adata: AnnData, key1: str, key2: str) -> float:
    """
    Calculate stability between two clustering results using Normalized Mutual Information (NMI).
    """
    if key1 not in adata.obs or key2 not in adata.obs:
        return np.nan

    labels1 = adata.obs[key1]
    labels2 = adata.obs[key2]

    return metrics.normalized_mutual_info_score(
        labels1, labels2, average_method="arithmetic"
    )


# ====================== Main Functions ======================


def find_resolution(
    adata: AnnData,
    config: Optional[ResolutionSearchConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Grid search clustering resolutions and provide comprehensive metrics to guide selection.

    This function is a decision support tool. It does NOT set a final "optimal" resolution.
    Instead, it returns a DataFrame and a plot showing how various metrics change
    with resolution, empowering the user to make an informed choice.

    Metrics:
    - n_clusters: The number of clusters found.
    - silhouette: Measures cluster separation (higher is better, but biased towards fewer clusters).
    - marker_abundance: Average number of significant markers per cluster (higher is better).
    - stability (NMI): How much the clustering changes from the previous resolution step.
      A sharp drop may indicate instability.

    Args:
        adata: AnnData object with a PCA representation.
        config: A ResolutionSearchConfig object.
        **kwargs: Keyword arguments to override config parameters.

    Returns:
        A pandas DataFrame containing the evaluation metrics for each resolution.
    """
    if config is None:
        active_config = ResolutionSearchConfig()
    else:
        active_config = dataclasses.replace(config)
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    # Unpack config
    start, end, steps = active_config.resolution_range
    resolutions = np.linspace(start, end, steps)
    method = active_config.method
    use_rep = active_config.use_rep

    if "neighbors" not in adata.uns:
        log.info(f"Neighbors graph not found. Computing on '{use_rep}'.")
        sc.pp.neighbors(adata, use_rep=use_rep)

    eval_results = []
    previous_key = None

    log.info(f"Searching for optimal resolution in range [{start}, {end}]...")
    for res in resolutions:
        current_key = f"{method}_res_{res:.3f}"
        log.info(f"Testing resolution: {res:.3f}")

        if method == "leiden":
            sc.tl.leiden(adata, resolution=res, key_added=current_key)
        else:  # louvain
            sc.tl.louvain(adata, resolution=res, key_added=current_key)

        n_clusters = adata.obs[current_key].nunique()
        result = {"resolution": res, "n_clusters": n_clusters}

        # Calculate requested metrics
        if active_config.compute_silhouette:
            result["silhouette"] = metrics.silhouette_score(
                adata.obsm[use_rep],
                adata.obs[current_key],
                sample_size=min(20000, adata.n_obs),
            )

        if active_config.compute_marker_abundance:
            result["marker_abundance"] = _get_marker_abundance(
                adata,
                current_key,
                active_config.de_method_for_markers,
                active_config.min_log2fc_for_markers,
                active_config.min_pct_for_markers,
            )

        if active_config.compute_stability:
            result["stability"] = (
                _get_clustering_stability(adata, previous_key, current_key)
                if previous_key
                else np.nan
            )

        eval_results.append(result)
        previous_key = current_key

    eval_df = pd.DataFrame(eval_results)

    # Store results in .uns for traceability
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "clustering", {}
    )
    adata.uns["sclucid"]["analysis"]["clustering"]["resolution_search"] = {
        "results_df": eval_df,
        "parameters": active_config.to_dict(),
    }

    if active_config.plot:
        metrics_to_plot = [
            col
            for col in ["silhouette", "marker_abundance", "stability"]
            if col in eval_df.columns
        ]
        n_plots = len(metrics_to_plot) + 1  # +1 for n_clusters

        fig, axes = plt.subplots(n_plots, 1, figsize=(10, 4 * n_plots), sharex=True)
        fig.suptitle("Clustering Resolution Optimization Guide", fontsize=16)

        # Plot n_clusters
        ax = axes[0]
        sns.lineplot(data=eval_df, x="resolution", y="n_clusters", marker="o", ax=ax)
        ax.set_title("Number of Clusters vs. Resolution")
        ax.set_ylabel("Count")
        ax.grid(True, linestyle="--")

        # Plot other metrics
        for i, metric in enumerate(metrics_to_plot):
            ax = axes[i + 1]
            sns.lineplot(data=eval_df, x="resolution", y=metric, marker="o", ax=ax)
            ax.set_title(f"{metric.replace('_', ' ').title()} vs. Resolution")
            ax.set_ylabel("Score")
            ax.grid(True, linestyle="--")
            if metric == "silhouette":
                ax.text(
                    0.02,
                    0.05,
                    "Note: Higher is 'better' but often favors fewer clusters",
                    transform=ax.transAxes,
                    fontsize=9,
                    style="italic",
                )

        axes[-1].set_xlabel("Resolution")
        plt.tight_layout(rect=[0, 0.03, 1, 0.96])

        if active_config.save_dir:
            save_path = Path(active_config.save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            figure_path = save_path / "resolution_search_guide.png"
            plt.savefig(figure_path, dpi=300)
            log.info(f"Saved resolution search plot to {figure_path}")

        plt.show()

    return eval_df


def cluster_cells(
    adata: AnnData,
    config: Optional[ClusteringConfig] = None,
    **kwargs,
) -> AnnData:
    """
    Perform clustering using Leiden, Louvain, KMeans, or HDBSCAN.

    Results and parameters are saved to adata.uns['sclucid']['analysis']['clustering'].
    """
    if config is None:
        active_config = ClusteringConfig()
    else:
        active_config = dataclasses.replace(config)
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    method = active_config.method
    use_rep = active_config.use_rep
    if use_rep not in adata.obsm:
        raise ValueError(
            f"Representation '{use_rep}' not found in adata.obsm. Please run PCA or another dimensionality reduction first."
        )
    key_added = active_config.key_added or f"{method}_clusters"

    if method in ["leiden", "louvain"] and "neighbors" not in adata.uns:
        log.info(
            f"Neighbors graph not found for Leiden/Louvain. Computing on '{use_rep}'."
        )
        sc.pp.neighbors(adata, use_rep=use_rep, random_state=active_config.random_state)

    log.info(f"Running {method} clustering...")
    if method == "leiden":
        sc.tl.leiden(
            adata,
            resolution=active_config.resolution,
            key_added=key_added,
            random_state=active_config.random_state,
            **active_config.extra_params,
        )
    elif method == "louvain":
        sc.tl.louvain(
            adata,
            resolution=active_config.resolution,
            key_added=key_added,
            random_state=active_config.random_state,
            **active_config.extra_params,
        )
    elif method == "kmeans":
        from sklearn.cluster import KMeans

        if active_config.n_clusters is None:
            raise ValueError("n_clusters must be specified for KMeans.")
        X = adata.obsm[use_rep]
        kmeans = KMeans(
            n_clusters=active_config.n_clusters,
            random_state=active_config.random_state,
            n_init=10,
            **active_config.extra_params,
        )
        labels = kmeans.fit_predict(X)
        adata.obs[key_added] = pd.Categorical(labels.astype(str))
    elif method == "hdbscan":
        try:
            import hdbscan
        except ImportError:
            raise ImportError("Please install hdbscan: pip install hdbscan")
        X = adata.obsm[use_rep]
        clusterer = hdbscan.HDBSCAN(**active_config.extra_params)
        labels = clusterer.fit_predict(X)
        # Convert -1 (noise) to a dedicated category
        adata.obs[key_added] = pd.Categorical(
            [str(l) if l != -1 else "Noise" for l in labels]
        )
    else:
        raise ValueError(f"Unknown clustering method: {method}")

    n_clusters = adata.obs[key_added].nunique()
    log.info(
        f"Clustering ({method}) finished: {n_clusters} clusters in obs['{key_added}']"
    )

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "clustering", {}
    )
    adata.uns["sclucid"]["analysis"]["clustering"][key_added] = {
        "config": active_config.to_dict(),
        "n_clusters": n_clusters,
    }

    if active_config.plot:
        if "X_umap" not in adata.obsm:
            sc.tl.umap(adata)

        sc.pl.umap(
            adata,
            color=key_added,
            legend_loc="on data",
            title=f"{method.capitalize()} Clustering (n={n_clusters})",
            show=False,
        )

        if active_config.save_dir:
            save_path = Path(active_config.save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            figure_path = save_path / f"{key_added}_umap.png"
            plt.savefig(figure_path, dpi=300, bbox_inches="tight")
            log.info(f"Saved UMAP to {figure_path}")

        plt.show()

    return adata


def merge_clusters(
    adata: AnnData,
    config: Optional[MergeClustersConfig] = None,
    **kwargs,
) -> AnnData:
    """
    Merge similar clusters based on marker overlap or expression correlation.
    """
    if config is None:
        active_config = MergeClustersConfig()
    else:
        active_config = dataclasses.replace(config)
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    cluster_key = active_config.cluster_key
    threshold = active_config.similarity_threshold
    method = active_config.method
    key_added = active_config.key_added or f"{cluster_key}_merged"

    if cluster_key not in adata.obs:
        raise ValueError(f"cluster_key '{cluster_key}' not found in adata.obs")

    original_clusters = adata.obs[cluster_key].cat.categories
    n_original = len(original_clusters)
    sim_matrix = pd.DataFrame(
        np.eye(n_original), index=original_clusters, columns=original_clusters
    )

    log.info(f"Calculating similarity matrix using '{method}' method...")
    if method == "marker_overlap":
        rank_key = f"rank_genes_{cluster_key}_for_merge"
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=active_config.de_method_for_markers,
            key_added=rank_key,
        )
        markers = sc.get.rank_genes_groups_df(adata, key=rank_key, group=None)
        top_markers = {
            g: set(
                df.query("pvals_adj < 0.05 & logfoldchanges > 0.5").head(50)["names"]
            )
            for g, df in markers.groupby("group")
        }
        for i, c1 in enumerate(original_clusters):
            for j, c2 in enumerate(original_clusters[i + 1 :], i + 1):
                set1, set2 = top_markers.get(c1, set()), top_markers.get(c2, set())
                if not set1 or not set2:
                    continue
                sim = len(set1 & set2) / len(set1 | set2)
                sim_matrix.loc[c1, c2] = sim_matrix.loc[c2, c1] = sim
    elif method == "expression_correlation":
        mean_profiles = adata.obs.groupby(cluster_key).apply(
            lambda x: adata[x.index].X.mean(axis=0)
        )
        mean_profiles = np.vstack(mean_profiles)
        corr_matrix = np.corrcoef(mean_profiles)
        sim_matrix = pd.DataFrame(
            corr_matrix, index=original_clusters, columns=original_clusters
        )

    G = nx.from_pandas_adjacency(sim_matrix > threshold)
    components = list(nx.connected_components(G))

    mapping = {}
    new_cluster_names = []
    for i, comp_nodes in enumerate(components):
        # Name merged clusters by concatenating original names, e.g., "0_3_5"
        new_name = "_".join(
            sorted(list(comp_nodes), key=lambda x: int(x) if x.isdigit() else x)
        )
        new_cluster_names.append(new_name)
        for node in comp_nodes:
            mapping[node] = new_name

    adata.obs[key_added] = adata.obs[cluster_key].map(mapping).astype("category")

    n_merged = len(new_cluster_names)
    log.info(
        f"Clusters merged: {n_original} -> {n_merged} (method={method}, threshold={threshold})"
    )

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "clustering", {}
    )
    adata.uns["sclucid"]["analysis"]["clustering"][f"{key_added}_params"] = {
        "source_clusters": cluster_key,
        "method": method,
        "similarity_threshold": threshold,
        "original_clusters": n_original,
        "merged_clusters": n_merged,
        "mapping": mapping,
        "config": active_config.to_dict(),
    }

    return adata
