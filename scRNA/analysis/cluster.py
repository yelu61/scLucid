"""
Clustering functions for single-cell RNA-seq data.

This module provides functions for marker-guided clustering,
optimal resolution selection, and cluster merging.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy
from sklearn import metrics
from typing import Optional, Tuple, Literal
import matplotlib.pyplot as plt
import networkx as nx

from .manager import Manager

# --- Helper Evaluation Functions ---


def _evaluate_marker_separation(adata, cluster_key, marker_genes):
    """Calculates how well clusters separate known marker gene sets."""
    if "log1p_norm" not in adata.layers:
        print(
            "Warning: 'log1p_norm' layer not found. Using adata.X for marker separation score."
        )
        X = adata.X
    else:
        X = adata.layers["log1p_norm"]

    clusters = adata.obs[cluster_key].cat.categories
    n_clusters = len(clusters)
    if n_clusters <= 1:
        return 0.0

    scores = []
    for cell_type, markers in marker_genes.items():
        marker_indices = [
            i for i, gene in enumerate(adata.var_names) if gene in markers
        ]
        if not marker_indices:
            continue

        # Mean expression of marker set in each cluster
        cluster_means = np.array(
            [X[adata.obs[cluster_key] == c][:, marker_indices].mean() for c in clusters]
        )

        if np.sum(cluster_means) == 0:
            continue

        # Normalize to a probability distribution
        cluster_probs = cluster_means / np.sum(cluster_means)

        # Calculate score as 1 - normalized_entropy
        marker_entropy = entropy(cluster_probs)
        max_entropy = np.log(n_clusters)
        score = 1 - (marker_entropy / max_entropy) if max_entropy > 0 else 0.0
        scores.append(score)

    return np.mean(scores) if scores else 0.0


def _evaluate_silhouette(adata, cluster_key, use_rep="X_pca"):
    """Calculates silhouette score."""
    labels = adata.obs[cluster_key].cat.codes
    if len(np.unique(labels)) <= 1:
        return 0.0
    return metrics.silhouette_score(adata.obsm[use_rep], labels)


# --- Main Clustering Functions ---


def find_resolution(
    adata: sc.AnnData,
    resolution_range: Tuple[float, float, int] = (0.1, 2.0, 10),
    clustering_method: Literal["leiden", "louvain"] = "leiden",
    metric: Literal["marker_separation", "silhouette"] = "marker_separation",
    marker_config: Optional[str] = None,
    use_rep: str = "X_pca",
    neighbors_key: Optional[str] = None,
    plot: bool = True,
    random_state: int = 42,
) -> sc.AnnData:
    """
    Perform clustering across a range of resolutions and find the optimal one.

    This function evaluates clustering quality using either marker gene separation
    or silhouette score to guide the selection of the best resolution.

    Args:
        adata: AnnData object.
        resolution_range: Tuple of (start, end, steps) for resolution search.
        clustering_method: 'leiden' or 'louvain'.
        metric: Metric for evaluation ('marker_separation' or 'silhouette').
        marker_config: Path to marker TOML file, required for 'marker_separation' metric.
        use_rep: Representation to use for silhouette score (e.g., 'X_pca').
        neighbors_key: Key in `adata.obsp` for the neighbors graph.
        plot: Whether to plot the evaluation metrics vs. resolution.
        random_state: Seed for clustering reproducibility.

    Returns:
        AnnData object with the optimal clustering stored in `adata.obs`.
    """
    if metric == "marker_separation" and marker_config is None:
        raise ValueError(
            "`marker_config` must be provided for 'marker_separation' metric."
        )

    mgr = Manager(marker_config) if marker_config else None
    marker_genes = (
        {cell.name: cell.markers for cell in mgr.CELLS.values()} if mgr else {}
    )

    start, end, steps = resolution_range
    resolutions = np.linspace(start, end, steps)
    eval_results = []

    print(f"Searching for optimal resolution using '{metric}' metric...")
    for res in resolutions:
        key = f"{clustering_method}_res_{res:.2f}"
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

        if metric == "marker_separation":
            score = _evaluate_marker_separation(adata, key, marker_genes)
        elif metric == "silhouette":
            score = _evaluate_silhouette(adata, key, use_rep)
        else:
            raise ValueError(f"Unknown metric: {metric}")

        n_clusters = adata.obs[key].nunique()
        eval_results.append(
            {"resolution": res, "n_clusters": n_clusters, "score": score}
        )
        print(
            f"  Resolution {res:.2f}: {n_clusters} clusters, {metric} score = {score:.4f}"
        )

    eval_df = pd.DataFrame(eval_results).dropna()
    if eval_df.empty:
        raise RuntimeError("Clustering evaluation failed for all resolutions.")

    optimal_idx = eval_df["score"].idxmax()
    optimal_res = eval_df.loc[optimal_idx, "resolution"]

    # Final clustering with optimal resolution
    key_added = f"{clustering_method}_optimal"
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

    adata.uns[f"{key_added}_evaluation"] = eval_df.to_dict("list")
    print(
        f"\nOptimal resolution found: {optimal_res:.2f}. Final clustering stored in `adata.obs['{key_added}']`."
    )

    if plot:
        fig, ax1 = plt.subplots(figsize=(8, 5))
        ax2 = ax1.twinx()
        ax1.plot(
            eval_df["resolution"],
            eval_df["score"],
            "o-",
            color="b",
            label=f"{metric} score",
        )
        ax2.plot(
            eval_df["resolution"],
            eval_df["n_clusters"],
            "o-",
            color="g",
            label="# Clusters",
        )
        ax1.set_xlabel("Resolution")
        ax1.set_ylabel(f"{metric.replace('_', ' ').title()} Score", color="b")
        ax2.set_ylabel("Number of Clusters", color="g")
        ax1.axvline(
            x=optimal_res,
            color="r",
            linestyle="--",
            label=f"Optimal Res: {optimal_res:.2f}",
        )
        fig.legend(
            loc="upper right", bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes
        )
        plt.title("Clustering Resolution Optimization")
        plt.show()

    return adata


def merge_clusters(
    adata: sc.AnnData,
    cluster_key: str,
    similarity_threshold: float = 0.8,
    method: Literal["marker_overlap", "expression_correlation"] = "marker_overlap",
    key_added: Optional[str] = None,
) -> sc.AnnData:
    """
    Merge similar clusters based on marker overlap or expression correlation.
    """
    if key_added is None:
        key_added = f"{cluster_key}_merged"

    clusters = adata.obs[cluster_key].cat.categories
    n_clusters = len(clusters)

    if method == "marker_overlap":
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            key_added=f"rank_genes_{cluster_key}",
        )
        markers_df = sc.get.rank_genes_groups_df(
            adata, key=f"rank_genes_{cluster_key}", group=None
        )
        top_markers = {
            group: markers_df[markers_df["group"] == group]["names"].head(50).tolist()
            for group in clusters
        }

        sim_matrix = np.zeros((n_clusters, n_clusters))
        for i, c1 in enumerate(clusters):
            for j, c2 in enumerate(clusters):
                if i >= j:
                    continue
                set1, set2 = set(top_markers[c1]), set(top_markers[c2])
                sim = (
                    len(set1.intersection(set2)) / len(set1.union(set2))
                    if set1.union(set2)
                    else 0
                )
                sim_matrix[i, j] = sim_matrix[j, i] = sim

    elif method == "expression_correlation":
        mean_profiles = pd.DataFrame(
            index=adata.var_names,
            columns=clusters,
            data={
                c: adata[adata.obs[cluster_key] == c, :].X.mean(axis=0).A1
                for c in clusters
            },
        )
        sim_matrix = np.corrcoef(mean_profiles.T)
    else:
        raise ValueError(f"Unknown method: {method}")

    G = nx.from_numpy_array(sim_matrix > similarity_threshold)
    components = list(nx.connected_components(G))

    mapping = {
        clusters[node]: f"M{i + 1}"
        for i, comp in enumerate(components)
        for node in comp
    }
    adata.obs[key_added] = adata.obs[cluster_key].map(mapping).astype("category")

    print(
        f"Merged {n_clusters} clusters into {len(components)} super-clusters based on '{method}'. Result in `adata.obs['{key_added}']`."
    )
    return adata
