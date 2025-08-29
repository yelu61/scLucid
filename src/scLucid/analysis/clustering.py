"""
Clustering functions for single-cell RNA-seq analysis.

This module provides:
- Unsupervised clustering (Leiden, Louvain, K-means, HDBSCAN)
- Resolution optimization with marker-separation or silhouette metrics
- Cluster merging based on marker overlap or expression correlation
- Standardized config dataclasses, logging, and results traceability
"""

import dataclasses
import logging
from pathlib import Path
from typing import Dict, List, Optional

import entropy
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from sklearn import metrics

from ..utils.marker_manager import Manager
from .config import ClusteringConfig, MergeClustersConfig, ResolutionSearchConfig

log = logging.getLogger(__name__)

__all__ = [
    "find_resolution",
    "cluster_cells",
    "merge_clusters",
]


# ====================== Clustering Evaluation ======================
def _evaluate_marker_separation(
    adata: sc.AnnData,
    cluster_key: str,
    marker_genes: Dict[str, List[str]],
    use_raw: bool = False,
) -> float:
    """
    Evaluate clustering by marker gene separation.
    Returns an average entropy-based separation score (0-1, higher is better).
    """
    if use_raw:
        if adata.raw is None:
            log.error("adata.raw is None but use_raw=True")
            raise ValueError("adata.raw must be set for use_raw=True")
        X = adata.raw.X
        var_names = adata.raw.var_names
    else:
        X = adata.layers.get("log1p_norm", adata.X)
        var_names = adata.var_names

    clusters = adata.obs[cluster_key].cat.categories
    n_clusters = len(clusters)
    if n_clusters <= 1:
        log.warning("Only one cluster found. Marker separation score is 0.")
        return 0.0

    scores = []
    for cell_type, markers in marker_genes.items():
        marker_indices = [i for i, gene in enumerate(var_names) if gene in markers]
        if not marker_indices:
            continue
        try:
            cluster_means = np.array(
                [
                    X[(adata.obs[cluster_key] == c).values][:, marker_indices].mean()
                    for c in clusters
                ]
            )
            if np.sum(cluster_means) == 0:
                continue
            cluster_probs = cluster_means / np.sum(cluster_means)
            marker_entropy = entropy(cluster_probs)
            max_entropy = np.log(n_clusters)
            score = 1 - (marker_entropy / max_entropy) if max_entropy > 0 else 0.0
            scores.append(score)
        except Exception as e:
            log.warning(f"Error for cell type '{cell_type}': {e}")
    if not scores:
        return 0.0
    final_score = float(np.mean(scores))
    log.info(
        f"Marker separation score: {final_score:.4f} from {len(scores)} cell types"
    )
    return final_score


def _evaluate_silhouette(
    adata: sc.AnnData, cluster_key: str, use_rep: str = "X_pca"
) -> float:
    """
    Evaluate clustering by silhouette score (-1~1, higher is better).
    """
    if use_rep not in adata.obsm:
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")
    labels = adata.obs[cluster_key].cat.codes
    n_labels = len(np.unique(labels))
    if n_labels <= 1:
        log.warning("Only one cluster found. Silhouette score is 0.")
        return 0.0
    try:
        score = metrics.silhouette_score(
            adata.obsm[use_rep],
            labels,
            sample_size=min(20000, adata.n_obs),
            random_state=42,
        )
        log.info(f"Silhouette score: {score:.4f}")
        return float(score)
    except Exception as e:
        log.warning(f"Silhouette error: {e}")
        return 0.0


# ====================== Main Functions ======================
def find_resolution(
    adata: sc.AnnData,
    config: Optional[ResolutionSearchConfig] = None,
    **kwargs,
) -> sc.AnnData:
    """
    Grid search clustering resolution and select optimal by marker separation or silhouette.
    Stores evaluation results and final clusters in adata.uns['sclucid']['analysis']['clustering'].
    """
    if config is None:
        active_config = ResolutionSearchConfig()
    else:
        active_config = dataclasses.replace(config)
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    # Parse parameters
    start, end, steps = active_config.resolution_range
    resolutions = np.linspace(start, end, steps)
    method = "leiden"  # Only leiden/louvain supported here

    # Prepare marker_manager if needed
    marker_genes = {}
    if active_config.metric == "marker_separation":
        marker_mgr = None
        if isinstance(active_config.marker_config, Manager):
            marker_mgr = active_config.marker_config
        elif isinstance(active_config.marker_config, str):
            marker_mgr = Manager(active_config.marker_config)
        else:
            raise ValueError(
                "marker_config must be str or Manager for marker_separation"
            )
        marker_mgr.intersect_with(
            adata.raw if active_config.use_raw_for_markers else adata
        )
        marker_genes = {cell.name: cell.markers for cell in marker_mgr.CELLS.values()}

    eval_results = []
    for res in resolutions:
        key = f"{method}_res_{res:.2f}"
        sc.pp.neighbors(adata, use_rep=active_config.use_rep)
        if method == "leiden":
            sc.tl.leiden(adata, resolution=res, key_added=key)
        else:
            sc.tl.louvain(adata, resolution=res, key_added=key)
        if active_config.metric == "marker_separation":
            score = _evaluate_marker_separation(
                adata, key, marker_genes, use_raw=active_config.use_raw_for_markers
            )
        else:
            score = _evaluate_silhouette(adata, key, active_config.use_rep)
        n_clusters = adata.obs[key].nunique()
        eval_results.append(
            {"resolution": res, "n_clusters": n_clusters, "score": score}
        )
        log.info(
            f"Res={res:.2f}  clusters={n_clusters}  {active_config.metric}={score:.4f}"
        )

    eval_df = pd.DataFrame(eval_results)
    if eval_df.empty:
        raise RuntimeError("No valid clustering results obtained.")

    best_idx = eval_df["score"].idxmax()
    best_res = eval_df.loc[best_idx, "resolution"]
    best_clusters = eval_df.loc[best_idx, "n_clusters"]
    best_score = eval_df.loc[best_idx, "score"]
    log.info(
        f"Optimal resolution: {best_res:.2f} clusters={best_clusters}, {active_config.metric}={best_score:.4f}"
    )

    # Final clustering
    key_added = f"{method}_optimal"
    sc.pp.neighbors(adata, use_rep=active_config.use_rep)
    if method == "leiden":
        sc.tl.leiden(adata, resolution=best_res, key_added=key_added)
    else:
        sc.tl.louvain(adata, resolution=best_res, key_added=key_added)

    # Save results
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "clustering", {}
    )
    adata.uns["sclucid"]["analysis"]["clustering"][f"{key_added}_evaluation"] = {
        "metric": active_config.metric,
        "optimal_resolution": float(best_res),
        "optimal_score": float(best_score),
        "resolutions": eval_df["resolution"].tolist(),
        "scores": eval_df["score"].tolist(),
        "n_clusters": eval_df["n_clusters"].tolist(),
        "parameters": active_config.to_dict(),
    }

    if active_config.plot:
        plt.figure(figsize=(8, 6))
        plt.plot(
            eval_df["resolution"],
            eval_df["score"],
            "o-",
            label=f"{active_config.metric} score",
        )
        plt.plot(
            eval_df["resolution"],
            eval_df["n_clusters"],
            "s--",
            label="#Clusters",
            color="orange",
        )
        plt.axvline(
            best_res, color="red", linestyle="--", label=f"Optimal: {best_res:.2f}"
        )
        plt.xlabel("Resolution")
        plt.title("Clustering Resolution Optimization")
        plt.legend()
        plt.tight_layout()
        if active_config.save_dir:
            save_path = Path(active_config.save_dir)
            # This line creates the directory if it doesn't exist
            save_path.mkdir(parents=True, exist_ok=True)
            # Now, save the file to the guaranteed-to-exist directory
            figure_path = save_path / "resolution_search.png"
            plt.savefig(figure_path, dpi=300)
            log.info(f"Saved resolution search plot to {figure_path}")

        plt.show()
    return adata


def cluster_cells(
    adata: sc.AnnData,
    config: Optional[ClusteringConfig] = None,
    **kwargs,
) -> sc.AnnData:
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
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm. Please run PCA or another dimensionality reduction first.")
    key_added = active_config.key_added or f"{method}_cluster"
    
    # Compute neighbors if needed
    if method in ["leiden", "louvain"] and "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, use_rep=use_rep, random_state=active_config.random_state)

    # Run clustering
    if method == "leiden":
        sc.tl.leiden(
            adata,
            resolution=active_config.resolution,
            key_added=key_added,
            random_state=active_config.random_state,
        )
    elif method == "louvain":
        sc.tl.louvain(
            adata,
            resolution=active_config.resolution,
            key_added=key_added,
            random_state=active_config.random_state,
        )
    elif method == "kmeans":
        from sklearn.cluster import KMeans

        X = adata.obsm[use_rep]
        kmeans_params = dict(n_init=10, max_iter=300, **config.extra_params)
        kmeans = KMeans(
            n_clusters=active_config.n_clusters,
            random_state=active_config.random_state,
            **kmeans_params,
        )
        labels = kmeans.fit_predict(X)
        adata.obs[key_added] = labels.astype(str)
        adata.obs[key_added] = adata.obs[key_added].astype("category")
    elif method == "hdbscan":
        try:
            import hdbscan
        except ImportError:
            raise ImportError("Please install hdbscan: pip install hdbscan")
        X = adata.obsm[use_rep]
        hdb_params = dict(min_cluster_size=50, min_samples=10, **config.extra_params)
        clusterer = hdbscan.HDBSCAN(**hdb_params)
        labels = clusterer.fit_predict(X)
        adata.obs[key_added] = labels.astype(str)
        adata.obs[key_added] = adata.obs[key_added].astype("category")
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

    # Visualization
    if config.plot:
        if "X_umap" not in adata.obsm:
            sc.tl.umap(adata)
        plt.figure(figsize=(8, 7))
        sc.pl.umap(
            adata,
            color=key_added,
            legend_loc="on data",
            title=f"{method.capitalize()} Clustering",
            show=False,
        )
        if config.save_dir:
            save_path = Path(config.save_dir)
            # This line creates the directory if it doesn't exist
            save_path.mkdir(parents=True, exist_ok=True)
            # Now, save the file to the guaranteed-to-exist directory
            figure_path = save_path / f"{key_added}_umap.png"
            plt.savefig(figure_path, dpi=300)
            log.info(f"Saved UMAP to {figure_path}")

        plt.show()

    return adata


def merge_clusters(
    adata: sc.AnnData,
    config: Optional[MergeClustersConfig] = None,
    **kwargs,
) -> sc.AnnData:
    """
    Merge similar clusters based on marker overlap or expression correlation.

    Stores merged cluster assignments in adata.obs[key_added] and parameters in .uns.
    """
    if config is None:
        active_config = MergeClustersConfig()
    else:
        active_config = dataclasses.replace(config)
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    # Extract parameters from the final config
    cluster_key = active_config.cluster_key
    similarity_threshold = active_config.similarity_threshold
    method = active_config.method
    key_added = active_config.key_added or f"{cluster_key}_merged"

    if cluster_key not in adata.obs:
        raise ValueError(f"cluster_key '{cluster_key}' not found in adata.obs")

    clusters = adata.obs[cluster_key].cat.categories
    n_clusters = len(clusters)

    # Compute similarity matrix
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
            group: markers_df[markers_df["group"] == group]
            .query("pvals_adj < 0.05 and logfoldchanges > 0.5")
            .head(50)["names"]
            .tolist()
            for group in clusters
        }
        sim_matrix = np.zeros((n_clusters, n_clusters))
        for i, c1 in enumerate(clusters):
            for j, c2 in enumerate(clusters):
                if i >= j:
                    continue
                set1, set2 = set(top_markers[c1]), set(top_markers[c2])
                if not set1 or not set2:
                    continue
                sim = len(set1 & set2) / len(set1 | set2)
                sim_matrix[i, j] = sim_matrix[j, i] = sim
    elif method == "expression_correlation":
        if scipy.sparse.issparse(adata.X):
            mean_profiles = {
                c: np.array(adata[adata.obs[cluster_key] == c].X.mean(axis=0)).flatten()
                for c in clusters
            }
        else:
            mean_profiles = {
                c: adata[adata.obs[cluster_key] == c].X.mean(axis=0) for c in clusters
            }
        mean_df = pd.DataFrame(mean_profiles, index=adata.var_names)
        sim_matrix = np.corrcoef(mean_df.values.T)
    else:
        raise ValueError(f"Unknown merge method: {method}")

    # Build graph and find connected components
    G = nx.from_numpy_array(sim_matrix > similarity_threshold)
    components = list(nx.connected_components(G))

    # Create the mapping from old cluster to new merged cluster
    mapping = {}
    for i, comp in enumerate(components):
        new_name = f"M{i + 1}"  # Merged cluster names
        for node_idx in comp:
            original_cluster_name = clusters[node_idx]
            mapping[original_cluster_name] = new_name

    adata.obs[key_added] = adata.obs[cluster_key].map(mapping).astype("category")

    # Store results and parameters in .uns for traceability
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "clustering", {}
    )[f"{key_added}_params"] = {
        "source_clusters": cluster_key,
        "method": method,
        "similarity_threshold": similarity_threshold,
        "original_clusters": n_clusters,
        "merged_clusters": len(components),
        "mapping": {str(k): v for k, v in mapping.items()},
        "config": active_config.to_dict(),
    }
    log.info(
        f"Clusters merged: {n_clusters} -> {len(components)} (method={method}, threshold={similarity_threshold})"
    )
    return adata
