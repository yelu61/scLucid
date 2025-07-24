"""
Clustering functions for single-cell RNA-seq data.

This module provides functions for marker-guided clustering,
optimal resolution selection, and cluster merging.
"""

from typing import Literal, Optional, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy
from sklearn import metrics

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
            [
                X[(adata.obs[cluster_key] == c).values][:, marker_indices].mean()
                for c in clusters
            ]
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
def optimize_neighbors_pcs(
    adata,
    n_neighbors_list,
    n_pcs_list,
    use_rep="X_scvi",
    clustering_method="leiden",
    progress=True,
    save_path=None,
    n_jobs=-1,  # Use all available cores by default
    compute_umap=False,  # Option to skip UMAP for faster processing
    subsample=None,  # Optional subsampling for very large datasets
):
    """
    Memory-efficient grid search for optimal n_neighbors and n_pcs parameters.

    Args:
        adata: AnnData object (will not be modified).
        n_neighbors_list: List of n_neighbors values to evaluate.
        n_pcs_list: List of n_pcs values to evaluate.
        use_rep: Dimensionality reduction to use (e.g., 'X_scvi').
        clustering_method: Clustering method ('leiden' or 'louvain').
        progress: Whether to show progress bar.
        save_path: If specified, save results to CSV.
        n_jobs: Number of parallel jobs for silhouette score calculation.
        compute_umap: Whether to compute UMAP (can be skipped for speed).
        subsample: Number of cells to subsample for large datasets (None=use all).

    Returns:
        pandas.DataFrame with clustering results and silhouette scores.
    """
    import gc
    import os

    import numpy as np
    from joblib import Parallel, delayed
    from sklearn.metrics import silhouette_score

    # If dataset is very large, subsample for parameter tuning
    if subsample is not None and subsample < adata.n_obs:
        print(f"Subsampling {subsample} cells from {adata.n_obs} for optimization")
        # Use deterministic sampling with fixed seed
        import numpy as np

        np.random.seed(42)
        indices = np.random.choice(adata.n_obs, subsample, replace=False)
        adata_opt = adata[indices].copy()
    else:
        adata_opt = adata.copy()  # Make just one copy instead of many

    results = []

    # Pre-extract the dimensional reduction embedding to avoid repeatedly accessing it
    X_embed = adata_opt.obsm[use_rep].copy() if use_rep in adata_opt.obsm else None

    # Define silhouette calculation function for parallel processing
    def compute_silhouette_for_params(n_neighbors, n_pcs):
        nonlocal np
        # Set up neighbors graph with specific parameters
        sc.pp.neighbors(
            adata_opt,
            use_rep=use_rep,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            key_added=f"neighbors_{n_neighbors}_{n_pcs}",  # Use unique key to avoid overwriting
        )

        # Compute clustering
        if clustering_method == "leiden":
            sc.tl.leiden(
                adata_opt,
                neighbors_key=f"neighbors_{n_neighbors}_{n_pcs}",
                key_added=f"leiden_{n_neighbors}_{n_pcs}",
                resolution=0.5,  # Use moderate resolution
            )
            cluster_key = f"leiden_{n_neighbors}_{n_pcs}"
        else:
            sc.tl.louvain(
                adata_opt,
                neighbors_key=f"neighbors_{n_neighbors}_{n_pcs}",
                key_added=f"louvain_{n_neighbors}_{n_pcs}",
                resolution=0.5,
            )
            cluster_key = f"louvain_{n_neighbors}_{n_pcs}"

        n_clusters = adata_opt.obs[cluster_key].nunique()

        # Compute silhouette score
        sil_score = np.nan
        if n_clusters > 1:
            try:
                if compute_umap:
                    # Only compute UMAP if specifically requested
                    sc.tl.umap(
                        adata_opt, neighbors_key=f"neighbors_{n_neighbors}_{n_pcs}"
                    )
                    embedding = adata_opt.obsm["X_umap"]
                else:
                    # Use the original embedding for silhouette calculation
                    # This is much faster but less visually interpretable
                    embedding = X_embed

                # Only use silhouette score if more than one cluster
                labels = adata_opt.obs[cluster_key].cat.codes
                sil_score = silhouette_score(
                    embedding, labels, sample_size=min(10000, len(labels))
                )

                # Clean up to save memory
                if compute_umap and "X_umap" in adata_opt.obsm:
                    del adata_opt.obsm["X_umap"]
            except Exception as e:
                print(
                    f"Error computing silhouette for n_neighbors={n_neighbors}, n_pcs={n_pcs}: {e}"
                )

        return {
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "n_clusters": n_clusters,
            "silhouette_score": sil_score,
        }

    # Generate parameter combinations
    param_combinations = [(n, p) for n in n_neighbors_list for p in n_pcs_list]

    # Run parameter search with progress bar
    if progress:
        from tqdm import tqdm

        results = [
            compute_silhouette_for_params(n, p)
            for n, p in tqdm(param_combinations, desc="Parameter optimization")
        ]
    else:
        # Optional: Use parallel processing for faster computation on multiple cores
        results = Parallel(n_jobs=n_jobs)(
            delayed(compute_silhouette_for_params)(n, p) for n, p in param_combinations
        )

    # Clean up temporary neighbor and cluster annotations
    for key in list(adata_opt.obs.keys()):
        if key.startswith(("leiden_", "louvain_")):
            del adata_opt.obs[key]

    for key in list(adata_opt.uns.keys()):
        if key.startswith("neighbors_"):
            del adata_opt.uns[key]

    for key in list(adata_opt.obsp.keys()):
        if "_connectivities" in key or "_distances" in key:
            if not key.startswith(
                ("connectivities", "distances")
            ):  # Keep the original ones
                del adata_opt.obsp[key]

    # Free memory
    del adata_opt
    gc.collect()

    # Create and return results DataFrame
    df_results = pd.DataFrame(results)
    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        df_results.to_csv(save_path, index=False)

    return df_results


# In scRNA/analysis/cluster.py


def find_resolution(
    adata: sc.AnnData,
    resolution_range: Tuple[float, float, int] = (0.1, 2.0, 10),
    clustering_method: Literal["leiden", "louvain"] = "leiden",
    metric: Literal["marker_separation", "silhouette"] = "marker_separation",
    marker_config: Optional[str | Manager] = None,
    use_rep: str = "X_pca",
    neighbors_key: Optional[str] = None,
    plot: bool = True,
    random_state: int = 42,
) -> sc.AnnData:
    """
    Perform clustering across a range of resolutions and find the optimal one.

    This function evaluates clustering quality. If a marker_config is provided,
    it defaults to using marker gene separation. Otherwise, it automatically
    falls back to using the silhouette score.

    Args:
        adata: AnnData object.
        resolution_range: Tuple of (start, end, steps) for resolution search.
        clustering_method: 'leiden' or 'louvain'.
        metric: Metric for evaluation ('marker_separation' or 'silhouette').
        marker_config: (Optional) A Manager instance or path to a marker TOML file.
        use_rep: Representation to use for silhouette score (e.g., 'X_pca').
        neighbors_key: Key in `adata.obsp` for the neighbors graph.
        plot: Whether to plot the evaluation metrics vs. resolution.
        random_state: Seed for clustering reproducibility.

    Returns:
        AnnData object with the optimal clustering stored in `adata.obs`.
    """
    # --- THIS IS THE NEW, SMARTER LOGIC ---
    effective_metric = metric
    if marker_config is None and metric == "marker_separation":
        effective_metric = "silhouette"
        print(
            "Warning: `marker_config` not provided. "
            f"Falling back to '{effective_metric}' metric for resolution optimization."
        )

    if effective_metric == "marker_separation" and marker_config is None:
        raise ValueError(
            "`marker_config` must be provided when explicitly using 'marker_separation' metric."
        )
    # --- END OF NEW LOGIC ---

    if isinstance(marker_config, Manager):
        mgr = marker_config
    elif isinstance(marker_config, str):
        mgr = Manager(marker_config)
    else:
        mgr = None

    if mgr:
        mgr.intersect_with(adata)

    marker_genes = (
        {cell.name: cell.markers for cell in mgr.CELLS.values()} if mgr else {}
    )

    start, end, steps = resolution_range
    resolutions = np.linspace(start, end, steps)
    eval_results = []

    print(f"Searching for optimal resolution using '{effective_metric}' metric...")
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

        if effective_metric == "marker_separation":
            score = _evaluate_marker_separation(adata, key, marker_genes)
        else:  # Silhouette
            score = _evaluate_silhouette(adata, key, use_rep)

        n_clusters = adata.obs[key].nunique()
        eval_results.append(
            {"resolution": res, "n_clusters": n_clusters, "score": score}
        )
        print(
            f"  Resolution {res:.2f}: {n_clusters} clusters, {effective_metric} score = {score:.4f}"
        )

    eval_df = pd.DataFrame(eval_results).dropna()
    if eval_df.empty:
        raise RuntimeError("Clustering evaluation failed for all resolutions.")

    optimal_idx = eval_df["score"].idxmax()
    optimal_res = eval_df.loc[optimal_idx, "resolution"]

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
            label=f"{effective_metric} score",
        )
        ax2.plot(
            eval_df["resolution"],
            eval_df["n_clusters"],
            "o-",
            color="g",
            label="# Clusters",
        )
        ax1.set_xlabel("Resolution")
        ax1.set_ylabel(f"{effective_metric.replace('_', ' ').title()} Score", color="b")
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

    Args:
        adata: AnnData object.
        cluster_key: Key in adata.obs containing cluster assignments.
        similarity_threshold: Threshold for similarity to merge clusters.
        method: Method to calculate similarity ('marker_overlap' or 'expression_correlation').
        key_added: Key in adata.obs to store merged cluster assignments.

    Returns:
        AnnData object with merged clusters in adata.obs[key_added].
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
