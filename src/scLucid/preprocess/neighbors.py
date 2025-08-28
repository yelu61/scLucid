"""
Functions for building and optimizing the nearest neighbor graph and evaluating parameters
with silhouette score and visualization for single-cell analysis.
"""

import gc
import logging
import os
from typing import Dict, List, Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

log = logging.getLogger(__name__)

__all__ = [
    "optimize_neighbors_pcs",
]


# --- Helper Functions ---#
def _compute_silhouette_for_params(
    adata_opt: AnnData,
    use_rep: str,
    clustering_method: Literal["leiden", "louvain"] = "leiden",
    resolution: float = 1.0,
    n_neighbors: int = 15,
    n_pcs: int = 30,
    compute_umap: bool = False,
) -> Dict:
    """
    Compute clustering and silhouette score for a given parameter set.

    Returns a dict of parameter values and silhouette score.
    """
    from sklearn.metrics import silhouette_score

    key_suffix = f"{n_neighbors}_{n_pcs}"
    # Pre-extract embedding for silhouette
    if use_rep not in adata_opt.obsm:
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")
    X_embed = adata_opt.obsm[use_rep].copy()

    try:
        sc.pp.neighbors(
            adata_opt,
            use_rep=use_rep,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            key_added=f"neighbors_{key_suffix}",
        )
        # Clustering
        if clustering_method == "leiden":
            sc.tl.leiden(
                adata_opt,
                neighbors_key=f"neighbors_{key_suffix}",
                key_added=f"leiden_{key_suffix}",
                resolution=resolution,
            )
            cluster_key = f"leiden_{key_suffix}"
        elif clustering_method == "louvain":
            sc.tl.louvain(
                adata_opt,
                neighbors_key=f"neighbors_{key_suffix}",
                key_added=f"louvain_{key_suffix}",
                resolution=resolution,
            )
            cluster_key = f"louvain_{key_suffix}"
        else:
            raise ValueError(f"Unknown clustering method: {clustering_method}")

        n_clusters = adata_opt.obs[cluster_key].nunique()
        if n_clusters <= 1:
            log.warning(f"Only one cluster: n_neighbors={n_neighbors}, n_pcs={n_pcs}")
            return {
                "n_neighbors": n_neighbors,
                "n_pcs": n_pcs,
                "n_clusters": n_clusters,
                "silhouette_score": np.nan,
            }

        # Silhouette
        sil_score = np.nan
        try:
            if compute_umap:
                sc.tl.umap(adata_opt, neighbors_key=f"neighbors_{key_suffix}")
                embedding = adata_opt.obsm["X_umap"]
            else:
                embedding = X_embed
            labels = adata_opt.obs[cluster_key].cat.codes
            sil_score = silhouette_score(
                embedding, labels, sample_size=min(10000, len(labels))
            )
            if compute_umap and "X_umap" in adata_opt.obsm:
                del adata_opt.obsm["X_umap"]
        except Exception as e:
            log.warning(
                f"Silhouette error: n_neighbors={n_neighbors}, n_pcs={n_pcs}: {str(e)}"
            )

        return {
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "n_clusters": n_clusters,
            "silhouette_score": sil_score,
        }
    except Exception as e:
        log.error(
            f"Grid search error for n_neighbors={n_neighbors}, n_pcs={n_pcs}: {str(e)}"
        )
        return {
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "n_clusters": np.nan,
            "silhouette_score": np.nan,
        }


def _plot_neighbors_grid_search(
    param_df: pd.DataFrame,
    title: str = "Silhouette Score Grid Search",
    save_path: Optional[str] = None,
    annot: bool = True,
    cmap: str = "YlGnBu",
) -> plt.Figure:
    """
    Visualize parameter grid search (n_neighbors vs n_pcs) as a heatmap.

    Args:
        param_df: DataFrame from optimize_neighbors_pcs.
        title: Plot title.
        save_path: Path to save figure.
        annot: Annotate heatmap cells.
        cmap: Matplotlib colormap.

    Returns:
        Matplotlib Figure object.
    """
    required_cols = {"n_neighbors", "n_pcs", "silhouette_score"}
    if not required_cols.issubset(param_df.columns):
        raise ValueError(f"param_df must contain columns: {required_cols}")
    pivot = param_df.pivot(
        index="n_neighbors", columns="n_pcs", values="silhouette_score"
    )
    fig, ax = plt.subplots(
        figsize=(max(8, 1.2 * len(pivot.columns)), max(6, 0.8 * len(pivot)))
    )
    sns.heatmap(pivot, annot=annot, fmt=".3f", cmap=cmap, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("n_pcs")
    ax.set_ylabel("n_neighbors")
    plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
        plt.savefig(save_path, dpi=200)
    plt.show()
    return fig


# --- Main Functions ---#
def optimize_neighbors_pcs(
    adata: AnnData,
    n_neighbors_list: List[int],
    n_pcs_list: List[int],
    use_rep: str = "X_pca",
    clustering_method: Literal["leiden", "louvain"] = "leiden",
    resolution: float = 1.0,
    compute_umap: bool = False,
    progress: bool = True,
    save_path: Optional[str] = None,
    n_jobs: int = 1,
    subsample: Optional[int] = None,
    copy: bool = False,
    plot: bool = True,
    plot_save_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Grid search for optimal n_neighbors and n_pcs parameters using silhouette score.
    Optionally visualize the results as a heatmap.

    Args:
        adata: AnnData object (not modified unless copy=True)
        n_neighbors_list: List of n_neighbors to try
        n_pcs_list: List of n_pcs to try
        use_rep: Dimensionality reduction (e.g., 'X_pca')
        clustering_method: 'leiden' or 'louvain'
        resolution: Clustering resolution parameter
        compute_umap: Whether to use UMAP for silhouette (slower)
        progress: Show progress bar
        save_path: Path to save CSV results
        n_jobs: Parallel jobs for grid search (set >1 for parallel)
        subsample: Subsample cells for speed (None=all)
        copy: Work on a copy of AnnData
        plot: Show heatmap of results
        plot_save_path: Path to save heatmap

    Returns:
        DataFrame with columns: n_neighbors, n_pcs, n_clusters, silhouette_score

    Example:
        >>> results = optimize_neighbors_pcs(
        ...     adata,
        ...     n_neighbors_list=[10, 20, 30],
        ...     n_pcs_list=[30, 50, 70],
        ...     use_rep="X_pca",
        ...     plot=True
        ... )
    """
    from joblib import Parallel, delayed

    log.info("Starting grid search for n_neighbors and n_pcs")
    log.info(
        f"Testing {len(n_neighbors_list)} neighbor values and {len(n_pcs_list)} PC values"
    )

    # Validation
    if use_rep not in adata.obsm:
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")
    if not n_neighbors_list or not n_pcs_list:
        raise ValueError("n_neighbors_list and n_pcs_list cannot be empty")
    if min(n_neighbors_list) <= 0 or min(n_pcs_list) <= 0:
        raise ValueError("n_neighbors and n_pcs values must be positive")

    # Copy or subsample
    adata_work = adata.copy() if copy else adata
    if subsample is not None and subsample < adata_work.n_obs:
        log.info(
            f"Subsampling {subsample} cells from {adata_work.n_obs} for grid search"
        )
        np.random.seed(42)
        indices = np.random.choice(adata_work.n_obs, subsample, replace=False)
        adata_opt = adata_work[indices].copy()
    else:
        adata_opt = adata_work.copy()

    # Prepare parameter grid
    param_combinations = [(n, p) for n in n_neighbors_list for p in n_pcs_list]
    log.info(f"Total parameter combinations: {len(param_combinations)}")

    # Compute results
    results = []
    if progress:
        try:
            from tqdm import tqdm

            iterator = tqdm(param_combinations, desc="Grid search")
        except ImportError:
            log.warning("tqdm not installed. Progress bar disabled.")
            iterator = param_combinations
    else:
        iterator = param_combinations

    if n_jobs == 1:
        for n, p in iterator:
            param_result = _compute_silhouette_for_params(
                adata_opt,
                use_rep,
                clustering_method,
                resolution,
                n_neighbors=n,
                n_pcs=p,
                compute_umap=compute_umap,
            )
            results.append(param_result)
    else:
        results = Parallel(n_jobs=n_jobs)(
            delayed(_compute_silhouette_for_params)(
                adata_opt, use_rep, clustering_method, resolution, n, p, compute_umap
            )
            for n, p in iterator
        )

    # Clean up temp clustering keys
    for key in list(adata_opt.obs.keys()):
        if key.startswith(("leiden_", "louvain_")):
            del adata_opt.obs[key]
    for key in list(adata_opt.uns.keys()):
        if key.startswith("neighbors_"):
            del adata_opt.uns[key]
    for key in list(adata_opt.obsp.keys()):
        if "_connectivities" in key or "_distances" in key:
            if not key.startswith(("connectivities", "distances")):
                del adata_opt.obsp[key]
    del adata_opt
    gc.collect()

    # Results DataFrame
    df_results = pd.DataFrame(results)
    valid_results = df_results.dropna(subset=["silhouette_score"])
    if len(valid_results) > 0:
        best_idx = valid_results["silhouette_score"].idxmax()
        best_params = valid_results.loc[best_idx]
        log.info(
            f"Best params: n_neighbors={int(best_params['n_neighbors'])}, "
            f"n_pcs={int(best_params['n_pcs'])}, silhouette={best_params['silhouette_score']:.4f}"
        )
    else:
        log.warning("No valid parameter combinations found")
    log.info("Grid search complete.")

    # Save results
    if save_path is not None:
        try:
            os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
            df_results.to_csv(save_path, index=False)
            log.info(f"Saved parameter grid search results to {save_path}")
        except Exception as e:
            log.error(f"Error saving results to {save_path}: {str(e)}")

    # Plot heatmap if requested
    if plot:
        try:
            _plot_neighbors_grid_search(
                df_results,
                title="Silhouette Score Grid Search",
                save_path=plot_save_path,
            )
        except Exception as e:
            log.warning(f"Failed to plot grid search heatmap: {str(e)}")

    return df_results
