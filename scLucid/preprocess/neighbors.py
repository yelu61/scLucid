"""
Functions for building and optimizing the nearest neighbor graph.
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
    "optimize_neighbors_pcs",
]


# --- Helper functions ---
def _compute_silhouette_for_params(adata_opt, use_rep, clustering_method, resolution, n_neighbors, n_pcs, compute_umap):
        # Pre-extract the dimensional reduction embedding to avoid repeatedly accessing it
        if use_rep in adata_opt.obsm:
            X_embed = adata_opt.obsm[use_rep].copy()
        else:
            log.error(f"Representation '{use_rep}' not found in adata.obsm")
            raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")

        # Define silhouette calculation function for parallel processing
        log.debug(f"Testing n_neighbors={n_neighbors}, n_pcs={n_pcs}")
        from sklearn.metrics import silhouette_score
        # Create a unique key for this parameter combination
        key_suffix = f"{n_neighbors}_{n_pcs}"

        try:
            # Set up neighbors graph with specific parameters
            sc.pp.neighbors(
                adata_opt,
                use_rep=use_rep,
                n_neighbors=n_neighbors,
                n_pcs=n_pcs,
                key_added=f"neighbors_{key_suffix}",
            )

            # Compute clustering
            if clustering_method == "leiden":
                sc.tl.leiden(
                    adata_opt,
                    neighbors_key=f"neighbors_{key_suffix}",
                    key_added=f"leiden_{key_suffix}",
                    resolution=resolution,
                )
                cluster_key = f"leiden_{key_suffix}"
            else:
                sc.tl.louvain(
                    adata_opt,
                    neighbors_key=f"neighbors_{key_suffix}",
                    key_added=f"louvain_{key_suffix}",
                    resolution=resolution,
                )
                cluster_key = f"louvain_{key_suffix}"

            n_clusters = adata_opt.obs[cluster_key].nunique()
            log.debug(
                f"Found {n_clusters} clusters with n_neighbors={n_neighbors}, n_pcs={n_pcs}"
            )

            # Skip silhouette calculation if only one cluster is found
            if n_clusters <= 1:
                log.warning(
                    f"Only one cluster found with n_neighbors={n_neighbors}, n_pcs={n_pcs}"
                )
                return {
                    "n_neighbors": n_neighbors,
                    "n_pcs": n_pcs,
                    "n_clusters": n_clusters,
                    "silhouette_score": np.nan,
                }

            # Compute silhouette score
            sil_score = np.nan
            try:
                if compute_umap:
                    # Only compute UMAP if specifically requested
                    sc.tl.umap(adata_opt, neighbors_key=f"neighbors_{key_suffix}")
                    embedding = adata_opt.obsm["X_umap"]
                else:
                    # Use the original embedding for silhouette calculation
                    embedding = X_embed

                # Calculate silhouette score
                labels = adata_opt.obs[cluster_key].cat.codes
                sil_score = silhouette_score(
                    embedding, labels, sample_size=min(10000, len(labels))
                )
                log.debug(f"Silhouette score: {sil_score:.4f}")

                # Clean up to save memory
                if compute_umap and "X_umap" in adata_opt.obsm:
                    del adata_opt.obsm["X_umap"]

            except Exception as e:
                log.warning(
                    f"Error computing silhouette for n_neighbors={n_neighbors}, n_pcs={n_pcs}: {str(e)}"
                )

            return {
                "n_neighbors": n_neighbors,
                "n_pcs": n_pcs,
                "n_clusters": n_clusters,
                "silhouette_score": sil_score,
            }

        except Exception as e:
            log.error(
                f"Error in parameter testing for n_neighbors={n_neighbors}, n_pcs={n_pcs}: {str(e)}"
            )
            return {
                "n_neighbors": n_neighbors,
                "n_pcs": n_pcs,
                "n_clusters": np.nan,
                "silhouette_score": np.nan,
            }

# --- Main Functions ---
def optimize_neighbors_pcs(
    adata: sc.AnnData,
    n_neighbors_list: List[int],
    n_pcs_list: List[int],
    use_rep: str = "X_pca",
    progress: bool = True,
    save_path: Optional[str] = None,
    n_jobs: int = -1,
    subsample: Optional[int] = None,
    copy: bool = False,
) -> pd.DataFrame:
    """
    Performs grid search to find optimal n_neighbors and n_pcs parameters.

    This function evaluates different combinations of nearest neighbors and principal
    components to find the best parameters for clustering. It uses silhouette score
    as the evaluation metric.

    Args:
        adata: AnnData object (will not be modified unless copy=True)
        n_neighbors_list: List of n_neighbors values to evaluate
        n_pcs_list: List of n_pcs values to evaluate
        use_rep: Dimensionality reduction to use (e.g., 'X_pca')
        clustering_method: Clustering method ('leiden' or 'louvain')
        resolution: Resolution parameter for clustering
        progress: Whether to show progress bar
        save_path: If specified, save results to CSV
        n_jobs: Number of parallel jobs for silhouette score calculation
        compute_umap: Whether to compute UMAP for evaluation (slower but often better)
        subsample: Number of cells to subsample for large datasets (None=use all)
        copy: If True, work on a copy of the AnnData object

    Returns:
        pandas.DataFrame with clustering results and silhouette scores

    Examples:
        >>> results = optimize_neighbors_pcs(
        ...     adata,
        ...     n_neighbors_list=[10, 20, 30],
        ...     n_pcs_list=[30, 50, 70],
        ...     use_rep="X_pca"
        ... )
        >>> print(results.sort_values('silhouette_score', ascending=False).head())
    """
    from joblib import Parallel, delayed
    

    log.info("Starting grid search for optimal n_neighbors and n_pcs parameters")
    log.info(
        f"Testing {len(n_neighbors_list)} neighbor values and {len(n_pcs_list)} PC values"
    )

    # Parameter validation
    if use_rep not in adata.obsm:
        log.error(f"Representation '{use_rep}' not found in adata.obsm")
        raise ValueError(f"Representation '{use_rep}' not found in adata.obsm")

    if not n_neighbors_list or not n_pcs_list:
        log.error("Empty parameter lists provided")
        raise ValueError("n_neighbors_list and n_pcs_list cannot be empty")

    if min(n_neighbors_list) <= 0 or min(n_pcs_list) <= 0:
        log.error("Parameters must be positive")
        raise ValueError("n_neighbors and n_pcs values must be positive")

    # Create a copy of the data for optimization if requested
    if copy:
        adata_work = adata.copy()
    else:
        adata_work = adata

    # If dataset is very large, subsample for parameter tuning
    if subsample is not None and subsample < adata.n_obs:
        log.info(f"Subsampling {subsample} cells from {adata.n_obs} for optimization")
        # Use deterministic sampling with fixed seed
        np.random.seed(42)
        indices = np.random.choice(adata.n_obs, subsample, replace=False)
        adata_opt = adata_work[indices].copy()
    else:
        adata_opt = adata_work.copy()  # Make just one copy instead of many

    # Generate parameter combinations
    param_combinations = [(n, p) for n in n_neighbors_list for p in n_pcs_list]
    log.info(f"Testing {len(param_combinations)} parameter combinations")

    # Run parameter search with progress bar
    results = []
    if progress:
        try:
            from tqdm import tqdm

            results = [
                _compute_silhouette_for_params(n, p)
                for n, p in tqdm(param_combinations, desc="Parameter optimization")
            ]
        except ImportError:
            log.warning("tqdm not installed. Progress bar disabled.")
            results = [
                _compute_silhouette_for_params(n, p) for n, p in param_combinations
            ]
    else:
        # Use parallel processing for faster computation on multiple cores
        log.info(f"Running in parallel with {n_jobs} jobs")
        results = Parallel(n_jobs=n_jobs)(
            delayed(_compute_silhouette_for_params)(n, p) for n, p in param_combinations
        )

    # Clean up temporary neighbor and cluster annotations
    log.info("Cleaning up temporary data")
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

    # Create results DataFrame and remove NaN rows
    df_results = pd.DataFrame(results)
    valid_results = df_results.dropna(subset=["silhouette_score"])

    # Log results summary
    if len(valid_results) > 0:
        best_idx = valid_results["silhouette_score"].idxmax()
        best_params = valid_results.loc[best_idx]
        log.info(
            f"Best parameters: n_neighbors={int(best_params['n_neighbors'])}, "
            f"n_pcs={int(best_params['n_pcs'])}"
        )
        log.info(f"Best silhouette score: {best_params['silhouette_score']:.4f}")
    else:
        log.warning("No valid parameter combinations found")

    log.info(
        f"Found {len(valid_results)}/{len(param_combinations)} valid parameter combinations"
    )

    # Save results if requested
    if save_path is not None:
        try:
            os.makedirs(os.path.dirname(os.path.abspath(save_path)), exist_ok=True)
            df_results.to_csv(save_path, index=False)
            log.info(f"Saved results to {save_path}")
        except Exception as e:
            log.error(f"Error saving results to {save_path}: {str(e)}")

    return df_results


