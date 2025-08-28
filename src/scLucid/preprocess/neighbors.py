"""
Functions for building and optimizing the nearest neighbor graph and evaluating parameters
with silhouette score and visualization for single-cell analysis.
"""

import dataclasses
import logging
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from joblib import Parallel, delayed

from .config import NeighborsConfig

log = logging.getLogger(__name__)

__all__ = ["optimize_neighbors_pcs"]


# --- Helper Functions ---#
def _compute_silhouette_for_params(
    adata_view: AnnData, config: NeighborsConfig, n_neighbors: int, n_pcs: int
) -> Dict:
    """
    Compute clustering and silhouette score for a given parameter set.
    Now driven by the config object for consistency.
    """
    from sklearn.metrics import silhouette_score

    # Use a copy to avoid race conditions in parallel processing
    adata_local = adata_view.copy()

    try:
        sc.pp.neighbors(
            adata_local,
            use_rep=config.use_rep,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
        )
        if config.clustering_method == "leiden":
            sc.tl.leiden(adata_local, resolution=config.resolution, key_added="cluster")
        else:  # louvain
            sc.tl.louvain(
                adata_local, resolution=config.resolution, key_added="cluster"
            )

        n_clusters = adata_local.obs["cluster"].nunique()
        if n_clusters <= 1:
            return {
                "n_neighbors": n_neighbors,
                "n_pcs": n_pcs,
                "n_clusters": n_clusters,
                "silhouette_score": np.nan,
            }

        labels = adata_local.obs["cluster"].cat.codes
        sil_score = silhouette_score(
            adata_local.obsm[config.use_rep],
            labels,
            sample_size=min(10000, len(labels)),  # Sample for speed
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
    save_path: Optional[Path] = None,  # Use Path object
) -> plt.Figure:
    """
    Visualize parameter grid search results as a heatmap.
    """
    pivot = param_df.pivot(
        index="n_neighbors", columns="n_pcs", values="silhouette_score"
    )
    fig, ax = plt.subplots(
        figsize=(max(8, 1.2 * len(pivot.columns)), max(6, 0.8 * len(pivot)))
    )
    sns.heatmap(pivot, annot=True, fmt=".3f", cmap="viridis", ax=ax)
    ax.set_title("Silhouette Score for Neighbor and PC Optimization")
    ax.set_xlabel("Number of Principal Components (n_pcs)")
    ax.set_ylabel("Number of Neighbors (n_neighbors)")
    plt.tight_layout()

    # CHANGED: Use pathlib for robust path handling
    if save_path:
        save_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved optimization heatmap to {save_path}")

    return fig


# --- Main Functions ---#
def optimize_neighbors_pcs(
    adata: AnnData,
    config: Optional[NeighborsConfig] = None,
    **kwargs,
) -> pd.DataFrame:
    """
    Grid search for optimal n_neighbors and n_pcs using silhouette score.

    This function is driven by a config object and supports kwargs for overrides.

    Args:
        adata: AnnData object.
        config: A NeighborsConfig object. If None, a default config is created.
        **kwargs: Keyword arguments to override config parameters
                  (e.g., `n_pcs_list=[20, 30]`, `subsample=None`).

    Returns:
        DataFrame with grid search results.
    """
    # --- 1. Establish the final configuration ---
    if config is None:
        active_config = NeighborsConfig()
    else:
        active_config = dataclasses.replace(config)

    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)
        else:
            log.warning(f"Ignoring unknown neighbors parameter: '{key}'")

    # --- 2. Extract parameters and validate inputs ---
    save_dir = Path(active_config.save_dir) if active_config.save_dir else None

    if active_config.use_rep not in adata.obsm:
        raise ValueError(
            f"Representation '{active_config.use_rep}' not found in adata.obsm."
        )
    if not active_config.n_neighbors_list or not active_config.n_pcs_list:
        raise ValueError("n_neighbors_list and n_pcs_list cannot be empty.")

    log.info("Starting grid search for n_neighbors and n_pcs...")
    log.info(
        f"Testing {len(active_config.n_neighbors_list)} neighbor values and {len(active_config.n_pcs_list)} PC values."
    )

    # --- 3. Prepare data (subsampling) ---
    if active_config.subsample and active_config.subsample < adata.n_obs:
        log.info(
            f"Subsampling to {active_config.subsample} cells for optimization speed."
        )
        np.random.seed(42)
        indices = np.random.choice(adata.n_obs, active_config.subsample, replace=False)
        adata_opt = adata[indices].copy()
    else:
        adata_opt = adata.copy()

    # --- 4. Run grid search in parallel ---
    param_combinations = [
        (n, p) for n in active_config.n_neighbors_list for p in active_config.n_pcs_list
    ]

    try:
        from tqdm import tqdm

        iterator = tqdm(param_combinations, desc="Optimizing Neighbors/PCs")
    except ImportError:
        iterator = param_combinations

    results = Parallel(n_jobs=active_config.n_jobs)(
        delayed(_compute_silhouette_for_params)(adata_opt, active_config, n, p)
        for n, p in iterator
    )

    # --- 5. Process and display results ---
    df_results = pd.DataFrame(results).dropna(subset=["silhouette_score"])

    if not df_results.empty:
        best_params = df_results.loc[df_results["silhouette_score"].idxmax()]
        log.info("=" * 30)
        log.info("Optimization Complete")
        log.info("Optimal parameters found:")
        log.info(f"  n_neighbors: {int(best_params['n_neighbors'])}")
        log.info(f"  n_pcs: {int(best_params['n_pcs'])}")
        log.info(f"  Best Silhouette Score: {best_params['silhouette_score']:.4f}")
        log.info("=" * 30)
    else:
        log.warning("Grid search yielded no valid results.")

    # --- 6. Save results and plot ---
    if save_dir:
        results_path = save_dir / "neighbors_optimization_results.csv"
        df_results.to_csv(results_path, index=False)
        log.info(f"Saved optimization results to {results_path}")

    if active_config.plot:
        fig = _plot_neighbors_grid_search(
            df_results,
            save_path=save_dir / "neighbors_optimization_heatmap.png"
            if save_dir
            else None,
        )
        plt.show()

    return df_results
