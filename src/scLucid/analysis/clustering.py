"""
Clustering functions for single-cell RNA-seq analysis.

This module provides:
- Unsupervised clustering (Leiden, Louvain, K-means, HDBSCAN)
- Resolution optimization providing comprehensive metrics to guide user choice.
- Cluster merging based on marker overlap or expression correlation.
- Standardized config dataclasses, logging, and results traceability.
"""

import logging
from pathlib import Path
from typing import Optional, Literal

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from sklearn import metrics

from ..utils import sanitize_for_hdf5
from .config import ClusteringConfig, MergeClustersConfig, ResolutionSearchConfig

log = logging.getLogger(__name__)

__all__ = [
    "find_resolution",
    "cluster_cells",
    "merge_clusters",
]

# ====================== Clustering Evaluation Helpers ======================
def _auto_select_resolution(
    eval_df: pd.DataFrame,
    strategy: str = "balanced",
) -> float:
    """
    Automatically select the best resolution based on computed metrics.
    
    Strategies:
    1. 'elbow': Find elbow point in n_clusters vs resolution curve
    2. 'peak': Maximum marker abundance
    3. 'balanced': Composite score balancing all metrics
    """
    if strategy == "elbow":
        # Detect elbow using kneedle algorithm
        from kneed import KneeLocator
        
        kn = KneeLocator(
            eval_df['resolution'],
            eval_df['n_clusters'],
            curve='convex',
            direction='increasing'
        )
        return kn.elbow if kn.elbow is not None else eval_df['resolution'].median()
    
    elif strategy == "peak":
        # Simply pick max marker abundance
        if 'marker_abundance' in eval_df.columns:
            idx = eval_df['marker_abundance'].idxmax()
            return eval_df.loc[idx, 'resolution']
        else:
            return eval_df['resolution'].median()
    
    elif strategy == "balanced":
        # Composite score: normalize each metric and compute weighted sum
        score = pd.Series(0.0, index=eval_df.index)
        
        # Normalize each metric to [0, 1]
        if 'silhouette' in eval_df.columns:
            sil_norm = (eval_df['silhouette'] - eval_df['silhouette'].min()) / \
                       (eval_df['silhouette'].max() - eval_df['silhouette'].min() + 1e-8)
            score += 0.3 * sil_norm
        
        if 'marker_abundance' in eval_df.columns:
            ma_norm = (eval_df['marker_abundance'] - eval_df['marker_abundance'].min()) / \
                      (eval_df['marker_abundance'].max() - eval_df['marker_abundance'].min() + 1e-8)
            score += 0.5 * ma_norm  # Higher weight for biological signal
        
        if 'stability' in eval_df.columns:
            # Penalize low stability
            stab_norm = (eval_df['stability'] - eval_df['stability'].min()) / \
                        (eval_df['stability'].max() - eval_df['stability'].min() + 1e-8)
            score += 0.2 * stab_norm
        
        # Pick resolution with max composite score
        idx = score.idxmax()
        return eval_df.loc[idx, 'resolution']
    
    else:
        raise ValueError(f"Unknown selection strategy: {strategy}")
    
    

def _get_marker_abundance(
    adata: AnnData, cluster_key: str, de_method: str, min_log2fc: float, min_pct: float
) -> float:
    """
    Calculate the average number of significant markers per cluster.
    A higher number suggests better biological separability.

    Notes:
    - Uses a temporary rank_genes_groups key to avoid interfering with user's results.
    - Interprets pct columns in their native scale (expects 0–1 fraction here).
    """
    try:
        rank_key = f"rank_genes_{cluster_key}_temp"
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method=de_method,
            key_added=rank_key,
            n_genes=100,  # enough to judge presence
            pts=True,
            use_raw=True,
        )

        markers_df = sc.get.rank_genes_groups_df(adata, key=rank_key, group=None)
        # Harmonize pct column if needed
        if "pct_nz_group" not in markers_df.columns and "pct_nz" in markers_df.columns:
            markers_df = markers_df.rename(columns={"pct_nz": "pct_nz_group"})

        sig_markers = markers_df[
            (markers_df.get("logfoldchanges", 0) > float(min_log2fc))
            & (markers_df.get("pvals_adj", 1.0) < 0.05)
            & (markers_df.get("pct_nz_group", 0) > float(min_pct))
        ]

        # Average marker count per cluster (fill missing with 0)
        if (
            cluster_key not in adata.obs.columns
            or not pd.api.types.is_categorical_dtype(adata.obs[cluster_key])
        ):
            # ensure categorical for stable category order
            adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

        marker_counts = (
            sig_markers.groupby("group").size()
            if not sig_markers.empty
            else pd.Series(dtype=int)
        )
        avg_markers = (
            marker_counts.reindex(
                adata.obs[cluster_key].cat.categories, fill_value=0
            ).mean()
            if len(adata.obs[cluster_key].cat.categories) > 0
            else 0.0
        )

        # Cleanup temporary results
        if rank_key in adata.uns:
            del adata.uns[rank_key]

        return float(avg_markers)

    except Exception as e:
        log.warning(f"Could not compute marker abundance for '{cluster_key}': {e}")
        # Cleanup best effort
        try:
            if rank_key in adata.uns:
                del adata.uns[rank_key]
        except Exception:
            pass
        return 0.0


def _get_clustering_stability(adata: AnnData, key1: str, key2: str) -> float:
    """
    Calculate stability between two clustering results using Normalized Mutual Information (NMI).
    """
    if key1 not in adata.obs or key2 not in adata.obs:
        return np.nan
    labels1 = adata.obs[key1]
    labels2 = adata.obs[key2]
    try:
        return metrics.normalized_mutual_info_score(
            labels1, labels2, average_method="arithmetic"
        )
    except Exception:
        return np.nan


def _print_resolution_guidance(eval_df: pd.DataFrame) -> None:
    """
    Print guidance based on peak silhouette and marker abundance.
    """
    if eval_df.empty:
        return
    # Guard against columns missing
    if "silhouette" in eval_df.columns and not eval_df["silhouette"].isna().all():
        best_silhouette = eval_df.loc[eval_df["silhouette"].idxmax()]
        log.info("--- Resolution Selection Guidance ---")
        log.info(
            f"Peak silhouette at resolution {best_silhouette['resolution']:.2f} "
            f"({int(best_silhouette['n_clusters'])} clusters): better separation."
        )
    if (
        "marker_abundance" in eval_df.columns
        and not eval_df["marker_abundance"].isna().all()
    ):
        best_markers = eval_df.loc[eval_df["marker_abundance"].idxmax()]
        log.info(
            f"Peak marker abundance at resolution {best_markers['resolution']:.2f} "
            f"({int(best_markers['n_clusters'])} clusters): higher biological separability."
        )
    log.info(
        "Suggestion: Choose a resolution around an 'elbow/plateau' where marker abundance is high, "
        "the silhouette is reasonable, and stability does not drop abruptly."
    )


# ====================== Main Functions ======================


def find_resolution(
    adata: AnnData,
    config: Optional[ResolutionSearchConfig] = None,
    auto_select: bool = True,
    selection_strategy: Literal["elbow", "peak", "balanced"] = "balanced",
    **kwargs,
) -> pd.DataFrame:
    """
    Grid search clustering resolutions and provide comprehensive metrics to guide selection.

    Enhancements:
    - Robust neighbors pre-check and informative logging.
    - Silhouette computation with sampling and exception safety.
    - Marker abundance calculation guarded against missing columns and empty DE results.
    - Stable plotting even when some metrics are missing.
    - Full trace saved under .uns['sclucid']['analysis']['clustering']['resolution_search'].

    Args:
        auto_select: If True, automatically recommend a resolution
        selection_strategy: How to pick the best resolution:
            - 'elbow': Elbow point in n_clusters curve
            - 'peak': Maximum marker abundance
            - 'balanced': Balance between silhouette, markers, and stability
    
    Returns:
        Tuple of (evaluation DataFrame, recommended resolution or None)
    """
    if config is None:
        active_config = ResolutionSearchConfig()
    else:
        active_config = config.model_copy()
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    start, end, steps = active_config.resolution_range
    resolutions = np.linspace(start, end, steps).tolist()
    method = active_config.method
    use_rep = active_config.use_rep

    # Check use_rep availability
    if use_rep not in adata.obsm:
        raise ValueError(
            f"Representation '{use_rep}' not found in adata.obsm. "
            "Please compute PCA or the selected embedding first."
        )

    # Ensure neighbors graph exists
    if "neighbors" not in adata.uns:
        log.info(f"Neighbors graph not found. Computing neighbors on '{use_rep}'.")
        sc.pp.neighbors(adata, use_rep=use_rep)

    eval_results = []
    previous_key = None

    log.info(
        f"Searching clustering resolution in [{start:.2f}, {end:.2f}] with {steps} steps..."
    )
    for res in resolutions:
        current_key = f"{method}_res_{res:.3f}"
        log.info(f"Testing resolution: {res:.3f}")

        # Clustering at this resolution
        try:
            if method == "leiden":
                sc.tl.leiden(adata, resolution=float(res), key_added=current_key)
            else:
                sc.tl.louvain(adata, resolution=float(res), key_added=current_key)
        except Exception as e:
            log.error(f"Clustering failed at resolution {res:.3f}: {e}")
            eval_results.append({"resolution": res, "n_clusters": np.nan})
            previous_key = None
            continue

        # Ensure categorical for consistency
        if not pd.api.types.is_categorical_dtype(adata.obs[current_key]):
            adata.obs[current_key] = adata.obs[current_key].astype("category")

        n_clusters = int(adata.obs[current_key].nunique())
        result = {"resolution": res, "n_clusters": n_clusters}

        # Silhouette score with optional sampling
        if active_config.compute_silhouette:
            try:
                from sklearn.metrics import silhouette_score

                labels = adata.obs[current_key].cat.codes.values
                X = adata.obsm[use_rep]
                max_cells = min(20000, adata.n_obs)
                if adata.n_obs > max_cells:
                    rng = np.random.RandomState(42)
                    idx = rng.choice(adata.n_obs, size=max_cells, replace=False)
                    sil = silhouette_score(X[idx], labels[idx], sample_size=None)
                else:
                    sil = silhouette_score(X, labels, sample_size=None)
                result["silhouette"] = float(sil)
            except Exception as e:
                log.warning(f"Silhouette failed at resolution {res:.3f}: {e}")
                result["silhouette"] = np.nan

        # Marker abundance (fix: use min_pct_for_markers from config)
        if active_config.compute_marker_abundance:
            try:
                min_pct = float(getattr(active_config, "min_pct_for_markers", 0.25))
                score = _get_marker_abundance(
                    adata=adata,
                    cluster_key=current_key,
                    de_method=active_config.de_method_for_markers,
                    min_log2fc=float(active_config.min_log2fc_for_markers),
                    min_pct=min_pct,
                )
            except Exception as e:
                log.warning(f"Marker abundance metric failed at {res:.3f}: {e}")
                score = np.nan
            result["marker_abundance"] = float(score) if score is not None else np.nan

        # Stability (NMI vs previous)
        if active_config.compute_stability:
            try:
                result["stability"] = (
                    _get_clustering_stability(adata, previous_key, current_key)
                    if previous_key is not None
                    else np.nan
                )
            except Exception as e:
                log.warning(f"Stability metric failed at {res:.3f}: {e}")
                result["stability"] = np.nan

        eval_results.append(result)
        previous_key = current_key

    eval_df = pd.DataFrame(eval_results)

    recommended_res = None
    if auto_select and not eval_df.empty:
        recommended_res = _auto_select_resolution(eval_df, strategy=selection_strategy)
        log.info(f"🎯 Auto-selected resolution: {recommended_res:.3f}")
    
    # Store recommendation
    adata.uns["sclucid"]["analysis"]["clustering"]["resolution_search"].update({
        "recommended_resolution": recommended_res,
        "selection_strategy": selection_strategy,
    })

    # Plotting
    if active_config.plot and not eval_df.empty:
        metrics_to_plot = [
            m
            for m in ["silhouette", "marker_abundance", "stability"]
            if m in eval_df.columns
        ]
        n_plots = 1 + len(metrics_to_plot)

        fig, axes = plt.subplots(n_plots, 1, figsize=(10, 4 * n_plots), sharex=True)
        if n_plots == 1:
            axes = [axes]
        fig.suptitle("Clustering Resolution Optimization Guide", fontsize=16)

        # n_clusters
        ax = axes[0]
        try:
            sns.lineplot(
                data=eval_df, x="resolution", y="n_clusters", marker="o", ax=ax
            )
        except Exception as e:
            ax.text(0.5, 0.5, f"Plot failed: {e}", ha="center", va="center")
        ax.set_title("Number of Clusters vs. Resolution")
        ax.set_ylabel("Count")
        ax.grid(True, linestyle="--")

        # Other metrics
        for i, metric in enumerate(metrics_to_plot):
            ax = axes[i + 1]
            try:
                sns.lineplot(data=eval_df, x="resolution", y=metric, marker="o", ax=ax)
                ax.set_title(f"{metric.replace('_', ' ').title()} vs. Resolution")
                ax.set_ylabel("Score")
                if metric == "silhouette":
                    ax.text(
                        0.02,
                        0.05,
                        "Note: Higher values may favor fewer clusters.",
                        transform=ax.transAxes,
                        fontsize=9,
                        style="italic",
                    )
            except Exception as e:
                ax.text(0.5, 0.5, f"Plot failed: {e}", ha="center", va="center")
            ax.grid(True, linestyle="--")

        axes[-1].set_xlabel("Resolution")
        plt.tight_layout(rect=[0, 0.03, 1, 0.96])

        if active_config.save_dir:
            save_path = Path(active_config.save_dir)
            save_path.mkdir(parents=True, exist_ok=True)
            figure_path = save_path / "resolution_search_guide.png"
            plt.savefig(figure_path, dpi=300)
            log.info(f"Saved resolution search plot to {figure_path}")

        plt.show()

    # Guidance message
    if not eval_df.empty:
        _print_resolution_guidance(eval_df)

    return eval_df, recommended_res


def cluster_cells(
    adata: AnnData,
    config: Optional[ClusteringConfig] = None,
    **kwargs,
) -> AnnData:
    """
    Perform clustering using Leiden, Louvain, KMeans, or HDBSCAN.

    Enhancements:
    - Safer neighbor graph handling for Leiden/Louvain with explicit random_state.
    - KMeans/HDBSCAN robust behavior and informative errors.
    - Result trace saved with full config and scanpy/sklearn/hdbscan versions.
    """
    if config is None:
        active_config = ClusteringConfig()
    else:
        active_config = config.model_copy()
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    method = active_config.method
    use_rep = active_config.use_rep

    # Representation check
    if use_rep not in adata.obsm:
        raise ValueError(
            f"Representation '{use_rep}' not found in adata.obsm. Please compute PCA or your representation first."
        )

    # Ensure neighbors for graph-based methods
    if method in ["leiden", "louvain"]:
        if "neighbors" not in adata.uns:
            log.info(
                f"Neighbors graph not found for {method}. Computing on '{use_rep}'."
            )
            sc.pp.neighbors(
                adata, use_rep=use_rep, random_state=active_config.random_state
            )

    key_added = active_config.key_added or f"{method}_clusters"
    log.info(f"Running {method} clustering...")

    try:
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

            if active_config.n_clusters is None or int(active_config.n_clusters) < 2:
                raise ValueError("n_clusters must be specified and >= 2 for KMeans.")
            X = adata.obsm[use_rep]
            kmeans = KMeans(
                n_clusters=int(active_config.n_clusters),
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
            adata.obs[key_added] = pd.Categorical(
                [str(l) if l != -1 else "Noise" for l in labels]
            )
        else:
            raise ValueError(f"Unknown clustering method: {method}")
    except Exception as e:
        log.error(f"Clustering failed: {e}")
        raise

    # Ensure categorical type
    if not pd.api.types.is_categorical_dtype(adata.obs[key_added]):
        adata.obs[key_added] = adata.obs[key_added].astype("category")

    n_clusters = int(adata.obs[key_added].nunique())
    log.info(
        f"Clustering ({method}) finished: {n_clusters} clusters in obs['{key_added}']"
    )

    # Trace
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "clustering", {}
    )
    trace = {
        "config": active_config.to_dict(),
        "n_clusters": n_clusters,
        "scanpy_version": getattr(sc, "__version__", "unknown"),
    }
    try:
        import sklearn

        trace["sklearn_version"] = getattr(sklearn, "__version__", "unknown")
    except Exception:
        pass
    try:
        import hdbscan as _hb

        trace["hdbscan_version"] = getattr(_hb, "__version__", "unknown")
    except Exception:
        pass
    adata.uns["sclucid"]["analysis"]["clustering"][key_added] = sanitize_for_hdf5(trace)

    # Optional UMAP plot
    if active_config.plot:
        if "X_umap" not in adata.obsm:
            try:
                sc.tl.umap(adata)
            except Exception as e:
                log.warning(f"UMAP computation failed: {e}")
        if "X_umap" in adata.obsm:
            color_key = f"{key_added}_colors"
            if color_key in adata.uns:
                del adata.uns[color_key]
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

    Enhancements:
    - Safer handling when DE results are empty or columns missing.
    - Proper cleanup of temporary rank_genes keys.
    - Robust expression correlation with dense guard.
    """
    if config is None:
        active_config = MergeClustersConfig()
    else:
        active_config = config.model_copy()
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)

    cluster_key = active_config.cluster_key
    threshold = active_config.similarity_threshold
    method = active_config.method
    key_added = active_config.key_added or f"{cluster_key}_merged"

    if cluster_key not in adata.obs:
        raise ValueError(f"cluster_key '{cluster_key}' not found in adata.obs")

    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
    original_clusters = list(adata.obs[cluster_key].cat.categories)
    n_original = len(original_clusters)

    sim_matrix = pd.DataFrame(
        np.eye(n_original, dtype=float),
        index=original_clusters,
        columns=original_clusters,
    )

    log.info(f"Calculating similarity matrix using '{method}' method...")
    if method == "marker_overlap":
        rank_key = f"rank_genes_{cluster_key}_for_merge"
        try:
            sc.tl.rank_genes_groups(
                adata,
                groupby=cluster_key,
                method=active_config.de_method_for_markers,
                key_added=rank_key,
            )
            markers = sc.get.rank_genes_groups_df(adata, key=rank_key, group=None)
            if (
                markers.empty
                or "names" not in markers.columns
                or "group" not in markers.columns
            ):
                log.warning(
                    "Rank genes results appear empty or malformed; cannot compute marker overlap."
                )
            else:
                top_markers = {}
                for g, df in markers.groupby("group"):
                    df2 = df.copy()
                    if "pvals_adj" in df2.columns and "logfoldchanges" in df2.columns:
                        df2 = df2.query("pvals_adj < 0.05 & logfoldchanges > 0.5")
                    top_markers[str(g)] = set(
                        df2["names"].head(50).astype(str).tolist()
                    )
                for i, c1 in enumerate(original_clusters):
                    s1 = top_markers.get(str(c1), set())
                    for j, c2 in enumerate(original_clusters[i + 1 :], i + 1):
                        s2 = top_markers.get(str(c2), set())
                        if not s1 or not s2:
                            continue
                        sim = len(s1 & s2) / len(s1 | s2)
                        sim_matrix.loc[c1, c2] = sim_matrix.loc[c2, c1] = float(sim)
        finally:
            if rank_key in adata.uns:
                del adata.uns[rank_key]
    elif method == "expression_correlation":
        # Compute mean expression per cluster, then Pearson correlation
        means = []
        for c in original_clusters:
            idx = (adata.obs[cluster_key] == c).values
            if not idx.any():
                means.append(np.zeros((adata.n_vars,), dtype=float))
                continue
            Xc = adata[idx].X
            if hasattr(Xc, "toarray"):
                Xc = Xc.toarray()
            vec = np.asarray(Xc.mean(axis=0)).ravel()
            means.append(vec)
        mean_profiles = (
            np.vstack(means) if len(means) > 0 else np.zeros((0, adata.n_vars))
        )
        if mean_profiles.shape[0] >= 2:
            corr_matrix = np.corrcoef(mean_profiles)
            sim_matrix = pd.DataFrame(
                corr_matrix, index=original_clusters, columns=original_clusters
            )
        else:
            log.warning(
                "Not enough clusters to compute correlation; fallback to identity."
            )
    else:
        raise ValueError(f"Unknown merge method: {method}")

    # Graph-based merging
    G = nx.from_pandas_adjacency(sim_matrix > threshold)
    components = list(nx.connected_components(G))

    mapping = {}
    new_cluster_names = []
    for comp_nodes in components:
        comp_nodes = sorted(list(comp_nodes), key=lambda x: (str(x).isdigit(), str(x)))
        new_name = "_".join(comp_nodes)
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
    adata.uns["sclucid"]["analysis"]["clustering"][f"{key_added}_params"] = (
        sanitize_for_hdf5(
            {
                "source_clusters": cluster_key,
                "method": method,
                "similarity_threshold": threshold,
                "original_clusters": n_original,
                "merged_clusters": n_merged,
                "mapping": mapping,
                "config": active_config.to_dict(),
                "scanpy_version": getattr(sc, "__version__", "unknown"),
            }
        )
    )
    return adata
