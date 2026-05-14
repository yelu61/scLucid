"""
Clustering functions for single-cell RNA-seq analysis.

This module provides:
- Unsupervised clustering (Leiden, Louvain, K-means, HDBSCAN)
- Practical resolution review with marker and composition evidence.
- Cluster merging based on marker overlap or expression correlation.
- Standardized config dataclasses, logging, and results traceability.
"""

import logging
from importlib.metadata import version
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Sequence

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ..base_config import apply_config_overrides
from ..utils import sanitize_for_hdf5
from .config import ClusteringConfig, MergeClustersConfig

log = logging.getLogger(__name__)

__all__ = [
    "run_clustering_review",
    "cluster_cells",
    "merge_clusters",
]


_CLUSTERING_NOISE_EXACT_GENES = {
    "MALAT1": "housekeeping",
    "NEAT1": "housekeeping",
    "XIST": "housekeeping",
}

_CLUSTERING_NOISE_PREFIXES = {
    "ribosomal": ("RPL", "RPS", "MRPL", "MRPS"),
    "mitochondrial": ("MT-",),
    "stress": ("HSP", "HSPA", "HSPB", "DNAJ", "FOS", "JUN", "ATF3", "IER", "DDIT3"),
}


def _classify_marker_noise(gene: Any) -> Optional[str]:
    """Classify common low-information marker genes for clustering review."""
    if pd.isna(gene):
        return None
    gene_upper = str(gene).upper()
    if gene_upper in _CLUSTERING_NOISE_EXACT_GENES:
        return _CLUSTERING_NOISE_EXACT_GENES[gene_upper]
    for category, prefixes in _CLUSTERING_NOISE_PREFIXES.items():
        if any(gene_upper.startswith(prefix) for prefix in prefixes):
            return category
    return None


def _format_top_distribution(values: pd.Series, n: int = 3) -> str:
    """Return a compact top-category distribution string."""
    if values.empty:
        return ""
    counts = values.astype(str).value_counts(normalize=True).head(n)
    return ", ".join(f"{name}:{frac:.2f}" for name, frac in counts.items())


# ====================== Clustering Evaluation Helpers ======================
def run_clustering_review(
    adata: AnnData,
    resolutions: Optional[Sequence[float]] = None,
    *,
    method: Literal["leiden", "louvain"] = "leiden",
    use_rep: str = "X_pca",
    key_prefix: Optional[str] = None,
    random_state: int = 42,
    de_method: str = "wilcoxon",
    use_raw: bool = True,
    n_top_markers: int = 15,
    min_cluster_fraction: float = 0.005,
    min_cluster_cells: int = 20,
    min_informative_markers: int = 3,
    marker_quality_weight: float = 0.55,
    size_balance_weight: float = 0.25,
    sample_mixing_weight: float = 0.20,
    sample_col: Optional[str] = None,
    batch_col: Optional[str] = None,
    copy: bool = False,
) -> pd.DataFrame:
    """
    Review practical clustering resolutions with marker and composition evidence.

    This lightweight path is intended for day-to-day annotation work. It does
    not claim a mathematically optimal resolution; it summarizes whether each
    candidate resolution yields interpretable clusters.
    """
    if copy:
        adata = adata.copy()
    if resolutions is None:
        if adata.n_obs < 5_000:
            resolutions = [0.3, 0.5, 0.8]
        elif adata.n_obs < 50_000:
            resolutions = [0.5, 0.6, 0.8, 1.0]
        else:
            resolutions = [0.6, 0.8, 1.0, 1.2]
    resolutions = [float(r) for r in resolutions]
    if not resolutions:
        raise ValueError("resolutions must contain at least one value.")
    if use_rep not in adata.obsm:
        raise ValueError(
            f"Representation '{use_rep}' not found in adata.obsm. "
            "Please compute PCA or the selected embedding first."
        )
    if "neighbors" not in adata.uns:
        log.info(f"Neighbors graph not found. Computing neighbors on '{use_rep}'.")
        sc.pp.neighbors(adata, use_rep=use_rep, random_state=random_state)

    key_prefix = key_prefix or f"{method}_review"
    review_rows: List[Dict[str, Any]] = []
    cluster_rows: List[Dict[str, Any]] = []
    marker_frames: List[pd.DataFrame] = []

    for resolution in resolutions:
        cluster_key = f"{key_prefix}_{resolution:g}"
        log.info(f"Reviewing {method} resolution {resolution:g} -> obs['{cluster_key}']")
        if method == "leiden":
            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=cluster_key,
                random_state=random_state,
            )
        elif method == "louvain":
            sc.tl.louvain(
                adata,
                resolution=resolution,
                key_added=cluster_key,
                random_state=random_state,
            )
        else:
            raise ValueError("method must be 'leiden' or 'louvain'.")

        if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
            adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

        cluster_series = adata.obs[cluster_key].astype(str)
        counts = cluster_series.value_counts()
        n_clusters = int(counts.shape[0])
        small_threshold = max(
            int(min_cluster_cells), int(np.ceil(min_cluster_fraction * adata.n_obs))
        )
        n_small_clusters = int((counts < small_threshold).sum())
        min_cluster_size = int(counts.min()) if not counts.empty else 0
        median_cluster_size = float(counts.median()) if not counts.empty else 0.0

        rank_key = f"rank_genes_{cluster_key}"
        markers_df = pd.DataFrame()
        try:
            sc.tl.rank_genes_groups(
                adata,
                groupby=cluster_key,
                method=de_method,
                use_raw=use_raw and adata.raw is not None,
                pts=True,
                key_added=rank_key,
            )
            markers_df = sc.get.rank_genes_groups_df(adata, key=rank_key, group=None)
            if not markers_df.empty:
                markers_df["resolution"] = resolution
                markers_df["cluster_key"] = cluster_key
                markers_df["noise_category"] = markers_df["names"].map(_classify_marker_noise)
                marker_frames.append(markers_df)
        except Exception as exc:
            log.warning(f"Marker review failed for resolution {resolution:g}: {exc}")

        marker_quality_scores: List[float] = []
        noise_fractions: List[float] = []
        informative_counts: List[int] = []
        sample_dominance: List[float] = []
        batch_dominance: List[float] = []

        for cluster, n_cells in counts.items():
            if not markers_df.empty and "group" in markers_df.columns:
                cluster_marker_df = markers_df[
                    markers_df["group"].astype(str).to_numpy() == str(cluster)
                ].copy()
            else:
                cluster_marker_df = pd.DataFrame()
            if not cluster_marker_df.empty and "scores" in cluster_marker_df.columns:
                cluster_marker_df = cluster_marker_df.sort_values("scores", ascending=False)
            top_markers = cluster_marker_df.head(n_top_markers)
            noise_fraction = (
                float(top_markers["noise_category"].notna().mean())
                if "noise_category" in top_markers.columns and not top_markers.empty
                else np.nan
            )
            informative_count = (
                int(top_markers["noise_category"].isna().sum())
                if "noise_category" in top_markers.columns and not top_markers.empty
                else 0
            )
            informative_counts.append(informative_count)
            if pd.notna(noise_fraction):
                noise_fractions.append(float(noise_fraction))
            marker_quality = min(1.0, informative_count / max(1, min_informative_markers))
            marker_quality_scores.append(marker_quality)

            mask = cluster_series == str(cluster)
            obs_subset = adata.obs.loc[mask]
            sample_top = ""
            batch_top = ""
            sample_dom = np.nan
            batch_dom = np.nan
            if sample_col and sample_col in obs_subset.columns:
                sample_counts = obs_subset[sample_col].astype(str).value_counts(normalize=True)
                sample_dom = float(sample_counts.iloc[0]) if not sample_counts.empty else np.nan
                sample_dominance.append(sample_dom)
                sample_top = _format_top_distribution(obs_subset[sample_col])
            if batch_col and batch_col in obs_subset.columns:
                batch_counts = obs_subset[batch_col].astype(str).value_counts(normalize=True)
                batch_dom = float(batch_counts.iloc[0]) if not batch_counts.empty else np.nan
                batch_dominance.append(batch_dom)
                batch_top = _format_top_distribution(obs_subset[batch_col])

            cluster_rows.append(
                {
                    "resolution": resolution,
                    "cluster_key": cluster_key,
                    "cluster": str(cluster),
                    "n_cells": int(n_cells),
                    "pct_cells": float(n_cells / max(1, adata.n_obs)),
                    "n_informative_top_markers": informative_count,
                    "noise_marker_fraction": noise_fraction,
                    "top_markers": ", ".join(
                        top_markers.get("names", pd.Series(dtype=str)).astype(str).head(8)
                    ),
                    "top_samples": sample_top,
                    "top_batches": batch_top,
                    "sample_dominance": sample_dom,
                    "batch_dominance": batch_dom,
                }
            )

        mean_marker_quality = (
            float(np.mean(marker_quality_scores)) if marker_quality_scores else 0.0
        )
        mean_noise_fraction = float(np.nanmean(noise_fractions)) if noise_fractions else np.nan
        marker_poor_fraction = (
            float(np.mean([x < min_informative_markers for x in informative_counts]))
            if informative_counts
            else np.nan
        )
        small_cluster_fraction = float(n_small_clusters / max(1, n_clusters))
        size_balance_score = float(np.clip(1.0 - small_cluster_fraction, 0.0, 1.0))

        sample_mixing_score = 1.0
        if sample_dominance:
            sample_mixing_score = float(np.clip(1.0 - np.nanmean(sample_dominance), 0.0, 1.0))
        batch_mixing_score = np.nan
        if batch_dominance:
            batch_mixing_score = float(np.clip(1.0 - np.nanmean(batch_dominance), 0.0, 1.0))

        interpretability_score = (
            marker_quality_weight * mean_marker_quality
            + size_balance_weight * size_balance_score
            + sample_mixing_weight * sample_mixing_score
        )
        warnings: List[str] = []
        if n_small_clusters:
            warnings.append("small_clusters")
        if pd.notna(marker_poor_fraction) and marker_poor_fraction > 0.25:
            warnings.append("marker_poor_clusters")
        if pd.notna(mean_noise_fraction) and mean_noise_fraction > 0.4:
            warnings.append("noise_marker_dominance")
        if sample_dominance and float(np.nanmean(sample_dominance)) > 0.8:
            warnings.append("sample_specific_clusters")

        review_rows.append(
            {
                "resolution": resolution,
                "cluster_key": cluster_key,
                "n_clusters": n_clusters,
                "min_cluster_size": min_cluster_size,
                "median_cluster_size": median_cluster_size,
                "n_small_clusters": n_small_clusters,
                "small_cluster_fraction": small_cluster_fraction,
                "mean_informative_top_markers": (
                    float(np.mean(informative_counts)) if informative_counts else 0.0
                ),
                "marker_poor_cluster_fraction": marker_poor_fraction,
                "mean_noise_marker_fraction": mean_noise_fraction,
                "marker_quality_score": mean_marker_quality,
                "size_balance_score": size_balance_score,
                "sample_mixing_score": sample_mixing_score,
                "batch_mixing_score": batch_mixing_score,
                "interpretability_score": float(interpretability_score),
                "warnings": ",".join(warnings),
            }
        )

    review_df = pd.DataFrame(review_rows)
    cluster_review_df = pd.DataFrame(cluster_rows)
    marker_review_df = (
        pd.concat(marker_frames, ignore_index=True) if marker_frames else pd.DataFrame()
    )

    recommended_resolution = None
    recommended_key = None
    rationale = ""
    if not review_df.empty:
        ranked = review_df.sort_values(
            ["interpretability_score", "marker_quality_score", "n_small_clusters"],
            ascending=[False, False, True],
        )
        best = ranked.iloc[0]
        recommended_resolution = float(best["resolution"])
        recommended_key = str(best["cluster_key"])
        rationale = (
            f"Selected {recommended_resolution:g}: "
            f"{int(best['n_clusters'])} clusters, "
            f"marker quality {best['marker_quality_score']:.2f}, "
            f"small-cluster fraction {best['small_cluster_fraction']:.2f}."
        )

    clustering_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("clustering", {})
    )
    clustering_ns["clustering_review"] = review_df
    clustering_ns["clustering_review_by_cluster"] = cluster_review_df
    clustering_ns["clustering_review_markers"] = marker_review_df
    clustering_ns["clustering_review_summary"] = sanitize_for_hdf5(
        {
            "method": method,
            "use_rep": use_rep,
            "resolutions": resolutions,
            "recommended_resolution": recommended_resolution,
            "recommended_cluster_key": recommended_key,
            "rationale": rationale,
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
        }
    )
    return review_df


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
        active_config = apply_config_overrides(config, **kwargs)

    method = active_config.method
    use_rep = active_config.use_rep

    # Representation check
    if use_rep not in adata.obsm:
        raise ValueError(
            f"[analysis] Clustering failed: representation '{use_rep}' not found in adata.obsm. "
            "Please compute PCA or the selected embedding first."
        )

    # Ensure neighbors for graph-based methods
    if method in ["leiden", "louvain"]:
        if "neighbors" not in adata.uns:
            log.info(f"Neighbors graph not found for {method}. Computing on '{use_rep}'.")
            sc.pp.neighbors(adata, use_rep=use_rep, random_state=active_config.random_state)

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
                raise ValueError(
                    "[analysis] KMeans clustering failed: n_clusters must be specified and >= 2."
                )
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
            adata.obs[key_added] = pd.Categorical([str(l) if l != -1 else "Noise" for l in labels])
        else:
            raise ValueError(
                f"[analysis] Unknown clustering method '{method}'. "
                "Expected one of: leiden, louvain, kmeans, hdbscan."
            )
    except Exception as e:
        raise RuntimeError(
            f"[analysis] Clustering failed: {e}. Check input data and configuration."
        ) from e

    # Ensure categorical type
    if not pd.api.types.is_categorical_dtype(adata.obs[key_added]):
        adata.obs[key_added] = adata.obs[key_added].astype("category")

    n_clusters = int(adata.obs[key_added].nunique())
    log.info(f"Clustering ({method}) finished: {n_clusters} clusters in obs['{key_added}']")

    # Trace
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("clustering", {})
    trace = {
        "config": active_config.to_dict(),
        "n_clusters": n_clusters,
        "scanpy_version": version("scanpy"),
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
        active_config = apply_config_overrides(config, **kwargs)

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
            if markers.empty or "names" not in markers.columns or "group" not in markers.columns:
                log.warning(
                    "Rank genes results appear empty or malformed; cannot compute marker overlap."
                )
            else:
                top_markers = {}
                for g, df in markers.groupby("group"):
                    df2 = df.copy()
                    if "pvals_adj" in df2.columns and "logfoldchanges" in df2.columns:
                        df2 = df2.query("pvals_adj < 0.05 & logfoldchanges > 0.5")
                    top_markers[str(g)] = set(df2["names"].head(50).astype(str).tolist())
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
        mean_profiles = np.vstack(means) if len(means) > 0 else np.zeros((0, adata.n_vars))
        if mean_profiles.shape[0] >= 2:
            corr_matrix = np.corrcoef(mean_profiles)
            sim_matrix = pd.DataFrame(
                corr_matrix, index=original_clusters, columns=original_clusters
            )
        else:
            log.warning("Not enough clusters to compute correlation; fallback to identity.")
    else:
        raise ValueError(
            f"[analysis] Unknown merge method '{method}'. "
            "Expected one of: marker_overlap, expression_correlation."
        )

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

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("clustering", {})
    adata.uns["sclucid"]["analysis"]["clustering"][f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "source_clusters": cluster_key,
            "method": method,
            "similarity_threshold": threshold,
            "original_clusters": n_original,
            "merged_clusters": n_merged,
            "mapping": mapping,
            "config": active_config.to_dict(),
            "scanpy_version": version("scanpy"),
        }
    )
    return adata
