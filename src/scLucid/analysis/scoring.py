import logging
from typing import Dict, List, Literal, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from scipy.stats import mannwhitneyu, ttest_ind

from ..utils import sanitize_for_hdf5

log = logging.getLogger(__name__)

__all__ = [
    "score_by_gene_sets",
    "compare_scores",
    "plot_score_comparison",
    "batch_compare_scores",
    "batch_plot_score_comparison",
]


def _ensure_scoring_namespace(adata: AnnData) -> dict:
    """
    Ensure the scoring namespace exists in adata.uns and return it.
    """
    return (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("scoring", {})
    )


def _cohens_d(x: np.ndarray, y: np.ndarray) -> Optional[float]:
    """
    Compute Cohen's d for two independent samples (unequal n).
    Returns None if not computable.
    """
    try:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        nx, ny = len(x), len(y)
        if nx < 2 or ny < 2:
            return None
        # Pooled standard deviation
        sx, sy = x.std(ddof=1), y.std(ddof=1)
        # Guard against zero variance
        if sx == 0 and sy == 0:
            return 0.0
        sp = np.sqrt(((nx - 1) * sx**2 + (ny - 1) * sy**2) / (nx + ny - 2))
        if sp == 0:
            return None
        return (x.mean() - y.mean()) / sp
    except Exception:
        return None


def score_by_gene_sets(
    adata: AnnData,
    gene_sets: Dict[str, List[str]],
    layer: Optional[str] = "log1p_norm",
    use_raw: bool = False,
    ctrl_size: int = 50,
    score_name_suffix: str = "_score",
    preserve_missing: bool = True,
    min_genes_required: int = 1,
    **kwargs,
) -> AnnData:
    """
    Score cells for multiple gene sets (e.g., pathways, cell types) and add columns to adata.obs.

    Parameters
    - gene_sets: dict of {set_name: [genes]}
    - layer: expression layer when use_raw=False; ignored if use_raw=True
    - use_raw: if True, score on adata.raw (must exist)
    - ctrl_size: control gene set size for scanpy.tl.score_genes
    - score_name_suffix: suffix for new obs columns
    - preserve_missing: if True, allow partial matches (use only found genes);
                        if False, skip sets where not all genes are found
    - min_genes_required: minimum number of found genes to attempt scoring
    - kwargs: forwarded to scanpy.tl.score_genes

    Behavior
    - Records summary and parameters under adata.uns['sclucid']['analysis']['scoring']['gene_set_scoring'].
    - Logs per-set found/total genes and whether it was scored.
    """
    ns = _ensure_scoring_namespace(adata)

    if use_raw:
        if adata.raw is None:
            raise ValueError("adata.raw is not set, but use_raw=True.")
        source_adata = adata.raw
        target_layer = None
    else:
        source_adata = adata
        target_layer = layer
        if target_layer is None:
            log.warning("layer=None and use_raw=False; using adata.X for scoring.")
        elif target_layer not in adata.layers:
            raise ValueError(f"Layer '{target_layer}' not found in adata.layers.")

    total_sets = len(gene_sets)
    scored_count = 0
    skipped_sets: List[str] = []
    per_set_stats: Dict[str, Dict[str, int]] = {}

    log.info(f"Scoring cells for {total_sets} gene sets (use_raw={use_raw}, layer={target_layer})...")

    for set_name, genes in gene_sets.items():
        genes = [g for g in genes if isinstance(g, str) and len(g) > 0]
        if len(genes) == 0:
            log.warning(f"Gene set '{set_name}' is empty. Skipping.")
            skipped_sets.append(set_name)
            per_set_stats[set_name] = {"n_input": 0, "n_found": 0, "scored": 0}
            continue

        found = [g for g in genes if g in source_adata.var_names]
        n_input, n_found = len(genes), len(found)
        per_set_stats[set_name] = {"n_input": n_input, "n_found": n_found, "scored": 0}

        if n_found < min_genes_required:
            log.warning(
                f"Gene set '{set_name}': found {n_found}/{n_input} (<{min_genes_required}). Skipping."
            )
            skipped_sets.append(set_name)
            continue

        if not preserve_missing and n_found < n_input:
            log.warning(
                f"Gene set '{set_name}': partial match {n_found}/{n_input} and preserve_missing=False. Skipping."
            )
            skipped_sets.append(set_name)
            continue

        score_name = f"{set_name}{score_name_suffix}"
        try:
            sc.tl.score_genes(
                adata,
                found,
                score_name=score_name,
                use_raw=use_raw,
                layer=target_layer,
                ctrl_size=ctrl_size,
                **kwargs,
            )
            per_set_stats[set_name]["scored"] = 1
            scored_count += 1
        except Exception as e:
            log.warning(f"Failed to score set '{set_name}': {e}")
            skipped_sets.append(set_name)

    ns["gene_set_scoring"] = sanitize_for_hdf5({
        "n_sets_input": total_sets,
        "n_sets_scored": scored_count,
        "n_sets_skipped": len(skipped_sets),
        "skipped_sets": skipped_sets,
        "per_set_stats": per_set_stats,
        "params": {
            "use_raw": use_raw,
            "layer": layer,
            "ctrl_size": ctrl_size,
            "score_name_suffix": score_name_suffix,
            "preserve_missing": preserve_missing,
            "min_genes_required": min_genes_required,
            "scanpy_version": getattr(sc, "__version__", "unknown"),
        },
    })
    log.info(f"Completed scoring: {scored_count}/{total_sets} sets scored, {len(skipped_sets)} skipped.")
    return adata


def compare_scores(
    adata: AnnData,
    score_key: str,
    groupby: str,
    group1: str,
    group2: Optional[str] = "rest",
    method: Literal["ttest", "wilcoxon"] = "wilcoxon",
) -> pd.DataFrame:
    """
    Differential comparison of a continuous score across two groups.

    Parameters
    - score_key: column in adata.obs with continuous scores (e.g., 'Tcell_score')
    - groupby: grouping column in adata.obs (categorical recommended)
    - group1: target group name
    - group2: comparison group name, or 'rest' to compare group1 vs all others
    - method: 'wilcoxon' (Mann–Whitney U) or 'ttest' (Welch's t-test)

    Returns
    - Single-row DataFrame with statistics. Also stored under
      adata.uns['sclucid']['analysis']['scoring'][f'compare_{score_key}_{group1}_vs_{group2}'].
    """
    ns = _ensure_scoring_namespace(adata)

    if score_key not in adata.obs:
        raise ValueError(f"Score key '{score_key}' not found in adata.obs.")
    if groupby not in adata.obs:
        raise ValueError(f"groupby '{groupby}' not found in adata.obs.")

    if group2 != "rest":
        # Ensure both groups exist
        groups_available = set(adata.obs[groupby].astype(str).unique())
        if str(group1) not in groups_available or str(group2) not in groups_available:
            raise ValueError(
                f"Groups not found in '{groupby}'. Available: {sorted(groups_available)}"
            )

    s = adata.obs[[groupby, score_key]].dropna()
    if s.empty:
        log.warning("No data available after dropping NaN.")
        return pd.DataFrame()

    mask1 = s[groupby].astype(str) == str(group1)
    if group2 == "rest":
        mask2 = s[groupby].astype(str) != str(group1)
    else:
        mask2 = s[groupby].astype(str) == str(group2)

    scores1 = s.loc[mask1, score_key].astype(float)
    scores2 = s.loc[mask2, score_key].astype(float)

    n1, n2 = len(scores1), len(scores2)
    if n1 < 3 or n2 < 3:
        log.warning("Not enough data points in one or both groups for statistical test (n<3).")
        return pd.DataFrame()

    if method == "wilcoxon":
        stat, pval = mannwhitneyu(scores1, scores2, alternative="two-sided")
    elif method == "ttest":
        stat, pval = ttest_ind(scores1, scores2, equal_var=False)
    else:
        raise ValueError("Method must be 'ttest' or 'wilcoxon'.")

    mean1, mean2 = float(scores1.mean()), float(scores2.mean())
    effect_size = mean1 - mean2
    direction = "group1>group2" if effect_size > 0 else ("group1<group2" if effect_size < 0 else "equal")
    d_cohen = _cohens_d(scores1.values, scores2.values)

    results = pd.DataFrame(
        {
            "score": [score_key],
            "groupby": [groupby],
            "group1": [group1],
            "group2": [group2],
            "method": [method],
            "statistic": [float(stat)],
            "pvalue": [float(pval)],
            "mean_group1": [mean1],
            "mean_group2": [mean2],
            "effect_size (mean_diff)": [float(effect_size)],
            "direction": [direction],
            "cohens_d": [d_cohen if d_cohen is None else float(d_cohen)],
            "n_cells_group1": [int(n1)],
            "n_cells_group2": [int(n2)],
        }
    )

    ns[f"compare_{score_key}_{group1}_vs_{group2}"] = sanitize_for_hdf5(results)
    return results


def plot_score_comparison(
    adata: AnnData,
    score_key: str,
    groupby: str,
    groups_to_compare: Optional[Sequence[str]] = None,
    plot_type: Literal["violin", "boxplot"] = "violin",
    order: Optional[Sequence[str]] = None,
    palette: Optional[str] = None,
    figsize: Tuple[float, float] = (8, 6),
    show: bool = True,
    save_path: Optional[str] = None,
    **kwargs,
) -> plt.Axes:
    """
    Plot comparison of a score across groups.

    Parameters
    - groups_to_compare: optional subset of groups to display (keeps given order)
    - order: explicit order of x-axis categories; overrides groups_to_compare order if provided
    - palette: seaborn palette name
    - figsize: figure size
    - show: display the plot
    - save_path: optional path to save figure
    - kwargs: forwarded to seaborn plotting function
    """
    ns = _ensure_scoring_namespace(adata)

    if score_key not in adata.obs or groupby not in adata.obs:
        raise ValueError(f"'{score_key}' or '{groupby}' not found in adata.obs.")

    plot_df = adata.obs[[groupby, score_key]].copy().dropna()
    if plot_df.empty:
        raise ValueError("No data available to plot (after dropping NaN).")

    # Filter to groups_to_compare
    if groups_to_compare:
        groups_set = set(map(str, groups_to_compare))
        plot_df = plot_df[plot_df[groupby].astype(str).isin(groups_set)]
        # Preserve the given order if not overridden by 'order'
        if order is None:
            order = list(map(str, groups_to_compare))

    # If order is not given, try to respect categorical order
    if order is None and pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        order = list(adata.obs[groupby].cat.categories)

    plt.figure(figsize=figsize)
    common_kwargs = dict(x=groupby, y=score_key, data=plot_df, order=order, palette=palette)
    if plot_type == "violin":
        ax = sns.violinplot(**common_kwargs, cut=0, inner="quartile", **kwargs)
    elif plot_type == "boxplot":
        ax = sns.boxplot(**common_kwargs, **kwargs)
        try:
            sns.swarmplot(x=groupby, y=score_key, data=plot_df, color=".25", size=3, ax=ax)
        except Exception:
            # Swarm may fail with many points; silently ignore
            pass
    else:
        raise ValueError("plot_type must be 'violin' or 'boxplot'.")

    ax.set_title(f"Comparison of '{score_key}' across '{groupby}'")
    ax.set_ylabel(score_key)
    ax.set_xlabel(groupby)
    ax.tick_params(axis="x", rotation=45)
    plt.tight_layout()

    ns[f"{score_key}_{groupby}_{plot_type}_plot"] = sanitize_for_hdf5({
        "groups": list(map(str, groups_to_compare)) if groups_to_compare else None,
        "order": list(order) if order is not None else None,
        "plot_type": plot_type,
        "palette": palette,
        "figsize": figsize,
    })

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved score comparison plot to {save_path}")
    if show:
        plt.show()
    return ax


def batch_compare_scores(
    adata: AnnData,
    score_keys: List[str],
    groupby: str,
    group_pairs: Optional[List[Tuple[str, str]]] = None,
    method: Literal["ttest", "wilcoxon"] = "wilcoxon",
) -> pd.DataFrame:
    """
    Batch differential comparison for multiple scores and group pairs.

    Parameters
    - score_keys: list of obs columns containing scores
    - groupby: grouping column
    - group_pairs: list of (group1, group2) pairs; if None, default to all groups vs 'rest'
    - method: 'wilcoxon' or 'ttest'

    Returns
    - Concatenated DataFrame of results. Also stored under
      adata.uns['sclucid']['analysis']['scoring']['batch_compare_results'].
    """
    ns = _ensure_scoring_namespace(adata)

    available_scores = [s for s in score_keys if s in adata.obs.columns]
    missing_scores = [s for s in score_keys if s not in adata.obs.columns]
    if missing_scores:
        log.warning(f"Missing score columns skipped: {missing_scores}")
    if len(available_scores) == 0:
        log.warning("No valid score columns to compare. Returning empty DataFrame.")
        return pd.DataFrame()

    if group_pairs is None:
        groups = list(map(str, adata.obs[groupby].astype(str).unique()))
        group_pairs = [(g, "rest") for g in groups]

    results = []
    for score in available_scores:
        for g1, g2 in group_pairs:
            df = compare_scores(adata, score, groupby, g1, g2, method)
            if not df.empty:
                df = df.copy()
                df["source_score"] = score
                results.append(df)

    if results:
        all_results = pd.concat(results, ignore_index=True)
        ns["batch_compare_results"] = sanitize_for_hdf5(all_results)
        return all_results
    else:
        log.info("No valid comparisons were produced.")
        return pd.DataFrame()


def batch_plot_score_comparison(
    adata: AnnData,
    score_keys: List[str],
    groupby: str,
    groups_to_compare: Optional[List[str]] = None,
    plot_type: Literal["violin", "boxplot"] = "violin",
    ncols: int = 2,
    figsize_per_ax: Tuple[float, float] = (6, 4),
    palette: Optional[str] = None,
    show: bool = True,
    save_path: Optional[str] = None,
    **kwargs,
):
    """
    Batch plot multiple score comparisons.

    Parameters
    - score_keys: list of obs score columns
    - groupby: grouping column in obs
    - groups_to_compare: optional subset of groups to show
    - ncols: number of columns in subplot grid
    - figsize_per_ax: size of each subplot
    - palette: seaborn palette name
    - show/save_path: display or save the figure

    Returns
    - List of matplotlib Axes for the plotted scores (length == n_plotted).
    """
    available_scores = [s for s in score_keys if s in adata.obs.columns]
    missing_scores = [s for s in score_keys if s not in adata.obs.columns]
    if missing_scores:
        log.warning(f"Missing score columns skipped: {missing_scores}")
    if len(available_scores) == 0:
        log.warning("No valid score columns to plot.")
        return []

    n = len(available_scores)
    ncols = max(1, ncols)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(figsize_per_ax[0] * ncols, figsize_per_ax[1] * nrows)
    )
    # Normalize axes array
    if isinstance(axes, np.ndarray):
        axes = axes.flatten()
    else:
        axes = np.array([axes])

    plotted_axes = []
    for i, score in enumerate(available_scores):
        ax = axes[i]
        plot_df = adata.obs[[groupby, score]].copy().dropna()
        if groups_to_compare:
            groups_set = set(map(str, groups_to_compare))
            plot_df = plot_df[plot_df[groupby].astype(str).isin(groups_set)]
            order = list(map(str, groups_to_compare))
        else:
            order = (
                list(adata.obs[groupby].cat.categories)
                if pd.api.types.is_categorical_dtype(adata.obs[groupby])
                else None
            )

        common_kwargs = dict(x=groupby, y=score, data=plot_df, order=order, palette=palette)
        if plot_type == "violin":
            sns.violinplot(**common_kwargs, cut=0, inner="quartile", ax=ax, **kwargs)
        elif plot_type == "boxplot":
            sns.boxplot(**common_kwargs, ax=ax, **kwargs)
            try:
                sns.swarmplot(x=groupby, y=score, data=plot_df, color=".25", size=3, ax=ax)
            except Exception:
                pass
        else:
            raise ValueError("plot_type must be 'violin' or 'boxplot'.")

        ax.set_title(score)
        ax.set_ylabel(score)
        ax.set_xlabel(groupby)
        ax.tick_params(axis="x", rotation=45)
        plotted_axes.append(ax)

    # Hide unused axes if any
    for j in range(len(plotted_axes), len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved batch score comparison plots to {save_path}")
    if show:
        plt.show()
    return plotted_axes
