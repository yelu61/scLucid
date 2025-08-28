import logging
from typing import Dict, List, Optional, Literal, Union
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from anndata import AnnData
from scipy.stats import ttest_ind, mannwhitneyu

log = logging.getLogger(__name__)

__all__ = [
    "score_by_gene_sets",
    "compare_scores",
    "plot_score_comparison",
    "batch_compare_scores",
    "batch_plot_score_comparison",
]

def score_by_gene_sets(
    adata: AnnData,
    gene_sets: Dict[str, List[str]],
    layer: Optional[str] = "log1p_norm",
    use_raw: bool = False,
    ctrl_size: int = 50,
    score_name_suffix: str = "_score",
    **kwargs,
) -> AnnData:
    """
    Scores cells for multiple gene sets (e.g., pathways, cell types).
    Adds columns to adata.obs.
    Stores scoring info in adata.uns['scrnatk']['analysis']['scoring'].
    """
    log.info(f"Scoring cells for {len(gene_sets)} gene sets...")
    source_adata = adata.raw if use_raw else adata
    if use_raw and adata.raw is None:
        raise ValueError("adata.raw is not set, but use_raw=True.")
    scored_count = 0
    for set_name, genes in gene_sets.items():
        genes_found = [g for g in genes if g in source_adata.var_names]
        if not genes_found:
            log.warning(f"No genes from '{set_name}' found in data. Skipping.")
            continue
        score_name = f"{set_name}{score_name_suffix}"
        sc.tl.score_genes(
            adata, genes_found, score_name=score_name,
            use_raw=use_raw, layer=layer if not use_raw else None,
            ctrl_size=ctrl_size, **kwargs
        )
        scored_count += 1
    # Store parameters
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('scoring', {})
    adata.uns['scrnatk']['analysis']['scoring']['gene_set_scoring'] = {
        'n_sets_scored': scored_count,
        'params': {'use_raw': use_raw, 'layer': layer, 'ctrl_size': ctrl_size}
    }
    log.info(f"Completed scoring for {scored_count} gene sets.")
    return adata

def compare_scores(
    adata: AnnData,
    score_key: str,
    groupby: str,
    group1: str,
    group2: Optional[str] = "rest",
    method: Literal["ttest", "wilcoxon"] = "wilcoxon"
) -> pd.DataFrame:
    """
    Performs differential analysis on a continuous score in adata.obs.
    Returns a DataFrame with results; stores in adata.uns['scrnatk']['analysis']['scoring'].
    """
    log.info(f"Comparing score '{score_key}' for '{group1}' vs '{group2}' in '{groupby}'")
    if score_key not in adata.obs:
        raise ValueError(f"Score key '{score_key}' not found in adata.obs.")
    scores1 = adata.obs.loc[adata.obs[groupby] == group1, score_key].dropna()
    scores2 = (
        adata.obs.loc[adata.obs[groupby] != group1, score_key].dropna()
        if group2 == "rest" else
        adata.obs.loc[adata.obs[groupby] == group2, score_key].dropna()
    )
    n1, n2 = len(scores1), len(scores2)
    if n1 < 3 or n2 < 3:
        log.warning("Not enough data points in one or both groups for statistical test.")
        return pd.DataFrame()
    if method == "wilcoxon":
        stat, pval = mannwhitneyu(scores1, scores2, alternative='two-sided')
    elif method == "ttest":
        stat, pval = ttest_ind(scores1, scores2, equal_var=False)
    else:
        raise ValueError("Method must be 'ttest' or 'wilcoxon'")
    mean1, mean2 = scores1.mean(), scores2.mean()
    effect_size = mean1 - mean2
    results = pd.DataFrame({
        "score": [score_key],
        "group1": [group1],
        "group2": [group2],
        "method": [method],
        "statistic": [stat],
        "pvalue": [pval],
        "mean_group1": [mean1],
        "mean_group2": [mean2],
        "effect_size (mean_diff)": [effect_size],
        "n_cells_group1": [n1],
        "n_cells_group2": [n2]
    })
    # Store result
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('scoring', {})
    adata.uns['scrnatk']['analysis']['scoring'][f"compare_{score_key}_{group1}_vs_{group2}"] = results
    return results

def plot_score_comparison(
    adata: AnnData,
    score_key: str,
    groupby: str,
    groups_to_compare: Optional[List[str]] = None,
    plot_type: Literal["violin", "boxplot"] = "violin",
    **kwargs
) -> plt.Axes:
    """
    Creates a plot comparing a score across different groups.
    Returns matplotlib Axes.
    """
    plot_df = adata.obs[[groupby, score_key]].copy()
    if groups_to_compare:
        plot_df = plot_df[plot_df[groupby].isin(groups_to_compare)]
    plt.figure(figsize=kwargs.pop("figsize", (8, 6)))
    if plot_type == "violin":
        ax = sns.violinplot(data=plot_df, x=groupby, y=score_key, **kwargs)
    elif plot_type == "boxplot":
        ax = sns.boxplot(data=plot_df, x=groupby, y=score_key, **kwargs)
        sns.swarmplot(data=plot_df, x=groupby, y=score_key, color=".25", size=3, ax=ax)
    else:
        raise ValueError("plot_type must be 'violin' or 'boxplot'")
    ax.set_title(f"Comparison of '{score_key}' across '{groupby}'")
    ax.set_ylabel(score_key)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    # Store plotting info
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {}).setdefault('scoring', {})
    adata.uns['scrnatk']['analysis']['scoring'][f"{score_key}_{groupby}_{plot_type}_plot"] = {
        "groups": groups_to_compare, "plot_type": plot_type
    }
    plt.show()
    return ax

def batch_compare_scores(
    adata: AnnData,
    score_keys: List[str],
    groupby: str,
    group_pairs: Optional[List[tuple]] = None,
    method: Literal["ttest", "wilcoxon"] = "wilcoxon"
) -> pd.DataFrame:
    """
    Batch differential comparison for multiple scores and group pairs.
    Returns a concatenated DataFrame.
    """
    results = []
    if group_pairs is None:
        groups = adata.obs[groupby].unique().tolist()
        group_pairs = [(g, "rest") for g in groups]
    for score in score_keys:
        for g1, g2 in group_pairs:
            df = compare_scores(adata, score, groupby, g1, g2, method)
            if not df.empty:
                results.append(df)
    if results:
        all_results = pd.concat(results, ignore_index=True)
        adata.uns['scrnatk']['analysis']['scoring']['batch_compare_results'] = all_results
        return all_results
    else:
        return pd.DataFrame()

def batch_plot_score_comparison(
    adata: AnnData,
    score_keys: List[str],
    groupby: str,
    groups_to_compare: Optional[List[str]] = None,
    plot_type: Literal["violin", "boxplot"] = "violin",
    ncols: int = 2,
    **kwargs
):
    """
    Batch plot for multiple scores. Returns list of Axes.
    """
    n = len(score_keys)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(8 * ncols, 5 * nrows))
    axes = axes.flatten()
    for i, score in enumerate(score_keys):
        ax = axes[i]
        plot_df = adata.obs[[groupby, score]].copy()
        if groups_to_compare:
            plot_df = plot_df[plot_df[groupby].isin(groups_to_compare)]
        if plot_type == "violin":
            sns.violinplot(data=plot_df, x=groupby, y=score, ax=ax, **kwargs)
        elif plot_type == "boxplot":
            sns.boxplot(data=plot_df, x=groupby, y=score, ax=ax, **kwargs)
            sns.swarmplot(data=plot_df, x=groupby, y=score, color=".25", size=3, ax=ax)
        ax.set_title(score)
        ax.set_ylabel(score)
        ax.set_xlabel(groupby)
        ax.tick_params(axis='x', rotation=45)
    plt.tight_layout()
    plt.show()
    return axes[:n]