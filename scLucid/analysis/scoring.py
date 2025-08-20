
import logging
from typing import Dict, List, Optional, Literal
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
from anndata import AnnData
from scipy.stats import ttest_ind, mannwhitneyu
# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = [
    "score_by_gene_sets",
    "compare_scores",
    "plot_score_comparison",
]

# --- Main Functions ---
def score_by_gene_sets(
    adata: AnnData,
    gene_sets: Dict[str, List[str]], # MODIFIED: Takes a dictionary directly
    layer: Optional[str] = "log1p_norm",
    use_raw: bool = False,
    ctrl_size: int = 50,
    score_name_suffix: str = "_score",
    **kwargs,
) -> AnnData:
    """
    Scores cells for multiple gene sets (e.g., pathways, cell types).

    This function calculates enrichment scores for each gene set by comparing the
    expression of genes in the set to randomly selected control genes.

    Args:
        adata: AnnData object.
        gene_sets: Dictionary where keys are set names (e.g., 'HALLMARK_APOPTOSIS')
                   and values are lists of gene symbols.
        layer: Layer to use for scoring.
        use_raw: If True, use `adata.raw` for scoring.
        ctrl_size: Number of control genes.
        score_name_suffix: Suffix to add to the gene set name for the .obs column.
        **kwargs: Additional arguments passed to sc.tl.score_genes.
    """
    log.info(f"Scoring cells for {len(gene_sets)} gene sets...")
    
    source_adata = adata.raw if use_raw else adata
    if source_adata is None:
        raise ValueError("adata.raw is not set, but use_raw=True.")

    scored_count = 0
    for set_name, genes in gene_sets.items():
        # Filter for genes present in the data
        genes_found = [g for g in genes if g in source_adata.var_names]
        if not genes_found:
            log.warning(f"No genes from '{set_name}' found in data. Skipping.")
            continue
            
        score_name = f"{set_name}{score_name_suffix}"
        sc.tl.score_genes(
            adata,
            gene_list=genes_found,
            score_name=score_name,
            use_raw=use_raw,
            layer=layer if not use_raw else None,
            ctrl_size=ctrl_size,
            **kwargs
        )
        scored_count += 1

    log.info(f"Completed scoring for {scored_count} gene sets.")
    # Store parameters in a unified namespace
    adata.uns.setdefault('scrnatk', {}).setdefault('analysis', {})['gene_set_scoring'] = {
        'n_sets_scored': scored_count,
        'params': {'use_raw': use_raw, 'layer': layer, 'ctrl_size': ctrl_size}
    }
    
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

    Args:
        adata: AnnData object.
        score_key: The column in adata.obs containing the score to test.
        groupby: The key in adata.obs for group assignments.
        group1: The first group for comparison.
        group2: The second group for comparison, or 'rest' for one-vs-rest.
        method: The statistical test to use ('ttest' or 'wilcoxon').

    Returns:
        A DataFrame with the results of the comparison.
    """
    log.info(f"Comparing score '{score_key}' for '{group1}' vs '{group2}' in '{groupby}'")
    if score_key not in adata.obs:
        raise ValueError(f"Score key '{score_key}' not found in adata.obs.")

    # Get data for group 1
    scores1 = adata.obs.loc[adata.obs[groupby] == group1, score_key].dropna()

    # Get data for group 2
    if group2 == "rest":
        scores2 = adata.obs.loc[adata.obs[groupby] != group1, score_key].dropna()
    else:
        scores2 = adata.obs.loc[adata.obs[groupby] == group2, score_key].dropna()

    if len(scores1) < 3 or len(scores2) < 3:
        log.warning("Not enough data points in one or both groups for statistical test.")
        return pd.DataFrame()

    # Perform statistical test
    if method == "wilcoxon":
        stat, pval = mannwhitneyu(scores1, scores2, alternative='two-sided')
    elif method == "ttest":
        stat, pval = ttest_ind(scores1, scores2, equal_var=False) # Welch's t-test
    else:
        raise ValueError("Method must be 'ttest' or 'wilcoxon'")

    # Calculate effect size (mean difference)
    mean1 = scores1.mean()
    mean2 = scores2.mean()
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
        "effect_size (mean_diff)": [effect_size]
    })
    
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
    Creates a publication-quality plot comparing a score across different groups.

    Args:
        adata: AnnData object.
        score_key: The column in adata.obs containing the score to plot.
        groupby: The key in adata.obs for group assignments.
        groups_to_compare: Optional list of specific groups to include in the plot.
        plot_type: The type of plot to generate.
        **kwargs: Additional arguments passed to the seaborn plotting function.
    """
    plot_df = adata.obs[[groupby, score_key]].copy()
    if groups_to_compare:
        plot_df = plot_df[plot_df[groupby].isin(groups_to_compare)]
    
    plt.figure(figsize=kwargs.pop("figsize", (8, 6)))
    
    if plot_type == "violin":
        ax = sns.violinplot(data=plot_df, x=groupby, y=score_key, **kwargs)
    elif plot_type == "boxplot":
        ax = sns.boxplot(data=plot_df, x=groupby, y=score_key, **kwargs)
        # Overlay swarmplot for better visualization of data points
        sns.swarmplot(data=plot_df, x=groupby, y=score_key, color=".25", size=3, ax=ax)
    else:
        raise ValueError("plot_type must be 'violin' or 'boxplot'")
        
    ax.set_title(f"Comparison of '{score_key}' across '{groupby}'")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
    return ax