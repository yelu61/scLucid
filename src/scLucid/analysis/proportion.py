"""
Cell type proportion statistics and visualization.

This module provides comprehensive tools for analyzing and visualizing cell type
proportions in single-cell RNA-seq data, including:
- Proportion & count computation
- Multiple statistical tests (DESeq2, t-test, Wilcoxon, ANOVA, paired tests)
- Effect size calculation (Cohen's d, Cliff's Delta)
- Rich visualization suite (15+ plot types)
- Batch effect detection
- Time-series analysis support

Quick Start
-----------
>>> from scrnaseq_pipeline.proportion import celltype_proportion_analysis
>>> from scrnaseq_pipeline.config import ProportionConfig
>>>
>>> # Basic configuration
>>> config = ProportionConfig(
...     celltype_col='cell_type',
...     sample_col='sample_id',
...     condition_col='condition',
...     test_method='wilcoxon',
...     plot_types=['bar', 'box', 'diff', 'heatmap', 'volcano'],
...     out_dir='./proportion_results'
... )
>>>
>>> # Run analysis
>>> prop_df, stat_df = celltype_proportion_analysis(adata, config)
>>>
>>> # Access results
>>> print(stat_df[stat_df['padj'] < 0.05])  # Significant cell types

Advanced Usage
--------------
**Multi-group ANOVA:**
>>> config = ProportionConfig(
...     test_method='anova',
...     condition_col='treatment',  # 3+ groups
...     plot_types=['box', 'heatmap']
... )

**Paired analysis (e.g., pre/post treatment):**
>>> config = ProportionConfig(
...     test_method='paired-wilcoxon',
...     pairing_col='patient_id',  # Subject identifier
...     condition_col='timepoint'
... )

**Time series tracking:**
>>> from scrnaseq_pipeline.proportion import plot_proportion_timeseries
>>> plot_proportion_timeseries(
...     prop_df,
...     timepoint_col=adata.obs['timepoint'],
...     celltype='CD8_T',
...     group_col=adata.obs['treatment']
... )

Plot Types Guide
----------------
Core plots:
  - 'counts': Total cell counts per sample/group
  - 'bar': Stacked proportion bars
  - 'bar_composition': Condition contribution per cell type
  - 'diff': Barplot with significance brackets
  - 'box': Boxplot with strip overlay
  - 'alluvial': Sankey-style flow diagram
  - 'heatmap': Proportion heatmap across samples
  - 'ci': Confidence interval plot
  - 'correlation': Cell type correlation matrix

Advanced plots:
  - 'volcano': Effect size volcano plot
  - 'variability': Coefficient of variation (CV) plot
  - 'batch_pca': Batch effect visualization
  - 'timeseries': Longitudinal tracking

Statistical Methods
-------------------
- 'deseq2': DESeq2 differential abundance (requires pydeseq2)
- 't-test': Independent samples t-test
- 'wilcoxon': Mann-Whitney U test (non-parametric)
- 'anova': One-way ANOVA (>2 groups)
- 'paired-t-test': Paired t-test
- 'paired-wilcoxon': Paired Wilcoxon signed-rank test
"""

import logging
from functools import wraps
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

# Try importing pydeseq2
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    HAS_DESEQ2 = True
except ImportError:
    HAS_DESEQ2 = False

from .config import ProportionConfig

log = logging.getLogger(__name__)

__all__ = [
    "celltype_proportion_analysis",
    "compute_celltype_proportion",
    "run_statistical_test",
    "export_analysis_data",
    "plot_cell_counts",
    "plot_proportion_bar",
    "plot_diff_stats",
    "plot_composition",
    "plot_box_summary",
    "plot_proportion_shifts",
    "plot_individual_boxplots",
    "plot_proportion_heatmap",
    "plot_proportion_with_ci",
    "plot_celltype_correlation",
    "plot_effect_size_volcano",
    "plot_celltype_variability",
    "plot_batch_effect",
    "plot_proportion_timeseries",
]

# ================= Publication-ready Style =================


def _set_publication_style():
    """Set publication-quality plotting style."""
    plt.style.use("seaborn-v0_8-paper")

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 11,
            "axes.labelsize": 12,
            "axes.titlesize": 14,
            "axes.titleweight": "bold",
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.linewidth": 1.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


# ================= Decorator for Plot Auto-Save =================


def save_and_close(plot_name: str):
    """
    Decorator to automatically save and close plots.

    Parameters
    ----------
    plot_name : str
        Base name for the saved file (without extension)

    Returns
    -------
    function
        Decorated plotting function
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, out_dir=None, **kwargs):
            fig = func(*args, **kwargs)

            if fig is not None:
                plt.tight_layout()

                if out_dir:
                    out_path = Path(out_dir) / f"{plot_name}.pdf"
                    plt.savefig(out_path, dpi=300, bbox_inches="tight")
                    log.debug(f"Saved plot to {out_path}")

                plt.close(fig)

            return fig

        return wrapper

    return decorator


# ================= Helper Functions =================


def _get_sig_stars(p_val: float) -> str:
    """
    Convert p-value to significance stars.

    Parameters
    ----------
    p_val : float
        P-value or adjusted p-value

    Returns
    -------
    str
        Significance annotation ('***', '**', '*', 'ns')
    """
    if pd.isna(p_val):
        return "ns"
    if p_val < 0.001:
        return "***"
    if p_val < 0.01:
        return "**"
    if p_val < 0.05:
        return "*"
    return "ns"


def _ensure_palette(
    palette: Optional[Dict], keys: pd.Index, default_cmap: str = "husl"
) -> Dict:
    """
    Ensure a color palette exists for the given keys.

    Parameters
    ----------
    palette : Optional[Dict]
        User-provided palette mapping
    keys : pd.Index
        Categories to assign colors to
    default_cmap : str
        Default colormap name

    Returns
    -------
    Dict
        Complete color palette
    """
    if palette is None:
        # Sort string keys for consistent color assignment
        sorted_keys = sorted(keys) if all(isinstance(k, str) for k in keys) else keys
        colors = sns.color_palette(default_cmap, len(sorted_keys)).as_hex()
        return dict(zip(sorted_keys, colors))
    return palette


def _calculate_bracket_height(
    ax: plt.Axes, y_data: np.ndarray, num_brackets: int = 1, base_gap: float = 0.03
) -> float:
    """
    Dynamically calculate the height for statistical annotation brackets.

    Parameters
    ----------
    ax : plt.Axes
        Matplotlib axes object
    y_data : np.ndarray
        Y-axis data values
    num_brackets : int
        Number of brackets to stack
    base_gap : float
        Base gap as fraction of y-range

    Returns
    -------
    float
        Calculated gap height
    """
    ylim = ax.get_ylim()
    y_range = ylim[1] - ylim[0]

    # Base gap (percentage of total height)
    gap = y_range * base_gap

    # If all data values are very small, use absolute gap
    data_max = np.max(y_data)
    if data_max < y_range * 0.1:
        gap = max(gap, data_max * 0.1)

    return gap


def _natural_sort_key(text):
    """
    Natural sorting key that handles numbers correctly.
    Example: ['Day1', 'Day2', 'Day10'] instead of ['Day1', 'Day10', 'Day2']
    """
    import re
    return [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', str(text))]


# ================= Core Logic: Computation =================


def compute_celltype_proportion(
    adata: AnnData,
    celltype_col: str,
    sample_col: str,
    condition_col: Optional[str] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[pd.Series]]:
    """
    Compute raw counts and proportions for cell types across samples.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell data object
    celltype_col : str
        Column name in .obs containing cell type annotations
    sample_col : str
        Column name in .obs containing sample identifiers
    condition_col : Optional[str]
        Column name in .obs containing experimental conditions

    Returns
    -------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : Optional[pd.Series]
        Mapping from sample to condition

    Raises
    ------
    KeyError
        If required columns are missing from adata.obs

    Examples
    --------
    >>> prop_df, count_df, mapping = compute_celltype_proportion(
    ...     adata, 'cell_type', 'sample_id', 'condition'
    ... )
    """
    if celltype_col not in adata.obs.columns or sample_col not in adata.obs.columns:
        raise KeyError(f"Missing columns: {celltype_col} or {sample_col}")

    # Compute count matrix
    count_df = pd.crosstab(adata.obs[sample_col], adata.obs[celltype_col])

    # Compute proportion matrix
    totals = count_df.sum(axis=1).replace(0, np.nan)
    prop_df = count_df.div(totals, axis=0).fillna(0.0)

    sample_to_cond = None
    if condition_col and condition_col in adata.obs.columns:
        # Create sample-to-condition mapping
        sample_to_cond = (
            adata.obs[[sample_col, condition_col]]
            .drop_duplicates(subset=[sample_col])
            .set_index(sample_col)[condition_col]
        )

        # Align indices
        common = prop_df.index.intersection(sample_to_cond.index)
        prop_df = prop_df.loc[common]
        count_df = count_df.loc[common]
        sample_to_cond = sample_to_cond.loc[common]

    log.info(
        f"Computed proportions for {prop_df.shape[1]} cell types across {prop_df.shape[0]} samples"
    )

    return prop_df, count_df, sample_to_cond


# ================= Statistical Testing =================


def _run_deseq2(
    count_df: pd.DataFrame, sample_to_cond: pd.Series, design_factor: str = "condition"
) -> pd.DataFrame:
    """
    Run DESeq2 differential abundance analysis.

    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    design_factor : str
        Column name for the design factor

    Returns
    -------
    pd.DataFrame
        Results with columns: celltype, log2FoldChange, pvalue, padj, comparison
    """
    if not HAS_DESEQ2:
        log.warning("pydeseq2 not installed. Skipping DESeq2 analysis.")
        return pd.DataFrame()

    # Filter low-count cell types (at least 2 samples with >0 cells)
    keep = (count_df > 0).sum(axis=0) >= 2
    count_df_filtered = count_df.loc[:, keep]

    if count_df_filtered.shape[1] < 1:
        log.warning("No cell types passed the filtering criteria for DESeq2.")
        return pd.DataFrame()

    # Prepare metadata
    metadata = sample_to_cond.loc[count_df_filtered.index].to_frame(name=design_factor)
    conditions = metadata[design_factor].unique()

    if len(conditions) != 2:
        log.warning(
            f"DESeq2 requires exactly 2 conditions, found {len(conditions)}. Skipping."
        )
        return pd.DataFrame()

    try:
        # Initialize DESeq2 dataset
        dds = DeseqDataSet(
            counts=count_df_filtered,
            metadata=metadata,
            design_factors=design_factor,
            quiet=True,
        )

        # Calculate size factors based on total cell counts
        total_cells = count_df_filtered.sum(axis=1)
        size_factors = total_cells / np.exp(np.mean(np.log(total_cells + 1)))
        dds.obsm["size_factors"] = size_factors.values

        # Run DESeq2
        dds.deseq2()

        # Perform statistical testing (Condition2 vs Condition1)
        res = DeseqStats(
            dds, contrast=(design_factor, conditions[1], conditions[0]), quiet=True
        )
        res.summary()

        res_df = res.results_df.reset_index().rename(columns={"index": "celltype"})
        res_df["comparison"] = f"{conditions[1]}_vs_{conditions[0]}"

        log.info(f"DESeq2 completed: {len(res_df)} cell types analyzed.")
        return res_df

    except Exception as e:
        log.error(f"DESeq2 analysis failed: {str(e)}", exc_info=True)
        return pd.DataFrame()


def _run_ttest(count_df: pd.DataFrame, sample_to_cond: pd.Series) -> pd.DataFrame:
    """
    Run independent t-test for each cell type.

    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping

    Returns
    -------
    pd.DataFrame
        Results with columns: celltype, statistic, pvalue, padj, mean_group1, mean_group2, log2fc
    """
    # Convert counts to proportions
    prop_df = count_df.div(count_df.sum(axis=1), axis=0)

    conditions = sample_to_cond.unique()
    if len(conditions) != 2:
        log.warning(f"t-test requires 2 conditions, found {len(conditions)}.")
        return pd.DataFrame()

    results = []
    for ct in count_df.columns:
        groups = prop_df.groupby(sample_to_cond)[ct].apply(list)

        if len(groups) == 2:
            g1, g2 = groups.iloc[0], groups.iloc[1]

            # Skip if either group has no variance
            if len(g1) < 2 or len(g2) < 2:
                continue

            try:
                stat, pval = stats.ttest_ind(g1, g2)

                results.append(
                    {
                        "celltype": ct,
                        "statistic": stat,
                        "pvalue": pval,
                        "mean_group1": np.mean(g1),
                        "mean_group2": np.mean(g2),
                        "log2fc": np.log2(
                            (np.mean(g2) + 1e-10) / (np.mean(g1) + 1e-10)
                        ),
                    }
                )
            except Exception as e:
                log.warning(f"t-test failed for {ct}: {e}")
                continue

    if not results:
        return pd.DataFrame()

    res_df = pd.DataFrame(results)
    res_df["padj"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]
    res_df["comparison"] = f"{conditions[1]}_vs_{conditions[0]}"

    log.info(f"t-test completed: {len(res_df)} cell types analyzed.")
    return res_df


def _run_wilcoxon(count_df: pd.DataFrame, sample_to_cond: pd.Series) -> pd.DataFrame:
    """
    Run Mann-Whitney U test (non-parametric alternative to t-test).

    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping

    Returns
    -------
    pd.DataFrame
        Results with columns: celltype, statistic, pvalue, padj, median_group1, median_group2
    """
    # Convert counts to proportions
    prop_df = count_df.div(count_df.sum(axis=1), axis=0)

    conditions = sample_to_cond.unique()
    if len(conditions) != 2:
        log.warning(f"Wilcoxon test requires 2 conditions, found {len(conditions)}.")
        return pd.DataFrame()

    results = []
    for ct in count_df.columns:
        groups = prop_df.groupby(sample_to_cond)[ct].apply(list)

        if len(groups) == 2:
            g1, g2 = groups.iloc[0], groups.iloc[1]

            if len(g1) < 2 or len(g2) < 2:
                continue

            try:
                stat, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")

                results.append(
                    {
                        "celltype": ct,
                        "statistic": stat,
                        "pvalue": pval,
                        "median_group1": np.median(g1),
                        "median_group2": np.median(g2),
                        "log2fc": np.log2(
                            (np.median(g2) + 1e-10) / (np.median(g1) + 1e-10)
                        ),
                    }
                )
            except Exception as e:
                log.warning(f"Wilcoxon test failed for {ct}: {e}")
                continue

    if not results:
        return pd.DataFrame()

    res_df = pd.DataFrame(results)
    res_df["padj"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]
    res_df["comparison"] = f"{conditions[1]}_vs_{conditions[0]}"

    log.info(f"Wilcoxon test completed: {len(res_df)} cell types analyzed.")
    return res_df


def _run_anova(count_df: pd.DataFrame, sample_to_cond: pd.Series) -> pd.DataFrame:
    """
    Run one-way ANOVA for multiple group comparison.

    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping

    Returns
    -------
    pd.DataFrame
        Results with columns: celltype, statistic, pvalue, padj
    """
    from scipy.stats import f_oneway

    # Convert counts to proportions
    prop_df = count_df.div(count_df.sum(axis=1), axis=0)

    results = []
    for ct in count_df.columns:
        groups = [g.values for _, g in prop_df.groupby(sample_to_cond)[ct]]

        if len(groups) < 2:
            continue

        # Filter out groups with insufficient data
        groups = [g for g in groups if len(g) >= 2]

        if len(groups) < 2:
            continue

        try:
            stat, pval = f_oneway(*groups)
            results.append({"celltype": ct, "statistic": stat, "pvalue": pval})
        except Exception as e:
            log.warning(f"ANOVA failed for {ct}: {e}")
            continue

    if not results:
        return pd.DataFrame()

    res_df = pd.DataFrame(results)
    res_df["padj"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]

    log.info(f"ANOVA completed: {len(res_df)} cell types analyzed.")
    return res_df


def _run_contingency_test(
    count_df: pd.DataFrame, 
    sample_to_cond: pd.Series, 
    method: str = "chi-square"
) -> pd.DataFrame:
    """
    Run Chi-square or Fisher's exact test for N=1 vs N=1 comparison.
    
    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples x cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    method : str
        'chi-square' or 'fisher'
        
    Returns
    -------
    pd.DataFrame
        Results with p-values and odds ratios (for Fisher)
    """
    from scipy.stats import chi2_contingency, fisher_exact
    
    conditions = sample_to_cond.unique()
    if len(conditions) != 2:
        log.warning(f"{method} requires exactly 2 conditions, found {len(conditions)}.")
        return pd.DataFrame()
        
    cond1, cond2 = conditions[0], conditions[1]
    
    samples_per_cond = count_df.groupby(sample_to_cond).size()

    if (samples_per_cond > 1).any():
        log.warning(
            f"Chi-square test detected multiple samples per condition. "
            f"This test is designed for N=1 comparisons. Consider using t-test/Wilcoxon instead."
        )
    
    # Aggregate counts by condition (Since N=1, this just grabs the single row)
    # If N>1 by accident, this sums them up, treating them as one pool
    group_counts = count_df.groupby(sample_to_cond).sum()
    
    total_cond1 = group_counts.loc[cond1].sum()
    total_cond2 = group_counts.loc[cond2].sum()
    
    results = []
    
    for ct in count_df.columns:
        # Construct 2x2 Contingency Table
        #           Cond1    Cond2
        # Target      a        b
        # Others      c        d
        
        a = group_counts.loc[cond1, ct]
        b = group_counts.loc[cond2, ct]
        c = total_cond1 - a
        d = total_cond2 - b
        
        table = [[a, b], [c, d]]
        
        try:
            if method == "fisher":
                oddsr, pval = fisher_exact(table, alternative='two-sided')
                stat = oddsr # Odds ratio as statistic
            else: # chi-square
                # correction=False is often used for large N in scRNA-seq
                stat, pval, _, _ = chi2_contingency(table, correction=True)
            
            # Calculate Log2 Fold Change of Proportions
            prop1 = a / total_cond1
            prop2 = b / total_cond2
            log2fc = np.log2((prop2 + 1e-6) / (prop1 + 1e-6))
            
            results.append({
                "celltype": ct,
                "statistic": stat,
                "pvalue": pval,
                "log2fc": log2fc,
                "mean_group1": prop1,
                "mean_group2": prop2
            })
            
        except Exception as e:
            log.warning(f"{method} failed for {ct}: {e}")
            continue
            
    if not results:
        return pd.DataFrame()
        
    res_df = pd.DataFrame(results)
    # Still perform FDR correction, though interpretations should be cautious
    res_df["padj"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]
    res_df["comparison"] = f"{cond2}_vs_{cond1}_(cell-level)"
    
    log.info(f"{method} test completed: {len(res_df)} cell types analyzed.")
    return res_df


def _run_paired_test(
    count_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    pairing_col: pd.Series,
    method: str = "wilcoxon",
) -> pd.DataFrame:
    """
    Run paired statistical tests (🆕 NEW).

    For paired experimental designs (e.g., same patient pre/post treatment).

    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    pairing_col : pd.Series
        Mapping from sample to subject/patient ID
    method : str
        'wilcoxon' (paired Wilcoxon signed-rank) or 't-test' (paired t-test)

    Returns
    -------
    pd.DataFrame
        Results with columns: celltype, statistic, pvalue, padj, median_diff, n_pairs

    Examples
    --------
    >>> pairing = pd.Series(
    ...     ['Patient1', 'Patient1', 'Patient2', 'Patient2'],
    ...     index=['Sample1', 'Sample2', 'Sample3', 'Sample4']
    ... )
    >>> stat_df = _run_paired_test(count_df, sample_to_cond, pairing, method='wilcoxon')
    """
    common_samples = count_df.index.intersection(sample_to_cond.index).intersection(pairing_col.index)
    if len(common_samples) < 3:
        log.warning(f"Insufficient aligned samples (n={len(common_samples)}) for paired test.")
        return pd.DataFrame()
    
    count_df = count_df.loc[common_samples]
    sample_to_cond = sample_to_cond.loc[common_samples]
    pairing_col = pairing_col.loc[common_samples]
    
    # Convert counts to proportions
    prop_df = count_df.div(count_df.sum(axis=1), axis=0)

    # Align pairing information
    pairing_col = pairing_col.loc[prop_df.index]

    # Identify valid pairs
    paired_data = []
    conditions = sample_to_cond.unique()

    if len(conditions) != 2:
        log.warning(
            f"Paired test requires exactly 2 conditions, found {len(conditions)}."
        )
        return pd.DataFrame()

    cond1, cond2 = conditions[0], conditions[1]

    for subject in pairing_col.unique():
        subject_samples = pairing_col[pairing_col == subject].index
        subject_conds = sample_to_cond.loc[subject_samples]

        # Must have exactly 2 conditions per subject
        if len(subject_conds.unique()) == 2:
            sample1 = subject_samples[subject_conds == cond1]
            sample2 = subject_samples[subject_conds == cond2]

            if len(sample1) == 1 and len(sample2) == 1:
                paired_data.append(
                    {"subject": subject, "sample1": sample1[0], "sample2": sample2[0]}
                )

    if len(paired_data) < 3:
        log.warning(
            f"Insufficient paired samples (n={len(paired_data)}). Need at least 3 pairs."
        )
        return pd.DataFrame()

    log.info(f"Found {len(paired_data)} valid paired samples")

    results = []
    for ct in count_df.columns:
        before_vals, after_vals = [], []

        for pair in paired_data:
            before_vals.append(prop_df.loc[pair["sample1"], ct])
            after_vals.append(prop_df.loc[pair["sample2"], ct])

        if method == "wilcoxon":
            try:
                stat, pval = stats.wilcoxon(before_vals, after_vals)
            except Exception as e:
                log.warning(f"Paired Wilcoxon failed for {ct}: {e}")
                continue

        elif method == "t-test":
            try:
                stat, pval = stats.ttest_rel(before_vals, after_vals)
            except Exception as e:
                log.warning(f"Paired t-test failed for {ct}: {e}")
                continue
        else:
            raise ValueError(f"Unsupported paired test: {method}")

        results.append(
            {
                "celltype": ct,
                "statistic": stat,
                "pvalue": pval,
                "median_diff": np.median(np.array(after_vals) - np.array(before_vals)),
                "mean_diff": np.mean(np.array(after_vals) - np.array(before_vals)),
                "n_pairs": len(before_vals),
            }
        )

    if not results:
        return pd.DataFrame()

    res_df = pd.DataFrame(results)
    res_df["padj"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]
    res_df["comparison"] = f"{cond2}_vs_{cond1}_paired"

    log.info(f"Paired {method} completed: {len(res_df)} cell types analyzed.")
    return res_df


def run_statistical_test(
    count_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    method: str = "deseq2",
    design_factor: str = "condition",
    pairing_col: Optional[pd.Series] = None,
    correction_scope: str = "per_test",
) -> pd.DataFrame:
    """
    Unified interface for running statistical tests.

    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    method : str
        Statistical test method:
        - 'deseq2': DESeq2 differential abundance
        - 't-test': Independent samples t-test
        - 'wilcoxon': Mann-Whitney U test
        - 'anova': One-way ANOVA (>2 groups)
        - 'paired-t-test': Paired t-test (🆕 NEW)
        - 'paired-wilcoxon': Paired Wilcoxon signed-rank (🆕 NEW)
    design_factor : str
        Column name for the design factor (for DESeq2)
    pairing_col : Optional[pd.Series]
        Subject/patient IDs for paired tests (🆕 NEW)
    correction_scope : str
        'per_test' (default): Correct p-values within each test
        'global': Correct across all tests (more conservative) (🆕 NEW)

    Returns
    -------
    pd.DataFrame
        Statistical test results

    Raises
    ------
    ValueError
        If an unsupported method is specified

    Examples
    --------
    >>> # Independent test
    >>> stat_df = run_statistical_test(count_df, sample_to_cond, method='wilcoxon')
    >>>
    >>> # Paired test
    >>> stat_df = run_statistical_test(
    ...     count_df, sample_to_cond,
    ...     method='paired-wilcoxon',
    ...     pairing_col=patient_ids
    ... )
    """
    method = method.lower()

    log.info(f"Running {method} statistical test...")

    # Run appropriate test
    if method == "deseq2":
        res_df = _run_deseq2(count_df, sample_to_cond, design_factor)

    elif method == "t-test":
        res_df = _run_ttest(count_df, sample_to_cond)

    elif method == "wilcoxon":
        res_df = _run_wilcoxon(count_df, sample_to_cond)

    elif method == "anova":
        res_df = _run_anova(count_df, sample_to_cond)
    
    elif method in ["chi-square", "fisher"]:
        res_df = _run_contingency_test(count_df, sample_to_cond, method)

    elif method in ["paired-t-test", "paired-wilcoxon"]:
        if pairing_col is None:
            raise ValueError(f"{method} requires pairing_col parameter")

        paired_method = "t-test" if "t-test" in method else "wilcoxon"
        res_df = _run_paired_test(count_df, sample_to_cond, pairing_col, paired_method)

    else:
        raise ValueError(
            f"Unsupported method: {method}. Choose from 'deseq2', 't-test', "
            f"'wilcoxon', 'anova', 'chi-square', 'fisher', 'paired-t-test', 'paired-wilcoxon'."
        )

    # Apply global correction if requested
    if correction_scope == "global" and not res_df.empty:
        res_df["padj"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]
        log.info("Applied global multiple testing correction")

    return res_df


# ================= Effect Size Calculation =================


def _calculate_effect_size(
    group1: np.ndarray, group2: np.ndarray, method: str = "cohen_d"
) -> float:
    """
    Calculate effect size between two groups.

    Parameters
    ----------
    group1 : np.ndarray
        First group values
    group2 : np.ndarray
        Second group values
    method : str
        Effect size metric ('cohen_d' or 'cliff_delta')

    Returns
    -------
    float
        Calculated effect size

    Notes
    -----
    - Cohen's d: Standardized mean difference
      - Small: 0.2, Medium: 0.5, Large: 0.8
    - Cliff's Delta: Non-parametric effect size
      - Negligible: < 0.147, Small: 0.147-0.33, Medium: 0.33-0.474, Large: > 0.474
    """
    
    if len(group1) < 2 or len(group2) < 2:
        return np.nan
    
    if method == "cohen_d":
        # Pooled standard deviation
        std1 = np.std(group1, ddof=1)
        std2 = np.std(group2, ddof=1)
        
        pooled_std = np.sqrt((std1 ** 2 + std2 ** 2) / 2)
        
        if pooled_std < 1e-10:  # ✅ Better than == 0
            return 0.0
        
        return (np.mean(group2) - np.mean(group1)) / pooled_std

    elif method == "cliff_delta":
        # Cliff's Delta (non-parametric effect size)
        n1, n2 = len(group1), len(group2)
        if n1 == 0 or n2 == 0:
            return 0.0

        # Count dominance: how many times group2 > group1
        dominance = sum(1 for x in group2 for y in group1 if x > y)
        return (dominance - (n1 * n2 / 2)) / (n1 * n2)

    else:
        raise ValueError(
            f"Unsupported method: {method}. Choose 'cohen_d' or 'cliff_delta'."
        )


def _add_effect_sizes(
    stat_df: pd.DataFrame, prop_df: pd.DataFrame, sample_to_cond: pd.Series
) -> pd.DataFrame:
    """
    Add effect size columns to statistical test results.

    Parameters
    ----------
    stat_df : pd.DataFrame
        Statistical test results
    prop_df : pd.DataFrame
        Proportion matrix
    sample_to_cond : pd.Series
        Sample to condition mapping

    Returns
    -------
    pd.DataFrame
        Results with added 'cohen_d' and 'cliff_delta' columns
    """
    if stat_df.empty:
        return stat_df

    conditions = sample_to_cond.unique()
    if len(conditions) != 2:
        log.warning("Effect size calculation requires exactly 2 conditions.")
        return stat_df

    cohen_d_list = []
    cliff_delta_list = []

    for ct in stat_df["celltype"]:
        groups = prop_df.groupby(sample_to_cond)[ct].apply(lambda x: x.values)

        if len(groups) == 2:
            g1, g2 = groups.iloc[0], groups.iloc[1]

            cohen_d = _calculate_effect_size(g1, g2, method="cohen_d")
            cliff_delta = _calculate_effect_size(g1, g2, method="cliff_delta")

            cohen_d_list.append(cohen_d)
            cliff_delta_list.append(cliff_delta)
        else:
            cohen_d_list.append(np.nan)
            cliff_delta_list.append(np.nan)

    stat_df["cohen_d"] = cohen_d_list
    stat_df["cliff_delta"] = cliff_delta_list

    log.info("Effect sizes (Cohen's d and Cliff's Delta) calculated.")
    return stat_df


# ================= Data Export =================


def export_analysis_data(
    prop_df: pd.DataFrame,
    count_df: pd.DataFrame,
    stat_df: pd.DataFrame,
    out_dir: Path,
    format: str = "csv",
):
    """
    Export all analysis data for downstream use (🆕 NEW).

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix
    count_df : pd.DataFrame
        Count matrix
    stat_df : pd.DataFrame
        Statistical results
    out_dir : Path
        Output directory
    format : str
        Export format ('csv', 'excel', 'parquet')

    Examples
    --------
    >>> export_analysis_data(prop_df, count_df, stat_df, './results', format='excel')
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if format == "csv":
        prop_df.to_csv(out_dir / "proportions.csv")
        count_df.to_csv(out_dir / "counts.csv")
        if not stat_df.empty:
            stat_df.to_csv(out_dir / "statistics.csv", index=False)

    elif format == "excel":
        with pd.ExcelWriter(out_dir / "proportion_analysis.xlsx") as writer:
            prop_df.to_excel(writer, sheet_name="Proportions")
            count_df.to_excel(writer, sheet_name="Counts")
            if not stat_df.empty:
                stat_df.to_excel(writer, sheet_name="Statistics", index=False)

    elif format == "parquet":
        prop_df.to_parquet(out_dir / "proportions.parquet")
        count_df.to_parquet(out_dir / "counts.parquet")
        if not stat_df.empty:
            stat_df.to_parquet(out_dir / "statistics.parquet")

    else:
        raise ValueError(f"Unsupported format: {format}")

    log.info(f"Data exported to {out_dir} in {format} format")


# ================= Plotting Functions =================


def plot_cell_counts(
    count_df: pd.DataFrame,
    sample_to_cond: Optional[pd.Series],
    mode: str = "group",
    sample_style: str = "grouped",
    ct_palette: Optional[Dict] = None,
    sample_palette: Optional[Dict] = None,
    cond_palette: Optional[Dict] = None,
    ct_order: Optional[List] = None,
    cond_order: Optional[List] = None,
    out_dir: Optional[Path] = None,
    figsize: Tuple[float, float] = (10, 6),
):
    """
    Plot total cell counts per sample or group.

    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : Optional[pd.Series]
        Sample to condition mapping
    mode : str
        'group' (aggregated) or 'sample' (individual samples)
    sample_style : str
        'stacked' or 'grouped' (only for mode='sample')
    ct_palette : Optional[Dict]
        Color palette for cell types
    sample_palette : Optional[Dict]
        Color palette for samples
    cond_palette : Optional[Dict]
        Color palette for conditions
    ct_order : Optional[List]
        Order of cell types
    cond_order : Optional[List]
        Order of conditions
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size
    """
    fig, ax = plt.subplots(figsize=figsize)

    # --- Mode: Group (Aggregated) ---
    if mode == "group" and sample_to_cond is not None:
        grouped = count_df.groupby(sample_to_cond).sum().T
        df_long = grouped.reset_index().melt(
            id_vars="celltype", var_name="condition", value_name="count"
        )

        sns.barplot(
            data=df_long,
            x="celltype",
            y="count",
            hue="condition",
            palette=cond_palette,
            order=ct_order,
            hue_order=cond_order,
            ax=ax,
        )

        for c in ax.containers:
            try:
                ax.bar_label(c, fmt="%d", padding=2, fontsize=8)
            except:
                pass

        ax.set_title("Total Cell Counts (Aggregated by Condition)", fontweight="bold")
        ax.set_xlabel("Cell Type")
        ax.set_ylabel("Total Count")

    # --- Mode: Sample ---
    elif mode == "sample":
        plot_df = count_df.copy()

        # Sort by condition
        if sample_to_cond is not None:
            plot_df["condition"] = plot_df.index.map(sample_to_cond)
            if cond_order:
                plot_df["condition"] = pd.Categorical(
                    plot_df["condition"], categories=cond_order, ordered=True
                )
            plot_df = plot_df.sort_values("condition").drop(columns=["condition"])

        # Reorder columns
        if ct_order:
            avail_cols = [c for c in ct_order if c in plot_df.columns]
            plot_df = plot_df[avail_cols]

        if sample_style == "stacked":
            colors = (
                [ct_palette.get(c, "gray") for c in plot_df.columns]
                if ct_palette
                else None
            )
            plot_df.plot(kind="bar", stacked=True, color=colors, ax=ax, width=0.9)
            ax.set_title("Total Cell Counts per Sample (Stacked)", fontweight="bold")
            ax.set_xlabel("Sample")
            ax.set_ylabel("Total Count")

        elif sample_style == "grouped":
            df_long = plot_df.reset_index().melt(
                id_vars=plot_df.index.name, var_name="celltype", value_name="count"
            )
            sns.barplot(
                data=df_long,
                x="celltype",
                y="count",
                hue=plot_df.index.name,
                order=ct_order,
                palette=sample_palette,
                ax=ax,
            )

            for c in ax.containers:
                try:
                    labels = [f"{int(v)}" if v > 0 else "" for v in c.datavalues]
                    ax.bar_label(c, labels=labels, padding=2, fontsize=6, rotation=90)
                except:
                    pass

            ax.set_title(
                "Cell Counts per Sample (Grouped by Cell Type)", fontweight="bold"
            )
            ax.set_xlabel("Cell Type")
            ax.set_ylabel("Count")

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    if ax.get_legend():
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False)

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / f"cell_counts_{mode}.pdf", bbox_inches="tight")

    plt.close(fig)
    return fig


def plot_proportion_bar(
    prop_df: pd.DataFrame,
    sample_to_cond: Optional[pd.Series],
    level: str,
    ct_palette: Optional[Dict],
    ct_order: Optional[List],
    cond_order: Optional[List],
    out_dir: Optional[Path],
    figsize: Tuple[float, float] = (5, 9),
):
    """
    Plot stacked bar chart of cell type proportions.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : Optional[pd.Series]
        Sample to condition mapping
    level : str
        'sample' or 'group'
    ct_palette : Optional[Dict]
        Color palette for cell types
    ct_order : Optional[List]
        Order of cell types
    cond_order : Optional[List]
        Order of conditions
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size
    """
    fig, ax = plt.subplots(figsize=figsize)
    plot_df = prop_df.copy()

    # Reorder columns
    if ct_order:
        available_cts = [c for c in ct_order if c in plot_df.columns]
        plot_df = plot_df[available_cts]

    colors = (
        [ct_palette.get(c, "gray") for c in plot_df.columns] if ct_palette else None
    )

    if level == "sample":
        if sample_to_cond is not None:
            plot_df["grp"] = plot_df.index.map(sample_to_cond)
            if cond_order:
                plot_df["grp"] = pd.Categorical(
                    plot_df["grp"], categories=cond_order, ordered=True
                )
            plot_df = plot_df.sort_values("grp").drop(columns=["grp"])

        plot_df.plot(kind="bar", stacked=True, color=colors, width=0.9, ax=ax)
        ax.set_xlabel("Sample")
        ax.set_title("Cell Type Proportions per Sample", fontweight="bold")

    elif level == "group" and sample_to_cond is not None:
        grp_mean = plot_df.groupby(sample_to_cond).mean()
        if cond_order:
            grp_mean = grp_mean.reindex(cond_order)

        grp_mean.plot(kind="bar", stacked=True, color=colors, width=0.8, ax=ax)
        ax.set_xlabel("Condition")
        ax.set_title("Cell Type Proportions per Condition", fontweight="bold")

    ax.set_ylabel("Proportion")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc="upper left", title="Cell Type", frameon=False
    )

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / f"proportion_bar_{level}.pdf", bbox_inches="tight")

    plt.close(fig)
    return fig


def plot_diff_stats(
    prop_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    stat_df: Optional[pd.DataFrame],
    cond_palette: Optional[Dict],
    ct_order: Optional[List],
    cond_order: Optional[List],
    out_dir: Optional[Path],
    figsize: Tuple[float, float] = (8, 6),
):
    """
    Plot barplot with error bars and statistical significance brackets.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    stat_df : Optional[pd.DataFrame]
        Statistical test results
    cond_palette : Optional[Dict]
        Color palette for conditions
    ct_order : Optional[List]
        Order of cell types
    cond_order : Optional[List]
        Order of conditions
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size
    """
    df_long = prop_df.reset_index().melt(
        id_vars=prop_df.index.name, var_name="celltype", value_name="proportion"
    )
    df_long["condition"] = df_long[prop_df.index.name].map(sample_to_cond)

    fig, ax = plt.subplots(figsize=figsize)

    # Plot barplot with error bars
    bar_width = 0.8
    sns.barplot(
        data=df_long,
        x="celltype",
        y="proportion",
        hue="condition",
        order=ct_order,
        hue_order=cond_order,
        palette=cond_palette,
        errorbar=("se", 1),
        capsize=0.1,
        width=bar_width,
        ax=ax,
    )

    if stat_df is not None and not stat_df.empty:
        pval_map = (
            stat_df.set_index("celltype")["padj"].to_dict()
            if "padj" in stat_df.columns
            else {}
        )

        # Calculate error bar tops (Mean + SEM)
        stats_summary = df_long.groupby(["celltype", "condition"])["proportion"].agg(
            ["mean", "sem"]
        )
        stats_summary["top"] = stats_summary["mean"] + stats_summary["sem"]

        # Find the max height per celltype
        max_heights = stats_summary.groupby("celltype")["top"].max()

        # Determine layout parameters
        plotted_cts = ct_order if ct_order else sorted(df_long["celltype"].unique())
        num_conditions = df_long["condition"].nunique()
        y_limit = df_long["proportion"].max()

        # Only draw brackets if we have exactly 2 conditions
        if num_conditions == 2:
            offset_center = bar_width / num_conditions / 2

            for i, ct in enumerate(plotted_cts):
                if ct not in pval_map:
                    continue

                p_val = pval_map[ct]
                sig = _get_sig_stars(p_val)

                # Determine Y positions
                y_bar_top = max_heights.get(ct, 0)
                if y_bar_top == 0:
                    y_bar_top = y_limit * 0.05

                gap_line = _calculate_bracket_height(ax, df_long["proportion"].values)
                gap_text = gap_line * 0.4

                y_line = y_bar_top + gap_line
                y_text = y_line + gap_text

                # X positions for bracket
                x_left = i - offset_center
                x_right = i + offset_center

                # Draw bracket line
                ax.plot([x_left, x_right], [y_line, y_line], color="black", lw=1.5)

                # Draw text
                font_weight = "bold" if sig != "ns" else "normal"
                font_size = 12 if sig != "ns" else 10

                ax.text(
                    i,
                    y_text,
                    sig,
                    ha="center",
                    va="bottom",
                    color="black",
                    fontsize=font_size,
                    fontweight=font_weight,
                )
        else:
            # Fallback for > 2 groups: floating text
            for i, ct in enumerate(plotted_cts):
                if ct not in pval_map:
                    continue
                p_val = pval_map[ct]
                sig = _get_sig_stars(p_val)
                y_pos = max_heights.get(ct, 0) + (y_limit * 0.05)
                ax.text(
                    i, y_pos, sig, ha="center", va="bottom", color="black", fontsize=10
                )

    ax.set_title("Cell Type Proportion Differences (Mean ± SEM)", fontweight="bold")
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Proportion")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc="upper left", title="Condition", frameon=False
    )

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "proportion_diff_barplot.pdf", bbox_inches="tight")

    plt.close(fig)


def plot_composition(
    count_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    cond_palette: Optional[Dict],
    ct_order: Optional[List],
    cond_order: Optional[List],
    out_dir: Optional[Path],
    figsize: Tuple[float, float] = (8, 6),
):
    """
    Plot composition: fraction of each cell type contributed by each condition.

    Parameters
    ----------
    count_df : pd.DataFrame
        Count matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    cond_palette : Optional[Dict]
        Color palette for conditions
    ct_order : Optional[List]
        Order of cell types
    cond_order : Optional[List]
        Order of conditions
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size
    """
    group_counts = count_df.groupby(sample_to_cond).sum()
    if cond_order:
        group_counts = group_counts.reindex(cond_order)

    # Normalize by column (cell type)
    comp_df = group_counts.div(group_counts.sum(axis=0), axis=1).T

    if ct_order:
        comp_df = comp_df.reindex(ct_order)

    fig, ax = plt.subplots(figsize=figsize)
    colors = (
        [cond_palette.get(c, "gray") for c in comp_df.columns] if cond_palette else None
    )

    comp_df.plot(kind="bar", stacked=True, width=0.8, color=colors, ax=ax)

    ax.set_ylabel("Fraction of Total Cells")
    ax.set_xlabel("Cell Type")
    ax.set_title("Condition Composition per Cell Type", fontweight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(
        bbox_to_anchor=(1.05, 1), loc="upper left", title="Condition", frameon=False
    )

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "composition_per_celltype.pdf", bbox_inches="tight")

    plt.close(fig)


def plot_box_summary(
    prop_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    cond_palette: Optional[Dict],
    ct_order: Optional[List],
    cond_order: Optional[List],
    out_dir: Optional[Path],
    figsize: Tuple[float, float] = (8, 6),
):
    """
    Plot boxplot summary of proportion distribution across conditions.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    cond_palette : Optional[Dict]
        Color palette for conditions
    ct_order : Optional[List]
        Order of cell types
    cond_order : Optional[List]
        Order of conditions
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size
    """
    df_long = prop_df.reset_index().melt(
        id_vars=prop_df.index.name, var_name="celltype", value_name="proportion"
    )

    hue = None
    if sample_to_cond is not None:
        df_long["condition"] = df_long[prop_df.index.name].map(sample_to_cond)
        hue = "condition"

    fig, ax = plt.subplots(figsize=figsize)

    sns.boxplot(
        data=df_long,
        x="celltype",
        y="proportion",
        hue=hue,
        order=ct_order,
        hue_order=cond_order,
        palette=cond_palette,
        ax=ax,
        fliersize=0,
    )

    sns.stripplot(
        data=df_long,
        x="celltype",
        y="proportion",
        hue=hue,
        order=ct_order,
        hue_order=cond_order,
        dodge=True,
        color="black",
        alpha=0.6,
        size=3,
        ax=ax,
        legend=False,
    )

    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Proportion")
    ax.set_title("Proportion Distribution per Condition", fontweight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "proportion_box_summary.pdf", bbox_inches="tight")

    plt.close(fig)


def plot_individual_boxplots(
    prop_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    stat_df: Optional[pd.DataFrame],
    cond_palette: Optional[Dict],
    cond_order: Optional[List],
    out_dir: Optional[Path],
):
    """
    Plot individual boxplots per cell type with significance annotations.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    stat_df : Optional[pd.DataFrame]
        Statistical test results
    cond_palette : Optional[Dict]
        Color palette for conditions
    cond_order : Optional[List]
        Order of conditions
    out_dir : Optional[Path]
        Directory to save the plots
    """
    if out_dir is None:
        return

    indiv_dir = out_dir / "individual_boxplots"
    indiv_dir.mkdir(parents=True, exist_ok=True)

    df_long = prop_df.reset_index().melt(
        id_vars=prop_df.index.name, var_name="celltype", value_name="proportion"
    )
    df_long["condition"] = df_long[prop_df.index.name].map(sample_to_cond)

    pval_map = {}
    if stat_df is not None and not stat_df.empty:
        pval_map = (
            stat_df.set_index("celltype")["padj"].to_dict()
            if "padj" in stat_df.columns
            else {}
        )

    for ct in prop_df.columns:
        fig, ax = plt.subplots(figsize=(4, 5))
        sub_df = df_long[df_long["celltype"] == ct]

        # Draw boxplot & stripplot
        sns.boxplot(
            data=sub_df,
            x="condition",
            y="proportion",
            order=cond_order,
            palette=cond_palette,
            ax=ax,
            width=0.5,
            fliersize=0,
        )

        sns.stripplot(
            data=sub_df,
            x="condition",
            y="proportion",
            order=cond_order,
            color="black",
            size=6,
            jitter=True,
            alpha=0.7,
            ax=ax,
        )

        # Annotate significance
        p_val = pval_map.get(ct, 1.0)
        sig_text = _get_sig_stars(p_val)

        # Dynamic Y calculation
        y_max = sub_df["proportion"].max()

        ylim_current = ax.get_ylim()
        y_range = ylim_current[1] - ylim_current[0]
        y_offset = y_range * 0.05

        if y_max == 0:
            y_max = y_offset * 2

        bar_y = y_max + y_offset
        text_y = bar_y + y_offset * 0.5

        # Only draw line if 2 groups
        conditions = sub_df["condition"].unique()
        if len(conditions) == 2:
            ax.plot([0, 1], [bar_y, bar_y], lw=1.5, color="black")

            if sig_text == "ns":
                ax.text(
                    0.5,
                    text_y,
                    f"p-adj = {p_val:.2e}",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )
            else:
                ax.text(
                    0.5,
                    text_y,
                    sig_text,
                    ha="center",
                    va="bottom",
                    fontsize=12,
                    weight="bold",
                )

            # Ensure plot is tall enough
            ax.set_ylim(top=text_y + y_offset)
        else:
            # Fallback for >2 groups
            ax.text(
                0.95,
                0.95,
                f"p={p_val:.2e}",
                transform=ax.transAxes,
                ha="right",
                va="top",
            )

        ax.set_title(ct, fontweight="bold")
        ax.set_ylabel("Proportion")
        ax.set_xlabel("Condition")

        plt.tight_layout()

        safe_ct = str(ct).replace("/", "_").replace(" ", "_")
        plt.savefig(indiv_dir / f"boxplot_{safe_ct}.pdf", bbox_inches="tight")
        plt.close(fig)
        
    log.info(f"Generated {len(prop_df.columns)} individual boxplots in {indiv_dir}")


def plot_proportion_shifts(
    prop_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    ct_palette: Optional[Dict],
    ct_order: Optional[List],
    cond_order: Optional[List],
    out_dir: Optional[Path],
    figsize: Tuple[float, float] = (8, 10),
):
    """
    Plot Alluvial/Sankey-style proportion shifts between two conditions.

    Uses smooth Bézier curves for better visual flow.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    ct_palette : Optional[Dict]
        Color palette for cell types
    ct_order : Optional[List]
        Order of cell types (vertical)
    cond_order : Optional[List]
        Order of conditions (horizontal)
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size
    """
    if sample_to_cond is None:
        log.warning("Alluvial plot requires condition mapping.")
        return

    # Calculate mean proportions
    df_grouped = prop_df.groupby(sample_to_cond).mean().T  # (CellType x Condition)

    # Determine conditions (left/right)
    if cond_order and len(cond_order) >= 2:
        cond1, cond2 = cond_order[0], cond_order[1]
    else:
        conditions = df_grouped.columns.unique()
        if len(conditions) < 2:
            log.warning("Alluvial plot requires at least 2 conditions.")
            return
        cond1, cond2 = conditions[0], conditions[1]

    # Filter to these two conditions
    df_grouped = df_grouped[[cond1, cond2]]

    # Determine cell type order (vertical)
    if ct_order:
        valid_order = [c for c in ct_order if c in df_grouped.index]
        df_grouped = df_grouped.reindex(valid_order)

    # Normalize to 0-1 for consistent spacing
    max_prop = max(df_grouped[cond1].sum(), df_grouped[cond2].sum())

    fig, ax = plt.subplots(figsize=figsize)

    bottom1, bottom2 = 0, 0
    x1, x2 = 0.2, 0.8
    bar_width = 0.15

    for ct in df_grouped.index:
        h1 = df_grouped.loc[ct, cond1] / max_prop
        h2 = df_grouped.loc[ct, cond2] / max_prop
        c = ct_palette.get(ct, "gray") if ct_palette else "gray"

        # Draw left and right bars
        ax.barh(
            bottom1 + h1 / 2,
            bar_width,
            height=h1,
            left=x1 - bar_width / 2,
            color=c,
            edgecolor="white",
            linewidth=1.5,
        )
        ax.barh(
            bottom2 + h2 / 2,
            bar_width,
            height=h2,
            left=x2 - bar_width / 2,
            color=c,
            edgecolor="white",
            linewidth=1.5,
        )

        # Use Bézier curves for smooth connections
        y1_start, y1_end = bottom1, bottom1 + h1
        y2_start, y2_end = bottom2, bottom2 + h2

        ctrl_x = (x1 + x2) / 2

        vertices = []
        codes = []

        # Top edge curve
        vertices.extend(
            [
                (x1 + bar_width / 2, y1_end),
                (ctrl_x, y1_end),
                (ctrl_x, y2_end),
                (x2 - bar_width / 2, y2_end),
            ]
        )
        codes.extend([MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])

        # Right edge
        vertices.append((x2 - bar_width / 2, y2_start))
        codes.append(MplPath.LINETO)

        # Bottom edge curve
        vertices.extend(
            [
                (ctrl_x, y2_start),
                (ctrl_x, y1_start),
                (x1 + bar_width / 2, y1_start),
            ]
        )
        codes.extend([MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])

        # Close path
        vertices.append((x1 + bar_width / 2, y1_end))
        codes.append(MplPath.CLOSEPOLY)

        path = MplPath(vertices, codes)
        patch = PathPatch(path, facecolor=c, alpha=0.4, edgecolor="none")
        ax.add_patch(patch)

        bottom1 += h1
        bottom2 += h2

    # Set axis labels
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([x1, x2])
    ax.set_xticklabels([cond1, cond2], fontsize=14, fontweight="bold")
    ax.set_yticks([])
    ax.set_ylabel("Proportion", fontsize=12)
    ax.set_title("Cell Type Distribution Shift", fontsize=14, fontweight="bold", pad=20)

    # Add legend
    handles = [
        plt.Rectangle((0, 0), 1, 1, fc=ct_palette.get(ct, "gray"))
        for ct in df_grouped.index
    ]
    ax.legend(
        handles,
        df_grouped.index,
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        frameon=False,
        fontsize=10,
        title="Cell Type",
    )

    # Clean spines
    for spine in ["top", "right", "bottom", "left"]:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "proportion_alluvial.pdf", dpi=300, bbox_inches="tight")

    plt.close(fig)
    return fig


def plot_proportion_heatmap(
    prop_df: pd.DataFrame,
    sample_to_cond: Optional[pd.Series],
    ct_order: Optional[List],
    cond_order: Optional[List],
    out_dir: Optional[Path],
    figsize: Tuple[float, float] = (12, 8),
):
    """
    Plot heatmap of cell type proportions across samples.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : Optional[pd.Series]
        Sample to condition mapping
    ct_order : Optional[List]
        Order of cell types
    cond_order : Optional[List]
        Order of conditions
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size
    """
    plot_df = prop_df.copy()

    # Sort samples by condition
    if sample_to_cond is not None:
        plot_df["condition"] = plot_df.index.map(sample_to_cond)
        if cond_order:
            plot_df["condition"] = pd.Categorical(
                plot_df["condition"], categories=cond_order, ordered=True
            )
        plot_df = plot_df.sort_values("condition").drop(columns=["condition"])

    # Reorder cell types
    if ct_order:
        plot_df = plot_df[[c for c in ct_order if c in plot_df.columns]]

    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        plot_df.T,
        cmap="RdYlBu_r",
        cbar_kws={"label": "Proportion"},
        linewidths=0.5,
        linecolor="lightgray",
        ax=ax,
    )

    ax.set_xlabel("Sample", fontsize=12)
    ax.set_ylabel("Cell Type", fontsize=12)
    ax.set_title("Cell Type Proportion Heatmap", fontsize=14, fontweight="bold")

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "proportion_heatmap.pdf", dpi=300, bbox_inches="tight")

    plt.close(fig)
    
    return fig


def plot_proportion_with_ci(
    prop_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    ct_order: Optional[List],
    cond_order: Optional[List],
    cond_palette: Optional[Dict],
    out_dir: Optional[Path],
    ci_level: float = 0.95,
    figsize: Tuple[float, float] = (10, 6),
):
    """
    Plot cell type proportions with confidence intervals.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    ct_order : Optional[List]
        Order of cell types
    cond_order : Optional[List]
        Order of conditions
    cond_palette : Optional[Dict]
        Color palette for conditions
    out_dir : Optional[Path]
        Directory to save the plot
    ci_level : float
        Confidence interval level (default: 0.95 for 95% CI)
    figsize : Tuple[float, float]
        Figure size
    """
    df_long = prop_df.reset_index().melt(
        id_vars=prop_df.index.name, var_name="celltype", value_name="proportion"
    )
    df_long["condition"] = df_long[prop_df.index.name].map(sample_to_cond)

    fig, ax = plt.subplots(figsize=figsize)

    # Calculate mean, SEM, and confidence intervals
    summary = (
        df_long.groupby(["celltype", "condition"])["proportion"]
        .agg(["mean", "sem", "count"])
        .reset_index()
    )

    # t-distribution confidence intervals
    alpha = 1 - ci_level
    summary["ci"] = summary.apply(
        lambda row: stats.t.ppf(1 - alpha / 2, row["count"] - 1) * row["sem"]
        if row["count"] > 1
        else 0,
        axis=1,
    )

    # Prepare positions
    celltypes = ct_order if ct_order else sorted(summary["celltype"].unique())
    conditions = cond_order if cond_order else sorted(summary["condition"].unique())

    x_positions = {ct: i for i, ct in enumerate(celltypes)}
    width = 0.8 / len(conditions)

    # Plot bars with CI
    for i, cond in enumerate(conditions):
        data = summary[summary["condition"] == cond]

        x = [
            x_positions[ct] + (i - len(conditions) / 2 + 0.5) * width
            for ct in data["celltype"]
        ]

        color = cond_palette.get(cond, None) if cond_palette else None

        ax.bar(
            x,
            data["mean"],
            width,
            yerr=data["ci"],
            label=cond,
            capsize=5,
            alpha=0.8,
            color=color,
        )

    ax.set_xticks(list(x_positions.values()))
    ax.set_xticklabels(list(x_positions.keys()), rotation=45, ha="right")
    ax.set_ylabel(f"Proportion (Mean ± {int(ci_level * 100)}% CI)")
    ax.set_xlabel("Cell Type")
    ax.set_title("Cell Type Proportions with Confidence Intervals", fontweight="bold")
    ax.legend(title="Condition", frameon=False)

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "proportion_with_ci.pdf", dpi=300, bbox_inches="tight")

    plt.close(fig)
    
    return fig


def plot_celltype_correlation(
    prop_df: pd.DataFrame,
    ct_palette: Optional[Dict],
    out_dir: Optional[Path],
    figsize: Tuple[float, float] = (10, 8),
):
    """
    Plot correlation heatmap of cell type proportions.

    Uses hierarchical clustering for ordering.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    ct_palette : Optional[Dict]
        Color palette for cell types (not used directly)
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size
    """
    corr_matrix = prop_df.corr()

    # Hierarchical clustering for ordering
    linkage_matrix = linkage(corr_matrix, method="ward")
    dendro = dendrogram(linkage_matrix, no_plot=True)
    ordered_idx = dendro["leaves"]

    corr_matrix = corr_matrix.iloc[ordered_idx, ordered_idx]

    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        corr_matrix,
        cmap="coolwarm",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.5,
        linecolor="lightgray",
        cbar_kws={"label": "Correlation"},
        annot=True,
        fmt=".2f",
        ax=ax,
    )

    ax.set_title(
        "Cell Type Proportion Correlation", fontsize=14, fontweight="bold", pad=20
    )

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "celltype_correlation.pdf", dpi=300, bbox_inches="tight")

    plt.close(fig)

    return fig


def plot_effect_size_volcano(
    stat_df: pd.DataFrame,
    effect_col: str = "cohen_d",
    pval_col: str = "padj",
    threshold_effect: float = 0.5,
    threshold_pval: float = 0.05,
    out_dir: Optional[Path] = None,
    figsize: Tuple[float, float] = (8, 6),
):
    """
    Plot volcano plot: Effect Size vs -log10(p-value) (🆕 NEW).

    Highlights significant cell types with large effect sizes.

    Parameters
    ----------
    stat_df : pd.DataFrame
        Statistical test results (must contain effect_col and pval_col)
    effect_col : str
        Column name for effect size ('cohen_d' or 'cliff_delta')
    pval_col : str
        Column name for p-values ('pvalue' or 'padj')
    threshold_effect : float
        Effect size threshold for significance (default: 0.5 for Cohen's d)
    threshold_pval : float
        P-value threshold (default: 0.05)
    ct_palette : Optional[Dict]
        Color palette for cell types (not used, reserved for future)
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size

    Examples
    --------
    >>> plot_effect_size_volcano(
    ...     stat_df,
    ...     effect_col='cohen_d',
    ...     threshold_effect=0.5,
    ...     out_dir='./results'
    ... )
    """
    fig, ax = plt.subplots(figsize=figsize)
    
    if stat_df.empty or effect_col not in stat_df.columns:
        log.warning("Cannot create volcano plot: missing data or effect size column")
        return
    
    # Comprehensive validation
    if stat_df.empty:
        log.warning("Cannot create volcano plot: empty DataFrame")
        return None

    if effect_col not in stat_df.columns:
        log.warning(f"Effect size column '{effect_col}' not found")
        return None

    if pval_col not in stat_df.columns:
        log.warning(f"P-value column '{pval_col}' not found")
        return None

    
    # Handle NaN/inf in p-values
    stat_df = stat_df.copy()
    stat_df[pval_col] = stat_df[pval_col].replace([0, np.inf], 1e-300)
    stat_df = stat_df.dropna(subset=[pval_col, effect_col])

    if stat_df.empty:
        log.warning("No valid data after filtering NaN/inf values")
        return None
    
    stat_df["-log10_padj"] = -np.log10(stat_df[pval_col])
    
    # Categorize cell types
    stat_df["category"] = "Not Significant"
    stat_df.loc[
        (stat_df[pval_col] < threshold_pval)
        & (np.abs(stat_df[effect_col]) > threshold_effect),
        "category",
    ] = "Significant + Large Effect"
    stat_df.loc[
        (stat_df[pval_col] < threshold_pval)
        & (np.abs(stat_df[effect_col]) <= threshold_effect),
        "category",
    ] = "Significant + Small Effect"
    stat_df.loc[
        (stat_df[pval_col] >= threshold_pval)
        & (np.abs(stat_df[effect_col]) > threshold_effect),
        "category",
    ] = "Large Effect Only"

    # Define colors
    colors = {
        "Significant + Large Effect": "#E74C3C",
        "Significant + Small Effect": "#F39C12",
        "Large Effect Only": "#3498DB",
        "Not Significant": "#95A5A6",
    }

    # Plot
    for cat, color in colors.items():
        subset = stat_df[stat_df["category"] == cat]
        ax.scatter(
            subset[effect_col],
            subset["-log10_padj"],
            c=color,
            label=cat,
            alpha=0.7,
            s=100,
            edgecolors="black",
            linewidths=0.5,
        )

    # Add threshold lines
    ax.axhline(
        -np.log10(threshold_pval),
        ls="--",
        color="gray",
        alpha=0.5,
        label=f"p-adj = {threshold_pval}",
    )
    ax.axvline(threshold_effect, ls="--", color="gray", alpha=0.5)
    ax.axvline(-threshold_effect, ls="--", color="gray", alpha=0.5)

    # Annotate significant cell types
    sig_df = stat_df[stat_df["category"] == "Significant + Large Effect"]
    for _, row in sig_df.iterrows():
        ax.annotate(
            row["celltype"],
            (row[effect_col], row["-log10_padj"]),
            fontsize=9,
            ha="center",
            va="bottom",
            fontweight="bold",
        )

    ax.set_xlabel(f"Effect Size ({effect_col})", fontsize=12)
    ax.set_ylabel("-log10(Adjusted p-value)", fontsize=12)
    ax.set_title("Effect Size Volcano Plot", fontsize=14, fontweight="bold")
    ax.legend(loc="upper right", frameon=False, fontsize=9)

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "volcano_effect_size.pdf", dpi=300, bbox_inches="tight")

    plt.close(fig)


def plot_celltype_variability(
    prop_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    ct_order: Optional[List] = None,
    cond_palette: Optional[Dict] = None,
    out_dir: Optional[Path] = None,
    figsize: Tuple[float, float] = (10, 6),
):
    """
    Plot coefficient of variation (CV) for each cell type (🆕 NEW).

    CV = (std / mean) × 100%

    Helps identify which cell types have high inter-sample variability.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_to_cond : pd.Series
        Sample to condition mapping
    ct_order : Optional[List]
        Order of cell types
    cond_palette : Optional[Dict]
        Color palette for conditions
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size

    Examples
    --------
    >>> plot_celltype_variability(
    ...     prop_df,
    ...     sample_to_cond,
    ...     ct_order=['CD8_T', 'CD4_T', 'B', 'NK'],
    ...     out_dir='./results'
    ... )
    """
    cv_data = []

    for cond in sample_to_cond.unique():
        cond_samples = sample_to_cond[sample_to_cond == cond].index
        cond_props = prop_df.loc[cond_samples]

        for ct in prop_df.columns:
            vals = cond_props[ct]
            mean_val = vals.mean()
            std_val = vals.std()

            cv = (std_val / mean_val * 100) if mean_val > 0 else 0

            cv_data.append(
                {
                    "celltype": ct,
                    "condition": cond,
                    "cv": cv,
                    "mean": mean_val,
                    "std": std_val,
                }
            )

    cv_df = pd.DataFrame(cv_data)

    fig, ax = plt.subplots(figsize=figsize)

    sns.barplot(
        data=cv_df,
        x="celltype",
        y="cv",
        hue="condition",
        order=ct_order,
        palette=cond_palette,
        ax=ax,
    )

    ax.set_ylabel("Coefficient of Variation (%)", fontsize=12)
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_title(
        "Inter-Sample Variability (CV) per Cell Type", fontsize=14, fontweight="bold"
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(title="Condition", frameon=False)

    plt.tight_layout()

    if out_dir:
        plt.savefig(
            out_dir / "celltype_variability_cv.pdf", dpi=300, bbox_inches="tight"
        )

    plt.close(fig)
    return fig


def plot_batch_effect(
    prop_df: pd.DataFrame,
    batch_col: pd.Series,
    ct_order: Optional[List] = None,
    out_dir: Optional[Path] = None,
    figsize: Tuple[float, float] = (12, 6),
):
    """
    Visualize batch effects on cell type proportions using PCA (🆕 NEW).

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    batch_col : pd.Series
        Batch labels for each sample
    ct_order : Optional[List]
        Order of cell types (not used directly)
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size

    Examples
    --------
    >>> batch_labels = pd.Series(
    ...     ['Batch1', 'Batch1', 'Batch2', 'Batch2'],
    ...     index=['Sample1', 'Sample2', 'Sample3', 'Sample4']
    ... )
    >>> plot_batch_effect(prop_df, batch_labels, out_dir='./results')
    """
    # PCA on proportion matrix
    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(prop_df)

    pca_df = pd.DataFrame(
        {
            "PC1": pca_coords[:, 0],
            "PC2": pca_coords[:, 1],
            "batch": prop_df.index.map(batch_col),
        }
    )

    fig, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="batch", s=100, alpha=0.7, ax=ax)

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})", fontsize=12)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})", fontsize=12)
    ax.set_title(
        "Batch Effect Visualization (PCA on Proportions)",
        fontsize=14,
        fontweight="bold",
    )
    ax.legend(title="Batch", frameon=False)

    plt.tight_layout()

    if out_dir:
        plt.savefig(out_dir / "batch_effect_pca.pdf", dpi=300, bbox_inches="tight")

    plt.close(fig)
    return fig


def plot_proportion_timeseries(
    prop_df: pd.DataFrame,
    timepoint_col: pd.Series,
    celltype: str,
    group_col: Optional[pd.Series] = None,
    timepoint_order: Optional[List[str]] = None,
    out_dir: Optional[Path] = None,
    figsize: Tuple[float, float] = (8, 5),
):
    """
    Plot cell type proportion over time (longitudinal studies) (🆕 NEW).

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    timepoint_col : pd.Series
        Sample to timepoint mapping (e.g., 'Day 0', 'Day 7')
    celltype : str
        Which cell type to track over time
    group_col : Optional[pd.Series]
        Optional grouping (e.g., treatment vs control)
    out_dir : Optional[Path]
        Directory to save the plot
    figsize : Tuple[float, float]
        Figure size

    Examples
    --------
    >>> timepoints = pd.Series(
    ...     ['Day0', 'Day7', 'Day14', 'Day0', 'Day7', 'Day14'],
    ...     index=['S1', 'S2', 'S3', 'S4', 'S5', 'S6']
    ... )
    >>> plot_proportion_timeseries(
    ...     prop_df,
    ...     timepoints,
    ...     celltype='CD8_T',
    ...     group_col=treatment_groups
    ... )
    """
    if celltype not in prop_df.columns:
        raise ValueError(f"Cell type '{celltype}' not found in proportion matrix")

    plot_data = pd.DataFrame(
        {"proportion": prop_df[celltype], "timepoint": prop_df.index.map(timepoint_col)}
    )

    if group_col is not None:
        plot_data["group"] = prop_df.index.map(group_col)
    
    # Convert to categorical with order
    if timepoint_order is None:
        try:
            timepoint_order = sorted(plot_data['timepoint'].unique(), key=_natural_sort_key)
        except:
            timepoint_order = sorted(plot_data['timepoint'].unique())
    
    plot_data['timepoint'] = pd.Categorical(
        plot_data['timepoint'], 
        categories=timepoint_order, 
        ordered=True
    )

    fig, ax = plt.subplots(figsize=figsize)

    if group_col is not None:
        sns.lineplot(
            data=plot_data,
            x="timepoint",
            y="proportion",
            hue="group",
            marker="o",
            markersize=8,
            ax=ax,
        )
    else:
        sns.lineplot(
            data=plot_data,
            x="timepoint",
            y="proportion",
            marker="o",
            markersize=8,
            ax=ax,
        )

    ax.set_xlabel("Timepoint", fontsize=12)
    ax.set_ylabel(f"{celltype} Proportion", fontsize=12)
    ax.set_title(f"{celltype} Dynamics Over Time", fontsize=14, fontweight="bold")

    if group_col is not None:
        ax.legend(title="Group", frameon=False)

    plt.tight_layout()

    if out_dir:
        safe_ct = celltype.replace("/", "_").replace(" ", "_")
        plt.savefig(out_dir / f"timeseries_{safe_ct}.pdf", dpi=300, bbox_inches="tight")

    plt.close(fig)


# ================= Main Controller =================


def _auto_configure_analysis(
    adata: AnnData, config: ProportionConfig
) -> ProportionConfig:
    """
    Automatically configure test method and plot types based on data characteristics.
    
    Logic:
    - N=1 per group: Force chi-square/fisher; disable boxplots.
    - N=2 per group: Prefer DESeq2; enable basic plots.
    - N>=3 per group: Prefer Wilcoxon/DESeq2; enable boxplots, volcano.
    - Paired data: Prefer paired tests.
    - Multi-group (>2): Prefer ANOVA.
    """
    from copy import deepcopy
    # Create a copy to avoid mutating the original
    config = deepcopy(config)
    
    # 1. Extract Metadata info
    if config.condition_col not in adata.obs:
        log.warning(f"Condition column '{config.condition_col}' not found. Skipping auto-config.")
        return config

    sample_meta = adata.obs[[config.sample_col, config.condition_col]].drop_duplicates()
    condition_counts = sample_meta[config.condition_col].value_counts()
    n_groups = len(condition_counts)
    min_reps = condition_counts.min()
    max_reps = condition_counts.max()
    
    is_paired = False
    if config.pairing_col and config.pairing_col in adata.obs.columns:
        # Check if pairing is valid (samples per subject > 1)
        pair_counts = adata.obs[[config.sample_col, config.pairing_col]].drop_duplicates()[config.pairing_col].value_counts()
        if (pair_counts > 1).all():
            is_paired = True

    log.info(f"Auto-config detected: {n_groups} groups, min reps={min_reps}, max reps={max_reps}, paired={is_paired}")

    # 2. Auto-select Test Method (if default/not specific)
    # We only override if the user hasn't explicitly set a very specific method 
    # or if the current method is invalid for the data (e.g. deseq2 on N=1).
    
    suggested_method = config.test_method
    
    if n_groups > 2:
        suggested_method = "anova"
    elif is_paired:
        suggested_method = "paired-wilcoxon" if min_reps >= 5 else "paired-t-test"
    elif min_reps == 1:
        # Force change because other methods fail mathematically
        log.warning("Detected N=1 in at least one group. Forcing statistical test to 'chi-square'.")
        suggested_method = "chi-square"
    elif min_reps == 2:
        # DESeq2 is best for low N, t-test is okay. Wilcoxon has no power at N=2.
        if config.test_method == "wilcoxon":
             log.info("N=2 is too small for Wilcoxon power. Suggesting 'deseq2' or 't-test'.")
             suggested_method = "deseq2" if HAS_DESEQ2 else "t-test"
    
    # Update method
    if getattr(config, 'auto_configure', True):  # Default to True
        if suggested_method != config.test_method:
            log.warning(
                f"Auto-config suggests '{suggested_method}' instead of '{config.test_method}' "
                f"based on data characteristics (n_groups={n_groups}, min_reps={min_reps}). "
                f"Set config.auto_configure=False to disable."
            )
            config.test_method = suggested_method

    # 3. Auto-select Plot Types
    # Start with user defaults, then filter/add
    current_plots = set(config.plot_types)
    
    # N=1 specific adjustments
    if min_reps == 1:
        if "box" in current_plots:
            log.info("Removing 'box' plot (N=1 per group makes boxplots trivial).")
            current_plots.remove("box")
        if "volcano" in current_plots:
            # Volcano needs p-values, chi-square gives them, but effect sizes might be noisy
            pass 
        # Ensure we have barplots to show the data
        current_plots.add("bar")
        current_plots.add("diff")
    
    # High N adjustments
    if min_reps >= 5:
        current_plots.add("box")
        current_plots.add("volcano")
        current_plots.add("variability") # CV plot is useful here

    # Multi-group adjustments
    if n_groups > 2:
        current_plots.add("heatmap")
        if "diff" in current_plots:
            # Diff barplot with brackets is messy for >2 groups
            # We might want to keep it but it will fallback to floating stars
            pass

    config.plot_types = list(current_plots)
    return config


def celltype_proportion_analysis(
    adata: AnnData, config: ProportionConfig
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for cell type proportion analysis.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell data object
    config : ProportionConfig
        Configuration object containing analysis parameters

    Returns
    -------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    stat_df : pd.DataFrame
        Statistical test results

    Examples
    --------
    >>> from .config import ProportionConfig
    >>> config = ProportionConfig(
    ...     celltype_col='cell_type',
    ...     sample_col='sample_id',
    ...     condition_col='condition',
    ...     test_method='wilcoxon',
    ...     plot_types=['bar', 'box', 'heatmap', 'volcano'],
    ...     out_dir='./proportion_analysis'
    ... )
    >>> prop_df, stat_df = celltype_proportion_analysis(adata, config)
    """
    # Set publication style
    _set_publication_style()
    
    config = _auto_configure_analysis(adata, config)

    out_dir = Path(config.out_dir) if config.out_dir else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    log.info(f"Starting proportion analysis with {len(config.plot_types)} plot types")
    log.info(f"Analysis level: {config.analysis_level}")

    try:
        # 1. Compute proportions
        prop_df, count_df, sample_map = compute_celltype_proportion(
            adata, config.celltype_col, config.sample_col, config.condition_col
        )

        # 2. Statistical testing
        stat_df = pd.DataFrame()
        if config.condition_col and config.test_method:
            log.info(f"Running {config.test_method} statistical test")

            # Get pairing column if needed
            pairing_col = None
            if "paired" in config.test_method.lower() and hasattr(
                config, "pairing_col"
            ):
                if config.pairing_col and config.pairing_col in adata.obs.columns:
                    pairing_col = (
                        adata.obs[[config.sample_col, config.pairing_col]]
                        .drop_duplicates(subset=[config.sample_col])
                        .set_index(config.sample_col)[config.pairing_col]
                    )

            stat_df = run_statistical_test(
                count_df,
                sample_map,
                method=config.test_method,
                pairing_col=pairing_col,
                correction_scope=getattr(config, "correction_scope", "per_test"),
            )

            if not stat_df.empty:
                # Add effect sizes (only for 2-group comparisons)
                if len(sample_map.unique()) == 2:
                    stat_df = _add_effect_sizes(stat_df, prop_df, sample_map)

                sig_count = (stat_df["padj"] < 0.05).sum()
                log.info(f"Found {sig_count} significant differences (padj < 0.05)")

                # Save results
                if out_dir:
                    stat_df.to_csv(out_dir / "statistical_results.csv", index=False)
                    log.info(
                        f"Statistical results saved to {out_dir / 'statistical_results.csv'}"
                    )

        # 3. Prepare palettes
        ct_palette = _ensure_palette(config.ct_palette, prop_df.columns, "husl")

        cond_palette = None
        if sample_map is not None:
            conditions = sample_map.unique()
            cond_palette = _ensure_palette(config.condition_palette, conditions, "Set2")

        sample_palette = None
        if config.analysis_level == "sample" or "counts" in config.plot_types:
            samples = prop_df.index
            sample_palette = _ensure_palette(config.sample_palette, samples, "tab20")

        # 4. Generate plots
        log.info(f"Generating {len(config.plot_types)} plot types...")

        for ptype in config.plot_types:
            try:
                if ptype == "counts":
                    plot_cell_counts(
                        count_df,
                        sample_map,
                        mode=config.analysis_level,
                        sample_style=config.sample_plot_style,
                        ct_palette=ct_palette,
                        sample_palette=sample_palette,
                        cond_palette=cond_palette,
                        ct_order=config.celltype_order,
                        cond_order=config.condition_order,
                        out_dir=out_dir,
                        figsize=config.figsize,
                    )

                elif ptype == "bar":
                    plot_proportion_bar(
                        prop_df,
                        sample_map,
                        level=config.analysis_level,
                        ct_palette=ct_palette,
                        ct_order=config.celltype_order,
                        cond_order=config.condition_order,
                        out_dir=out_dir,
                        figsize=config.figsize,
                    )

                elif ptype == "bar_composition":
                    if sample_map is not None:
                        plot_composition(
                            count_df,
                            sample_map,
                            cond_palette,
                            config.celltype_order,
                            config.condition_order,
                            out_dir,
                            config.figsize,
                        )

                elif ptype == "diff" and config.analysis_level == "group":
                    if sample_map is not None:
                        plot_diff_stats(
                            prop_df,
                            sample_map,
                            stat_df,
                            cond_palette,
                            config.celltype_order,
                            config.condition_order,
                            out_dir,
                            config.figsize,
                        )

                elif ptype == "box":
                    if sample_map is not None:
                        plot_box_summary(
                            prop_df,
                            sample_map,
                            cond_palette,
                            config.celltype_order,
                            config.condition_order,
                            out_dir,
                            config.figsize,
                        )

                        # Auto-call individual boxplots if grouping exists
                        if config.analysis_level == "group":
                            plot_individual_boxplots(
                                prop_df,
                                sample_map,
                                stat_df,
                                cond_palette,
                                config.condition_order,
                                out_dir,
                            )

                elif ptype == "alluvial":
                    if sample_map is not None:
                        plot_proportion_shifts(
                            prop_df,
                            sample_map,
                            ct_palette=ct_palette,
                            ct_order=config.celltype_order,
                            cond_order=config.condition_order,
                            out_dir=out_dir,
                            figsize=config.figsize,
                        )

                elif ptype == "heatmap":
                    plot_proportion_heatmap(
                        prop_df,
                        sample_map,
                        ct_order=config.celltype_order,
                        cond_order=config.condition_order,
                        out_dir=out_dir,
                        figsize=(12, 8),
                    )

                elif ptype == "ci":
                    if sample_map is not None:
                        plot_proportion_with_ci(
                            prop_df,
                            sample_map,
                            ct_order=config.celltype_order,
                            cond_order=config.condition_order,
                            cond_palette=cond_palette,
                            out_dir=out_dir,
                            ci_level=0.95,
                            figsize=config.figsize,
                        )

                elif ptype == "correlation":
                    plot_celltype_correlation(
                        prop_df, ct_palette=ct_palette, out_dir=out_dir, figsize=(10, 8)
                    )

                elif ptype == "volcano":
                    if not stat_df.empty and "cohen_d" in stat_df.columns:
                        plot_effect_size_volcano(
                            stat_df,
                            effect_col="cohen_d",
                            pval_col="padj",
                            threshold_effect=0.5,
                            threshold_pval=0.05,
                            out_dir=out_dir,
                            figsize=(8, 6),
                        )

                elif ptype == "variability":
                    if sample_map is not None:
                        plot_celltype_variability(
                            prop_df,
                            sample_map,
                            ct_order=config.celltype_order,
                            cond_palette=cond_palette,
                            out_dir=out_dir,
                            figsize=(10, 6),
                        )

                elif ptype == "batch_pca":
                    if hasattr(config, "batch_col") and config.batch_col:
                        batch_labels = (
                            adata.obs[[config.sample_col, config.batch_col]]
                            .drop_duplicates(subset=[config.sample_col])
                            .set_index(config.sample_col)[config.batch_col]
                        )
                        plot_batch_effect(
                            prop_df,
                            batch_labels,
                            ct_order=config.celltype_order,
                            out_dir=out_dir,
                            figsize=(12, 6),
                        )

                elif ptype == "timeseries":
                    if hasattr(config, "timepoint_col") and config.timepoint_col:
                        timepoint_labels = (
                            adata.obs[[config.sample_col, config.timepoint_col]]
                            .drop_duplicates(subset=[config.sample_col])
                            .set_index(config.sample_col)[config.timepoint_col]
                        )

                        # Plot for each cell type
                        for ct in prop_df.columns:
                            plot_proportion_timeseries(
                                prop_df,
                                timepoint_labels,
                                celltype=ct,
                                group_col=sample_map,
                                out_dir=out_dir,
                                figsize=(8, 5),
                            )

                else:
                    log.warning(f"Unknown plot type: {ptype}")

            except Exception as e:
                log.error(f"Failed to generate {ptype} plot: {str(e)}", exc_info=True)
                continue

        # 5. Export data (🆕 NEW)
        if out_dir and getattr(config, "export_data", True):
            export_format = getattr(config, "export_format", "csv")
            export_analysis_data(
                prop_df, count_df, stat_df, out_dir, format=export_format
            )

        log.info("Proportion analysis completed successfully")
        return prop_df, stat_df

    except Exception as e:
        log.error(f"Proportion analysis failed: {str(e)}", exc_info=True)
        raise
