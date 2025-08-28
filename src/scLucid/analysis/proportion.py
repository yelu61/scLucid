"""
Cell type proportion statistics and visualization.

For each sample/group/condition, computes cell type proportions,
visualizes (barplot, boxplot, dotplot), and performs group-wise statistical tests.
"""

import logging
from typing import Optional, List, Literal, Tuple, Union
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind, mannwhitneyu, f_oneway

log = logging.getLogger(__name__)

def compute_celltype_proportion(
    adata,
    celltype_col: str = "celltype",
    group_col: str = "sample",
    condition_col: Optional[str] = None,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, pd.Series]]:
    """
    Computes cell type proportions per group.

    Args:
        adata: AnnData object.
        celltype_col: Column in .obs specifying cell type annotation.
        group_col: Column in .obs specifying group/sample.
        condition_col: Optional column in .obs specifying condition/grouping for tests.

    Returns:
        If condition_col is None:
            proportion_df: DataFrame indexed by group_col, columns = celltypes.
        If condition_col is provided:
            (proportion_df, sample_to_condition): 
                sample_to_condition is a Series mapping group_col to condition.
    """
    # Count cells in each group/celltype
    count_df = adata.obs.groupby([group_col, celltype_col]).size().unstack(fill_value=0)
    prop_df = count_df.div(count_df.sum(axis=1), axis=0)
    if condition_col:
        sample_to_cond = adata.obs.drop_duplicates(group_col)[[group_col, condition_col]].set_index(group_col)[condition_col]
        return prop_df, sample_to_cond
    else:
        return prop_df

def plot_celltype_proportion(
    prop_df: pd.DataFrame,
    sample_to_cond: Optional[pd.Series] = None,
    plot_type: Literal["bar", "box", "dot"] = "box",
    out_dir: Optional[Path] = None,
    group_col: str = "sample",
    figsize=(10, 6),
):
    """
    Generates main proportion plots for all celltypes.

    Args:
        prop_df: Cell type proportion DataFrame (index=sample, columns=celltypes).
        sample_to_cond: Optional Series mapping sample to condition.
        plot_type: "bar", "box", or "dot".
        out_dir: Optional Path to save plots.
        group_col: Name of the sample/group column.
        figsize: Figure size tuple.
    """
    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
    if plot_type == "bar":
        plt.figure(figsize=figsize)
        prop_df.plot(kind="bar", stacked=True, figsize=figsize)
        plt.ylabel("Proportion")
        plt.xlabel("Sample")
        plt.tight_layout()
        if out_dir:
            plt.savefig(out_dir / "proportion_barplot.png")
        else:
            plt.show()
        plt.close()
    else:
        # Melt to long-form DataFrame
        df_long = prop_df.reset_index().melt(id_vars=prop_df.index.name, var_name="celltype", value_name="proportion")
        if sample_to_cond is not None:
            df_long["condition"] = df_long[group_col].map(sample_to_cond)
        plt.figure(figsize=figsize)
        if plot_type == "box":
            sns.boxplot(
                data=df_long,
                x="celltype",
                y="proportion",
                hue="condition" if sample_to_cond is not None else None
            )
            plt.xticks(rotation=45)
        elif plot_type == "dot":
            sns.stripplot(
                data=df_long,
                x="celltype",
                y="proportion",
                hue="condition" if sample_to_cond is not None else None,
                dodge=True,
                jitter=True
            )
            plt.xticks(rotation=45)
        plt.tight_layout()
        if out_dir:
            plt.savefig(out_dir / f"proportion_{plot_type}plot.png")
        else:
            plt.show()
        plt.close()

def celltype_proportion_test(
    prop_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    test: Literal["t", "wilcoxon", "anova"] = "wilcoxon"
) -> pd.DataFrame:
    """
    For each celltype, performs group-wise statistical tests between conditions.

    Args:
        prop_df: Proportion DataFrame, index = sample, columns = celltypes.
        sample_to_cond: Series mapping sample to condition.
        test: Statistical test to use ("t", "wilcoxon", or "anova").

    Returns:
        DataFrame with columns: celltype, stat, pvalue, mean_{condition} for each group.
    """
    results = []
    conditions = sample_to_cond.unique()
    for ct in prop_df.columns:
        group_vals = [prop_df.loc[sample_to_cond == cond, ct].dropna() for cond in conditions]
        if test == "t" and len(conditions) == 2:
            stat, p = ttest_ind(*group_vals)
        elif test == "wilcoxon" and len(conditions) == 2:
            stat, p = mannwhitneyu(*group_vals)
        elif test == "anova" and len(conditions) > 2:
            stat, p = f_oneway(*group_vals)
        else:
            stat, p = np.nan, np.nan
        mean_per_group = {f"mean_{cond}": vals.mean() for cond, vals in zip(conditions, group_vals)}
        results.append(dict(celltype=ct, stat=stat, pvalue=p, **mean_per_group))
    return pd.DataFrame(results).sort_values("pvalue")

def celltype_proportion_analysis(
    adata,
    celltype_col: str = "celltype",
    group_col: str = "sample",
    condition_col: Optional[str] = None,
    out_dir: Optional[Path] = None,
    plot_types: List[str] = ["bar", "box", "dot"],
    test: Literal["t", "wilcoxon", "anova"] = "wilcoxon",
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """
    One-stop cell type proportion analysis: computes proportions, generates plots, performs statistical tests.

    Args:
        adata: AnnData object.
        celltype_col: Cell type annotation column in .obs.
        group_col: Sample/group column in .obs.
        condition_col: Optional column for group-wise comparison/statistics.
        out_dir: Optional Path to save results.
        plot_types: List of plot types to generate.
        test: Statistical test for group comparison.

    Returns:
        prop_df: Proportion DataFrame (index=sample, columns=celltypes).
        stat_df: Statistical test results DataFrame (if condition_col is provided), else None.
    """
    log.info("Starting cell type proportion analysis")
    if condition_col:
        prop_df, sample_to_cond = compute_celltype_proportion(adata, celltype_col, group_col, condition_col)
    else:
        prop_df = compute_celltype_proportion(adata, celltype_col, group_col)
        sample_to_cond = None

    # Save proportion table
    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        prop_df.to_csv(out_dir / "proportion_table.csv")
        if sample_to_cond is not None:
            sample_to_cond.to_csv(out_dir / "sample_to_condition.csv")
    
    # Draw plots
    for ptype in plot_types:
        plot_celltype_proportion(
            prop_df,
            sample_to_cond,
            plot_type=ptype,
            out_dir=out_dir,
            group_col=group_col
        )
    
    # Statistical test
    if sample_to_cond is not None:
        stat_df = celltype_proportion_test(prop_df, sample_to_cond, test)
        if out_dir:
            stat_df.to_csv(out_dir / "proportion_stats.csv", index=False)
        return prop_df, stat_df
    else:
        return prop_df, None