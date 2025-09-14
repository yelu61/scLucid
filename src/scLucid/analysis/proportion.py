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

    Enhancements:
    - Consistent index/name handling; avoid division-by-zero.
    - Returns sample_to_condition as a Series if requested.
    """
    if celltype_col not in adata.obs.columns or group_col not in adata.obs.columns:
        raise KeyError(f"'{celltype_col}' or '{group_col}' not found in adata.obs.")
    count_df = adata.obs.groupby([group_col, celltype_col]).size().unstack(fill_value=0)
    totals = count_df.sum(axis=1).replace(0, np.nan)
    prop_df = count_df.div(totals, axis=0).fillna(0.0)
    prop_df.index.name = group_col
    prop_df.columns.name = celltype_col
    if condition_col:
        if condition_col not in adata.obs.columns:
            raise KeyError(f"'{condition_col}' not found in adata.obs.")
        sample_to_cond = (
            adata.obs[[group_col, condition_col]]
            .drop_duplicates(subset=[group_col])
            .set_index(group_col)[condition_col]
        )
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

    Enhancements:
    - Guards for two-group vs multi-group tests.
    - Handles NaN and small-sample cases gracefully.
    """
    results = []
    # Align indices
    common = prop_df.index.intersection(sample_to_cond.index)
    if len(common) == 0:
        raise ValueError("No overlapping samples between prop_df and sample_to_cond.")
    prop_df = prop_df.loc[common]
    sample_to_cond = sample_to_cond.loc[common]

    conditions = sample_to_cond.dropna().unique()
    if len(conditions) < 2:
        log.warning("Fewer than 2 conditions available; skipping tests.")
        return pd.DataFrame()

    for ct in prop_df.columns:
        group_vals = [prop_df.loc[sample_to_cond == cond, ct].dropna() for cond in conditions]
        # Require at least one observation per group
        if any(len(v) < 2 for v in group_vals):
            stat, p = np.nan, np.nan
        else:
            if test == "t" and len(conditions) == 2:
                stat, p = ttest_ind(group_vals[0], group_vals[1], equal_var=False)
            elif test == "wilcoxon" and len(conditions) == 2:
                stat, p = mannwhitneyu(group_vals[0], group_vals[1], alternative="two-sided")
            elif test == "anova" and len(conditions) >= 3:
                stat, p = f_oneway(*group_vals)
            else:
                stat, p = np.nan, np.nan
        mean_per_group = {f"mean_{cond}": float(vals.mean()) if len(vals) > 0 else np.nan for cond, vals in zip(conditions, group_vals)}
        results.append(dict(celltype=ct, stat=stat, pvalue=p, **mean_per_group))

    out = pd.DataFrame(results)
    if not out.empty:
        out = out.sort_values("pvalue", na_position="last")
    return out

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