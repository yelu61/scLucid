"""
Cell type proportion statistics and visualization.

For each sample/group/condition, computes cell type proportions,
visualizes (barplot, boxplot, dotplot), and performs group-wise statistical tests.

This module is config-driven and stores all results in adata.uns['sclucid']['analysis']['proportion'].
"""

import logging
from pathlib import Path
from typing import List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData
from scipy.stats import f_oneway, mannwhitneyu, ttest_ind

from ..utils import sanitize_for_hdf5
from .config import ProportionConfig

log = logging.getLogger(__name__)

__all__ = [
    "compute_celltype_proportion",
    "plot_celltype_proportion",
    "celltype_proportion_test",
    "celltype_proportion_analysis",
]


def compute_celltype_proportion(
    adata: AnnData,
    celltype_col: str,
    sample_col: str,
    condition_col: Optional[str] = None,
) -> Tuple[pd.DataFrame, Optional[pd.Series]]:
    """
    Computes cell type proportions per sample.

    Returns:
        A tuple of:
        - prop_df: DataFrame (index=sample, columns=celltypes)
        - sample_to_cond: Optional Series (index=sample, values=condition)
    """
    if celltype_col not in adata.obs.columns or sample_col not in adata.obs.columns:
        raise KeyError(f"'{celltype_col}' or '{sample_col}' not found in adata.obs.")
    
    # Create the contingency table
    count_df = adata.obs.groupby([sample_col, celltype_col]).size().unstack(fill_value=0)
    
    # Calculate totals, replacing 0 with NaN to avoid division by zero
    totals = count_df.sum(axis=1).replace(0, np.nan)
    
    # Calculate proportions, filling NaNs (from division by zero) with 0
    prop_df = count_df.div(totals, axis=0).fillna(0.0)
    
    prop_df.index.name = sample_col
    prop_df.columns.name = celltype_col
    
    sample_to_cond_map = None
    if condition_col:
        if condition_col not in adata.obs.columns:
            raise KeyError(f"'{condition_col}' not found in adata.obs.")
        sample_to_cond_map = (
            adata.obs[[sample_col, condition_col]]
            .drop_duplicates(subset=[sample_col])
            .set_index(sample_col)[condition_col]
        )
    
    return prop_df, sample_to_cond_map


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
    """
    if out_dir:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
    
    # --- Bar Plot (Sample-level) ---
    if plot_type == "bar":
        fig, ax = plt.subplots(figsize=figsize)
        prop_df.plot(kind="bar", stacked=True, ax=ax)
        ax.set_ylabel("Proportion")
        ax.set_xlabel(group_col)
        ax.legend(title=prop_df.columns.name, bbox_to_anchor=(1.05, 1), loc="upper left")
        plt.tight_layout()
        if out_dir:
            plt.savefig(out_dir / "proportion_barplot.png", dpi=300, bbox_inches="tight")
        else:
            plt.show()
        plt.close(fig)
        
    # --- Box/Dot Plots (Grouped-level) ---
    elif plot_type in ["box", "dot"]:
        # Melt to long-form DataFrame
        df_long = prop_df.reset_index().melt(
            id_vars=prop_df.index.name, var_name="celltype", value_name="proportion"
        )
        
        hue_col = None
        if sample_to_cond is not None:
            df_long["condition"] = df_long[group_col].map(sample_to_cond)
            hue_col = "condition"

        fig, ax = plt.subplots(figsize=figsize)
        
        if plot_type == "box":
            sns.boxplot(
                data=df_long,
                x="celltype",
                y="proportion",
                hue=hue_col,
                ax=ax
            )
            sns.stripplot( # Add points over the boxplot
                data=df_long,
                x="celltype",
                y="proportion",
                hue=hue_col,
                dodge=True,
                jitter=0.2,
                color="black",
                alpha=0.3,
                s=3,
                ax=ax,
                legend=False # Avoid duplicate legend
            )
        elif plot_type == "dot":
            sns.stripplot(
                data=df_long,
                x="celltype",
                y="proportion",
                hue=hue_col,
                dodge=True,
                jitter=True,
                ax=ax
            )
        
        ax.tick_params(axis='x', rotation=45)
        plt.tight_layout()
        if out_dir:
            plt.savefig(out_dir / f"proportion_{plot_type}plot.png", dpi=300, bbox_inches="tight")
        else:
            plt.show()
        plt.close(fig)


def celltype_proportion_test(
    prop_df: pd.DataFrame,
    sample_to_cond: pd.Series,
    config: ProportionConfig,
    sample_metadata: Optional[pd.DataFrame] = None, # From adata.obs
) -> pd.DataFrame:
    """
    Statistical testing of cell type proportions, driven by ProportionConfig.
    """
    results = []
    
    # Align indices
    common_idx = prop_df.index.intersection(sample_to_cond.index)
    if len(common_idx) == 0:
        raise ValueError("No overlapping samples between prop_df and sample_to_cond.")
    
    prop_df = prop_df.loc[common_idx]
    sample_to_cond = sample_to_cond.loc[common_idx]

    conditions = sample_to_cond.dropna().unique()
    if len(conditions) < 2:
        log.warning("Fewer than 2 conditions available; skipping tests.")
        return pd.DataFrame()

    for celltype in prop_df.columns:
        
        # --- Advanced: Linear Mixed Model (LMM) ---
        if config.test_method == "lmm":
            if config.lmm_batch_key is None:
                log.warning("LMM test selected but 'lmm_batch_key' not provided. Falling back to 'wilcoxon'.")
                config.test_method = "wilcoxon" # Fallback for this iteration
            else:
                try:
                    import statsmodels.formula.api as smf
                    
                    lmm_data = pd.DataFrame({
                        'proportion': prop_df[celltype],
                        'condition': sample_to_cond
                    })
                    
                    # Align and add batch/covariates from full metadata
                    if sample_metadata is None or config.lmm_batch_key not in sample_metadata.columns:
                        raise ValueError(f"batch_key '{config.lmm_batch_key}' not found in sample_metadata (adata.obs).")
                    
                    lmm_data = lmm_data.join(sample_metadata)
                    
                    # Build formula
                    formula = config.lmm_formula
                    if formula is None:
                        formula = f"proportion ~ C(condition, Treatment(reference='{conditions[0]}'))"
                    
                    model = smf.mixedlm(formula, lmm_data, groups=lmm_data[config.lmm_batch_key])
                    result_lmm = model.fit()
                    
                    # Extract p-value for the condition effect
                    pval_col = f"C(condition, Treatment(reference='{conditions[0]}'))[T.{conditions[1]}]"
                    pval = result_lmm.pvalues.get(pval_col, np.nan)
                    stat = result_lmm.tvalues.get(pval_col, np.nan)
                    
                    results.append({
                        'celltype': celltype,
                        'test': 'lmm',
                        'stat': stat,
                        'pvalue': pval,
                        **{f'mean_{c}': prop_df.loc[sample_to_cond == c, celltype].mean() for c in conditions}
                    })
                    continue # Move to next celltype

                except ImportError:
                    log.warning("statsmodels not installed. LMM failed. Falling back to 'wilcoxon'.")
                    config.test_method = "wilcoxon" # Fallback
                except Exception as e:
                    log.warning(f"LMM failed for {celltype}: {e}. Falling back to 'wilcoxon'.")
                    config.test_method = "wilcoxon" # Fallback

        # --- Standard Statistical Tests ---
        group_vals = [
            prop_df.loc[sample_to_cond == cond, celltype].dropna() for cond in conditions
        ]
        
        stat, p = np.nan, np.nan
        test_used = config.test_method

        if any(len(v) < 2 for v in group_vals):
            log.warning(f"Skipping {celltype}: not enough data points in one group.")
        else:
            try:
                if config.test_method == "t-test" and len(conditions) == 2:
                    stat, p = ttest_ind(group_vals[0], group_vals[1], equal_var=False)
                elif config.test_method == "wilcoxon" and len(conditions) == 2:
                    stat, p = mannwhitneyu(group_vals[0], group_vals[1], alternative="two-sided")
                elif config.test_method == "anova" and len(conditions) >= 3:
                    stat, p = f_oneway(*group_vals)
                else:
                    log.warning(f"Test '{config.test_method}' not applicable for {len(conditions)} groups. Skipping {celltype}.")
                    test_used = "skipped"
            except Exception as e:
                log.error(f"Statistical test failed for {celltype}: {e}")

        mean_per_group = {
            f"mean_{cond}": float(vals.mean()) if len(vals) > 0 else np.nan
            for cond, vals in zip(conditions, group_vals)
        }
        results.append(dict(celltype=celltype, test=test_used, stat=stat, pvalue=p, **mean_per_group))

    out_df = pd.DataFrame(results)
    if not out_df.empty:
        out_df = out_df.sort_values("pvalue", na_position="last")
    return out_df


def celltype_proportion_analysis(
    adata: AnnData,
    config: ProportionConfig
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """
    One-stop, config-driven cell type proportion analysis.
    
    Computes proportions, generates plots, performs statistical tests,
    and stores all results in adata.uns.

    Args:
        adata: AnnData object.
        config: ProportionConfig object with all parameters.

    Returns:
        A tuple of (prop_df, stat_df).
    """
    log.info("Starting config-driven cell type proportion analysis...")
    
    cfg = config
    out_dir = Path(cfg.out_dir) if cfg.out_dir else None

    # --- 1. Compute Proportions ---
    prop_df, sample_to_cond = compute_celltype_proportion(
        adata, cfg.celltype_col, cfg.sample_col, cfg.condition_col
    )

    # --- 2. Save Proportion Table ---
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)
        prop_df.to_csv(out_dir / "proportion_table.csv")
        if sample_to_cond is not None:
            sample_to_cond.to_csv(out_dir / "sample_to_condition_map.csv")

    # --- 3. Draw Plots ---
    for ptype in cfg.plot_types:
        plot_celltype_proportion(
            prop_df,
            sample_to_cond,
            plot_type=ptype,
            out_dir=out_dir,
            group_col=cfg.sample_col,
            figsize=cfg.figsize
        )
    log.info(f"Generated plots: {', '.join(cfg.plot_types)}")

    # --- 4. Statistical Test ---
    stat_df = None
    if cfg.condition_col:
        log.info(f"Running statistical test: '{cfg.test_method}'")
        stat_df = celltype_proportion_test(
            prop_df, 
            sample_to_cond, 
            config=cfg,
            sample_metadata=adata.obs # Pass full metadata for LMM
        )
        if out_dir:
            stat_df.to_csv(out_dir / "proportion_stats.csv", index=False)
        log.info(f"Statistical analysis complete. Top result:\n{stat_df.head(1)}")

    # --- 5. Store in .uns ---
    uns_key = "proportion"
    store_root = (
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault(uns_key, {})
    )
    store_root["prop_df"] = prop_df
    store_root["sample_map"] = sample_to_cond
    store_root["stat_df"] = stat_df
    store_root["config"] = sanitize_for_hdf5(cfg.to_dict())
    
    log.info(f"Proportion analysis complete. Results stored in .uns['sclucid']['analysis']['{uns_key}']")

    return prop_df, stat_df