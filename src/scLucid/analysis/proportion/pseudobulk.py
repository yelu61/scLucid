"""
Cell type proportion analysis workflow (Pseudo-bulk method).

This module orchestrates Pseudo-bulk proportion analysis by combining
statistical testing and visualization functions from submodules.

Main workflow:
- Compute proportions
- Run statistical tests
- Generate visualizations
- Export results

For detailed statistical and plotting functions, see:
- stats.py: Statistical tests and effect sizes
- plots.py: Visualization functions
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from anndata import AnnData

from .config import ProportionConfig
from .stats import (
    compute_celltype_proportion,
    run_statistical_test,
    export_analysis_data,
)
from .plots import (
    plot_cell_counts,
    plot_proportion_bar,
    plot_box_summary,
    plot_proportion_heatmap,
    plot_celltype_correlation,
    plot_effect_size_volcano,
    plot_proportion_timeseries,
    plot_batch_effect,
)

log = logging.getLogger(__name__)


# ================= Helper Functions =================


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


def _natural_sort_key(text):
    """
    Generate key for natural sorting (handles numbers in strings).

    Parameters
    ----------
    text : str
        Text string to sort

    Returns
    -------
    list
        Key for sorting
    """
    import re

    return [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', text)]


# ================= Main Workflow =================


def _auto_configure_analysis(
    adata: AnnData, config: ProportionConfig
) -> ProportionConfig:
    """
    Automatically configure test method and plot types based on data characteristics.

    Logic:
    - N=1 per group: Force chi-square; disable boxplots
    - N=2 per group: Prefer DESeq2/t-test; enable basic plots
    - N>=3 per group: Prefer Wilcoxon/DESeq2; enable boxplots, volcano
    - Paired data: Prefer paired tests
    - Multi-group (>2): Prefer ANOVA

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    config : ProportionConfig
        Configuration object

    Returns
    -------
    ProportionConfig
        Auto-configured settings
    """
    from copy import deepcopy

    config = deepcopy(config)

    # Extract metadata
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
        pair_counts = adata.obs[[config.sample_col, config.pairing_col]].drop_duplicates()[config.pairing_col].value_counts()
        if (pair_counts > 1).all():
            is_paired = True

    log.info(f"Auto-config detected: {n_groups} groups, min reps={min_reps}, max reps={max_reps}, paired={is_paired}")

    # Auto-select test method
    suggested_method = config.test_method

    if n_groups > 2:
        suggested_method = "anova"
    elif is_paired:
        suggested_method = "paired-wilcoxon" if min_reps >= 5 else "paired-t-test"
    elif min_reps == 1:
        log.warning("Detected N=1 in at least one group. Forcing statistical test to 'chi-square'.")
        suggested_method = "chi-square"
    elif min_reps == 2:
        if config.test_method == "wilcoxon":
            log.info("N=2 is too small for Wilcoxon power. Suggesting 'deseq2' or 't-test'.")
            suggested_method = "deseq2"

    # Update method
    if getattr(config, 'auto_configure', True):
        if suggested_method != config.test_method:
            log.warning(
                f"Auto-config suggests '{suggested_method}' instead of '{config.test_method}' "
                f"based on data characteristics (n_groups={n_groups}, min_reps={min_reps}). "
                f"Set config.auto_configure=False to disable."
            )
            config.test_method = suggested_method

    # Auto-select plot types
    current_plots = set(config.plot_types)

    # N=1 specific adjustments
    if min_reps == 1:
        if "box" in current_plots:
            log.info("Removing 'box' plot (N=1 per group makes boxplots trivial).")
            current_plots.remove("box")
        current_plots.add("bar")
        current_plots.add("diff")

    # High N adjustments
    if min_reps >= 5:
        current_plots.add("box")
        current_plots.add("volcano")

    # Multi-group adjustments
    if n_groups > 2:
        current_plots.add("heatmap")

    config.plot_types = list(current_plots)
    return config


def celltype_proportion_analysis(
    adata: AnnData,
    config: ProportionConfig
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Main entry point for cell type proportion analysis.

    This function orchestrates the complete analysis workflow:
    1. Compute cell type proportions
    2. Run statistical tests
    3. Generate visualizations
    4. Export results

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
    >>> from scLucid.analysis import ProportionConfig, celltype_proportion_analysis
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

    # Auto-configure
    config = _auto_configure_analysis(adata, config)

    out_dir = Path(config.out_dir) if config.out_dir else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Compute proportions
    log.info("Computing cell type proportions...")
    prop_df = compute_celltype_proportion(
        adata,
        celltype_col=config.celltype_col,
        sample_col=config.sample_col,
        normalize=True
    )

    # 2. Run statistical tests
    stat_df = pd.DataFrame()
    if config.condition_col:
        log.info(f"Running statistical tests ({config.test_method})...")

        sample_to_cond = adata.obs[config.condition_col]
        sample_to_pair = adata.obs[config.pairing_col] if config.pairing_col else None

        stat_df = run_statistical_test(
            prop_df,
            condition_col=config.condition_col,
            test_method=config.test_method,
            sample_to_cond=sample_to_cond,
            sample_to_pair=sample_to_pair
        )

        # Add effect sizes
        if config.test_method in ['t-test', 'wilcoxon', 'paired-t-test', 'paired-wilcoxon']:
            from .proportion_stats import _add_effect_sizes
            stat_df = _add_effect_sizes(
                stat_df,
                prop_df,
                sample_to_cond,
                method='cohens_d'
            )

    # 3. Generate plots
    if config.plot_types:
        log.info(f"Generating {len(config.plot_types)} plots...")

        # Prepare data
        condition = adata.obs[config.condition_col] if config.condition_col else None
        palette = config.palette if hasattr(config, 'palette') else None

        for plot_type in config.plot_types:
            try:
                if plot_type == 'counts':
                    plot_cell_counts(
                        adata,
                        celltype_col=config.celltype_col,
                        sample_col=config.sample_col,
                        group_col=condition.name if condition is not None else None,
                        palette=palette,
                        out_dir=out_dir
                    )

                elif plot_type == 'bar':
                    sample_order = sorted(prop_df.index, key=_natural_sort_key)
                    plot_proportion_bar(
                        prop_df,
                        sample_order=sample_order,
                        palette=palette,
                        out_dir=out_dir
                    )

                elif plot_type == 'box':
                    if condition is not None:
                        plot_box_summary(
                            prop_df,
                            condition=condition,
                            palette=palette,
                            out_dir=out_dir
                        )

                elif plot_type == 'heatmap':
                    celltype_order = stat_df['cell_type'].values if not stat_df.empty else None
                    plot_proportion_heatmap(
                        prop_df,
                        celltype_order=celltype_order,
                        cluster_samples=True,
                        out_dir=out_dir
                    )

                elif plot_type == 'correlation':
                    plot_celltype_correlation(
                        prop_df,
                        out_dir=out_dir
                    )

                elif plot_type == 'volcano':
                    if not stat_df.empty:
                        plot_effect_size_volcano(
                            stat_df,
                            out_dir=out_dir
                        )

                elif plot_type == 'timeseries':
                    if config.timepoint_col and config.timepoint_col in adata.obs:
                        timepoints = adata.obs[config.timepoint_col]

                        # Plot top varying cell types
                        celltype_var = prop_df.var(axis=0)
                        top_celltypes = celltype_var.nlargest(3).index.tolist()

                        for celltype in top_celltypes:
                            plot_proportion_timeseries(
                                prop_df,
                                timepoints=timepoints,
                                celltype=celltype,
                                group_col=condition,
                                palette=palette,
                                out_dir=out_dir
                            )

                elif plot_type == 'batch_pca':
                    if config.batch_col and config.batch_col in adata.obs:
                        batch = adata.obs[config.batch_col]
                        plot_batch_effect(
                            prop_df,
                            batch=batch,
                            method='pca',
                            palette=palette,
                            out_dir=out_dir
                        )

            except Exception as e:
                log.error(f"Failed to generate {plot_type} plot: {e}")

    # 4. Export data
    if out_dir:
        export_analysis_data(prop_df, stat_df, out_dir)
        log.info(f"Analysis complete. Results saved to {out_dir}")

    return prop_df, stat_df
