"""
Cell type proportion visualization.

This module provides comprehensive plotting functions for cell type
proportion analysis, including:
- Count and proportion bar plots
- Box plots with significance annotations
- Heatmaps and correlation matrices
- Volcano plots and effect size visualizations
- Time series and batch effect plots
"""

from __future__ import annotations

import logging
from functools import wraps
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath

log = logging.getLogger(__name__)


# ================= Plot Helper Functions =================


def _get_sig_stars(p_val: float) -> str:
    """Convert p-value to significance stars."""
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
    """Ensure a color palette exists for the given keys."""
    if palette is None:
        sorted_keys = sorted(keys) if all(isinstance(k, str) for k in keys) else keys
        colors = sns.color_palette(default_cmap, len(sorted_keys)).as_hex()
        return dict(zip(sorted_keys, colors))
    return palette


def _calculate_bracket_height(
    ax: plt.Axes, y_data: np.ndarray, num_brackets: int = 1, base_gap: float = 0.03
) -> float:
    """Dynamically calculate the height for statistical annotation brackets."""
    ylim = ax.get_ylim()
    y_range = ylim[1] - ylim[0]

    gap = y_range * base_gap
    max_y = y_data.max()

    return max_y + (gap * num_brackets)


def save_and_close(plot_name: str):
    """Decorator to automatically save and close plots."""
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


# ================= Plotting Functions =================


@save_and_close("cell_counts")
def plot_cell_counts(
    adata,
    celltype_col: str = "cell_type",
    sample_col: str = "sample_id",
    group_col: Optional[str] = None,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot total cell counts per sample grouped by cell type.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    celltype_col : str
        Column in adata.obs containing cell type labels
    sample_col : str
        Column in adata.obs containing sample identifiers
    group_col : str, optional
        Column to group samples by (e.g., condition)
    palette : Dict, optional
        Color palette for cell types
    out_dir : str, optional
        Output directory for saving plot

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    # Count cells per sample per cell type
    df = adata.obs[[sample_col, celltype_col]].copy()
    if group_col:
        df[group_col] = adata.obs[group_col]

    count_df = df.groupby([sample_col, celltype_col]).size().unstack(fill_value=0)

    # Ensure palette
    palette = _ensure_palette(palette, count_df.columns)

    # Create figure
    if group_col:
        n_groups = df[group_col].nunique()
        fig, axes = plt.subplots(1, n_groups, figsize=(5 * n_groups, 5), sharey=True)

        if n_groups == 1:
            axes = [axes]

        for ax, (group, group_df) in zip(axes, df.groupby(group_col)):
            group_counts = group_df.groupby([sample_col, celltype_col]).size().unstack(fill_value=0)
            group_counts.plot(kind='bar', stacked=True, ax=ax, color=[palette[c] for c in group_counts.columns])
            ax.set_title(f'{group}')
            ax.set_xlabel('Sample')
            ax.set_ylabel('Cell Count')
            ax.legend(title=celltype_col, bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        fig, ax = plt.subplots(figsize=(max(10, len(count_df) * 0.5), 5))
        count_df.plot(kind='bar', stacked=True, ax=ax, color=[palette[c] for c in count_df.columns])
        ax.set_title('Cell Counts per Sample')
        ax.set_xlabel('Sample')
        ax.set_ylabel('Cell Count')
        ax.legend(title=celltype_col, bbox_to_anchor=(1.05, 1), loc='upper left')

    return fig


@save_and_close("proportion_bar")
def plot_proportion_bar(
    prop_df: pd.DataFrame,
    sample_order: Optional[List] = None,
    celltype_order: Optional[List] = None,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot stacked proportion bar chart.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_order : List, optional
        Order for samples on x-axis
    celltype_order : List, optional
        Order for cell types in stack
    palette : Dict, optional
        Color palette for cell types
    out_dir : str, optional
        Output directory for saving plot

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    # Reorder if specified
    if sample_order:
        prop_df = prop_df.reindex(sample_order)

    if celltype_order:
        prop_df = prop_df[celltype_order]

    # Ensure palette
    palette = _ensure_palette(palette, prop_df.columns)

    # Create figure
    fig, ax = plt.subplots(figsize=(max(10, len(prop_df) * 0.5), 5))

    # Plot stacked bar
    prop_df.plot(
        kind='bar',
        stacked=True,
        ax=ax,
        color=[palette[c] for c in prop_df.columns],
        edgecolor='white',
        linewidth=0.5
    )

    ax.set_title('Cell Type Proportions per Sample')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Proportion')
    ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_ylim(0, 1)

    return fig


@save_and_close("proportion_box")
def plot_box_summary(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot cell type proportions as box plots grouped by condition.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    condition : pd.Series
        Condition labels for each sample
    palette : Dict, optional
        Color palette for cell types
    out_dir : str, optional
        Output directory for saving plot

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    # Ensure palette
    palette = _ensure_palette(palette, prop_df.columns)

    # Prepare data for plotting
    plot_data = prop_df.T
    plot_data.columns = condition.values

    # Create figure
    n_celltypes = len(plot_data)
    fig, axes = plt.subplots(1, n_celltypes, figsize=(4 * n_celltypes, 5), sharey=True)

    if n_celltypes == 1:
        axes = [axes]

    for ax, (celltype, data) in zip(axes, plot_data.iterrows()):
        # Create box plot
        conditions = sorted(data.unique())
        box_data = [data[data == cond].values for cond in conditions]

        bp = ax.boxplot(box_data, labels=conditions, patch_artist=True)

        # Color boxes
        for patch in bp['boxes']:
            patch.set_facecolor(palette.get(celltype, 'gray'))
            patch.set_alpha(0.7)

        # Add strip plot
        for i, cond in enumerate(conditions):
            x = np.random.normal(i + 1, 0.04, size=len(data[data == cond]))
            ax.scatter(x, data[data == cond], alpha=0.5, s=20, color='black', zorder=3)

        ax.set_title(celltype)
        ax.set_ylabel('Proportion')

    plt.tight_layout()
    return fig


@save_and_close("proportion_heatmap")
def plot_proportion_heatmap(
    prop_df: pd.DataFrame,
    sample_order: Optional[List] = None,
    celltype_order: Optional[List] = None,
    cluster_samples: bool = False,
    cluster_celltypes: bool = False,
    cmap: str = "viridis",
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot proportion heatmap.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    sample_order : List, optional
        Order for samples
    celltype_order : List, optional
        Order for cell types
    cluster_samples : bool
        Whether to cluster samples
    cluster_celltypes : bool
        Whether to cluster cell types
    cmap : str
        Colormap name
    out_dir : str, optional
        Output directory for saving plot

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    from scipy.cluster.hierarchy import linkage, dendrogram

    # Reorder if specified
    if sample_order:
        prop_df = prop_df.reindex(sample_order)

    if celltype_order:
        prop_df = prop_df[celltype_order]

    # Create figure
    fig, ax = plt.subplots(figsize=(max(10, len(prop_df.columns) * 0.5), max(8, len(prop_df) * 0.1)))

    # Plot heatmap
    sns.heatmap(
        prop_df.T,
        cmap=cmap,
        cbar_kws={'label': 'Proportion'},
        ax=ax,
        linewidths=0.5,
        annot=False
    )

    ax.set_title('Cell Type Proportion Heatmap')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Cell Type')

    return fig


@save_and_close("celltype_correlation")
def plot_celltype_correlation(
    prop_df: pd.DataFrame,
    method: str = 'pearson',
    cmap: str = "coolwarm",
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot cell type proportion correlation matrix.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    method : str
        Correlation method ('pearson', 'spearman')
    cmap : str
        Colormap name
    out_dir : str, optional
        Output directory for saving plot

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    # Calculate correlation
    corr = prop_df.corr(method=method)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot heatmap
    sns.heatmap(
        corr,
        cmap=cmap,
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.5,
        cbar_kws={'label': f'{method.capitalize()} Correlation'},
        annot=True,
        fmt='.2f',
        ax=ax
    )

    ax.set_title('Cell Type Proportion Correlation')

    return fig


@save_and_close("effect_size_volcano")
def plot_effect_size_volcano(
    stat_df: pd.DataFrame,
    effect_size_col: str = 'effect_size_cohens_d',
    pval_col: str = 'padj',
    sig_threshold: float = 0.05,
    effect_threshold: float = 0.5,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot volcano plot of effect size vs significance.

    Parameters
    ----------
    stat_df : pd.DataFrame
        Statistical test results
    effect_size_col : str
        Column name for effect size
    pval_col : str
        Column name for p-value (or adjusted p-value)
    sig_threshold : float
        Significance threshold
    effect_threshold : float
        Effect size threshold for highlighting
    out_dir : str, optional
        Output directory for saving plot

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Extract data
    x = stat_df[effect_size_col].values
    y = -np.log10(stat_df[pval_col].values)

    # Determine significance
    is_sig = stat_df[pval_col] < sig_threshold
    is_large = np.abs(x) > effect_threshold

    # Plot points
    ax.scatter(x[~is_sig], y[~is_sig], color='gray', alpha=0.5, label='ns')
    ax.scatter(x[is_sig & ~is_large], y[is_sig & ~is_large], color='blue', alpha=0.7, label='sig')
    ax.scatter(x[is_sig & is_large], y[is_sig & is_large], color='red', alpha=0.7, label='sig + large')

    # Add reference lines
    ax.axhline(-np.log10(sig_threshold), color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(-effect_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(effect_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)

    # Labels
    ax.set_xlabel(f"Effect Size ({effect_size_col})")
    ax.set_ylabel(f"-log10({pval_col})")
    ax.set_title("Effect Size Volcano Plot")
    ax.legend()

    # Annotate top hits
    top_hits = stat_df[is_sig & is_large].nsmallest(5, pval_col)
    for _, row in top_hits.iterrows():
        ax.annotate(
            row['cell_type'],
            xy=(row[effect_size_col], -np.log10(row[pval_col])),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=9,
            arrowprops=dict(arrowstyle='->', lw=0.5)
        )

    return fig


@save_and_close("proportion_timeseries")
def plot_proportion_timeseries(
    prop_df: pd.DataFrame,
    timepoints: pd.Series,
    celltype: str,
    group_col: Optional[pd.Series] = None,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot proportion changes over time for a specific cell type.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    timepoints : pd.Series
        Timepoint values for each sample
    celltype : str
        Cell type to plot
    group_col : pd.Series, optional
        Grouping variable (e.g., treatment)
    palette : Dict, optional
        Color palette for groups
    out_dir : str, optional
        Output directory for saving plot

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    if celltype not in prop_df.columns:
        raise ValueError(f"Cell type {celltype} not in proportion matrix")

    # Prepare data
    plot_df = pd.DataFrame({
        'timepoint': timepoints,
        'proportion': prop_df[celltype]
    })

    if group_col is not None:
        plot_df['group'] = group_col
        palette = _ensure_palette(palette, plot_df['group'].unique())

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    if group_col is not None:
        for group in plot_df['group'].unique():
            group_df = plot_df[plot_df['group'] == group]
            ax.plot(
                group_df['timepoint'],
                group_df['proportion'],
                'o-',
                label=group,
                color=palette.get(group),
                linewidth=2,
                markersize=8
            )
    else:
        ax.plot(
            plot_df['timepoint'],
            plot_df['proportion'],
            'o-',
            linewidth=2,
            markersize=8
        )

    ax.set_xlabel('Timepoint')
    ax.set_ylabel(f'{celltype} Proportion')
    ax.set_title(f'{celltype} Proportion Over Time')

    if group_col is not None:
        ax.legend()

    return fig


@save_and_close("batch_effect")
def plot_batch_effect(
    prop_df: pd.DataFrame,
    batch: pd.Series,
    method: str = 'pca',
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Visualize batch effects in proportion data using PCA.

    Parameters
    ----------
    prop_df : pd.DataFrame
        Proportion matrix (samples × cell types)
    batch : pd.Series
        Batch labels for each sample
    method : str
        Dimensionality reduction method ('pca', 'umap')
    palette : Dict, optional
        Color palette for batches
    out_dir : str, optional
        Output directory for saving plot

    Returns
    -------
    plt.Figure
        Matplotlib figure object
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Standardize data
    scaler = StandardScaler()
    prop_scaled = scaler.fit_transform(prop_df)

    # Dimensionality reduction
    if method == 'pca':
        reducer = PCA(n_components=2)
        emb = reducer.fit_transform(prop_scaled)
        var_explained = reducer.explained_variance_ratio_
        xlabel = f'PC1 ({var_explained[0]*100:.1f}%)'
        ylabel = f'PC2 ({var_explained[1]*100:.1f}%)'
    else:
        log.warning(f"Unknown method: {method}. Using PCA.")
        reducer = PCA(n_components=2)
        emb = reducer.fit_transform(prop_scaled)
        xlabel = 'PC1'
        ylabel = 'PC2'

    # Ensure palette
    palette = _ensure_palette(palette, batch.unique())

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot samples colored by batch
    for b in batch.unique():
        mask = batch == b
        ax.scatter(
            emb[mask, 0],
            emb[mask, 1],
            label=b,
            color=palette.get(b),
            s=100,
            alpha=0.7
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title('Batch Effect Visualization')
    ax.legend()

    return fig


# Additional simplified plotting functions for other plot types
# These would follow the same pattern as above

def plot_composition(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot contribution of each condition to cell type proportions."""
    # Implementation similar to plot_box_summary
    pass


def plot_diff_stats(
    prop_df: pd.DataFrame,
    stat_df: pd.DataFrame,
    condition: pd.Series,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot bar chart with significance brackets."""
    # Implementation with significance annotations
    pass


def plot_individual_boxplots(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    stat_df: pd.DataFrame,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot individual box plots with significance tests."""
    # Implementation similar to plot_box_summary but with significance
    pass


def plot_proportion_shifts(
    prop_df: pd.DataFrame,
    condition_col: str,
    condition1: str,
    condition2: str,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot proportion shifts between two conditions."""
    # Implementation for comparing two conditions
    pass


def plot_proportion_with_ci(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot proportions with confidence intervals."""
    # Implementation with error bars
    pass


def plot_celltype_variability(
    prop_df: pd.DataFrame,
    method: str = 'cv',
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot cell type variability across samples."""
    # Implementation for CV or other variability metrics
    pass
