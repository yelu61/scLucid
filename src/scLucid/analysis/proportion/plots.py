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
from typing import Dict, List, Literal, Optional, Tuple

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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


def _ensure_palette(palette: Optional[Dict], keys: pd.Index, default_cmap: str = "husl") -> Dict:
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


def _resolve_order(values: pd.Index, requested: Optional[List] = None) -> List:
    """Resolve display order while keeping only present values."""
    present = list(values)
    if requested is None:
        return present
    requested_present = [value for value in requested if value in present]
    remainder = [value for value in present if value not in requested_present]
    return requested_present + remainder


def _resolve_plot_colors(columns: List[str], palette: Optional[Dict] = None) -> List:
    """Resolve ordered colors for a list of labels."""
    palette = _ensure_palette(palette, pd.Index(columns))
    return [palette.get(col, "#808080") for col in columns]


def _align_series_to_index(series: pd.Series, index: pd.Index, name: str) -> pd.Series:
    """Align a metadata series to a proportion matrix index."""
    if not isinstance(series, pd.Series):
        series = pd.Series(series, name=name)

    if series.index.is_unique and index.isin(series.index).all():
        aligned = series.reindex(index)
    elif len(series) == len(index):
        aligned = pd.Series(series.to_numpy(), index=index, name=series.name)
    else:
        raise ValueError(f"{name} must align to prop_df.index or have the same length")

    aligned.name = series.name or name
    return aligned


def _proportion_long_frame(prop_df: pd.DataFrame, condition: pd.Series) -> pd.DataFrame:
    """Convert a sample x cell-type matrix plus condition labels to long format."""
    if prop_df.empty:
        raise ValueError("prop_df is empty")

    condition = _align_series_to_index(condition, prop_df.index, "condition")
    sample_name = prop_df.index.name or "sample"
    condition_name = condition.name or "condition"

    plot_df = prop_df.copy()
    plot_df.index = plot_df.index.astype(str)
    plot_df[sample_name] = plot_df.index
    plot_df[condition_name] = condition.astype(str).to_numpy()
    return plot_df.melt(
        id_vars=[sample_name, condition_name],
        var_name="cell_type",
        value_name="proportion",
    )


def _resolve_pvalue_col(stat_df: pd.DataFrame) -> Optional[str]:
    """Return the preferred p-value column available in a stats table."""
    for col in ("padj", "pvals_adj", "qval", "pval", "p_value", "p-value"):
        if col in stat_df.columns:
            return col
    return None


def _resolve_celltype_col(stat_df: pd.DataFrame) -> str:
    """Return the cell-type column name used by a stats table."""
    for col in ("cell_type", "celltype", "Cell Type", "CellType"):
        if col in stat_df.columns:
            return col
    raise KeyError("stat_df must contain a cell type column")


def _stat_lookup(stat_df: pd.DataFrame, value_col: str) -> Dict[str, float]:
    """Build a cell-type keyed lookup for one stats column."""
    if stat_df.empty or value_col not in stat_df.columns:
        return {}
    celltype_col = _resolve_celltype_col(stat_df)
    values = pd.to_numeric(stat_df[value_col], errors="coerce")
    return dict(zip(stat_df[celltype_col].astype(str), values))


def transform_composition(
    prop_df: pd.DataFrame,
    method: Literal["proportion", "clr", "alr", "logit"] = "clr",
    pseudocount: float = 1e-6,
    reference_celltype: Optional[str] = None,
) -> pd.DataFrame:
    """
    Transform cell-type proportions for compositional visualization.

    CLR/ALR transforms reduce closed-sum artifacts when comparing sample
    compositions. They are intended for visualization and exploratory summaries,
    not as a replacement for the statistical model chosen upstream.
    """
    if prop_df.empty:
        raise ValueError("prop_df is empty")

    numeric = prop_df.apply(pd.to_numeric, errors="coerce").fillna(0)
    comp = numeric.clip(lower=0) + pseudocount
    comp = comp.div(comp.sum(axis=1), axis=0)

    if method == "proportion":
        return comp
    if method == "clr":
        log_comp = np.log(comp)
        return log_comp.sub(log_comp.mean(axis=1), axis=0)
    if method == "alr":
        if reference_celltype is None:
            reference_celltype = str(comp.mean(axis=0).idxmax())
        if reference_celltype not in comp.columns:
            raise KeyError(f"reference_celltype '{reference_celltype}' not in prop_df columns")
        transformed = np.log(comp.drop(columns=[reference_celltype]).div(comp[reference_celltype], axis=0))
        transformed.columns = [f"{col}/{reference_celltype}" for col in transformed.columns]
        return transformed
    if method == "logit":
        clipped = comp.clip(lower=pseudocount, upper=1 - pseudocount)
        return np.log(clipped / (1 - clipped))

    raise ValueError("method must be one of: 'proportion', 'clr', 'alr', 'logit'")


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
        Color palette for conditions
    out_dir : str, optional
        Output directory for saving plot

    Returns:
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
            group_counts.plot(
                kind="bar", stacked=True, ax=ax, color=[palette[c] for c in group_counts.columns]
            )
            ax.set_title(f"{group}")
            ax.set_xlabel("Sample")
            ax.set_ylabel("Cell Count")
            ax.legend(title=celltype_col, bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        fig, ax = plt.subplots(figsize=(max(10, len(count_df) * 0.5), 5))
        count_df.plot(kind="bar", stacked=True, ax=ax, color=[palette[c] for c in count_df.columns])
        ax.set_title("Cell Counts per Sample")
        ax.set_xlabel("Sample")
        ax.set_ylabel("Cell Count")
        ax.legend(title=celltype_col, bbox_to_anchor=(1.05, 1), loc="upper left")

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

    Returns:
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
        kind="bar",
        stacked=True,
        ax=ax,
        color=[palette[c] for c in prop_df.columns],
        edgecolor="white",
        linewidth=0.5,
    )

    ax.set_title("Cell Type Proportions per Sample")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Proportion")
    ax.legend(title="Cell Type", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.set_ylim(0, 1)

    return fig


@save_and_close("grouped_celltype_counts")
def plot_grouped_celltype_counts(
    count_df: pd.DataFrame,
    group_col: str = "group",
    celltype_col: str = "cell_type",
    count_col: str = "count",
    group_order: Optional[List] = None,
    celltype_order: Optional[List] = None,
    palette: Optional[Dict] = None,
    annotate: bool = False,
    figsize: Tuple[float, float] = (12, 6),
    title: str = "Cell Counts by Group",
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot grouped cell-type counts from a long-format count table.

    Parameters
    ----------
    count_df : pd.DataFrame
        Long-format table with group, cell type, and count columns.
    group_col : str
        Grouping column on the x-axis.
    celltype_col : str
        Cell-type column used as hue.
    count_col : str
        Count column.
    group_order : list, optional
        Display order for groups.
    celltype_order : list, optional
        Display order for cell types in the legend.
    palette : dict, optional
        Color map keyed by cell type.
    annotate : bool
        If True, add count labels above non-zero bars.
    figsize : tuple
        Figure size.
    title : str
        Plot title.
    out_dir : str, optional
        Output directory for saving plot.
    """
    required = {group_col, celltype_col, count_col}
    missing = required - set(count_df.columns)
    if missing:
        raise KeyError(f"count_df missing required columns: {sorted(missing)}")

    plot_df = count_df.copy()
    plot_df[group_col] = plot_df[group_col].astype(str)
    plot_df[celltype_col] = plot_df[celltype_col].astype(str)
    plot_df[count_col] = pd.to_numeric(plot_df[count_col], errors="coerce").fillna(0)

    group_order = _resolve_order(pd.Index(plot_df[group_col].unique()), group_order)
    celltype_order = _resolve_order(pd.Index(plot_df[celltype_col].unique()), celltype_order)
    palette = _ensure_palette(palette, pd.Index(celltype_order))

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(
        data=plot_df,
        x=group_col,
        y=count_col,
        hue=celltype_col,
        order=group_order,
        hue_order=celltype_order,
        palette=palette,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Cell Count")
    ax.set_title(title)
    ax.legend(title=celltype_col, bbox_to_anchor=(1.02, 1), loc="upper left")

    if annotate:
        for patch in ax.patches:
            height = patch.get_height()
            if pd.notna(height) and height > 0:
                ax.text(
                    x=patch.get_x() + patch.get_width() / 2,
                    y=height,
                    s=f"{int(round(height))}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    rotation=90,
                )

    return fig


@save_and_close("grouped_proportion_bar")
def plot_grouped_proportion_bar(
    group_props: pd.DataFrame,
    group_order: Optional[List] = None,
    celltype_order: Optional[List] = None,
    palette: Optional[Dict] = None,
    figsize: Tuple[float, float] = (9, 6),
    title: str = "Cell Type Composition by Group",
    xlabel: str = "Group",
    ylabel: str = "Proportion",
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot stacked cell-type proportions from a group x cell-type matrix.
    """
    if group_props.empty:
        raise ValueError("group_props is empty")

    plot_df = group_props.copy()
    plot_df.index = plot_df.index.astype(str)
    plot_df.columns = plot_df.columns.astype(str)

    resolved_groups = _resolve_order(pd.Index(plot_df.index), group_order)
    resolved_celltypes = _resolve_order(pd.Index(plot_df.columns), celltype_order)
    plot_df = plot_df.loc[resolved_groups, resolved_celltypes]

    fig, ax = plt.subplots(figsize=figsize)
    plot_df.plot(
        kind="bar",
        stacked=True,
        color=_resolve_plot_colors(list(plot_df.columns), palette),
        edgecolor="white",
        linewidth=0.5,
        ax=ax,
    )
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(title="Cell Type", loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.tick_params(axis="x", rotation=45)
    return fig


@save_and_close("celltype_alluvial")
def plot_celltype_alluvial(
    group_props: pd.DataFrame,
    celltype_order: Optional[List] = None,
    palette: Optional[Dict] = None,
    figsize: Tuple[float, float] = (12, 7),
    title: str = "Cell Type Alluvial",
    bar_width: float = 0.35,
    band_alpha: float = 0.28,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """
    Plot an alluvial-style stacked composition chart from group proportions.

    Parameters
    ----------
    group_props : pd.DataFrame
        Matrix indexed by group with cell types as columns. Each row should sum
        approximately to 1.
    celltype_order : list, optional
        Display order for cell types.
    palette : dict, optional
        Color map keyed by cell type.
    figsize : tuple
        Figure size.
    title : str
        Plot title.
    bar_width : float
        Width of each stacked bar.
    band_alpha : float
        Alpha value for connecting ribbons.
    out_dir : str, optional
        Output directory for saving plot.
    """
    if group_props.empty:
        raise ValueError("group_props is empty")

    plot_df = group_props.copy()
    plot_df.index = plot_df.index.astype(str)
    plot_df.columns = plot_df.columns.astype(str)
    resolved_celltypes = _resolve_order(pd.Index(plot_df.columns), celltype_order)
    plot_df = plot_df[resolved_celltypes]

    groups = list(plot_df.index)
    x = np.arange(len(groups))
    palette = _ensure_palette(palette, pd.Index(resolved_celltypes))

    fig, ax = plt.subplots(figsize=figsize)
    bottoms = dict.fromkeys(groups, 0.0)
    yspans: Dict[str, Dict[str, Tuple[float, float]]] = {group: {} for group in groups}

    for celltype in resolved_celltypes:
        color = palette.get(celltype, "#808080")
        for idx, group in enumerate(groups):
            height = float(plot_df.loc[group, celltype]) if celltype in plot_df.columns else 0.0
            y0 = bottoms[group]
            y1 = y0 + height
            yspans[group][celltype] = (y0, y1)
            ax.bar(
                x[idx],
                height,
                bottom=y0,
                width=bar_width,
                color=color,
                edgecolor="white",
                linewidth=0.6,
            )
            bottoms[group] = y1

    for idx in range(len(groups) - 1):
        left_group, right_group = groups[idx], groups[idx + 1]
        x_left = x[idx] + bar_width / 2
        x_right = x[idx + 1] - bar_width / 2
        for celltype in resolved_celltypes:
            y0_left, y1_left = yspans[left_group][celltype]
            y0_right, y1_right = yspans[right_group][celltype]
            if (y1_left - y0_left) <= 0 and (y1_right - y0_right) <= 0:
                continue
            polygon = patches.Polygon(
                [
                    (x_left, y0_left),
                    (x_right, y0_right),
                    (x_right, y1_right),
                    (x_left, y1_left),
                ],
                closed=True,
                facecolor=palette.get(celltype, "#808080"),
                edgecolor="none",
                alpha=band_alpha,
            )
            ax.add_patch(polygon)

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=45, ha="right")
    ax.set_xlim(-0.6, len(groups) - 1 + 0.6)
    ax.set_ylim(0, max(1.0, float(plot_df.sum(axis=1).max())))
    ax.set_ylabel("Proportion")
    ax.set_title(title)
    handles = [
        patches.Patch(color=palette.get(ct, "#808080"), label=ct) for ct in resolved_celltypes
    ]
    ax.legend(handles=handles, title="Cell Type", loc="center left", bbox_to_anchor=(1.02, 0.5))
    return fig


@save_and_close("proportion_box")
def plot_box_summary(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    palette: Optional[Dict] = None,
    celltype_order: Optional[List] = None,
    condition_order: Optional[List] = None,
    show_points: bool = True,
    figsize: Optional[Tuple[float, float]] = None,
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

    Returns:
    -------
    plt.Figure
        Matplotlib figure object
    """
    plot_df = _proportion_long_frame(prop_df, condition)
    condition_col = (
        condition.name if isinstance(condition, pd.Series) and condition.name else "condition"
    )
    celltype_order = _resolve_order(pd.Index(plot_df["cell_type"].unique()), celltype_order)
    condition_order = _resolve_order(pd.Index(plot_df[condition_col].unique()), condition_order)
    palette = _ensure_palette(palette, pd.Index(condition_order))

    if figsize is None:
        figsize = (max(8, len(celltype_order) * 1.2), 5)

    fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(
        data=plot_df,
        x="cell_type",
        y="proportion",
        hue=condition_col,
        order=celltype_order,
        hue_order=condition_order,
        palette=palette,
        ax=ax,
    )
    if show_points:
        sns.stripplot(
            data=plot_df,
            x="cell_type",
            y="proportion",
            hue=condition_col,
            order=celltype_order,
            hue_order=condition_order,
            dodge=True,
            palette=dict.fromkeys(condition_order, "black"),
            alpha=0.45,
            size=3,
            legend=False,
            ax=ax,
        )
        ax.legend(title=condition_col)
    else:
        ax.legend(title=condition_col)

    ax.set_title("Cell Type Proportions by Condition")
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Proportion")
    ax.tick_params(axis="x", rotation=45)
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

    Returns:
    -------
    plt.Figure
        Matplotlib figure object
    """
    # Reorder if specified
    if sample_order:
        prop_df = prop_df.reindex(sample_order)

    if celltype_order:
        present_celltypes = [celltype for celltype in celltype_order if celltype in prop_df.columns]
        prop_df = prop_df[present_celltypes]

    if cluster_samples or cluster_celltypes:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import pdist

        if cluster_samples and len(prop_df) > 1:
            sample_dist = pdist(prop_df.fillna(0).to_numpy(), metric="euclidean")
            prop_df = prop_df.iloc[leaves_list(linkage(sample_dist, method="average"))]
        if cluster_celltypes and prop_df.shape[1] > 1:
            celltype_dist = pdist(prop_df.fillna(0).T.to_numpy(), metric="euclidean")
            prop_df = prop_df.iloc[:, leaves_list(linkage(celltype_dist, method="average"))]

    # Create figure
    fig, ax = plt.subplots(
        figsize=(max(10, len(prop_df.columns) * 0.5), max(8, len(prop_df) * 0.1))
    )

    # Plot heatmap
    sns.heatmap(
        prop_df.T, cmap=cmap, cbar_kws={"label": "Proportion"}, ax=ax, linewidths=0.5, annot=False
    )

    ax.set_title("Cell Type Proportion Heatmap")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Cell Type")

    return fig


@save_and_close("composition_transform_heatmap")
def plot_composition_transform_heatmap(
    prop_df: pd.DataFrame,
    transform: Literal["proportion", "clr", "alr", "logit"] = "clr",
    sample_order: Optional[List] = None,
    celltype_order: Optional[List] = None,
    cmap: str = "RdBu_r",
    center: Optional[float] = 0,
    pseudocount: float = 1e-6,
    reference_celltype: Optional[str] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot a transformed composition heatmap, usually CLR for group comparison."""
    transformed = transform_composition(
        prop_df,
        method=transform,
        pseudocount=pseudocount,
        reference_celltype=reference_celltype,
    )

    if sample_order:
        transformed = transformed.reindex(sample_order)
    if celltype_order:
        present = [celltype for celltype in celltype_order if celltype in transformed.columns]
        transformed = transformed[present]

    fig, ax = plt.subplots(
        figsize=(max(10, len(transformed.columns) * 0.6), max(6, len(transformed) * 0.18))
    )
    sns.heatmap(
        transformed.T,
        cmap=cmap,
        center=center,
        cbar_kws={"label": transform.upper()},
        linewidths=0.5,
        ax=ax,
    )
    ax.set_title(f"{transform.upper()} Cell Type Composition Heatmap")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Cell Type")
    return fig


@save_and_close("celltype_correlation")
def plot_celltype_correlation(
    prop_df: pd.DataFrame,
    method: str = "pearson",
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

    Returns:
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
        cbar_kws={"label": f"{method.capitalize()} Correlation"},
        annot=True,
        fmt=".2f",
        ax=ax,
    )

    ax.set_title("Cell Type Proportion Correlation")

    return fig


@save_and_close("effect_size_volcano")
def plot_effect_size_volcano(
    stat_df: pd.DataFrame,
    effect_size_col: str = "effect_size_cohens_d",
    pval_col: str = "padj",
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

    Returns:
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
    ax.scatter(x[~is_sig], y[~is_sig], color="gray", alpha=0.5, label="ns")
    ax.scatter(x[is_sig & ~is_large], y[is_sig & ~is_large], color="blue", alpha=0.7, label="sig")
    ax.scatter(
        x[is_sig & is_large], y[is_sig & is_large], color="red", alpha=0.7, label="sig + large"
    )

    # Add reference lines
    ax.axhline(-np.log10(sig_threshold), color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(-effect_threshold, color="black", linestyle="--", linewidth=1, alpha=0.5)
    ax.axvline(effect_threshold, color="black", linestyle="--", linewidth=1, alpha=0.5)

    # Labels
    ax.set_xlabel(f"Effect Size ({effect_size_col})")
    ax.set_ylabel(f"-log10({pval_col})")
    ax.set_title("Effect Size Volcano Plot")
    ax.legend()

    # Annotate top hits
    top_hits = stat_df[is_sig & is_large].nsmallest(5, pval_col)
    for _, row in top_hits.iterrows():
        ax.annotate(
            row["cell_type"],
            xy=(row[effect_size_col], -np.log10(row[pval_col])),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=0.5),
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

    Returns:
    -------
    plt.Figure
        Matplotlib figure object
    """
    if celltype not in prop_df.columns:
        raise ValueError(f"Cell type {celltype} not in proportion matrix")

    # Prepare data
    plot_df = pd.DataFrame({"timepoint": timepoints, "proportion": prop_df[celltype]})

    if group_col is not None:
        plot_df["group"] = group_col
        palette = _ensure_palette(palette, plot_df["group"].unique())

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))

    if group_col is not None:
        for group in plot_df["group"].unique():
            group_df = plot_df[plot_df["group"] == group]
            ax.plot(
                group_df["timepoint"],
                group_df["proportion"],
                "o-",
                label=group,
                color=palette.get(group),
                linewidth=2,
                markersize=8,
            )
    else:
        ax.plot(plot_df["timepoint"], plot_df["proportion"], "o-", linewidth=2, markersize=8)

    ax.set_xlabel("Timepoint")
    ax.set_ylabel(f"{celltype} Proportion")
    ax.set_title(f"{celltype} Proportion Over Time")

    if group_col is not None:
        ax.legend()

    return fig


@save_and_close("batch_effect")
def plot_batch_effect(
    prop_df: pd.DataFrame,
    batch: pd.Series,
    method: Literal["pca"] = "pca",
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
        Dimensionality reduction method. Currently only 'pca' is supported.
    palette : Dict, optional
        Color palette for batches
    out_dir : str, optional
        Output directory for saving plot

    Returns:
    -------
    plt.Figure
        Matplotlib figure object
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    # Standardize data
    scaler = StandardScaler()
    prop_scaled = scaler.fit_transform(prop_df)

    if method != "pca":
        raise ValueError("plot_batch_effect currently supports only method='pca'")

    reducer = PCA(n_components=2)
    emb = reducer.fit_transform(prop_scaled)
    var_explained = reducer.explained_variance_ratio_
    xlabel = f"PC1 ({var_explained[0]*100:.1f}%)"
    ylabel = f"PC2 ({var_explained[1]*100:.1f}%)"

    # Ensure palette
    palette = _ensure_palette(palette, batch.unique())

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot samples colored by batch
    for b in batch.unique():
        mask = batch == b
        ax.scatter(emb[mask, 0], emb[mask, 1], label=b, color=palette.get(b), s=100, alpha=0.7)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Batch Effect Visualization")
    ax.legend()

    return fig


@save_and_close("composition_pca")
def plot_composition_pca(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    transform: Literal["proportion", "clr", "alr", "logit"] = "clr",
    palette: Optional[Dict] = None,
    top_loadings: int = 5,
    pseudocount: float = 1e-6,
    reference_celltype: Optional[str] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot sample-level composition PCA with optional cell-type loading arrows."""
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    condition = _align_series_to_index(condition, prop_df.index, "condition")
    condition_col = condition.name or "condition"
    matrix = transform_composition(
        prop_df,
        method=transform,
        pseudocount=pseudocount,
        reference_celltype=reference_celltype,
    )
    scaled = StandardScaler().fit_transform(matrix)
    pca = PCA(n_components=2)
    emb = pca.fit_transform(scaled)
    explained = pca.explained_variance_ratio_ * 100
    palette = _ensure_palette(palette, pd.Index(condition.astype(str).unique()), default_cmap="Set2")

    fig, ax = plt.subplots(figsize=(8, 6))
    for group in condition.astype(str).unique():
        mask = condition.astype(str) == group
        ax.scatter(
            emb[mask, 0],
            emb[mask, 1],
            label=group,
            color=palette.get(group, "#808080"),
            s=85,
            alpha=0.82,
            edgecolor="white",
            linewidth=0.6,
        )

    if top_loadings > 0:
        loadings = pd.DataFrame(
            pca.components_.T,
            index=matrix.columns.astype(str),
            columns=["PC1", "PC2"],
        )
        loadings["magnitude"] = np.sqrt(loadings["PC1"] ** 2 + loadings["PC2"] ** 2)
        top = loadings.sort_values("magnitude", ascending=False).head(top_loadings)
        span_x = max(np.ptp(emb[:, 0]), 1e-6)
        span_y = max(np.ptp(emb[:, 1]), 1e-6)
        arrow_scale = 0.35 * min(span_x, span_y)
        for celltype, row in top.iterrows():
            dx = row["PC1"] * arrow_scale
            dy = row["PC2"] * arrow_scale
            ax.arrow(0, 0, dx, dy, color="#333333", width=0.003, head_width=0.06, alpha=0.75)
            ax.text(dx * 1.12, dy * 1.12, celltype, fontsize=9, ha="center", va="center")

    ax.axhline(0, color="#dddddd", linewidth=0.8)
    ax.axvline(0, color="#dddddd", linewidth=0.8)
    ax.set_xlabel(f"PC1 ({explained[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1]:.1f}%)")
    ax.set_title(f"Sample Composition PCA ({transform.upper()})")
    ax.legend(title=condition_col, loc="best")
    return fig


# Additional proportion diagnostics


@save_and_close("composition")
def plot_composition(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot contribution of each condition to cell type proportions."""
    condition = _align_series_to_index(condition, prop_df.index, "condition")
    group_props = prop_df.groupby(condition.astype(str)).mean()
    group_props = group_props.div(group_props.sum(axis=1), axis=0).fillna(0)
    group_props.index = group_props.index.astype(str)

    fig, ax = plt.subplots(figsize=(max(8, len(group_props) * 1.4), 5))
    group_props.plot(
        kind="bar",
        stacked=True,
        color=_resolve_plot_colors(list(group_props.columns), palette),
        edgecolor="white",
        linewidth=0.5,
        ax=ax,
    )
    ax.set_ylim(0, 1)
    ax.set_xlabel(condition.name or "Condition")
    ax.set_ylabel("Mean Proportion")
    ax.set_title("Mean Cell Type Composition by Condition")
    ax.legend(title="Cell Type", loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.tick_params(axis="x", rotation=45)
    return fig


@save_and_close("diff_stats")
def plot_diff_stats(
    prop_df: pd.DataFrame,
    stat_df: pd.DataFrame,
    condition: pd.Series,
    palette: Optional[Dict] = None,
    value_col: Optional[str] = None,
    pval_col: Optional[str] = None,
    sort_by: Literal["value", "abs", "pvalue"] = "abs",
    top_n: Optional[int] = None,
    horizontal: bool = True,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot bar chart with significance brackets."""
    if stat_df.empty:
        raise ValueError("stat_df is empty")

    celltype_col = _resolve_celltype_col(stat_df)
    pval_col = pval_col or _resolve_pvalue_col(stat_df)
    if value_col is None:
        value_col = next(
            (
                col
                for col in (
                    "mean_diff",
                    "effect_size_cohens_d",
                    "effect_size_cliffs_delta",
                    "log2FoldChange",
                    "log-fold change",
                    "statistic",
                )
                if col in stat_df.columns
            ),
            None,
        )
    if value_col is None:
        raise KeyError("stat_df must contain a plottable effect/statistic column")

    plot_df = stat_df.copy()
    plot_df[celltype_col] = plot_df[celltype_col].astype(str)
    plot_df[value_col] = pd.to_numeric(plot_df[value_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[value_col])
    if pval_col:
        plot_df[pval_col] = pd.to_numeric(plot_df[pval_col], errors="coerce")

    if sort_by == "abs":
        plot_df = plot_df.reindex(plot_df[value_col].abs().sort_values(ascending=False).index)
    elif sort_by == "pvalue" and pval_col:
        plot_df = plot_df.sort_values(pval_col)
    else:
        plot_df = plot_df.sort_values(value_col)
    if top_n is not None:
        plot_df = plot_df.head(top_n)

    colors = ["#1f77b4" if val < 0 else "#d62728" for val in plot_df[value_col]]
    if palette:
        colors = [
            palette.get("negative" if val < 0 else "positive", color)
            for val, color in zip(plot_df[value_col], colors)
        ]

    if horizontal:
        fig, ax = plt.subplots(figsize=(8, max(5, len(plot_df) * 0.35)))
        bars = ax.barh(plot_df[celltype_col], plot_df[value_col], color=colors, alpha=0.8)
        ax.axvline(0, color="black", linewidth=1)
        ax.set_xlabel(value_col)
        ax.set_ylabel("Cell Type")
    else:
        fig, ax = plt.subplots(figsize=(max(8, len(plot_df) * 0.55), 5))
        bars = ax.bar(plot_df[celltype_col], plot_df[value_col], color=colors, alpha=0.8)
        ax.axhline(0, color="black", linewidth=1)
        ax.set_xlabel("Cell Type")
        ax.set_ylabel(value_col)
        ax.tick_params(axis="x", rotation=45)
    ax.set_title("Differential Cell Type Proportion Statistics")

    if pval_col:
        pvals = plot_df[pval_col].to_numpy()
        for bar, pval in zip(bars, pvals):
            if pd.isna(pval):
                continue
            if horizontal:
                x = bar.get_width()
                x_range = ax.get_xlim()[1] - ax.get_xlim()[0]
                offset = 0.02 * x_range if x >= 0 else -0.02 * x_range
                ha = "left" if x >= 0 else "right"
                ax.text(
                    x + offset,
                    bar.get_y() + bar.get_height() / 2,
                    _get_sig_stars(float(pval)),
                    ha=ha,
                    va="center",
                    fontsize=10,
                )
            else:
                y = bar.get_height()
                y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
                offset = 0.03 * y_range if y >= 0 else -0.06 * y_range
                va = "bottom" if y >= 0 else "top"
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    y + offset,
                    _get_sig_stars(float(pval)),
                    ha="center",
                    va=va,
                    fontsize=10,
                )
    return fig


@save_and_close("individual_boxplots")
def plot_individual_boxplots(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    stat_df: pd.DataFrame,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot individual box plots with significance tests."""
    plot_df = _proportion_long_frame(prop_df, condition)
    condition_col = (
        condition.name if isinstance(condition, pd.Series) and condition.name else "condition"
    )
    celltypes = list(prop_df.columns.astype(str))
    conditions = _resolve_order(pd.Index(plot_df[condition_col].unique()))
    palette = _ensure_palette(palette, pd.Index(conditions))
    pval_col = _resolve_pvalue_col(stat_df)
    pvals = _stat_lookup(stat_df, pval_col) if pval_col else {}

    ncols = min(4, max(1, len(celltypes)))
    nrows = int(np.ceil(len(celltypes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False)

    for idx, celltype in enumerate(celltypes):
        ax = axes[idx // ncols, idx % ncols]
        sub = plot_df[plot_df["cell_type"] == celltype]
        sns.boxplot(
            data=sub,
            x=condition_col,
            y="proportion",
            hue=condition_col,
            order=conditions,
            palette=palette,
            legend=False,
            ax=ax,
        )
        sns.stripplot(
            data=sub,
            x=condition_col,
            y="proportion",
            order=conditions,
            color="black",
            alpha=0.45,
            size=3,
            jitter=0.12,
            ax=ax,
        )
        ax.set_title(celltype)
        ax.set_xlabel("")
        ax.set_ylabel("Proportion")
        ax.tick_params(axis="x", rotation=30)

        if celltype in pvals and len(conditions) >= 2:
            y_data = sub["proportion"].to_numpy()
            y = _calculate_bracket_height(ax, y_data)
            ax.plot([0, 0, 1, 1], [y * 0.98, y, y, y * 0.98], color="black", linewidth=1)
            ax.text(0.5, y, _get_sig_stars(float(pvals[celltype])), ha="center", va="bottom")

    for idx in range(len(celltypes), nrows * ncols):
        axes[idx // ncols, idx % ncols].axis("off")

    return fig


@save_and_close("proportion_shifts")
def plot_proportion_shifts(
    prop_df: pd.DataFrame,
    condition_col: str,
    condition1: str,
    condition2: str,
    palette: Optional[Dict] = None,
    sort_by: Literal["shift", "abs"] = "abs",
    top_n: Optional[int] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot proportion shifts between two conditions."""
    if condition_col not in prop_df.columns:
        raise KeyError(
            "prop_df must include condition_col for plot_proportion_shifts. "
            "Use plot_composition or plot_proportion_with_ci when condition is a separate Series."
        )

    if {"cell_type", "proportion"}.issubset(prop_df.columns):
        grouped = (
            prop_df[prop_df[condition_col].isin([condition1, condition2])]
            .groupby([condition_col, "cell_type"])["proportion"]
            .mean()
            .unstack(fill_value=0)
        )
    else:
        celltype_cols = [
            col
            for col in prop_df.columns
            if col != condition_col and pd.api.types.is_numeric_dtype(prop_df[col])
        ]
        grouped = (
            prop_df[prop_df[condition_col].isin([condition1, condition2])]
            .groupby(condition_col)[celltype_cols]
            .mean()
        )

    missing = [condition for condition in (condition1, condition2) if condition not in grouped.index]
    if missing:
        raise ValueError(f"Conditions not present in prop_df: {missing}")

    shift_df = grouped.loc[[condition1, condition2]].T
    shift_df["shift"] = shift_df[condition2] - shift_df[condition1]
    if sort_by == "abs":
        shift_df = shift_df.reindex(shift_df["shift"].abs().sort_values(ascending=True).index)
    else:
        shift_df = shift_df.sort_values("shift")
    if top_n is not None:
        shift_df = shift_df.tail(top_n)

    palette = _ensure_palette(palette, pd.Index(shift_df.index))
    fig, ax = plt.subplots(figsize=(8, max(5, len(shift_df) * 0.35)))
    for y_pos, (celltype, row) in enumerate(shift_df.iterrows()):
        color = palette.get(celltype, "#808080")
        ax.plot([row[condition1], row[condition2]], [y_pos, y_pos], color=color, linewidth=2)
        ax.scatter(row[condition1], y_pos, color="#1f77b4", s=45, zorder=3)
        ax.scatter(row[condition2], y_pos, color="#d62728", s=45, zorder=3)

    ax.set_yticks(np.arange(len(shift_df)))
    ax.set_yticklabels(shift_df.index)
    ax.set_xlabel("Mean Proportion")
    ax.set_title(f"Cell Type Proportion Shifts: {condition2} vs {condition1}")
    ax.legend(
        handles=[
            patches.Patch(color="#1f77b4", label=condition1),
            patches.Patch(color="#d62728", label=condition2),
        ],
        title=condition_col,
        loc="best",
    )
    return fig


@save_and_close("paired_proportion_shifts")
def plot_paired_proportion_shifts(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    pair: pd.Series,
    condition1: str,
    condition2: str,
    celltypes: Optional[List[str]] = None,
    palette: Optional[Dict] = None,
    max_celltypes: int = 12,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot paired sample-level shifts for each cell type between two conditions."""
    condition = _align_series_to_index(condition, prop_df.index, "condition")
    pair = _align_series_to_index(pair, prop_df.index, "pair")
    condition_col = condition.name or "condition"
    pair_col = pair.name or "pair"

    long_df = _proportion_long_frame(prop_df, condition)
    sample_col = prop_df.index.name or "sample"
    long_df[pair_col] = long_df[sample_col].map(pair.astype(str))
    long_df = long_df[long_df[condition_col].isin([condition1, condition2])]

    if celltypes is None:
        pivot = long_df.pivot_table(
            index=pair_col,
            columns=[condition_col, "cell_type"],
            values="proportion",
            aggfunc="mean",
        )
        shifts = {}
        for celltype in prop_df.columns.astype(str):
            if (condition1, celltype) in pivot.columns and (condition2, celltype) in pivot.columns:
                delta = pivot[(condition2, celltype)] - pivot[(condition1, celltype)]
                shifts[celltype] = delta.abs().median()
        celltypes = (
            pd.Series(shifts).sort_values(ascending=False).head(max_celltypes).index.tolist()
            if shifts
            else list(prop_df.columns.astype(str))[:max_celltypes]
        )

    ncols = min(4, max(1, len(celltypes)))
    nrows = int(np.ceil(len(celltypes) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.5 * nrows), squeeze=False)
    palette = _ensure_palette(palette, pd.Index([condition1, condition2]), default_cmap="Set2")

    for idx, celltype in enumerate(celltypes):
        ax = axes[idx // ncols, idx % ncols]
        sub = long_df[long_df["cell_type"] == celltype]
        paired = sub.pivot_table(
            index=pair_col,
            columns=condition_col,
            values="proportion",
            aggfunc="mean",
        ).dropna(subset=[condition1, condition2], how="any")

        for _, row in paired.iterrows():
            ax.plot([0, 1], [row[condition1], row[condition2]], color="#9a9a9a", alpha=0.55)
        ax.scatter(
            np.zeros(len(paired)),
            paired[condition1],
            color=palette.get(condition1, "#1f77b4"),
            s=25,
            zorder=3,
            label=condition1,
        )
        ax.scatter(
            np.ones(len(paired)),
            paired[condition2],
            color=palette.get(condition2, "#d62728"),
            s=25,
            zorder=3,
            label=condition2,
        )
        median_delta = (paired[condition2] - paired[condition1]).median() if not paired.empty else 0
        ax.set_title(f"{celltype}\nmedian shift={median_delta:.3f}", fontsize=10)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([condition1, condition2], rotation=25, ha="right")
        ax.set_ylabel("Proportion")

    for idx in range(len(celltypes), nrows * ncols):
        axes[idx // ncols, idx % ncols].axis("off")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles[:2], labels[:2], loc="upper right")
    return fig


@save_and_close("proportion_with_ci")
def plot_proportion_with_ci(
    prop_df: pd.DataFrame,
    condition: pd.Series,
    palette: Optional[Dict] = None,
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot proportions with confidence intervals."""
    long_df = _proportion_long_frame(prop_df, condition)
    condition_col = (
        condition.name if isinstance(condition, pd.Series) and condition.name else "condition"
    )
    summary = (
        long_df.groupby([condition_col, "cell_type"])["proportion"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    summary["ci95"] = 1.96 * summary["std"].fillna(0) / np.sqrt(summary["count"].clip(lower=1))

    conditions = _resolve_order(pd.Index(summary[condition_col].unique()))
    celltypes = _resolve_order(pd.Index(summary["cell_type"].unique()))
    palette = _ensure_palette(palette, pd.Index(conditions))

    x = np.arange(len(celltypes))
    width = min(0.8 / max(1, len(conditions)), 0.35)
    fig, ax = plt.subplots(figsize=(max(8, len(celltypes) * 0.7), 5))

    for idx, cond in enumerate(conditions):
        sub = summary[summary[condition_col] == cond].set_index("cell_type").reindex(celltypes)
        xpos = x + (idx - (len(conditions) - 1) / 2) * width
        ax.bar(
            xpos,
            sub["mean"].fillna(0),
            yerr=sub["ci95"].fillna(0),
            width=width,
            color=palette.get(cond, "#808080"),
            capsize=3,
            label=cond,
            edgecolor="white",
            linewidth=0.5,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(celltypes, rotation=45, ha="right")
    ax.set_ylabel("Mean Proportion")
    ax.set_title("Cell Type Proportions with 95% CI")
    ax.legend(title=condition_col)
    return fig


@save_and_close("celltype_variability")
def plot_celltype_variability(
    prop_df: pd.DataFrame,
    method: str = "cv",
    out_dir: Optional[str] = None,
) -> plt.Figure:
    """Plot cell type variability across samples."""
    if prop_df.empty:
        raise ValueError("prop_df is empty")

    if method == "cv":
        mean = prop_df.mean(axis=0)
        values = prop_df.std(axis=0) / mean.replace(0, np.nan)
        ylabel = "Coefficient of Variation"
    elif method in {"sd", "std"}:
        values = prop_df.std(axis=0)
        ylabel = "Standard Deviation"
    elif method == "variance":
        values = prop_df.var(axis=0)
        ylabel = "Variance"
    elif method == "iqr":
        values = prop_df.quantile(0.75, axis=0) - prop_df.quantile(0.25, axis=0)
        ylabel = "Interquartile Range"
    else:
        raise ValueError("method must be one of: 'cv', 'sd', 'std', 'variance', 'iqr'")

    values = values.replace([np.inf, -np.inf], np.nan).fillna(0).sort_values(ascending=False)
    fig, ax = plt.subplots(figsize=(max(8, len(values) * 0.55), 5))
    ax.bar(values.index.astype(str), values.to_numpy(), color="#4c72b0", alpha=0.85)
    ax.set_title("Cell Type Proportion Variability")
    ax.set_xlabel("Cell Type")
    ax.set_ylabel(ylabel)
    ax.tick_params(axis="x", rotation=45)
    return fig
