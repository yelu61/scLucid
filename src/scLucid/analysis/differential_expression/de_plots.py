"""
Visualization functions for differential expression results.

This module provides publication-quality plots:
- visualize_markers: Multi-panel marker visualization
- plot_volcano: Volcano plots for DE results
- plot_multi_cluster_deg: Heatmap of DE genes across clusters
""",

import logging
from typing import Dict, List, Literal, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from adjustText import adjust_text
from anndata import AnnData
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from ...base_config import SclucidBaseConfig

log = logging.getLogger(__name__)


def visualize_markers(
    adata: AnnData,
    markers: Union[pd.DataFrame, Dict[str, List[str]], List[str]],
    groupby: Optional[str] = None,
    n_genes_per_group: int = 5,
    plot_type: Literal[
        "dotplot", "heatmap", "stacked_violin", "violin", "matrixplot"
    ] = "dotplot",
    dendrogram: bool = False,
    standard_scale: Optional[Literal["var", "group"]] = "var",
    swap_axes: bool = False,
    layer: Optional[str] = None,
    use_raw: bool = False,
    save_path: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwargs,
) -> None:
    """
    Visualize marker genes with automatic formatting and error handling.

    Supports multiple input formats:
    - DataFrame with 'group' and 'names' columns
    - Dictionary mapping groups to gene lists
    - Simple list of genes

    Args:
        adata: AnnData object
        markers: Marker genes (DataFrame/dict/list)
        groupby: Grouping column (required for list input)
        n_genes_per_group: Top N genes per group (for DataFrame input)
        plot_type: Visualization type
        dendrogram: Add dendrogram
        standard_scale: Standardization ('var' or 'group')
        swap_axes: Swap x/y axes
        layer: Data layer to use
        use_raw: Use .raw
        save_path: Output file path
        figsize: Figure size (auto-calculated if None)
        **kwargs: Additional arguments to Scanpy plotting functions

    Example:
        >>> # From DataFrame
        >>> visualize_markers(
        ...     adata,
        ...     markers=filtered_markers,
        ...     groupby="leiden",
        ...     plot_type="dotplot"
        ... )
        >>>
        >>> # From dictionary
        >>> marker_dict = {
        ...     "T_cells": ["CD3D", "CD3E", "CD3G"],
        ...     "B_cells": ["CD19", "MS4A1", "CD79A"]
        ... }
        >>> visualize_markers(adata, markers=marker_dict)
    """
    gene_list: List[str] = []
    gene_dict: Dict[str, List[str]] = {}

    # Parse input markers
    if isinstance(markers, pd.DataFrame):
        # Standardize column names
        if "names" not in markers.columns:
            for alt in ("gene", "Gene", "feature", "symbol"):
                if alt in markers.columns:
                    markers = markers.rename(columns={alt: "names"})
                    break

        if "group" not in markers.columns:
            raise ValueError(
                "DataFrame must contain 'group' and 'names' columns "
                "for grouped visualization"
            )

        # Extract top genes per group
        for g in markers["group"].unique():
            group_markers = markers[markers["group"] == g]

            # Sort by logfoldchanges or scores
            if "logfoldchanges" in group_markers.columns:
                group_markers = group_markers.sort_values(
                    "logfoldchanges", ascending=False
                )
            elif "scores" in group_markers.columns:
                group_markers = group_markers.sort_values("scores", ascending=False)

            top_genes = group_markers["names"].head(n_genes_per_group).tolist()
            gene_list.extend(top_genes)
            gene_dict[str(g)] = top_genes

    elif isinstance(markers, dict):
        gene_dict = markers
        for genes in markers.values():
            gene_list.extend(list(genes))

    elif isinstance(markers, (list, tuple)):
        gene_list = list(markers)
        if groupby is None:
            raise ValueError("groupby must be specified when markers is a list")
        gene_dict = {"Selected Markers": gene_list}

    else:
        raise TypeError("markers must be a DataFrame, dictionary, or list")

    # Deduplicate and validate
    gene_list_unique = [g for g in dict.fromkeys(gene_list) if g in adata.var_names]
    if not gene_list_unique:
        raise ValueError("No valid genes found in adata.var_names")

    # Prepare gene_dict for grouped plots
    for g, glist in gene_dict.items():
        gene_dict[g] = [gene for gene in glist if gene in adata.var_names]

    # Auto-calculate figsize
    if figsize is None:
        if groupby and groupby in adata.obs:
            n_groups = (
                len(adata.obs[groupby].cat.categories)
                if pd.api.types.is_categorical_dtype(adata.obs[groupby])
                else len(adata.obs[groupby].unique())
            )
        else:
            n_groups = 1

        n_genes = len(gene_list_unique)

        if plot_type in ["heatmap", "dotplot", "matrixplot"]:
            width = max(6, min(16, n_groups * 0.5))
            height = max(4, min(25, n_genes * 0.3))
            if swap_axes:
                width, height = height, width
        elif plot_type == "stacked_violin":
            width = max(6, min(16, n_groups * 0.5))
            height = max(4, min(25, n_genes * 0.4))
        else:  # violin
            width = max(8, n_genes * 2)
            height = 6

        figsize = (width, height)

    # Plot
    plot_kwargs = {
        "groupby": groupby,
        "dendrogram": dendrogram,
        "standard_scale": standard_scale,
        "use_raw": use_raw,
        "layer": layer,
        "figsize": figsize,
        "show": False,
        **kwargs,
    }

    try:
        if plot_type == "dotplot":
            sc.pl.dotplot(
                adata, var_names=gene_dict, swap_axes=swap_axes, **plot_kwargs
            )
        elif plot_type == "heatmap":
            sc.pl.heatmap(
                adata, var_names=gene_dict, swap_axes=swap_axes, **plot_kwargs
            )
        elif plot_type == "stacked_violin":
            sc.pl.stacked_violin(
                adata, var_names=gene_list_unique, swap_axes=swap_axes, **plot_kwargs
            )
        elif plot_type == "matrixplot":
            sc.pl.matrixplot(
                adata, var_names=gene_dict, swap_axes=swap_axes, **plot_kwargs
            )
        elif plot_type == "violin":
            sc.pl.violin(adata, keys=gene_list_unique, **plot_kwargs)
        else:
            raise ValueError(f"Unknown plot type: {plot_type}")

    except Exception as e:
        log.error(f"Plotting failed: {e}")
        raise

    # Save or show
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"Saved visualization to {save_path}")
        plt.close()
    else:
        plt.show()


def plot_volcano(
    degs_df: pd.DataFrame,
    title: str,
    subtitle: Optional[str] = None,
    top_n_up: int = 15,
    top_n_down: int = 15,
    genes_to_highlight: Optional[List[str]] = None,
    lfc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
    palette: Optional[Dict[str, str]] = None,
    savepath: Optional[str] = None,
    figsize: tuple = (12, 12),
    dpi: int = 300,
) -> None:
    """
    Publication-quality volcano plot with intelligent label placement.

    Features:
    - Smart label selection based on ranking score (|LFC| * -log10(p))
    - adjustText for anti-overlap
    - Statistical summary box
    - Custom gene highlighting

    Args:
        degs_df: DE results DataFrame
        title: Main title
        subtitle: Subtitle (e.g., sample info)
        top_n_up: Number of top up-regulated genes to label
        top_n_down: Number of top down-regulated genes to label
        genes_to_highlight: Additional genes to highlight
        lfc_threshold: Log2 fold change threshold
        pval_threshold: Adjusted p-value threshold
        palette: Color scheme
        savepath: Output file path
        figsize: Figure size
        dpi: Resolution

    Example:
        >>> plot_volcano(
        ...     degs_df,
        ...     title="T cells: Treated vs Control",
        ...     subtitle="n=1234 cells",
        ...     top_n_up=20,
        ...     genes_to_highlight=["CD3D", "CD4"]
        ... )
    """
    df = degs_df.copy()

    # Calculate -log10(p-adj)
    df["-log10_pvals_adj"] = -np.log10(df["pvals_adj"].astype(float) + 1e-300)

    # Categorize genes
    df["status"] = "Not significant"
    df.loc[
        (df["logfoldchanges"] > lfc_threshold) & (df["pvals_adj"] < pval_threshold),
        "status",
    ] = "Up-regulated"
    df.loc[
        (df["logfoldchanges"] < -lfc_threshold) & (df["pvals_adj"] < pval_threshold),
        "status",
    ] = "Down-regulated"

    # Default palette
    if palette is None:
        palette = {
            "Up-regulated": "#d62728",
            "Down-regulated": "#1f77b4",
            "Not significant": "#cccccc",
        }

    # Create figure
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=figsize)

    # Plot points (layered for z-order)
    for status, color, alpha, size in [
        ("Not significant", palette["Not significant"], 0.4, 15),
        ("Up-regulated", palette["Up-regulated"], 0.8, 30),
        ("Down-regulated", palette["Down-regulated"], 0.8, 30),
    ]:
        mask = df["status"] == status
        ax.scatter(
            df.loc[mask, "logfoldchanges"],
            df.loc[mask, "-log10_pvals_adj"],
            s=size,
            alpha=alpha,
            c=color,
            label=status,
            ec="none",
            zorder=2 if status != "Not significant" else 1,
        )

    # Smart label selection
    df["ranking_score"] = np.abs(df["logfoldchanges"]) * df["-log10_pvals_adj"]

    up_genes = df[df["status"] == "Up-regulated"].nlargest(top_n_up, "ranking_score")
    down_genes = df[df["status"] == "Down-regulated"].nlargest(
        top_n_down, "ranking_score"
    )

    genes_to_label_df = pd.concat([up_genes, down_genes])

    # Add custom highlights
    if genes_to_highlight:
        highlight_df = df[df["names"].isin(genes_to_highlight)]
        genes_to_label_df = pd.concat(
            [genes_to_label_df, highlight_df]
        ).drop_duplicates(subset=["names"])

    # Add labels
    texts = []
    for _, row in genes_to_label_df.iterrows():
        txt = ax.text(
            row["logfoldchanges"],
            row["-log10_pvals_adj"],
            row["names"],
            fontsize=10,
            zorder=3,
        )
        texts.append(txt)

    # adjustText for anti-overlap
    if texts:
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.5, alpha=0.7),
            expand_points=(2.0, 2.0),
            expand_text=(1.3, 1.3),
            force_points=(0.3, 0.6),
            force_text=(0.5, 1.0),
            lim=1000,
            precision=0.01,
        )

    # Threshold lines
    ax.axhline(
        y=-np.log10(pval_threshold),
        color="grey",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
    )
    ax.axvline(x=lfc_threshold, color="grey", linestyle="--", linewidth=1, alpha=0.7)
    ax.axvline(x=-lfc_threshold, color="grey", linestyle="--", linewidth=1, alpha=0.7)

    # Statistical summary
    num_up = (df["status"] == "Up-regulated").sum()
    num_down = (df["status"] == "Down-regulated").sum()

    ax.text(
        0.02,
        0.98,
        f"Up: {num_up}\nDown: {num_down}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.8, ec="none"),
    )

    # Titles and labels
    fig.suptitle(title, fontsize=20, weight="bold", y=0.98)
    if subtitle:
        ax.set_title(subtitle, fontsize=14, pad=10)

    ax.set_xlabel("Log2 Fold Change", fontsize=14, weight="bold")
    ax.set_ylabel("-log10(Adjusted P-value)", fontsize=14, weight="bold")

    # Legend
    ax.legend(loc="upper right", frameon=True, fontsize=12, shadow=True)

    # Clean up
    sns.despine(ax=ax)
    plt.tight_layout()

    # Save or show
    if savepath:
        plt.savefig(savepath, dpi=dpi, bbox_inches="tight")
        log.info(f"Volcano plot saved: {savepath}")
        plt.close()
    else:
        plt.show()


def plot_multi_cluster_deg(
    df: pd.DataFrame,
    highlight_genes: Optional[List[str]] = None,
    pval_cutoff: float = 0.01,
    logfc_threshold: float = 1.0,
    top_n: int = 5,
    point_size_by_pval: bool = False,
    add_colored_bottom: bool = True,
    cluster_color_dict: Optional[Dict] = None,
    out_path: Optional[str] = None,
    figsize: Optional[tuple] = None,
    dpi: int = 300,
) -> None:
    """
    Multi-cluster differential expression overview plot.

    Shows:
    - All genes across all clusters in a single view
    - Significance-based coloring
    - Smart label placement for top genes
    - Optional colored bottom strip for cluster identification

    Args:
        df: DE DataFrame with columns: group, names, logfoldchanges, pvals_adj
        highlight_genes: Genes to highlight in green
        pval_cutoff: P-value cutoff for significance
        logfc_threshold: Log fold change threshold
        top_n: Top N genes to label per cluster (up and down)
        point_size_by_pval: Scale point size by -log10(p)
        add_colored_bottom: Add colored cluster strip at bottom
        cluster_color_dict: Custom cluster colors
        out_path: Output file path
        figsize: Figure size (auto-calculated if None)
        dpi: Resolution

    Example:
        >>> plot_multi_cluster_deg(
        ...     markers_df,
        ...     top_n=10,
        ...     highlight_genes=["CD3D", "CD19"],
        ...     out_path="cluster_overview.pdf"
        ... )
    """
    # Standardize column names
    if "names" in df.columns and "Gene" not in df.columns:
        df = df.rename(columns={"names": "Gene"})
    if "logfoldchanges" in df.columns and "avg_logFC" not in df.columns:
        df = df.rename(columns={"logfoldchanges": "avg_logFC"})
    if "group" in df.columns and "Cluster" not in df.columns:
        df = df.rename(columns={"group": "Cluster"})

    # Sort clusters
    try:
        clusters = sorted(df["Cluster"].unique(), key=int)
    except (ValueError, TypeError):
        clusters = sorted(df["Cluster"].unique())

    x_pos = np.arange(len(clusters))
    cluster_map = dict(zip(clusters, x_pos))

    # Colors
    if cluster_color_dict:
        color_map = cluster_color_dict
    else:
        cluster_colors = plt.cm.Spectral(np.linspace(0, 1, len(clusters)))
        color_map = dict(zip(clusters, cluster_colors))

    # Auto figsize
    if figsize is None:
        fig_width = max(16, len(clusters) * 1.8)
        fig_height = max(8, 8 + top_n * 0.2)
        figsize = (fig_width, fig_height)

    fig, ax = plt.subplots(figsize=figsize)

    texts = []
    points_coords = []

    for c in clusters:
        sub = df[df["Cluster"] == c].copy()
        idx = cluster_map[c]

        y = sub["avg_logFC"].values
        sig = sub["pvals_adj"].values < pval_cutoff

        up = (y > logfc_threshold) & sig
        down = (y < -logfc_threshold) & sig
        ns = ~sig

        # -log10(p) for sizing
        sub["neg_log_p"] = -np.log10(np.clip(sub["pvals_adj"], 1e-10, 1))

        # X jitter
        x = np.full(len(sub), idx)
        x_jitter = x + np.random.uniform(-0.45, 0.45, len(sub))

        # Point sizes
        base_size = 5
        if point_size_by_pval:
            sizes_ns = base_size * np.ones(sum(ns))
            sizes_up = base_size + 5 * sub.loc[up, "neg_log_p"]
            sizes_down = base_size + 5 * sub.loc[down, "neg_log_p"]
        else:
            sizes_ns = base_size
            sizes_up = base_size * 1.6
            sizes_down = base_size * 1.6

        # Plot points
        ax.scatter(x_jitter[ns], y[ns], c="#cccccc", s=sizes_ns, alpha=0.4, zorder=1)
        ax.scatter(x_jitter[up], y[up], c="#d62728", s=sizes_up, alpha=0.8, zorder=2)
        ax.scatter(
            x_jitter[down], y[down], c="#1f77b4", s=sizes_down, alpha=0.8, zorder=2
        )

        # Smart labeling
        sub["ranking_score"] = np.abs(sub["avg_logFC"]) * sub["neg_log_p"]

        # Top up
        top_up = (
            sub[up]
            .nlargest(top_n, "ranking_score")
            .sort_values("avg_logFC", ascending=False)
        )

        for j, (_, row) in enumerate(top_up.iterrows()):
            x_offset = [-0.25, 0, 0.25][j % 3]
            y_offset = 0.15

            txt = ax.text(
                idx + x_offset,
                row["avg_logFC"] + y_offset,
                row["Gene"],
                fontsize=6.5,
                ha="center",
                va="bottom",
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.7,
                ),
                zorder=3,
            )
            texts.append(txt)
            points_coords.append((idx, row["avg_logFC"]))

        # Top down
        top_down = (
            sub[down]
            .nlargest(top_n, "ranking_score")
            .sort_values("avg_logFC", ascending=True)
        )

        for j, (_, row) in enumerate(top_down.iterrows()):
            x_offset = [-0.25, 0, 0.25][j % 3]
            y_offset = -0.15

            txt = ax.text(
                idx + x_offset,
                row["avg_logFC"] + y_offset,
                row["Gene"],
                fontsize=6.5,
                ha="center",
                va="top",
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.7,
                ),
                zorder=3,
            )
            texts.append(txt)
            points_coords.append((idx, row["avg_logFC"]))

        # Highlight genes
        if highlight_genes:
            high_sub = sub[sub["Gene"].isin(highlight_genes)]
            for _, row in high_sub.iterrows():
                va = "bottom" if row["avg_logFC"] > 0 else "top"
                y_off = 0.2 if row["avg_logFC"] > 0 else -0.2

                txt = ax.text(
                    idx,
                    row["avg_logFC"] + y_off,
                    row["Gene"],
                    fontsize=7.5,
                    fontweight="bold",
                    color="green",
                    ha="center",
                    va=va,
                    bbox=dict(
                        boxstyle="round,pad=0.25",
                        facecolor="yellow",
                        edgecolor="green",
                        alpha=0.6,
                        linewidth=1.5,
                    ),
                    zorder=4,
                )
                texts.append(txt)
                points_coords.append((idx, row["avg_logFC"]))

    # adjustText
    if texts:
        adjust_text(
            texts,
            x=[p[0] for p in points_coords],
            y=[p[1] for p in points_coords],
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.4),
            expand_points=(2.5, 2.5),
            expand_text=(1.5, 1.8),
            expand_objects=(1.5, 1.5),
            force_points=(0.4, 0.8),
            force_text=(0.6, 1.2),
            force_objects=(0.4, 0.6),
            lim=2000,
            precision=0.001,
            only_move={"points": "xy", "text": "xy"},
            avoid_self=True,
            avoid_points=True,
            avoid_text=True,
            autoalign="xy",
        )

    # Threshold lines
    ax.axhline(logfc_threshold, ls="--", c="black", alpha=0.5, linewidth=1)
    ax.axhline(-logfc_threshold, ls="--", c="black", alpha=0.5, linewidth=1)
    ax.axhline(0, ls="--", c="gray", linewidth=0.8)

    # Colored bottom strip
    if add_colored_bottom:
        ylim = ax.get_ylim()
        dy = (ylim[1] - ylim[0]) * 0.035
        y_margin = (ylim[1] - ylim[0]) * 0.15
        ax.set_ylim(ylim[0] - dy, ylim[1] + y_margin)

        for i, c in enumerate(clusters):
            color = color_map.get(c, "gray")
            ax.add_patch(
                Rectangle(
                    (i - 0.5, ylim[0] - dy),
                    1,
                    dy,
                    color=color,
                    edgecolor="white",
                    linewidth=0.5,
                    clip_on=False,
                    zorder=0,
                )
            )

            # Auto text color
            if isinstance(color, str) and color.startswith("#"):
                rgb = [int(color.lstrip("#")[k : k + 2], 16) / 255 for k in (0, 2, 4)]
                text_color = "white" if np.mean(rgb) < 0.5 else "black"
            else:
                text_color = "black"

            ax.text(
                i,
                ylim[0] - dy / 2,
                str(c),
                ha="center",
                va="center",
                fontsize=9,
                color=text_color,
                weight="bold",
            )

        ax.set_xticks([])
        ax.set_xlabel("")
    else:
        ax.set_xticks(x_pos)
        ax.set_xticklabels(clusters, rotation=45, ha="right")
        ax.set_xlabel("Cluster")

    # Labels and title
    ax.set_ylabel("Average Log2 Fold Change", fontsize=12, weight="bold")
    ax.set_title(
        "Differential Expression per Cluster", fontsize=14, weight="bold", pad=20
    )

    # Grid
    ax.grid(True, ls="--", alpha=0.2, linewidth=0.5)

    # Spines
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    # Legend
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=f"Sig Up (P < {pval_cutoff})",
            markerfacecolor="#d62728",
            markersize=8,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=f"Sig Down (P < {pval_cutoff})",
            markerfacecolor="#1f77b4",
            markersize=8,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Non-Sig",
            markerfacecolor="#cccccc",
            markersize=8,
        ),
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper right",
        frameon=True,
        fancybox=True,
        shadow=True,
        fontsize=10,
    )

    plt.tight_layout()

    if out_path:
        plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
        log.info(f"Multi-cluster DEG plot saved: {out_path}")
        plt.close()
    else:
        plt.show()


# ==================== Result Management ====================


