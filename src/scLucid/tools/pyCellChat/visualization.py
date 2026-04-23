"""
Visualization functions for CellChat (R-free)
"""

import logging
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

log = logging.getLogger(__name__)


def plot_circle_network(
    cellchat_obj,
    sources_use: Optional[List[str]] = None,
    targets_use: Optional[List[str]] = None,
    signaling: Optional[List[str]] = None,
    remove_isolate: bool = True,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (8, 8),
    **kwargs,
):
    """
    Plot circle plot of cell-cell communication network

    Parameters
    ----------
    cellchat_obj : CellChat
        CellChat object
    sources_use : Optional[List[str]]
        Source cell groups to show
    targets_use : Optional[List[str]]
        Target cell groups to show
    signaling : Optional[List[str]]
        Specific signaling pathways to show
    remove_isolate : bool
        Remove isolated nodes
    thresh : float
        Threshold for edge display
    """
    # Get network data
    if signaling is not None:
        prob_matrix = np.zeros_like(list(cellchat_obj.netP["prob"].values())[0])
        for pathway in signaling:
            if pathway in cellchat_obj.netP["prob"]:
                prob_matrix += cellchat_obj.netP["prob"][pathway]
    else:
        prob_matrix = cellchat_obj.net["prob"].sum(axis=0)

    # Filter by threshold
    prob_matrix[prob_matrix < thresh] = 0

    # Filter sources and targets
    groups = cellchat_obj.unique_groups
    if sources_use is not None:
        source_idx = [i for i, g in enumerate(groups) if g in sources_use]
        prob_matrix = prob_matrix[source_idx, :]
        groups = [groups[i] for i in source_idx]

    if targets_use is not None:
        target_idx = [i for i, g in enumerate(groups) if g in targets_use]
        prob_matrix = prob_matrix[:, target_idx]

    # Remove isolated nodes
    if remove_isolate:
        active_nodes = (prob_matrix.sum(axis=0) > 0) | (prob_matrix.sum(axis=1) > 0)
        prob_matrix = prob_matrix[active_nodes][:, active_nodes]
        groups = [g for g, a in zip(groups, active_nodes) if a]

    # Create circular layout
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection="polar"))

    n_groups = len(groups)
    theta = np.linspace(0, 2 * np.pi, n_groups, endpoint=False)

    # Plot nodes
    ax.scatter(theta, np.ones(n_groups), s=1000, alpha=0.6)

    # Add labels
    for i, (t, g) in enumerate(zip(theta, groups)):
        ax.text(t, 1.2, g, ha="center", va="center")

    # Plot edges
    for i in range(n_groups):
        for j in range(n_groups):
            if prob_matrix[i, j] > thresh:
                ax.plot(
                    [theta[i], theta[j]],
                    [1, 1],
                    alpha=prob_matrix[i, j] / prob_matrix.max(),
                    linewidth=2,
                )

    ax.set_ylim(0, 1.5)
    ax.axis("off")
    plt.tight_layout()

    return fig, ax


def plot_chord_diagram(
    cellchat_obj,
    signaling: Optional[List[str]] = None,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (10, 10),
):
    """
    Plot chord diagram of communication network
    """
    from matplotlib.patches import Wedge

    # Get network data
    if signaling is not None:
        prob_matrix = np.zeros_like(list(cellchat_obj.netP["prob"].values())[0])
        for pathway in signaling:
            if pathway in cellchat_obj.netP["prob"]:
                prob_matrix += cellchat_obj.netP["prob"][pathway]
    else:
        prob_matrix = cellchat_obj.net["prob"].sum(axis=0)

    prob_matrix[prob_matrix < thresh] = 0
    groups = cellchat_obj.unique_groups
    n_groups = len(groups)

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")

    # Calculate angles for each group
    gap = 0.02
    total_angle = 2 * np.pi * (1 - gap * n_groups)
    angles = []
    current_angle = 0

    for i in range(n_groups):
        group_size = total_angle / n_groups
        angles.append((current_angle, current_angle + group_size))
        current_angle += group_size + 2 * np.pi * gap

    # Draw group arcs
    colors = plt.cm.tab20(np.linspace(0, 1, n_groups))

    for i, ((start, end), color, group) in enumerate(zip(angles, colors, groups)):
        wedge = Wedge(
            (0, 0),
            1,
            np.degrees(start),
            np.degrees(end),
            width=0.1,
            facecolor=color,
            edgecolor="white",
            linewidth=2,
        )
        ax.add_patch(wedge)

        # Add label
        mid_angle = (start + end) / 2
        x = 1.2 * np.cos(mid_angle)
        y = 1.2 * np.sin(mid_angle)
        ax.text(x, y, group, ha="center", va="center", fontsize=10)

    # Draw ribbons
    for i in range(n_groups):
        for j in range(n_groups):
            if prob_matrix[i, j] > thresh:
                start_angle = (angles[i][0] + angles[i][1]) / 2
                end_angle = (angles[j][0] + angles[j][1]) / 2

                t = np.linspace(0, 1, 100)
                x = (1 - t) * 0.9 * np.cos(start_angle) + t * 0.9 * np.cos(end_angle)
                y = (1 - t) * 0.9 * np.sin(start_angle) + t * 0.9 * np.sin(end_angle)

                alpha = prob_matrix[i, j] / prob_matrix.max()
                ax.plot(x, y, color=colors[i], alpha=alpha * 0.5, linewidth=2)

    plt.tight_layout()
    return fig, ax


def plot_heatmap(
    cellchat_obj,
    signaling: Optional[List[str]] = None,
    sources_use: Optional[List[str]] = None,
    targets_use: Optional[List[str]] = None,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (8, 6),
    cmap: str = "Reds",
    **kwargs,
):
    """
    Plot heatmap of communication probability
    """
    # Get network data
    if signaling is not None:
        prob_matrix = np.zeros_like(list(cellchat_obj.netP["prob"].values())[0])
        for pathway in signaling:
            if pathway in cellchat_obj.netP["prob"]:
                prob_matrix += cellchat_obj.netP["prob"][pathway]
    else:
        prob_matrix = cellchat_obj.net["prob"].sum(axis=0)

    # Create DataFrame
    groups = cellchat_obj.unique_groups
    df = pd.DataFrame(prob_matrix, index=groups, columns=groups)

    # Filter
    if sources_use is not None:
        df = df.loc[sources_use, :]
    if targets_use is not None:
        df = df.loc[:, targets_use]

    # Apply threshold
    df[df < thresh] = 0

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        df,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        square=True,
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "Communication Probability"},
    )

    ax.set_xlabel("Target")
    ax.set_ylabel("Source")
    ax.set_title("Cell-Cell Communication Heatmap")

    plt.tight_layout()
    return fig, ax


def plot_bubble(
    cellchat_obj,
    sources_use: Optional[List[str]] = None,
    targets_use: Optional[List[str]] = None,
    signaling: Optional[List[str]] = None,
    thresh: float = 0.05,
    figsize: Tuple[int, int] = (12, 8),
    **kwargs,
):
    """
    Plot bubble plot of L-R pairs
    """
    # Prepare data
    if signaling is not None:
        lr_pairs = cellchat_obj.LR[cellchat_obj.LR["pathway_name"].isin(signaling)]
    else:
        lr_pairs = cellchat_obj.LR

    # Get probabilities
    prob_data = []
    pval_data = []

    for idx in lr_pairs.index:
        lr_name = lr_pairs.loc[idx, "interaction_name"]
        prob = cellchat_obj.net["prob"][idx]
        pval = cellchat_obj.net["pval"][idx]

        for i, source in enumerate(cellchat_obj.unique_groups):
            for j, target in enumerate(cellchat_obj.unique_groups):
                if sources_use is None or source in sources_use:
                    if targets_use is None or target in targets_use:
                        if prob[i, j] > thresh:
                            prob_data.append(
                                {
                                    "source": source,
                                    "target": target,
                                    "lr_pair": lr_name,
                                    "prob": prob[i, j],
                                    "pval": pval[i, j],
                                }
                            )

    df = pd.DataFrame(prob_data)

    if df.empty:
        log.warning("No significant interactions to plot")
        return None, None

    # Create bubble plot
    fig, ax = plt.subplots(figsize=figsize)

    df["interaction"] = df["source"] + " -> " + df["target"]

    scatter = ax.scatter(
        range(len(df)),
        df["lr_pair"].astype("category").cat.codes,
        s=df["prob"] * 1000,
        c=-np.log10(df["pval"] + 1e-10),
        cmap="Reds",
        alpha=0.6,
        edgecolors="black",
        linewidth=0.5,
    )

    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(df["interaction"], rotation=90)
    ax.set_yticks(range(len(df["lr_pair"].unique())))
    ax.set_yticklabels(df["lr_pair"].unique())

    plt.colorbar(scatter, label="-log10(p-value)", ax=ax)
    ax.set_xlabel("Cell-Cell Interactions")
    ax.set_ylabel("L-R Pairs")
    ax.set_title("Communication Probability (bubble size)")

    plt.tight_layout()
    return fig, ax


def plot_contribution(
    cellchat_obj, signaling: str, thresh: float = 0.05, figsize: Tuple[int, int] = (10, 6), **kwargs
):
    """
    Plot contribution of each L-R pair to a pathway
    """
    # Get L-R pairs in this pathway
    lr_pairs = cellchat_obj.LR[cellchat_obj.LR["pathway_name"] == signaling]

    # Calculate contribution
    contributions = []

    for idx in lr_pairs.index:
        lr_name = lr_pairs.loc[idx, "interaction_name"]
        prob = cellchat_obj.net["prob"][idx].sum()
        contributions.append({"lr_pair": lr_name, "contribution": prob})

    df = pd.DataFrame(contributions).sort_values("contribution", ascending=False)

    if df.empty:
        log.warning(f"No L-R pairs found for pathway {signaling}")
        return None, None

    # Plot
    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(range(len(df)), df["contribution"])
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["lr_pair"])
    ax.set_xlabel("Contribution")
    ax.set_title(f"L-R Pair Contribution to {signaling}")

    plt.tight_layout()
    return fig, ax


def plot_signaling_gene_expression(
    cellchat_obj, signaling: str, figsize: Tuple[int, int] = (12, 8), **kwargs
):
    """
    Plot expression of signaling genes
    """
    # Get genes in this pathway
    lr_pairs = cellchat_obj.LR[cellchat_obj.LR["pathway_name"] == signaling]

    genes = set()
    for _, row in lr_pairs.iterrows():
        genes.add(row["ligand"])
        genes.add(row["receptor"])

    # Get expression
    gene_expr = []
    for gene in genes:
        expr = cellchat_obj._get_gene_expression(gene)
        if expr is not None:
            gene_expr.append(pd.Series(expr, name=gene, index=cellchat_obj.unique_groups))

    if not gene_expr:
        log.warning("No gene expression data available")
        return None, None

    df = pd.DataFrame(gene_expr).T

    # Plot heatmap
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        df.T,
        cmap="RdYlBu_r",
        center=0,
        annot=True,
        fmt=".2f",
        cbar_kws={"label": "Expression Level"},
        ax=ax,
    )

    ax.set_xlabel("Genes")
    ax.set_ylabel("Cell Groups")
    ax.set_title(f"Signaling Gene Expression: {signaling}")

    plt.tight_layout()
    return fig, ax
