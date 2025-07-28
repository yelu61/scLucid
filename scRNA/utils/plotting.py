"""
Dimensionality reduction and visualization functions for single-cell RNA-seq data.

This module provides functions for visualizing marker gene expression and
cell type compositions across clusters.
"""

import os
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from adjustText import adjust_text


def plot_embedding(
    adata: sc.AnnData,
    color_by: str,
    basis: str = "umap",
    title: Optional[str] = None,
    show_labels: bool = True,
    **kwargs,
) -> None:
    """
    Custom wrapper for `sc.pl.embedding` with improved labeling.

    Args:
        adata: AnnData object.
        color_by: Key in `adata.obs` to color by.
        basis: Basis for the embedding (e.g., 'umap').
        title: Plot title.
        show_labels: Whether to show cluster/group labels on the plot.
        **kwargs: Additional arguments passed to `sc.pl.embedding`.
    """
    fig, ax = plt.subplots(figsize=kwargs.pop("figsize", (8, 7)))

    sc.pl.embedding(
        adata,
        basis=basis,
        color=color_by,
        ax=ax,
        show=False,
        legend_loc="on data" if show_labels else "right margin",
        **kwargs,
    )

    if show_labels:
        ax.legend_.remove()  # Remove default legend to replace with adjust_text

        texts = []
        for label in adata.obs[color_by].cat.categories:
            x, y = np.median(
                adata[adata.obs[color_by] == label].obsm[f"X_{basis}"], axis=0
            )
            texts.append(ax.text(x, y, label, fontsize=10, fontweight="bold"))

        adjust_text(
            texts, ax=ax, arrowprops=dict(arrowstyle="-", color="black", lw=0.5)
        )

    ax.set_title(title if title else f"{basis.upper()} colored by {color_by}")
    plt.tight_layout()
    plt.show()


def plot_marker_heatmap(
    adata: sc.AnnData,
    markers_df: pd.DataFrame,
    groupby: str,
    n_genes: int = 5,
    **kwargs,
):
    """
    Plots a heatmap of top marker genes from a markers DataFrame.

    Args:
        adata: AnnData object.
        markers_df: DataFrame from `find_markers` or `filter_markers`.
        groupby: Key in `adata.obs` used for grouping.
        n_genes: Number of top genes to show per group.
        **kwargs: Additional arguments passed to `sc.pl.heatmap`.
    """
    top_markers = markers_df.groupby("group").head(n_genes)
    marker_dict = top_markers.groupby("group")["names"].apply(list).to_dict()

    sc.pl.heatmap(
        adata, marker_dict, groupby=groupby, dendrogram=True, swap_axes=True, **kwargs
    )


def plot_composition(
    adata: sc.AnnData, group_key: str, stack_key: str, normalize: bool = True, **kwargs
):
    """
    Generates a stacked bar plot showing the composition of groups.

    Args:
        adata: AnnData object.
        group_key: Key in `adata.obs` for the x-axis groups (e.g., 'leiden').
        stack_key: Key in `adata.obs` for the stacked bars (e.g., 'cell_type_annotated').
        normalize: If True, show proportions; otherwise, show raw counts.
    """
    composition = pd.crosstab(
        adata.obs[group_key],
        adata.obs[stack_key],
        normalize="index" if normalize else "all",
    )

    composition.plot.bar(stacked=True, figsize=kwargs.pop("figsize", (12, 7)), **kwargs)
    plt.ylabel("Proportion" if normalize else "Count")
    plt.legend(title=stack_key, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.show()


def plot_enrichment(
    adata: sc.AnnData,
    cluster: str,
    enrich_key: str = "enrichment",
    n_top_terms: int = 10,
    save_dir: Optional[str] = None,
    show: bool = True,
):
    """
    Visualize functional enrichment analysis results for a single cluster.

    Args:
        adata: AnnData object.
        cluster: Name of the cluster to visualize.
        enrich_key: Key in adata.uns storing enrichment analysis results.
        n_top_terms: Number of top terms to display.
        save_dir: Directory to save plot files. If None, don't save.
        show: Whether to display the image.
    """
    if enrich_key not in adata.uns or cluster not in adata.uns[enrich_key]:
        raise KeyError(
            f"Enrichment results for cluster '{cluster}' not found in `adata.uns['{enrich_key}']`. "
            "Please run `scRNA.analysis.run_enrichment()` first."
        )

    results_df = adata.uns[enrich_key][cluster]
    if results_df.empty:
        print(f"No enrichment results to plot for cluster '{cluster}'.")
        return

    top_terms = results_df.head(n_top_terms)

    # Create plot
    plt.figure(figsize=(8, max(5, n_top_terms * 0.5)))  # Dynamically adjust height
    plt.barh(
        top_terms["Term"], -np.log10(top_terms["Adjusted P-value"]), color="steelblue"
    )
    plt.title(f"Top GO Terms for {cluster}", fontsize=14)
    plt.xlabel("-log10(Adjusted P-value)", fontsize=12)
    plt.ylabel("GO Biological Process", fontsize=12)
    plt.gca().invert_yaxis()
    plt.grid(axis="x", linestyle="--", alpha=0.6)

    plt.tight_layout()

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(
            os.path.join(save_dir, f"{cluster}_GO_enrichment_plot.png"),
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )

    if show:
        plt.show()

    plt.close()

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def plot_slingshot_pseudotime(
    adata,
    embedding_key="X_umap",
    pseudotime_key="slingshot_pseudotime",
    lineage_key="slingshot_lineage",
    cluster_key=None,
    show_lineage=True,
    figsize=(7, 6),
    cmap="viridis",
    save=None
):
    """
    可视化 Slingshot 伪时序分析结果（Python版）

    参数:
        adata: AnnData对象
        embedding_key: obsm降维坐标key，通常是"X_umap"或"X_pca"
        pseudotime_key: obs中伪时序字段
        lineage_key: obs中分支/主干字段
        cluster_key: obs中聚类字段（可选，作为参考）
        show_lineage: 是否分面显示不同分支
        figsize: 图像大小
        cmap: 连续型配色
        save: 文件名（如不为None则保存图片）
    """
    X = adata.obsm[embedding_key]
    pt = adata.obs[pseudotime_key]
    lineage = adata.obs[lineage_key] if lineage_key in adata.obs else None

    plt.figure(figsize=figsize)
    if lineage is not None and show_lineage:
        uniq_lineage = np.unique(lineage.dropna())
        n_l = len(uniq_lineage)
        fig, axs = plt.subplots(1, n_l, figsize=(figsize[0]*n_l, figsize[1]))
        axs = np.array(axs).reshape(-1)
        for idx, l in enumerate(uniq_lineage):
            sel = (lineage == l)
            sc = axs[idx].scatter(
                X[sel, 0], X[sel, 1],
                c=pt[sel],
                cmap=cmap,
                s=8, alpha=0.8
            )
            axs[idx].set_title(f"Lineage {l}")
            axs[idx].set_xlabel("UMAP1")
            axs[idx].set_ylabel("UMAP2")
            plt.colorbar(sc, ax=axs[idx], label="Pseudotime")
    else:
        sc = plt.scatter(
            X[:, 0], X[:, 1],
            c=pt,
            cmap=cmap,
            s=8, alpha=0.8
        )
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title("Slingshot Pseudotime")
        plt.colorbar(sc, label="Pseudotime")

    if save is not None:
        plt.savefig(save, bbox_inches="tight", dpi=200)
    plt.show()

def plot_slingshot_lineage(
    adata,
    embedding_key="X_umap",
    lineage_key="slingshot_lineage",
    cluster_key=None,
    palette="tab10",
    figsize=(7,6),
    save=None
):
    """
    可视化 Slingshot 主分支分配
    """
    X = adata.obsm[embedding_key]
    lineage = adata.obs[lineage_key].astype(str)
    uniq_l = np.unique(lineage.dropna())
    color_dict = dict(zip(uniq_l, sns.color_palette(palette, len(uniq_l))))
    colors = lineage.map(color_dict)

    plt.figure(figsize=figsize)
    plt.scatter(X[:, 0], X[:, 1], c=colors, s=8, alpha=0.8)
    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")
    plt.title("Slingshot Lineage Assignment")
    # legend
    for l in uniq_l:
        plt.scatter([], [], c=[color_dict[l]], label=f"Lineage {l}")
    plt.legend(markerscale=2)
    if save is not None:
        plt.savefig(save, bbox_inches="tight", dpi=200)
    plt.show()