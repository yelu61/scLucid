"""HVG metrics visualization.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ..config import HVGConfig

log = logging.getLogger(__name__)

def plot_hvg_metrics(
    adata: AnnData,
    hvg_key: str = "highly_variable",
    n_top_genes: int = 20,
    metrics: Optional[List[str]] = None,
    show_gene_labels: bool = True,
    size_by_expr: bool = True,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Creates visualizations of HVG metrics to evaluate selection quality.
    """
    if hvg_key not in adata.var:
        raise KeyError(f"HVG key '{hvg_key}' not found in adata.var")
    available_metrics = {}
    for column in adata.var.columns:
        if any(x in column for x in ["mean", "disp", "var", "score"]):
            if "mean" in column:
                available_metrics["mean"] = column
            elif any(x in column for x in ["disp", "var"]):
                available_metrics["dispersion"] = column
            elif "norm" in column:
                available_metrics["norm_dispersion"] = column
            elif "score" in column:
                available_metrics["score"] = column
    if metrics is not None:
        for metric in metrics:
            if metric not in adata.var.columns:
                log.warning(f"Metric '{metric}' not found in adata.var columns")
    if "norm_dispersion" in available_metrics and "mean" in available_metrics:
        x = available_metrics["mean"]
        y = available_metrics["norm_dispersion"]
        plot_type = "dispersion_vs_mean"
    elif "dispersion" in available_metrics and "mean" in available_metrics:
        x = available_metrics["mean"]
        y = available_metrics["dispersion"]
        plot_type = "dispersion_vs_mean"
    elif "score" in available_metrics:
        x = available_metrics.get("mean")
        y = available_metrics["score"]
        plot_type = "score"
    else:
        method_specific_x = f"{hvg_key}_means"
        method_specific_y = f"{hvg_key}_dispersions_norm"
        if method_specific_x in adata.var and method_specific_y in adata.var:
            x = method_specific_x
            y = method_specific_y
            plot_type = "dispersion_vs_mean"
        else:
            log.warning("[HVG plot] Could not find appropriate metrics for plotting")
            adata.var["_temp_mean"] = np.array(adata.X.mean(axis=0)).flatten()
            x = "_temp_mean"
            if "hvg_selection_frequency" in adata.var:
                y = "hvg_selection_frequency"
                plot_type = "stability"
            else:
                raise ValueError("[HVG plot] Cannot create HVG plot: no appropriate metrics found")
    fig, ax = plt.subplots(figsize=(10, 8))
    if size_by_expr and "mean" in available_metrics:
        sizes = np.clip(adata.var[available_metrics["mean"]] * 20, 5, 200)
    else:
        sizes = 30
    scatter = ax.scatter(
        adata.var[x],
        adata.var[y],
        s=sizes,
        c=adata.var[hvg_key].astype(int),
        cmap="coolwarm",
        alpha=0.7,
    )
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_ticks([0, 1])
    cbar.set_ticklabels(["Not HVG", "HVG"])
    if plot_type == "dispersion_vs_mean":
        ax.set_xlabel("Mean Expression")
        ax.set_ylabel("Dispersion (normalized)")
        ax.set_title("Highly Variable Genes: Dispersion vs. Mean")
        ax.set_xscale("log")
    elif plot_type == "score":
        ax.set_xlabel("Mean Expression" if x else "Gene Index")
        ax.set_ylabel("HVG Score")
        ax.set_title("Highly Variable Genes: Score Distribution")
    elif plot_type == "stability":
        ax.set_xlabel("Mean Expression")
        ax.set_ylabel("Selection Frequency")
        ax.set_title("HVG Selection Stability")
    if show_gene_labels and n_top_genes > 0:
        hvg_mask = adata.var[hvg_key]
        if hvg_mask.sum() > 0:
            top_indices = adata.var.loc[hvg_mask, y].nlargest(n_top_genes).index
            for idx in top_indices:
                gene_name = idx
                ax.annotate(
                    gene_name,
                    (adata.var.loc[idx, x], adata.var.loc[idx, y]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
                )
    n_hvgs = adata.var[hvg_key].sum()
    total_genes = len(adata.var)
    ax.set_title(
        f"{ax.get_title()}\n{n_hvgs} HVGs selected ({n_hvgs / total_genes:.1%} of {total_genes} genes)"
    )
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        log.info(f"[HVG plot] Saved metrics plot to {save_path}")
    if "_temp_mean" in adata.var:
        del adata.var["_temp_mean"]
    return fig
