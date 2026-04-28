"""HVG stability evaluation via bootstrap resampling.

Extracted from core.py for maintainability.
"""

from __future__ import annotations

import logging
import random
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

from ..config import HVGConfig
from .core import find_hvgs

log = logging.getLogger(__name__)

def evaluate_hvg_stability(
    adata: AnnData,
    hvg_key: str = "highly_variable",
    n_bootstrap: int = 20,
    sample_fraction: float = 0.8,
    method: str = "scanpy",
    flavor: str = "seurat",
    n_top_genes: Optional[int] = 2000,
    layer: Optional[str] = None,
    random_state: Optional[int] = 42,
    plot: bool = True,
    save_path: Optional[str] = None,
) -> AnnData:
    """
    Evaluates the stability of HVG selection through bootstrap resampling.
    Adds stability info to .uns['sclucid']['preprocess']['hvg_stability'].
    """

    if hvg_key not in adata.var:
        raise KeyError(f"HVG key '{hvg_key}' not found in adata.var")
    if n_bootstrap <= 0:
        raise ValueError("n_bootstrap must be a positive integer")
    if not 0 < sample_fraction < 1:
        raise ValueError("sample_fraction must be between 0 and 1")
    current_hvgs = set(adata.var_names[adata.var[hvg_key]])
    log.info(f"[HVG stability] Evaluating stability of {len(current_hvgs)} HVGs")
    if random_state is not None:
        random.seed(random_state)
        np.random.seed(random_state)
    gene_selection_count = dict.fromkeys(adata.var_names, 0)
    n_cells_per_bootstrap = int(adata.n_obs * sample_fraction)
    report_interval = max(1, n_bootstrap // 10)
    for i in range(n_bootstrap):
        if i % report_interval == 0:
            log.info(f"[HVG stability] Bootstrap iteration {i + 1}/{n_bootstrap}")
        cell_indices = np.random.choice(adata.n_obs, size=n_cells_per_bootstrap, replace=False)
        bootstrap_adata_view = adata[cell_indices, :]
        bootstrap_adata = sc.AnnData(X=bootstrap_adata_view.X, var=bootstrap_adata_view.var)
        find_hvgs(
            bootstrap_adata,
            HVGConfig(method=method, n_top_genes=n_top_genes, flavor=flavor),
            force=True,
            plot=False,
            input_layer="X",
        )
        bootstrap_hvgs = set(
            bootstrap_adata.var_names[bootstrap_adata.var[f"highly_variable_{method}_{flavor}"]]
            if method == "scanpy"
            else bootstrap_adata.var_names[bootstrap_adata.var[f"highly_variable_{method}"]]
        )
        for gene in bootstrap_hvgs:
            if gene in gene_selection_count:
                gene_selection_count[gene] += 1
    selection_frequency = {
        gene: count / n_bootstrap for gene, count in gene_selection_count.items()
    }
    adata.var["hvg_selection_frequency"] = pd.Series(
        [selection_frequency.get(gene, 0) for gene in adata.var_names],
        index=adata.var_names,
    )
    stability_score = np.mean([selection_frequency.get(gene, 0) for gene in current_hvgs])
    top_quartile = np.quantile([selection_frequency.get(gene, 0) for gene in current_hvgs], 0.75)
    bottom_quartile = np.quantile([selection_frequency.get(gene, 0) for gene in current_hvgs], 0.25)
    log.info("[HVG stability] Stability metrics:")
    log.info(f"  - Overall stability score: {stability_score:.3f}")
    log.info(f"  - Top 25% of HVGs selected with frequency >= {top_quartile:.3f}")
    log.info(f"  - Bottom 25% of HVGs selected with frequency <= {bottom_quartile:.3f}")
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["hvg_stability"] = {
        "overall_score": stability_score,
        "top_quartile": top_quartile,
        "bottom_quartile": bottom_quartile,
        "n_bootstrap": n_bootstrap,
        "sample_fraction": sample_fraction,
        "method": method,
    }
    if plot:
        try:
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))
            sns.histplot(adata.var["hvg_selection_frequency"], bins=30, kde=True, ax=axes[0])
            axes[0].set_title("HVG Selection Frequency Distribution")
            axes[0].set_xlabel("Selection Frequency")
            axes[0].set_ylabel("Number of Genes")
            axes[0].axvline(
                stability_score,
                color="red",
                linestyle="--",
                label=f"Avg HVG Stability: {stability_score:.3f}",
            )
            axes[0].legend()
            if "means" in adata.var:
                x = "means"
            elif f"{hvg_key}_means" in adata.var:
                x = f"{hvg_key}_means"
            else:
                if layer is None:
                    adata.var["temp_means"] = np.array(adata.X.mean(axis=0)).flatten()
                else:
                    adata.var["temp_means"] = np.array(adata.layers[layer].mean(axis=0)).flatten()
                x = "temp_means"
            scatter = axes[1].scatter(
                adata.var[x],
                adata.var["hvg_selection_frequency"],
                c=adata.var[hvg_key].astype(int),
                alpha=0.6,
                cmap="coolwarm",
                s=10,
            )
            axes[1].set_title("HVG Stability vs. Mean Expression")
            axes[1].set_xlabel("Mean Expression")
            axes[1].set_ylabel("Selection Frequency")
            cbar = plt.colorbar(scatter, ax=axes[1])
            cbar.set_ticks([0, 1])
            cbar.set_ticklabels(["Not HVG", "HVG"])
            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                log.info(f"Saved stability plot to {save_path}")
            if "temp_means" in adata.var:
                del adata.var["temp_means"]
        except Exception as e:
            log.warning(f"[HVG stability] Failed to create stability plot: {str(e)}")
    return adata


