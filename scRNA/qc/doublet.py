"""
Doublet detection for single-cell RNA-seq data.

This module provides functions for identifying potential doublet cells
using algorithmic methods like Scrublet and custom filtering approaches.
"""

import gc
import logging
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scrublet as scr
from anndata import AnnData

from ..utils.utils import identify_outliers

log = logging.getLogger(__name__)

__all__ = [
    "generate_doublet_rates",
    "predict_doublets",
]

def _identify_doublet_expression_patterns(
    adata: AnnData,
    sample_key: str = "sampleID",
    metrics: Optional[List[Tuple[str, str, Optional[float]]]] = None,
    nmads: float = 5.0,
) -> pd.Series:
    """
    Identify cells with expression patterns characteristic of doublets using heuristics.

    This function uses the identify_outliers utility to detect cells with outlier
    expression patterns that may indicate they are doublets (e.g., unusually high
    gene counts or expression of mutually exclusive markers).

    Args:
        adata: AnnData object containing single-cell data.
        sample_key: Key in adata.obs for sample identification.
        metrics: List of tuples (metric, direction, threshold) for detecting patterns.
                 If None, uses default metrics: upper outliers in gene counts and total counts.
        nmads: Number of MADs for outlier detection when threshold is not fixed.

    Returns:
        Boolean pd.Series indicating potential doublets based on expression patterns.
    """
    if metrics is None:
        # Default metrics focus on cells with unusually high complexity, a classic doublet sign.
        log.info("Using default heuristic metrics for doublet detection: log1p_n_genes_by_counts (upper) and log1p_total_counts (upper).")
        metrics = [
            ("log1p_n_genes_by_counts", "upper", None),
            ("log1p_total_counts", "upper", None),
        ]

    log.info("Identifying potential doublets based on heuristic expression patterns...")
    
    heuristic_doublets = identify_outliers(
        adata, metrics=metrics, sample_key=sample_key, nmads=nmads
    )
    
    return heuristic_doublets

def generate_doublet_rates(
    adata: AnnData,
    sample_key: str = "sampleID",
    rate_per_1000_cells: float = 0.008,
) -> Dict[str, float]:
    """
    Automatically generate expected doublet rates dictionary based on cell count per sample.

    This function is based on 10x Genomics general guideline: for every 1000 cells,
    the multiplet rate increases by approximately 0.8% (0.008).

    Args:
        adata (AnnData): AnnData object containing cell count information.
        sample_key (str, optional): Column name in adata.obs used to distinguish samples.
                                   Defaults to "sampleID".
        rate_per_1000_cells (float, optional): Expected doublet rate per 1000 cells.
                                              0.008 corresponds to standard 3' v3.1 kit.
                                              For high-throughput (HT) kits, this value may be higher (e.g., 0.016).
                                              Defaults to 0.008.

    Returns:
        Dict[str, float]: A dictionary with sample IDs as keys and calculated doublet rates as values.

    Example:
        >>> doublet_rates = generate_doublet_rates(adata, sample_key="sampleID")
        >>> print(doublet_rates)
        {'sample_A': 0.04, 'sample_B': 0.08}
    """
    log.info(
        "Automatically generating doublet rates based on cell counts per sample..."
    )
    # Group and sum cell counts by sample ID
    cell_counts = adata.obs[sample_key].value_counts()

    doublet_rates = {}
    for sample, n_cells in cell_counts.items():
        # Apply 10x Genomics linear formula
        rate = (n_cells / 1000) * rate_per_1000_cells
        # Cap the rate at a reasonable maximum, e.g., not exceeding 20%
        rate = min(rate, 0.20)
        doublet_rates[sample] = rate
        log.info(
            f"  - Sample '{sample}': {n_cells} cells -> Calculated doublet rate: {rate:.4f}"
        )

    return doublet_rates


def predict_doublets(
    adata: AnnData,
    method: str = "scrublet",
    sample_key: str = "sampleID",
    rate: Union[float, Dict[str, float]] = 0.1,
    use_heuristics: bool = True,
    heuristic_metrics: Optional[List[Tuple[str, str, Optional[float]]]] = None,
    plot_umap: bool = True,
    save_dir: Optional[str] = None,
    show: bool = True,
    **kwargs,
) -> AnnData:
    """
    Identify potential doublet cells using multiple complementary methods.

    Args:
        adata: AnnData object.
        method: Algorithmic method to use. Currently supports 'scrublet'.
        sample_key: Key in adata.obs to identify different samples.
        rate: Expected doublet rate. Can be a single float or a dict mapping samples to rates.
        use_heuristics: Whether to use expression-based heuristics as a complementary filter.
        heuristic_metrics: Custom metrics for the heuristic filter.
        plot_umap: For 'scrublet', whether to plot UMAP embedding with doublet scores.
        save_dir: Directory to save plots.
        show: Whether to display plots.
        **kwargs: Additional arguments passed to the specific method (e.g., n_pcs for scrublet).

    Returns:
        AnnData object with doublet scores and predictions in .obs.
    """
    adata.obs[f"{method}_score"] = np.nan
    adata.obs[f"{method}_predicted"] = False

    samples = adata.obs[sample_key].unique()

    if method == "scrublet":
        n_pcs = kwargs.get("n_pcs", 30)
        for sample in samples:
            log.info(f"Processing sample '{sample}' with Scrublet...")
            data_view = adata[adata.obs[sample_key] == sample]
            
            if data_view.n_obs < 10:
                log.warning(f"Skipping {sample}: fewer than 10 cells.")
                continue

            current_rate = rate.get(sample, 0.1) if isinstance(rate, dict) else rate
            actual_n_pcs = min(n_pcs, data_view.n_obs - 1, data_view.n_vars - 1)

            try:
                scrub = scr.Scrublet(data_view.X, expected_doublet_rate=current_rate)
                scores, predicted = scrub.scrub_doublets(n_prin_comps=actual_n_pcs, verbose=False)
                final_doublets = scrub.call_doublets(verbose=False)

                adata.obs.loc[data_view.obs.index, "scrublet_score"] = scores
                adata.obs.loc[data_view.obs.index, "scrublet_predicted"] = final_doublets
                log.info(f"  Found {sum(final_doublets)} potential doublets via Scrublet.")

                if plot_umap:
                    try:
                        scrub.set_embedding(
                            "UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)
                        )
                        fig = scrub.plot_embedding("UMAP", order_points=True)

                        if save_dir:
                            import os

                            os.makedirs(save_dir, exist_ok=True)
                            plt.savefig(
                                os.path.join(save_dir, f"{sample}_doublets_umap.png"),
                                dpi=300,
                            )
                        if show:
                            plt.show()
                        plt.close()
                    except Exception as e:
                        print(f"  Warning: Could not generate UMAP for doublets: {e}")
                    pass
                
            except Exception as e:
                    log.error(f"  Scrublet failed for sample {sample}: {e}")
            finally:
                    gc.collect()
    # NEW: Add elif for other methods here in the future
    # elif method == "solo":
    #     ...
    
    else:
        raise ValueError(f"Method '{method}' is not supported. Choose from ['scrublet'].")


    # --- Heuristic-based Doublet Identification ---
    if use_heuristics:
        heuristic_doublets = _identify_doublet_expression_patterns(
            adata, sample_key=sample_key, metrics=heuristic_metrics
        )
        adata.obs["heuristic_predicted"] = heuristic_doublets
        log.info(f"Found {heuristic_doublets.sum()} potential doublets via heuristics.")
    
    # --- Final Combination ---
    adata.obs["predicted_doublet"] = adata.obs[f"{method}_predicted"]
    if use_heuristics:
        adata.obs["predicted_doublet"] |= adata.obs["heuristic_predicted"]

    log.info("\n--- Overall Doublet Detection Summary ---")
    total_cells = adata.n_obs
    algo_count = adata.obs[f"{method}_predicted"].sum()
    log.info(f"Algorithm ({method}) predicted doublets: {algo_count} ({algo_count/total_cells:.2%})")
    if use_heuristics:
        heuristic_count = adata.obs["heuristic_predicted"].sum()
        log.info(f"Heuristic predicted doublets: {heuristic_count} ({heuristic_count/total_cells:.2%})")
    
    final_count = adata.obs["predicted_doublet"].sum()
    log.info(f"Total cells marked as doublets: {final_count} ({final_count/total_cells:.2%})")

    return adata
