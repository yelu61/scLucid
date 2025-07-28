"""
Doublet detection for single-cell RNA-seq data.

This module provides functions for identifying potential doublet cells
using the Scrublet algorithm and custom filtering approaches.
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
    "is_doublet",
    "_identify_doublet_expression_patterns",
]


def _identify_doublet_expression_patterns(
    adata: AnnData,
    sample_key: str = "sampleID",
    metrics: List[Tuple[str, str, Optional[float]]] = None,
    nmads: float = 5.0,
) -> pd.Series:
    """
    Identify cells with expression patterns characteristic of doublets.

    This function uses the identify_outliers utility to detect cells with expression
    patterns that may indicate they are doublets.

    Args:
        adata: AnnData object containing single-cell data.
        sample_key: Key in adata.obs for sample identification.
        metrics: List of tuples (metric, direction, threshold) for detecting doublet-specific patterns.
                 If None, default metrics focused on doublet characteristics will be used.
        nmads: Number of MADs for outlier detection.

    Returns:
        Boolean pd.Series indicating potential doublets based on expression patterns.
    """
    if metrics is None:
        # Default metrics focused on doublet detection
        metrics = [
            ("pct_counts_in_top_20_genes", "upper", None),
            ("log1p_n_genes_by_counts", "upper", None),
        ]

    final_metrics = []
    for metric in metrics:
        if len(metric) == 3 and metric[2] is not None:
            metric_name, direction, threshold = metric

            if direction == "upper":
                adata.obs[f"outlier_{metric_name}"] = adata.obs[metric_name] > threshold
                log.info(f"Using fixed threshold {threshold} for {metric_name}")
            elif direction == "lower":
                adata.obs[f"outlier_{metric_name}"] = adata.obs[metric_name] < threshold
                log.info(f"Using fixed threshold {threshold} for {metric_name}")
            else:  # "both"
                raise ValueError("Fixed threshold not supported for direction 'both'")
        else:
            if len(metric) == 2:
                final_metrics.append((metric[0], metric[1], None))
            else:
                final_metrics.append(metric)

    # Pre-calculate log1p versions if needed
    for base_metric in ["total_counts", "n_genes_by_counts"]:
        log_metric = f"log1p_{base_metric}"
        if log_metric not in adata.obs.columns and base_metric in adata.obs.columns:
            adata.obs[log_metric] = np.log1p(adata.obs[base_metric])

    # Use the unified outlier detection function from utils
    log.info("Identifying cells with expression patterns characteristic of doublets...")
    if final_metrics:
        doublet_candidates_mad = identify_outliers(
            adata, metrics=final_metrics, sample_key=sample_key, nmads=nmads
        )
    else:
        doublet_candidates_mad = pd.Series(False, index=adata.obs_names)

    doublet_candidates = doublet_candidates_mad.copy()
    for metric in metrics:
        if len(metric) == 3 and metric[2] is not None:
            metric_name = metric[0]
            outlier_col = f"outlier_{metric_name}"
            if outlier_col in adata.obs:
                doublet_candidates = doublet_candidates | adata.obs[outlier_col]

    return doublet_candidates


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


def is_doublet(
    adata: AnnData,
    sample_key: str = "sampleID",
    rate: Union[float, Dict[str, float]] = 0.1,
    n_pcs: int = 30,
    threshold: Optional[float] = None,
    over_genes_q: float = 0.99,
    nmads: float = 5.0,
    check_expression_patterns: bool = True,
    doublet_expression_metrics: Optional[List[Tuple[str, str, Optional[float]]]] = None,
    top20_genes_threshold: Optional[float] = 80.0,
    use_fixed_top20_threshold: bool = False,
    plot_umap: bool = True,
    save_dir: Optional[str] = None,
    show: bool = True,
) -> AnnData:
    """
    Identify potential doublet cells using multiple complementary methods.

    This function combines algorithmic detection (Scrublet) with expression-based
    heuristics to comprehensively identify likely doublets.

    Args:
        adata: AnnData object containing single-cell data.
        sample_key: The key in adata.obs to identify different samples.
        rate: Expected doublet rate. Can be a single float for all samples,
              or a dictionary mapping sample IDs to specific rates.
        n_pcs: Number of principal components to use.
        threshold: Scrublet threshold for calling doublets. If None, scrublet auto-detects it.
        over_genes_q: Quantile threshold for identifying cells with an excessive
                      number of genes, as a secondary doublet indicator.
        nmads: Number of median absolute deviations (MADs) for defining outliers.
        check_expression_patterns: Whether to check for expression patterns
                                   characteristic of doublets.
        doublet_expression_metrics: List of metrics to use for identifying.
        top20_genes_threshold: Fixed threshold for pct_counts_in_top_20_genes if use_fixed_top20_threshold=True.
        use_fixed_top20_threshold: Whether to use fixed threshold instead of MAD for pct_counts_in_top_20_genes.
        plot_umap: Whether to plot UMAP embedding with doublet scores.
        save_dir: Directory to save plots. If None, plots are not saved.
        show: Whether to display the plots.

    Returns:
        AnnData object with doublet scores and predictions added to .obs.
    """
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet_scrublet"] = False

    if check_expression_patterns:
        adata.obs["doublet_expression_pattern"] = False

    samples = adata.obs[sample_key].unique()
    total_samples = len(samples)

    for i, sample in enumerate(samples):
        log.info(f"Processing sample {i + 1}/{total_samples}: {sample}")
        data_view = adata[adata.obs[sample_key] == sample]

        n_cells, n_features = data_view.shape
        if n_cells < 10:
            log.warning(
                f"Skipping doublet detection for sample {sample}: has fewer than 10 cells ({n_cells})."
            )
            continue

        if isinstance(rate, dict):
            current_rate = rate.get(sample, 0.1)
        else:
            current_rate = rate

        actual_n_pcs = min(n_pcs, n_cells - 1, n_features - 1)

        try:
            log.info(
                f"  Running Scrublet with n_pcs={actual_n_pcs} and expected_doublet_rate={current_rate:.3f}"
            )
            scrub = scr.Scrublet(data_view.X, expected_doublet_rate=current_rate)
            doublet_scores, predicted_doublets = scrub.scrub_doublets(
                n_prin_comps=actual_n_pcs, verbose=False
            )

            if threshold is None:
                final_doublets = scrub.call_doublets(verbose=False)
                log.info(f"  Auto-detected Scrublet threshold: {scrub.threshold_:.3f}")
            else:
                final_doublets = scrub.call_doublets(threshold=threshold, verbose=False)
                log.info(
                    f"  User-provided threshold: {threshold:.3f} (Scrublet auto-detected: {scrub.threshold_:.3f})"
                )

            log.info(f"  Found {sum(final_doublets)} potential doublets via Scrublet.")

            adata.obs.loc[data_view.obs.index, "doublet_score"] = doublet_scores
            adata.obs.loc[data_view.obs.index, "predicted_doublet_scrublet"] = (
                final_doublets
            )

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

        except Exception as e:
            log.error(f"  Scrublet failed for sample {sample}: {e}")
            adata.obs.loc[data_view.obs.index, "doublet_score"] = np.nan
            adata.obs.loc[data_view.obs.index, "predicted_doublet_scrublet"] = False
        finally:
            gc.collect()

    if check_expression_patterns:
        if doublet_expression_metrics is None:
            # Calculate quantile-based thresholds for high gene count
            gene_count_threshold = np.quantile(
                adata.obs["n_genes_by_counts"], over_genes_q
            )
            # Calculate MAD-based thresholds for other metrics
            if use_fixed_top20_threshold:
                doublet_expression_metrics = [
                    ("n_genes_by_counts", "upper", gene_count_threshold),
                    ("pct_counts_in_top_20_genes", "upper", None),
                    ("log1p_n_genes_by_counts", "upper", None),
                ]
                log.info(
                    f"Using fixed threshold {top20_genes_threshold}% for pct_counts_in_top_20_genes"
                )
            else:
                doublet_expression_metrics = [
                    ("n_genes_by_counts", "upper", gene_count_threshold),
                    ("pct_counts_in_top_20_genes", "upper", None),
                    ("log1p_n_genes_by_counts", "upper", None),
                ]
                log.info(
                    f"Using MAD-based threshold (nmads={nmads}) for pct_counts_in_top_20_genes"
                )

        log.info("Identifying potential doublets based on expression patterns...")
        expression_doublets = _identify_doublet_expression_patterns(
            adata,
            sample_key=sample_key,
            metrics=doublet_expression_metrics,
            nmads=nmads,
        )
        adata.obs["doublet_expression_pattern"] = expression_doublets

    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet_scrublet"]

    if check_expression_patterns:
        adata.obs["predicted_doublet"] |= adata.obs["doublet_expression_pattern"]

    log.info("\n--- Overall Doublet Detection Statistics ---")
    total_cells = adata.n_obs
    scrublet_count = adata.obs["predicted_doublet_scrublet"].sum()
    log.info(
        f"Scrublet-predicted doublets: {scrublet_count} ({scrublet_count / total_cells:.2%})"
    )
    custom_count = adata.obs["doublet_expression_pattern"].sum()
    log.info(
        f"Custom high-gene doublets (>{over_genes_q:.2f} quantile): {custom_count} ({custom_count / total_cells:.2%})"
    )
    combined = (
        adata.obs["predicted_doublet_scrublet"]
        | adata.obs["doublet_expression_pattern"]
    )
    combined_count = combined.sum()
    log.info(
        f"Total unique cells marked as potential doublets: {combined_count} ({combined_count / total_cells:.2%})"
    )

    return adata
