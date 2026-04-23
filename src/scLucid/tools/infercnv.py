"""
Copy number variation (CNV) analysis for single-cell RNA-seq data.

This module provides functions for inferring copy number variations from
scRNA-seq data, which is particularly useful for identifying tumor cells and
analyzing genomic aberrations in cancer samples.
"""

import logging
import os
from typing import List, Optional, Tuple, Union

import anndata as ad
import infercnvpy as cnv
import matplotlib.pyplot as plt
import numpy as np

# Configure logging
log = logging.getLogger(__name__)

# Export public functions
__all__ = [
    "run_cnv_analysis",
    "find_tumor",
]


# --- Main Functions ---
def find_tumor(
    adata: ad.AnnData,
    cnv_score_key: str = "cnv_score",
    alpha: float = 2.0,
    key_added: str = "tumor",
    percentile_threshold: Optional[float] = None,
    min_tumor_fraction: float = 0.05,
    max_tumor_fraction: float = 0.8,
    plot: bool = False,
    copy: bool = False,
) -> ad.AnnData:
    """
    Identify tumor cells based on CNV scores.

    This function classifies cells as tumor or non-tumor based on their copy number
    variation scores. It uses either a data-driven threshold based on the distribution
    of CNV scores or a user-specified percentile threshold.

    Args:
        adata: AnnData object with CNV scores in obs
        cnv_score_key: Column name in adata.obs containing CNV scores
        alpha: Number of standard deviations for identifying tumor threshold
        key_added: Name of the column to add with tumor classification (0/1)
        percentile_threshold: Optional percentile to use as threshold (e.g., 0.75)
        min_tumor_fraction: Minimum fraction of cells that should be classified as tumor
        max_tumor_fraction: Maximum fraction of cells that should be classified as tumor
        plot: Whether to plot the CNV score distribution with threshold
        copy: Whether to return a copy of adata or modify in place

    Returns:
        AnnData object with added tumor classification column

    Raises:
        ValueError: If CNV score column is not found in adata.obs
    """
    if copy:
        adata = adata.copy()

    # Check if CNV score exists
    if cnv_score_key not in adata.obs_keys():
        log.error(f"'{cnv_score_key}' not found in adata.obs")
        raise ValueError(f"'{cnv_score_key}' not found in adata.obs, please run infercnvpy first")

    log.info(f"Identifying tumor cells based on '{cnv_score_key}'")

    # Get CNV scores
    cnv_scores = adata.obs[cnv_score_key].values

    # Determine threshold
    if percentile_threshold is not None:
        # Use percentile-based threshold
        threshold = np.percentile(cnv_scores, percentile_threshold * 100)
        log.info(f"Using {percentile_threshold:.2f} percentile as threshold: {threshold:.4f}")
    else:
        # Use distribution-based threshold
        # First, find the largest gap in sorted CNV scores
        sorted_scores = np.sort(cnv_scores)
        gaps = np.diff(sorted_scores)

        if len(gaps) > 0:
            # Find the index of the largest gap
            gap_idx = np.argmax(gaps)
            candidate_threshold = sorted_scores[gap_idx + 1]

            # Now check if this threshold results in a reasonable tumor fraction
            tumor_fraction = np.mean(cnv_scores >= candidate_threshold)

            if tumor_fraction < min_tumor_fraction:
                # Too few tumor cells, use a lower threshold
                log.warning(
                    f"Gap-based threshold would classify only {tumor_fraction:.2%} of cells as tumor"
                )
                threshold = np.percentile(cnv_scores, 100 - min_tumor_fraction * 100)
                log.info(
                    f"Using lower threshold to ensure {min_tumor_fraction:.2%} tumor cells: {threshold:.4f}"
                )
            elif tumor_fraction > max_tumor_fraction:
                # Too many tumor cells, use a higher threshold
                log.warning(
                    f"Gap-based threshold would classify {tumor_fraction:.2%} of cells as tumor"
                )
                threshold = np.percentile(cnv_scores, 100 - max_tumor_fraction * 100)
                log.info(
                    f"Using higher threshold to limit to {max_tumor_fraction:.2%} tumor cells: {threshold:.4f}"
                )
            else:
                # Gap-based threshold is reasonable
                threshold = candidate_threshold
                log.info(
                    f"Using gap-based threshold: {threshold:.4f} ({tumor_fraction:.2%} tumor cells)"
                )
        else:
            # Fall back to a simple statistical threshold
            mean_score = np.mean(cnv_scores)
            std_score = np.std(cnv_scores)
            threshold = mean_score + alpha * std_score
            log.info(f"Using statistical threshold (mean + {alpha} SD): {threshold:.4f}")

    # Classify cells
    adata.obs[key_added] = (adata.obs[cnv_score_key] >= threshold).astype(int)
    tumor_fraction = adata.obs[key_added].mean()
    log.info(
        f"Classified {tumor_fraction:.2%} of cells ({np.sum(adata.obs[key_added])}/{len(adata)}) as tumor cells"
    )

    # Plot if requested
    if plot:
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot histogram of CNV scores
        ax.hist(cnv_scores, bins=50, alpha=0.7, color="lightblue", edgecolor="black")

        # Add vertical line for threshold
        ax.axvline(x=threshold, color="red", linestyle="--", linewidth=2)

        # Add text annotation
        ax.text(
            x=threshold + 0.05 * (max(cnv_scores) - min(cnv_scores)),
            y=0.9 * ax.get_ylim()[1],
            s=f"Threshold: {threshold:.4f}\nTumor cells: {tumor_fraction:.2%}",
            color="red",
            bbox=dict(facecolor="white", alpha=0.8),
        )

        # Set labels
        ax.set_xlabel("CNV Score", fontsize=12)
        ax.set_ylabel("Number of Cells", fontsize=12)
        ax.set_title("Distribution of CNV Scores", fontsize=14)

        # Show grid
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    return adata


# --- Main Functions ---
def run_cnv_analysis(
    adata: ad.AnnData,
    sample_key: Optional[str] = None,
    ref_obs: str = "cell_type",
    ref_keys: Union[str, List[str]] = "Immune",
    window_size: int = 250,
    step: int = 1,
    plot_heatmap: bool = True,
    heatmap_groupby: Optional[str] = None,
    plot_umap: bool = True,
    plot_tumor: bool = True,
    figsize: Tuple[float, float] = (15, 4),
    find_tumor_cells: bool = True,
    tumor_finding_alpha: float = 2.0,
    percentile_threshold: Optional[float] = None,
    save_dir: Optional[str] = None,
    key_added: str = "cnv",
    copy: bool = False,
) -> ad.AnnData:
    """
    Perform copy number variation (CNV) analysis on single-cell data.

    This function runs CNV inference, tumor cell identification, and generates
    visualizations to explore copy number variations in the dataset.

    Args:
        adata: AnnData object containing single-cell data
        sample_key: Optional key in adata.obs for sample stratification
        ref_obs: Column name in adata.obs containing cell type annotations
        ref_keys: Cell type(s) to use as reference (normal cells)
        window_size: Window size for CNV analysis (number of genes)
        step: Step size for CNV analysis
        plot_heatmap: Whether to plot chromosome heatmap
        heatmap_groupby: Column name for grouping cells in heatmap
        plot_umap: Whether to plot UMAP embedding with CNV scores and clusters
        plot_tumor: Whether to plot UMAP embedding with predicted tumor cells
        figsize: Figure size for UMAP plots
        find_tumor_cells: Whether to identify tumor cells
        tumor_finding_alpha: Number of standard deviations for identifying tumor threshold
        alpha: Number of standard deviations for identifying tumor threshold
        percentile_threshold: Optional percentile to use as threshold
        save_dir: Directory to save plots (None for no saving)
        key_added: Prefix for CNV analysis results
        copy: Whether to return a copy of adata

    Returns:
        AnnData object with CNV analysis results

    Notes:
        This function uses infercnvpy for CNV inference. For each sample (if sample_key
        is provided) or for the entire dataset, it:
        1. Infers CNVs using reference cell types
        2. Calculates CNV scores
        3. Identifies tumor cells
        4. Generates visualizations
    """
    if copy:
        adata = adata.copy()

    log.info("Starting CNV analysis")

    # Create save directory if needed
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        log.info(f"Created save directory: {save_dir}")

    # Set up heatmap groupby if not specified
    if heatmap_groupby is None:
        heatmap_groupby = ref_obs

    # Initialize columns for results
    adata.obs[f"{key_added}_score"] = 0.0
    adata.obs[f"{key_added}_tumor"] = 0

    # Convert ref_keys to list if it's a string
    if isinstance(ref_keys, str):
        ref_keys = [ref_keys]

    log.info(f"Using {', '.join(ref_keys)} from column '{ref_obs}' as reference cells")

    # Determine samples to process
    if sample_key is not None and sample_key in adata.obs:
        samples = adata.obs[sample_key].unique()
        log.info(f"Processing {len(samples)} samples from '{sample_key}'")
    else:
        samples = [None]
        log.info("Processing entire dataset as one sample")

    # Process each sample
    for sample in samples:
        # Select data for this sample
        if sample is not None:
            log.info(f"Processing sample: {sample}")
            data = adata[adata.obs[sample_key] == sample].copy()
            sample_name = str(sample)
        else:
            log.info("Processing entire dataset")
            data = adata.copy()
            sample_name = "all"

        # Run CNV inference
        try:
            log.info(f"Running infercnv with window_size={window_size}, step={step}")
            cnv.tl.infercnv(
                data,
                reference_key=ref_obs,
                reference_cat=ref_keys,
                window_size=window_size,
                key_added=key_added,
                step=step,
            )

            # Calculate CNV scores
            cnv.tl.cnv_score(data, key_added=key_added)
            log.info("CNV scores calculated")

            # Generate visualizations
            if plot_heatmap or plot_umap or plot_tumor:
                # Run dimensionality reduction if needed
                if not all(x in data.obsm for x in [f"X_{key_added}_pca", f"X_{key_added}_umap"]):
                    log.info("Running dimensionality reduction for visualizations")
                    cnv.tl.pca(data, key_added=key_added)
                    cnv.pp.neighbors(data, key_added=key_added)
                    cnv.tl.umap(data, key_added=key_added)
                    cnv.tl.leiden(data, key_added=key_added)

            # Plot chromosome heatmap
            if plot_heatmap:
                log.info(f"Plotting chromosome heatmap grouped by {heatmap_groupby}")

                # Save path for heatmap
                if save_dir is not None:
                    heatmap_path = os.path.join(save_dir, f"cnv_heatmap_{sample_name}.png")
                else:
                    heatmap_path = None

                # Plot heatmap by cell type
                cnv.pl.chromosome_heatmap(
                    data, groupby=heatmap_groupby, key_added=key_added, save=heatmap_path
                )

                # Plot heatmap by CNV clusters
                if save_dir is not None:
                    leiden_heatmap_path = os.path.join(
                        save_dir, f"cnv_heatmap_leiden_{sample_name}.png"
                    )
                else:
                    leiden_heatmap_path = None

                cnv.pl.chromosome_heatmap(
                    data,
                    groupby=f"{key_added}_leiden",
                    dendrogram=True,
                    key_added=key_added,
                    save=leiden_heatmap_path,
                )

            # Plot UMAP visualization
            if plot_umap:
                log.info("Creating UMAP visualization of CNV analysis")

                fig, axes = plt.subplots(nrows=1, ncols=3, figsize=figsize)

                # Plot Leiden clusters
                cnv.pl.umap(
                    data,
                    color=f"{key_added}_leiden",
                    legend_loc="on data",
                    legend_fontoutline=2,
                    key_added=key_added,
                    ax=axes[0],
                    show=False,
                )
                axes[0].set_title("UMAP (Leiden Clusters)")

                # Plot CNV scores
                cnv.pl.umap(
                    data, color=f"{key_added}_score", key_added=key_added, ax=axes[1], show=False
                )
                axes[1].set_title("UMAP (CNV Score)")

                # Plot cell types
                cnv.pl.umap(
                    data, color=heatmap_groupby, key_added=key_added, ax=axes[2], show=False
                )
                axes[2].set_title(f"UMAP ({heatmap_groupby})")

                plt.tight_layout()

                if save_dir is not None:
                    umap_path = os.path.join(save_dir, f"umap_cnv_{sample_name}.png")
                    plt.savefig(umap_path, bbox_inches="tight", dpi=300)
                    log.info(f"Saved UMAP plot to {umap_path}")
                else:
                    plt.show()

            # Identify tumor cells
            if find_tumor_cells:
                log.info("Identifying tumor cells within the sample...")
                data = find_tumor(
                    data,
                    cnv_score_key=f"{key_added}_score",
                    key_added=f"{key_added}_tumor",
                    alpha=tumor_finding_alpha,
                    percentile_threshold=percentile_threshold,
                    plot=False,
                )

            # Plot tumor cells
            if plot_tumor:
                log.info("Creating UMAP visualization of tumor cells")

                if save_dir is not None:
                    tumor_path = os.path.join(save_dir, f"umap_tumor_{sample_name}.png")
                else:
                    tumor_path = None

                cnv.pl.umap(data, color=f"{key_added}_tumor", key_added=key_added, save=tumor_path)

            # Copy results back to main AnnData object
            if sample is not None:
                adata.obs.loc[data.obs.index, f"{key_added}_score"] = data.obs[f"{key_added}_score"]
                adata.obs.loc[data.obs.index, f"{key_added}_tumor"] = data.obs[f"{key_added}_tumor"]
            else:
                adata.obs[f"{key_added}_score"] = data.obs[f"{key_added}_score"]
                adata.obs[f"{key_added}_tumor"] = data.obs[f"{key_added}_tumor"]

            # Save the CNV matrix if available
            if f"X_{key_added}" in data.layers:
                if f"X_{key_added}" not in adata.layers:
                    adata.layers[f"X_{key_added}"] = np.zeros((adata.n_obs, adata.n_vars))

                if sample is not None:
                    # We need to align the genes
                    mask = adata.obs[sample_key] == sample
                    adata.layers[f"X_{key_added}"][mask] = data.layers[f"X_{key_added}"]
                else:
                    adata.layers[f"X_{key_added}"] = data.layers[f"X_{key_added}"]

            log.info(f"Completed CNV analysis for sample: {sample_name}")

        except Exception as e:
            log.error(f"Error in CNV analysis for sample {sample_name}: {str(e)}")
            raise

    # Add summary statistics to uns
    adata.uns[f"{key_added}_stats"] = {
        "window_size": window_size,
        "step": step,
        "reference_cell_types": ref_keys,
        "tumor_cell_count": int(adata.obs[f"{key_added}_tumor"].sum()),
        "tumor_cell_fraction": float(adata.obs[f"{key_added}_tumor"].mean()),
        "mean_cnv_score": float(adata.obs[f"{key_added}_score"].mean()),
        "min_cnv_score": float(adata.obs[f"{key_added}_score"].min()),
        "max_cnv_score": float(adata.obs[f"{key_added}_score"].max()),
    }

    log.info(
        f"CNV analysis completed. Found {adata.uns[f'{key_added}_stats']['tumor_cell_count']} tumor cells "
        f"({adata.uns[f'{key_added}_stats']['tumor_cell_fraction']:.2%} of total)"
    )

    return adata
