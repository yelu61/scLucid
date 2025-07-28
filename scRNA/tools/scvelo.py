"""
RNA Velocity analysis using scVelo.

This module provides a streamlined workflow for running and visualizing RNA
velocity, which predicts the future state of cells based on spliced and
unspliced mRNA counts.
"""

import logging
import os
from typing import Optional

import anndata
import scvelo as scv

log = logging.getLogger(__name__)


def run_velocity_analysis(
    adata: anndata.AnnData,
    mode: str = "stochastic",
    min_shared_counts: int = 20,
    n_pcs: int = 30,
    n_neighbors: int = 30,
    copy: bool = False,
) -> anndata.AnnData:
    """
    Perform core RNA velocity analysis computations.

    This function runs the main steps of a scVelo analysis, including preprocessing,
    moment calculation, velocity estimation, and graph construction.

    Args:
        adata: AnnData object with 'spliced' and 'unspliced' layers.
        mode: Velocity model to use ('stochastic', 'deterministic', or 'dynamical').
        min_shared_counts: Minimum shared counts for gene filtering.
        n_pcs: Number of principal components for neighbor calculations.
        n_neighbors: Number of neighbors for velocity graph construction.
        copy: Whether to return a copy of adata.

    Returns:
        AnnData object with velocity results computed.
    """
    if copy:
        adata = adata.copy()

    log.info(f"Starting RNA velocity analysis using '{mode}' mode")
    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise ValueError(
            "AnnData object must contain 'spliced' and 'unspliced' layers."
        )

    # scVelo Preprocessing
    scv.pp.filter_and_normalize(adata, min_shared_counts=min_shared_counts)
    scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)

    # Velocity Estimation
    log.info("Estimating RNA velocity...")
    scv.tl.velocity(adata, mode=mode)
    if mode == "dynamical":
        log.info("Running dynamical model...")
        scv.tl.recover_dynamics(adata)
        scv.tl.latent_time(adata)

    # Velocity Graph and Embedding Projection
    log.info("Constructing velocity graph...")
    scv.tl.velocity_graph(adata)

    log.info("RNA velocity analysis complete.")
    return adata


def plot_velocity_results(
    adata: anndata.AnnData,
    basis: str = "umap",
    color: Optional[str] = None,
    stream: bool = True,
    plot_genes: bool = True,
    n_top_genes: int = 6,
    save_dir: Optional[str] = None,
):
    """
    Visualize RNA velocity results.

    This function generates key plots for interpreting RNA velocity, including
    the velocity stream on an embedding and phase portraits of top driver genes.

    Args:
        adata: AnnData object after running `run_velocity_analysis`.
        basis: Embedding to use for visualization (e.g., 'umap', 'tsne').
        color: Key in adata.obs to color cells by.
        stream: If True, plot velocity stream; otherwise, plot grid/arrows.
        plot_genes: Whether to plot phase portraits for top velocity genes.
        n_top_genes: Number of top velocity genes to plot.
        save_dir: Directory to save plots.
    """
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        log.info(f"Saving velocity plots to {save_dir}")

    # Plot velocity stream
    log.info(f"Plotting velocity on '{basis}' embedding, colored by '{color}'")
    plot_func = (
        scv.pl.velocity_embedding_stream if stream else scv.pl.velocity_embedding_grid
    )

    save_path = (
        os.path.join(save_dir, f"velocity_stream_{basis}.png") if save_dir else None
    )
    plot_func(
        adata, basis=basis, color=color, save=save_path, show=not save_dir, dpi=300
    )

    # Plot top velocity genes
    if plot_genes:
        log.info(f"Plotting phase portraits for top {n_top_genes} velocity genes")
        scv.tl.rank_velocity_genes(adata, n_genes=n_top_genes)
        top_genes = adata.var["velocity_genes"].index[:n_top_genes]

        save_path = (
            os.path.join(save_dir, "top_velocity_genes.png") if save_dir else None
        )
        scv.pl.velocity(
            adata,
            top_genes,
            ncols=min(3, n_top_genes),
            save=save_path,
            show=not save_dir,
            dpi=300,
        )
