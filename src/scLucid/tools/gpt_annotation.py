"""
Trajectory and dynamics analysis for single-cell RNA-seq data.

This module provides a unified interface to popular trajectory inference and
RNA velocity methods, including PAGA, scVelo, and Monocle3 (via rtools).
"""

import logging
import os
from typing import Literal, Optional, Union

import anndata
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import scvelo as scv

from .rtools import RTools

log = logging.getLogger(__name__)


def run_trajectory_analysis(
    adata: anndata.AnnData,
    method: Literal["paga", "velocity", "monocle3"] = "paga",
    # PAGA params
    paga_groupby: Optional[str] = None,
    paga_root: Optional[Union[str, int]] = None,
    # Velocity params
    velocity_mode: Literal["stochastic", "dynamical"] = "stochastic",
    # Monocle3 params
    monocle3_root_group_key: Optional[str] = None,
    monocle3_root_group_name: Optional[str] = None,
    r_tools_instance: Optional[RTools] = None,
    # General params
    copy: bool = False,
    **kwargs,
) -> anndata.AnnData:
    """
    Run trajectory or dynamics analysis using a specified method.

    This function is a high-level wrapper that performs the core computations
    for trajectory inference or RNA velocity analysis.

    Args:
        adata: AnnData object.
        method: The analysis method to use:
            - 'paga': Computes cluster connectivity and pseudotime (requires clusters).
            - 'velocity': Computes RNA velocity to predict cell dynamics (requires 'spliced'/'unspliced' layers).
            - 'monocle3': Runs the Monocle3 workflow via R (requires RTools instance).
        paga_groupby: For 'paga', the key in adata.obs with clustering results.
        paga_root: For 'paga', the root cluster/cell for pseudotime calculation.
        velocity_mode: For 'velocity', the scVelo model to use.
        r_tools_instance: For 'monocle3', a pre-initialized RTools instance.
        copy: Whether to return a copy of the AnnData object.
        **kwargs: Additional arguments passed to the underlying tool (e.g., scv.tl.velocity).

    Returns:
        AnnData object with trajectory results in .obs, .uns, and .obsm.
    """
    if copy:
        adata = adata.copy()

    log.info(f"Running trajectory analysis with method: '{method}'")

    if method == "paga":
        if paga_groupby is None or paga_groupby not in adata.obs:
            raise ValueError("For 'paga' method, a valid 'paga_groupby' key is required.")

        log.info(f"Computing PAGA graph based on '{paga_groupby}'")
        sc.tl.paga(adata, groups=paga_groupby)

        if paga_root is not None:
            log.info(f"Calculating diffusion pseudotime with root: {paga_root}")
            adata.uns["iroot"] = np.flatnonzero(adata.obs[paga_groupby] == paga_root)[0]
            sc.tl.diffmap(adata)
            sc.tl.dpt(adata)
            adata.obs["trajectory_pseudotime"] = adata.obs["dpt_pseudotime"]
            log.info("Stored results in adata.obs['trajectory_pseudotime']")

    elif method == "velocity":
        if "spliced" not in adata.layers or "unspliced" not in adata.layers:
            raise ValueError(
                "For 'velocity' method, 'spliced' and 'unspliced' layers are required."
            )

        log.info(f"Running RNA velocity analysis using '{velocity_mode}' model")
        scv.pp.filter_and_normalize(adata)
        scv.pp.moments(adata)
        scv.tl.velocity(adata, mode=velocity_mode, **kwargs)
        scv.tl.velocity_graph(adata)

        if velocity_mode == "dynamical":
            log.info("Running dynamical model recovery and latent time...")
            scv.tl.recover_dynamics(adata)
            scv.tl.latent_time(adata)
            adata.obs["trajectory_pseudotime"] = adata.obs["latent_time"]

    elif method == "monocle3":
        if r_tools_instance is None:
            raise ValueError(
                "For 'monocle3' method, an initialized 'r_tools_instance' must be provided."
            )
        if monocle3_root_group_key is None or monocle3_root_group_name is None:
            raise ValueError(
                "For 'monocle3' method, 'monocle3_root_group_key' and 'monocle3_root_group_name' must be provided to define the trajectory start."
            )

        log.info(
            f"Running Monocle3 workflow via rtools, with root defined by '{monocle3_root_group_name}' in '{monocle3_root_group_key}'"
        )
        adata = r_tools_instance.run_monocle3(
            adata,
            root_group_key=monocle3_root_group_key,
            root_group_name=monocle3_root_group_name,
            key_added="trajectory_pseudotime",
        )

    else:
        raise ValueError(f"Unknown method: '{method}'. Choose from 'paga', 'velocity', 'monocle3'.")

    log.info("Trajectory analysis complete.")
    return adata


def plot_trajectory(
    adata: anndata.AnnData,
    method: Literal["paga", "velocity", "monocle3"],
    basis: str = "umap",
    color: Optional[str] = None,
    save_dir: Optional[str] = None,
):
    """
    Visualize trajectory and dynamics analysis results.

    Args:
        adata: AnnData object after running `run_trajectory_analysis`.
        method: The analysis method that was used.
        basis: Embedding to use for visualization (e.g., 'umap').
        color: Key in adata.obs to color cells by.
        save_dir: Directory to save plots.
    """
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        log.info(f"Saving trajectory plots to {save_dir}")

    if method == "paga":
        log.info("Plotting PAGA graph and pseudotime")
        plt.figure()
        sc.pl.paga_compare(adata, basis=basis, show=False)
        if save_dir:
            plt.savefig(os.path.join(save_dir, "paga_compare.png"))
        if not save_dir:
            plt.show()
        plt.close()

        if "trajectory_pseudotime" in adata.obs:
            plt.figure()
            sc.pl.embedding(
                adata, basis=basis, color="trajectory_pseudotime", cmap="viridis", show=False
            )
            if save_dir:
                plt.savefig(os.path.join(save_dir, "pseudotime_umap.png"))
            if not save_dir:
                plt.show()
            plt.close()

    elif method == "velocity":
        log.info("Plotting RNA velocity stream")
        plt.figure()
        scv.pl.velocity_embedding_stream(adata, basis=basis, color=color, show=False)
        if save_dir:
            plt.savefig(os.path.join(save_dir, "velocity_stream.png"))
        if not save_dir:
            plt.show()
        plt.close()

    elif method == "monocle3":
        log.info("Plotting Monocle3 pseudotime")
        if "trajectory_pseudotime" in adata.obs:
            plt.figure()
            sc.pl.embedding(
                adata,
                basis=basis,
                color="trajectory_pseudotime",
                cmap="viridis",
                show=False,
                title="Monocle3 Pseudotime",
            )
            if save_dir:
                plt.savefig(os.path.join(save_dir, "monocle3_pseudotime.png"))
            if not save_dir:
                plt.show()
            plt.close()
