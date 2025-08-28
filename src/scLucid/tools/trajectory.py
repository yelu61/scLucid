"""
Unified trajectory and dynamics analysis for single-cell RNA-seq data.

Supports PAGA, scVelo, Monocle3, and Slingshot.
"""

import logging
from typing import Literal, Optional, Union
import os

import anndata
import numpy as np
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

from .rtools import RTools

log = logging.getLogger(__name__)

def run_trajectory_analysis(
    adata: anndata.AnnData,
    method: Literal["paga", "velocity", "monocle3", "slingshot"] = "paga",
    # PAGA
    paga_groupby: Optional[str] = None,
    paga_root: Optional[Union[str, int]] = None,
    paga_pseudotime_method: Literal["dpt", "paga_path"] = "dpt",
    # velocity
    velocity_mode: Literal["stochastic", "dynamical"] = "stochastic",
    min_shared_counts: int = 20,
    n_pcs: int = 30,
    n_neighbors: int = 30,
    # monocle3
    monocle3_root_group_key: Optional[str] = None,
    monocle3_root_group_name: Optional[str] = None,
    # slingshot
    slingshot_groupby: Optional[str] = None,
    slingshot_start: Optional[str] = None,
    r_tools_instance: Optional[RTools] = None,
    # General
    copy: bool = False,
    **kwargs,
) -> anndata.AnnData:
    """
    Unified interface for trajectory/dynamics inference.

    Results are stored in adata.uns['scrnatk']['trajectory'][method]
    """
    if copy:
        adata = adata.copy()
    adata.uns.setdefault('scrnatk', {}).setdefault('trajectory', {})
    traj_uns = adata.uns['scrnatk']['trajectory']

    log.info(f"Trajectory analysis: method = {method}")

    if method == "paga":
        # --- 自动聚类 ---
        if paga_groupby is None or paga_groupby not in adata.obs:
            log.info("No paga_groupby, running Leiden clustering...")
            sc.tl.leiden(adata, key_added="leiden")
            paga_groupby = "leiden"
        log.info(f"PAGA with groupby={paga_groupby}")

        sc.tl.paga(adata, groups=paga_groupby)
        traj_uns["paga"] = {"groupby": paga_groupby}

        # --- 伪时间 ---
        if paga_root is not None:
            log.info(f"Diffusion pseudotime, root={paga_root}")
            adata.uns["iroot"] = np.flatnonzero(adata.obs[paga_groupby] == paga_root)[0]
            sc.tl.diffmap(adata)
            if paga_pseudotime_method == "dpt":
                sc.tl.dpt(adata)
                adata.obs["trajectory_pseudotime"] = adata.obs["dpt_pseudotime"]
            # 可扩展paga_path等
            traj_uns["paga"]["root"] = paga_root

    elif method == "velocity":
        log.info(f"Running scVelo ({velocity_mode})")
        if "spliced" not in adata.layers or "unspliced" not in adata.layers:
            raise ValueError("AnnData must have 'spliced'/'unspliced' layers for velocity.")
        scv.pp.filter_and_normalize(adata, min_shared_counts=min_shared_counts)
        scv.pp.moments(adata, n_pcs=n_pcs, n_neighbors=n_neighbors)
        scv.tl.velocity(adata, mode=velocity_mode)
        if velocity_mode == "dynamical":
            scv.tl.recover_dynamics(adata)
            scv.tl.latent_time(adata)
            adata.obs["trajectory_pseudotime"] = adata.obs["latent_time"]
        scv.tl.velocity_graph(adata)
        traj_uns["velocity"] = {
            "mode": velocity_mode,
            "min_shared_counts": min_shared_counts,
            "n_pcs": n_pcs,
            "n_neighbors": n_neighbors,
        }

    elif method == "monocle3":
        if r_tools_instance is None:
            raise ValueError("r_tools_instance is required for 'monocle3'.")
        if monocle3_root_group_key is None or monocle3_root_group_name is None:
            raise ValueError("For 'monocle3', root_group_key and root_group_name must be given.")
        log.info(f"Running Monocle3: root={monocle3_root_name} in {monocle3_root_group_key}")
        adata = r_tools_instance.run_monocle3(
            adata,
            root_group_key=monocle3_root_group_key,
            root_group_name=monocle3_root_group_name,
            key_added="trajectory_pseudotime"
        )
        traj_uns["monocle3"] = {
            "root_group_key": monocle3_root_group_key,
            "root_group_name": monocle3_root_group_name
        }

    elif method == "slingshot":
        if r_tools_instance is None:
            raise ValueError("r_tools_instance is required for 'slingshot'.")
        if slingshot_groupby is None:
            raise ValueError("slingshot_groupby (cluster label) is required.")
        log.info(f"Running Slingshot: groupby={slingshot_groupby}, start={slingshot_start}")
        adata = r_tools_instance.run_slingshot(
            adata, groupby=slingshot_groupby, start=slingshot_start, key_added="trajectory_pseudotime"
        )
        traj_uns["slingshot"] = {
            "groupby": slingshot_groupby,
            "start": slingshot_start
        }

    else:
        raise ValueError(f"Unknown method: {method}")

    log.info("Trajectory analysis complete.")
    return adata


def plot_trajectory(
    adata: anndata.AnnData,
    method: Literal["paga", "velocity", "monocle3", "slingshot"],
    basis: str = "umap",
    color: Optional[str] = None,
    save_dir: Optional[str] = None,
    plot_types = ["embedding", "pseudotime", "tree"],
):
    """
    Visualize trajectory results for all supported methods.
    """
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        log.info(f"Saving plots to {save_dir}")

    # --- PAGA ---
    if method == "paga":
        if "embedding" in plot_types:
            plt.figure()
            sc.pl.paga_compare(adata, basis=basis, show=False)
            if save_dir: plt.savefig(os.path.join(save_dir, "paga_compare.png"))
            else: plt.show()
            plt.close()
        if "pseudotime" in plot_types and "trajectory_pseudotime" in adata.obs:
            plt.figure()
            sc.pl.embedding(
                adata, basis=basis, color="trajectory_pseudotime", cmap="viridis", show=False,
                title="PAGA pseudotime"
            )
            if save_dir: plt.savefig(os.path.join(save_dir, "pseudotime_umap.png"))
            else: plt.show()
            plt.close()
        if "tree" in plot_types:
            plt.figure()
            sc.pl.paga(adata, color=color, show=False)
            if save_dir: plt.savefig(os.path.join(save_dir, "paga_tree.png"))
            else: plt.show()
            plt.close()

    # --- Velocity ---
    elif method == "velocity":
        if "embedding" in plot_types:
            plt.figure()
            scv.pl.velocity_embedding_stream(
                adata, basis=basis, color=color, show=False
            )
            if save_dir: plt.savefig(os.path.join(save_dir, "velocity_stream.png"))
            else: plt.show()
            plt.close()
        if "pseudotime" in plot_types and "trajectory_pseudotime" in adata.obs:
            plt.figure()
            sc.pl.embedding(
                adata, basis=basis, color="trajectory_pseudotime", cmap="viridis", show=False,
                title="Velocity pseudotime"
            )
            if save_dir: plt.savefig(os.path.join(save_dir, "velocity_pseudotime_umap.png"))
            else: plt.show()
            plt.close()
        # 可扩展phase portrait等

    # --- Monocle3 ---
    elif method == "monocle3":
        if "embedding" in plot_types and "trajectory_pseudotime" in adata.obs:
            plt.figure()
            sc.pl.embedding(
                adata, basis=basis, color="trajectory_pseudotime", cmap="viridis", show=False,
                title="Monocle3 Pseudotime"
            )
            if save_dir: plt.savefig(os.path.join(save_dir, "monocle3_pseudotime.png"))
            else: plt.show()
            plt.close()

    # --- Slingshot ---
    elif method == "slingshot":
        if "embedding" in plot_types and "trajectory_pseudotime" in adata.obs:
            plt.figure()
            sc.pl.embedding(
                adata, basis=basis, color="trajectory_pseudotime", cmap="viridis", show=False,
                title="Slingshot Pseudotime"
            )
            if save_dir: plt.savefig(os.path.join(save_dir, "slingshot_pseudotime.png"))
            else: plt.show()
            plt.close()
        # 若需R主图导入，可用rpy2保存主图后os.system导入

    log.info("Trajectory plotting complete.")