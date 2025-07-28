"""
Gene regulatory network analysis using SCENIC.

This module provides a workflow to run pySCENIC for inferring regulons
(transcription factors and their target genes) and to analyze the resulting
regulon activity scores.
"""

import logging
import os
import subprocess
from typing import Optional

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

log = logging.getLogger(__name__)


def run_scenic(
    adata: anndata.AnnData,
    species: str,
    out_dir: str,
    scenic_db_dir: str,
    n_cpu: int = 8,
) -> anndata.AnnData:
    """
    Run the core pySCENIC workflow (GRN -> motifs -> AUCell).

    This function is a wrapper around the pySCENIC command-line interface,
    automating the three main steps of the analysis.

    Args:
        adata: AnnData object with raw counts in .X.
        species: Species for database selection ('hgnc' for human, 'mgi' for mouse).
        out_dir: Directory for all SCENIC output files.
        scenic_db_dir: Directory containing SCENIC databases (motif and TF files).
        n_cpu: Number of CPUs to use.

    Returns:
        AnnData object with SCENIC AUC matrix in .obsm['SCENIC_AUC'].
    """
    os.makedirs(out_dir, exist_ok=True)
    log.info(f"Starting pySCENIC workflow. Output will be in '{out_dir}'")

    # --- Step 0: Prepare input files ---
    raw_counts_file = os.path.join(out_dir, "raw_counts.loom")
    adata.write_loom(raw_counts_file)

    # Define database paths based on species
    if species == "hgnc":
        tf_names_file = os.path.join(scenic_db_dir, "allTFs_hgnc.txt")
        motif_db_file = os.path.join(
            scenic_db_dir, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
        )
    elif species == "mgi":
        tf_names_file = os.path.join(scenic_db_dir, "allTFs_mgi.txt")
        motif_db_file = os.path.join(
            scenic_db_dir, "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
        )
    else:
        raise ValueError("Species must be 'hgnc' (human) or 'mgi' (mouse)")

    # --- Step 1: GRN inference (arboreto) ---
    log.info("Step 1: Running arboreto for GRN inference...")
    adj_file = os.path.join(out_dir, "adj.tsv")
    cmd = (
        f"pyscenic grn {raw_counts_file} {tf_names_file} -o {adj_file} "
        f"--num_workers {n_cpu}"
    )
    subprocess.run(cmd, shell=True, check=True)

    # --- Step 2: Motif enrichment (cisTarget) ---
    log.info("Step 2: Running cisTarget for motif enrichment...")
    regulons_file = os.path.join(out_dir, "regulons.csv")
    cmd = (
        f"pyscenic ctx {adj_file} {motif_db_file} --annotations_fname {tf_names_file} "
        f"-o {regulons_file} --num_workers {n_cpu}"
    )
    subprocess.run(cmd, shell=True, check=True)

    # --- Step 3: AUCell scoring ---
    log.info("Step 3: Running AUCell to score regulon activity...")
    auc_matrix_file = os.path.join(out_dir, "auc_matrix.loom")
    cmd = (
        f"pyscenic aucell {raw_counts_file} {regulons_file} "
        f"-o {auc_matrix_file} --num_workers {n_cpu}"
    )
    subprocess.run(cmd, shell=True, check=True)

    # --- Step 4: Load results back into AnnData ---
    log.info("Loading results back into AnnData object")
    import loompy

    with loompy.connect(auc_matrix_file) as ds:
        auc_matrix = ds.ca.RegulonsAUC.T
        regulon_names = [r.split("(")[0] for r in ds.ra.Regulons]

    adata.obsm["SCENIC_AUC"] = auc_matrix
    adata.uns["SCENIC_regulons"] = pd.read_csv(regulons_file)
    adata.uns["SCENIC_regulon_names"] = regulon_names
    log.info("pySCENIC workflow complete.")

    return adata


def analyze_scenic_results(
    adata: anndata.AnnData,
    groupby: str,
    n_top_regulons: int = 12,
    save_dir: Optional[str] = None,
):
    """
    Analyze and visualize results from a SCENIC run.

    This function computes a SCENIC-based UMAP, performs clustering on the
    AUC matrix, and generates plots to visualize top regulon activities.
    """
    if "SCENIC_AUC" not in adata.obsm:
        raise ValueError("SCENIC results not found. Please run `run_scenic` first.")

    log.info("Analyzing SCENIC results and generating visualizations")
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    # --- 1. SCENIC-based UMAP and Clustering ---
    log.info("Computing SCENIC-based UMAP and Leiden clustering")
    sc.pp.neighbors(adata, use_rep="SCENIC_AUC", key_added="scenic_neighbors")
    sc.tl.umap(adata, neighbors_key="scenic_neighbors")
    sc.tl.leiden(adata, neighbors_key="scenic_neighbors", key_added="scenic_leiden")

    # --- 2. UMAP Visualization ---
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(adata, color=groupby, ax=axes[0], show=False, title=f"UMAP by {groupby}")
    sc.pl.umap(
        adata,
        color="scenic_leiden",
        ax=axes[1],
        show=False,
        title="UMAP by SCENIC clusters",
        legend_loc="on data",
    )

    if save_dir:
        plt.savefig(os.path.join(save_dir, "scenic_umap_clusters.png"))
        plt.close()
    else:
        plt.show()

    # --- 3. Plot Top Variable Regulons ---
    log.info(f"Plotting activity of top {n_top_regulons} variable regulons")
    regulon_names = adata.uns.get("SCENIC_regulon_names", [])
    auc_matrix = adata.obsm["SCENIC_AUC"]
    regulon_vars = np.var(auc_matrix, axis=0)
    top_regulon_idx = np.argsort(regulon_vars)[-n_top_regulons:]
    top_regulons = [regulon_names[i] for i in top_regulon_idx]

    save_path = (
        os.path.join(save_dir, "top_regulons_activity.png") if save_dir else None
    )
    sc.pl.umap(
        adata,
        color=top_regulons,
        title=top_regulons,
        save=save_path,
        show=not save_dir,
        cmap="viridis",
        ncols=4,
    )
