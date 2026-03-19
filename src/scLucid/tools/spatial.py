"""
Unified spatial transcriptomics analysis for single-cell and spatial data.

Supports clustering, marker detection, spatial statistics, and visualization.
Integrates with Squidpy, Scanpy, and outputs results in structured AnnData.
"""

import logging
import os
from typing import Optional, Literal, List, Dict, Any

import anndata
import scanpy as sc
import pandas as pd
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt

log = logging.getLogger(__name__)

def run_spatial_analysis(
    adata: anndata.AnnData,
    spatial_key: str = "spatial",
    method: Literal["squidpy", "scanpy"] = "squidpy",
    cluster_key: str = "spatial_leiden",
    n_clusters: Optional[int] = None,
    spot_size: Optional[float] = None,
    marker_n_top: int = 10,
    spatial_neighbors: int = 6,
    compute_moran: bool = True,
    compute_lisi: bool = False,
    copy: bool = False,
    save_dir: Optional[str] = None,
    **kwargs,
) -> anndata.AnnData:
    """
    Unified spatial transcriptomics analysis.
    Results are stored in adata.uns['sclucid']['spatial'].
    """
    if copy:
        adata = adata.copy()
    adata.uns.setdefault('sclucid', {}).setdefault('spatial', {})
    spat_uns = adata.uns['sclucid']['spatial']

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        log.info(f"Saving spatial results to {save_dir}")

    log.info("Spatial analysis: method = %s", method)

    # --- Spatial neighbors graph ---
    log.info("Calculating spatial neighbors...")
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neigh=spatial_neighbors, key_added="spatial_neighbors")

    # --- Clustering ---
    log.info("Running spatial clustering...")
    if method == "squidpy":
        # PCA, neighbors, leiden
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.leiden(adata, key_added=cluster_key)
    elif method == "scanpy":
        sc.pp.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.louvain(adata, key_added=cluster_key)
    else:
        raise ValueError(f"Unknown method: {method}")

    spat_uns["cluster_key"] = cluster_key

    # --- Marker detection ---
    log.info("Detecting spatial marker genes...")
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon", n_genes=marker_n_top)
    spat_uns["marker_genes"] = sc.get.rank_genes_groups_df(adata, group=None)
    # Save marker table
    if save_dir:
        spat_uns["marker_genes"].to_csv(os.path.join(save_dir, "markers.csv"), index=False)

    # --- Moran's I spatial autocorrelation ---
    if compute_moran:
        log.info("Computing Moran's I spatial autocorrelation...")
        sq.gr.spatial_autocorr(adata, mode="moran")
        moran_df = sq.gr.spatial_autocorr_results(adata, mode="moran")
        spat_uns["moran_top"] = moran_df.sort_values("I", ascending=False).head(marker_n_top)
        if save_dir:
            spat_uns["moran_top"].to_csv(os.path.join(save_dir, "moran_top.csv"), index=False)

    # --- LISI for spatial diversity (optional) ---
    if compute_lisi:
        log.info("Computing LISI spatial diversity...")
        try:
            import lisi
            spatial_lisi = lisi.compute_lisi(adata.obsm[spatial_key], adata.obs[cluster_key])
            spat_uns["lisi"] = spatial_lisi
            if save_dir:
                np.savetxt(os.path.join(save_dir, "lisi.txt"), spatial_lisi)
        except ImportError:
            log.warning("lisi package not installed. Skipping LISI computation.")

    # --- Save cluster assignment to obs ---
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

    log.info("Spatial analysis complete.")
    return adata

def plot_spatial(
    adata: anndata.AnnData,
    spatial_key: str = "spatial",
    cluster_key: str = "spatial_leiden",
    color: Optional[str] = None,
    markers: Optional[List[str]] = None,
    save_dir: Optional[str] = None,
    spot_size: Optional[float] = None,
):
    """
    Plot spatial clustering, marker expression, and spatial autocorrelation.
    """
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)

    log.info("Plotting spatial cluster map...")
    sq.pl.spatial_scatter(
        adata,
        color=color or cluster_key,
        library_id=None,
        size=spot_size or 1.5,
        show=False,
    )
    if save_dir:
        plt.savefig(os.path.join(save_dir, f"spatial_cluster.png"))
    else:
        plt.show()
    plt.close()

    # Plot marker expression
    if markers:
        for gene in markers:
            log.info(f"Plotting spatial expression for {gene}")
            sq.pl.spatial_scatter(
                adata,
                color=gene,
                size=spot_size or 1.5,
                cmap="viridis",
                show=False,
            )
            if save_dir:
                plt.savefig(os.path.join(save_dir, f"spatial_{gene}.png"))
            else:
                plt.show()
            plt.close()

    # Plot Moran's I top genes
    if "sclucid" in adata.uns and "spatial" in adata.uns["sclucid"]:
        moran_top = adata.uns["sclucid"]["spatial"].get("moran_top")
        if moran_top is not None:
            for gene in moran_top["gene"].head(5):
                log.info(f"Plotting spatial autocorr for {gene}")
                sq.pl.spatial_scatter(
                    adata,
                    color=gene,
                    size=spot_size or 1.5,
                    cmap="coolwarm",
                    show=False,
                )
                if save_dir:
                    plt.savefig(os.path.join(save_dir, f"spatial_moran_{gene}.png"))
                else:
                    plt.show()
                plt.close()

    log.info("Spatial plotting complete.")

def run_spatial_batch(
    adatas: List[anndata.AnnData],
    out_dir: str,
    cluster_key: str = "spatial_leiden",
    method: Literal["squidpy", "scanpy"] = "squidpy",
    **kwargs,
) -> Dict[str, Optional[anndata.AnnData]]:
    """
    Batch process spatial analysis for multiple AnnData objects.
    """
    results = {}
    os.makedirs(out_dir, exist_ok=True)
    for i, adata in enumerate(adatas):
        sample_id = getattr(adata, 'sample_id', f"sample{i+1}")
        sample_dir = os.path.join(out_dir, sample_id)
        try:
            results[sample_id] = run_spatial_analysis(
                adata, save_dir=sample_dir, cluster_key=cluster_key, method=method, **kwargs
            )
            plot_spatial(results[sample_id], save_dir=sample_dir, cluster_key=cluster_key)
            log.info(f"Spatial analysis completed for {sample_id}")
        except Exception as e:
            log.error(f"Spatial analysis failed for {sample_id}: {e}")
            results[sample_id] = None
    return results

def export_spatial_report(
    adata: anndata.AnnData,
    out_dir: str,
    top_n: int = 10
):
    """
    Export spatial analysis results and summary report.
    """
    os.makedirs(out_dir, exist_ok=True)
    spatial = adata.uns.get("sclucid", {}).get("spatial", {})
    # Export marker genes and Moran's I results
    if "marker_genes" in spatial:
        spatial["marker_genes"].head(top_n).to_csv(os.path.join(out_dir, "top_markers.csv"), index=False)
    if "moran_top" in spatial:
        spatial["moran_top"].head(top_n).to_csv(os.path.join(out_dir, "top_moran.csv"), index=False)
    # Export parameters
    with open(os.path.join(out_dir, "spatial_params.txt"), "w") as f:
        for k, v in spatial.items():
            if not isinstance(v, (pd.DataFrame, np.ndarray)):
                f.write(f"{k}: {v}\n")
    log.info(f"Spatial report exported to: {out_dir}")