"""
Gene regulatory network analysis using SCENIC.

This module provides a workflow to run pySCENIC for inferring regulons
(transcription factors and their target genes) and to analyze the resulting
regulon activity scores.
"""


import logging
import os
import subprocess
import datetime
from typing import Optional, List, Dict, Union

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
    Results are saved to adata.uns['sclucid']['scenic'].
    """
    os.makedirs(out_dir, exist_ok=True)
    log.info(f"Starting pySCENIC workflow. Output: {out_dir}")

    # --- Step 0: Prepare input files ---
    raw_counts_file = os.path.join(out_dir, "raw_counts.loom")
    adata.write_loom(raw_counts_file)

    # Define database paths by species
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
    cmd1 = (
        f"pyscenic grn {raw_counts_file} {tf_names_file} -o {adj_file} "
        f"--num_workers {n_cpu}"
    )
    subprocess.run(cmd1, shell=True, check=True)

    # --- Step 2: Motif enrichment (cisTarget) ---
    log.info("Step 2: Running cisTarget for motif enrichment...")
    regulons_file = os.path.join(out_dir, "regulons.csv")
    cmd2 = (
        f"pyscenic ctx {adj_file} {motif_db_file} --annotations_fname {tf_names_file} "
        f"-o {regulons_file} --num_workers {n_cpu}"
    )
    subprocess.run(cmd2, shell=True, check=True)

    # --- Step 3: AUCell scoring ---
    log.info("Step 3: Running AUCell to score regulon activity...")
    auc_matrix_file = os.path.join(out_dir, "auc_matrix.loom")
    cmd3 = (
        f"pyscenic aucell {raw_counts_file} {regulons_file} "
        f"-o {auc_matrix_file} --num_workers {n_cpu}"
    )
    subprocess.run(cmd3, shell=True, check=True)

    # --- Step 4: Load results back into AnnData ---
    log.info("Loading SCENIC results into AnnData")
    import loompy
    with loompy.connect(auc_matrix_file) as ds:
        auc_matrix = ds.ca.RegulonsAUC.T
        regulon_names = [r.split("(")[0] for r in ds.ra.Regulons]

    adata.obsm["SCENIC_AUC"] = auc_matrix
    adata.uns["SCENIC_regulons"] = pd.read_csv(regulons_file)
    adata.uns["SCENIC_regulon_names"] = regulon_names

    # --- Step 5: Structured results for downstream use ---
    adata.uns.setdefault("sclucid", {}).setdefault("scenic", {})
    adata.uns["sclucid"]["scenic"]["params"] = {
        "species": species,
        "scenic_db_dir": scenic_db_dir,
        "n_cpu": n_cpu,
        "out_dir": out_dir,
        "run_time": datetime.datetime.now().isoformat(),
    }
    adata.uns["sclucid"]["scenic"]["regulons"] = adata.uns["SCENIC_regulons"]
    adata.uns["sclucid"]["scenic"]["regulon_names"] = regulon_names
    adata.uns["sclucid"]["scenic"]["auc_matrix"] = auc_matrix

    log.info("pySCENIC workflow complete.")
    return adata

def run_scenic_batch(
    adatas: List[anndata.AnnData],
    species: str,
    out_dir: str,
    scenic_db_dir: str,
    n_cpu: int = 8,
    sample_ids: Optional[List[str]] = None,
) -> Dict[str, Optional[anndata.AnnData]]:
    """
    Batch run SCENIC for a list of AnnData objects.
    Returns dict[sample_id] = AnnData with SCENIC results or None if failed.
    """
    results = {}
    if sample_ids is None:
        sample_ids = [f"sample{i+1}" for i in range(len(adatas))]
    for adata, sid in zip(adatas, sample_ids):
        sample_dir = os.path.join(out_dir, sid)
        try:
            results[sid] = run_scenic(
                adata, species=species, out_dir=sample_dir,
                scenic_db_dir=scenic_db_dir, n_cpu=n_cpu
            )
            log.info(f"SCENIC completed for {sid}")
        except Exception as e:
            log.error(f"SCENIC failed for {sid}: {e}")
            results[sid] = None
    return results

def run_scenic_by_group(
    adata: anndata.AnnData,
    groupby: str,
    species: str,
    out_dir: str,
    scenic_db_dir: str,
    n_cpu: int = 8,
    min_cells: int = 100,
) -> Dict[str, Optional[anndata.AnnData]]:
    """
    Run SCENIC for each group (e.g., cluster/celltype) in AnnData.
    Returns dict[group] = AnnData with SCENIC results.
    """
    results = {}
    groups = adata.obs[groupby].unique()
    for group in groups:
        adata_sub = adata[adata.obs[groupby] == group].copy()
        if adata_sub.n_obs < min_cells:
            log.warning(f"Group {group} skipped (too few cells: {adata_sub.n_obs})")
            continue
        group_dir = os.path.join(out_dir, f"{group}")
        try:
            results[group] = run_scenic(
                adata_sub, species, group_dir, scenic_db_dir, n_cpu=n_cpu
            )
            log.info(f"SCENIC completed for group {group}")
        except Exception as e:
            log.error(f"SCENIC failed for group {group}: {e}")
            results[group] = None
    return results

def analyze_scenic_results(
    adata: anndata.AnnData,
    groupby: str,
    n_top_regulons: int = 12,
    save_dir: Optional[str] = None,
):
    """
    Analyze and visualize SCENIC results. UMAP, clusters, top regulons.
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

def export_scenic_report(
    adata: anndata.AnnData,
    out_dir: str,
    top_n: int = 10,
    groupby: Optional[str] = None,
):
    """
    Auto-generate SCENIC analysis report: parameter summary, main figures, top regulons.
    """
    os.makedirs(out_dir, exist_ok=True)
    params = adata.uns.get("sclucid", {}).get("scenic", {}).get("params", {})
    with open(os.path.join(out_dir, "scenic_params.txt"), "w") as f:
        for k, v in params.items():
            f.write(f"{k}: {v}\n")
    regulons = adata.uns.get("sclucid", {}).get("scenic", {}).get("regulon_names", [])
    pd.Series(regulons[:top_n]).to_csv(os.path.join(out_dir, "top_regulons.csv"))
    # Optionally, export main figures if available
    # ... call analyze_scenic_results() with save_dir=out_dir ...
    log.info(f"SCENIC report exported to: {out_dir}")

# 可选：AI辅助调控因子功能注释（需有AI接口/LLM支持）
def ai_annotate_regulons(regulon_names, ai_model):
    """
    Use AI (e.g., GPT API) to annotate regulons.
    """
    comments = []
    for reg in regulon_names:
        prompt = f"What is the main biological function of the transcription factor '{reg}' in mammals?"
        result = ai_model.ask(prompt)
        comments.append({"regulon": reg, "annotation": result})
    return pd.DataFrame(comments)