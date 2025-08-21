# scRNA/qc/workflow.py

import os
from anndata import AnnData
import scanpy as sc
from .metrics import calculate_qc_metric
from .filtering import mark_low_quality_cell, filter_cells, generate_qc_report, QCThresholds, FilterConfig, suggest_qc_thresholds
from .doublet import predict_doublets, DoubletConfig
from ..utils.marker_manager import get_marker_manager

def run_standard_qc(
    adata_in: AnnData,
    sample_key: str = "sampleID",
    results_dir: str = "./qc_results",
    species: str = "human"
) -> AnnData:
    """
    Run a standard single-cell RNA-seq QC workflow with sensible defaults.
    Returns filtered AnnData.
    """
    os.makedirs(results_dir, exist_ok=True)
    adata = adata_in.copy()
    adata = calculate_qc_metric(
        adata,
        sample_key=sample_key,
        save_dir=os.path.join(results_dir, "metrics"),
        calculate_cell_cycle=True,
        cell_cycle_species=species
    )
    adata = predict_doublets(
        adata,
        sample_key=sample_key,
        marker_species=species,
        save_dir=os.path.join(results_dir, "doublet")
    )
    adata = mark_low_quality_cell(
        adata,
        sample_key=sample_key,
        min_genes=200,
        pc_mt=20.0,
        nmads=5.0,
        save_dir=os.path.join(results_dir, "low_quality")
    )
    adata_filtered = filter_cells(adata, copy=True)
    generate_qc_report(
        adata_filtered,
        save_dir=os.path.join(results_dir, "report"),
        sample_key=sample_key,
        adata_before=adata
    )
    adata_filtered.uns.setdefault("qc", {})["workflow"] = "standard"
    return adata_filtered

def run_advanced_qc(
    adata_in: AnnData,
    sample_key: str = "sampleID",
    results_dir: str = "./qc_results",
    species: str = "human",
    tissue: str = "Lung"
) -> AnnData:
    """
    Run an advanced single-cell RNA-seq QC workflow with custom settings.
    Returns filtered AnnData.
    """
    os.makedirs(results_dir, exist_ok=True)
    adata = adata_in.copy()
    custom_sets = {
        "stress": ["HSPA1A", "HSPB1", "FOS", "JUN"],
        "hypoxia": r"^(HIF|EGLN|ADM)"
    }
    adata = calculate_qc_metric(
        adata,
        sample_key=sample_key,
        extra_gene_sets=custom_sets,
        save_dir=os.path.join(results_dir, "metrics"),
        calculate_cell_cycle=True,
        cell_cycle_species=species
    )
    marker_mgr = get_marker_manager(species=species, tissue=tissue)
    marker_mgr.intersect_with(adata)
    doublet_cfg = DoubletConfig(
        method="scrublet",
        merge_strategy="union",
        use_heuristics=True,
        marker_configs=marker_mgr.get_markers_by_level("major"),
        min_lineages_for_doublet=2,
        save_dir=os.path.join(results_dir, "doublet")
    )
    adata = predict_doublets(
        adata,
        config=doublet_cfg,
        sample_key=sample_key
    )
    suggested_thresholds = suggest_qc_thresholds(adata, method="mad")
    qc_thresholds = QCThresholds(
        min_genes=300,
        max_genes=7000,
        pc_mt=15.0,
        nmads=4.0
    )
    adata = mark_low_quality_cell(
        adata,
        sample_key=sample_key,
        thresholds=qc_thresholds,
        plot_outliers=True,
        save_dir=os.path.join(results_dir, "low_quality")
    )
    filter_cfg = FilterConfig(
        combination_logic="custom",
        custom_logic_expr="predicted_doublet | outlier_mt | outlier_min_genes"
    )
    adata_filtered = filter_cells(adata, config=filter_cfg, copy=True)
    generate_qc_report(
        adata_filtered,
        save_dir=os.path.join(results_dir, "report"),
        sample_key=sample_key,
        adata_before=adata
    )
    adata_filtered.uns.setdefault("qc", {})["workflow"] = "advanced"
    return adata_filtered