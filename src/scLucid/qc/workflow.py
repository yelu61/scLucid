"""
High-level QC workflow functions for single-cell RNA-seq data.

This module provides turn-key workflows for standard and advanced
quality control analysis using all components of the package.
"""

import logging
from pathlib import Path
from typing import Optional, Tuple

from anndata import AnnData

from .config import QCWorkflowConfig
from .doublet import predict_doublets
from .filtering import (
    filter_cells,
    generate_qc_report,
    mark_low_quality_cell,
)
from .metrics import calculate_qc_metric

log = logging.getLogger(__name__)


def _setup_workflow(
    adata_in: AnnData, results_dir: str, overwrite: bool
) -> Tuple[AnnData, Path]:
    """
    Prepares the AnnData object and results directory for a workflow.
    """
    results_path = Path(results_dir)
    if results_path.exists() and not overwrite:
        log.warning(
            f"{results_dir} already exists. Old results may be overwritten. "
            "Consider setting overwrite=True in the config."
        )
    results_path.mkdir(parents=True, exist_ok=True)
    adata = adata_in.copy()
    return adata, results_path


def run_standard_qc(
    adata_in: AnnData,
    config: Optional[QCWorkflowConfig] = None,
    overwrite: bool = False,
) -> AnnData:
    """
    Run a standard single-cell RNA-seq QC workflow driven by a configuration object.

    If no config is provided, sensible defaults are used. This workflow includes:
    1. QC metric calculation (with cell cycle scoring).
    2. Doublet detection (Scrublet + basic heuristics).
    3. Low-quality cell marking (standard thresholds).
    4. Filtering of marked cells.
    5. Generation of a final report.

    Args:
        adata_in: Input AnnData object (raw or pre-normalized).
        config: A QCWorkflowConfig object. If None, a default config is created.
        overwrite: If True, overwrite existing results directory.

    Returns:
        Filtered AnnData object after QC.
    """
    if config is None:
        log.info("No QCWorkflowConfig provided, using standard defaults.")
        config = QCWorkflowConfig()
    config.validate()

    adata, results_path = _setup_workflow(adata_in, config.results_dir, overwrite)

    # --- 1. QC metric calculation ---
    adata = calculate_qc_metric(
        adata,
        sample_key=config.sample_key,
        reporting_config=config.metrics_reporting_config,
        calculate_cell_cycle=True,
        cell_cycle_species=config.species,
    )

    # --- 2. Doublet detection ---
    config.doublet_config.save_dir = str(results_path / "doublet")
    adata = predict_doublets(
        adata,
        config=config.doublet_config,
        sample_key=config.sample_key,
    )

    # --- 4. Low-quality cell marking ---
    config.marking_config.save_dir = str(results_path / "low_quality")
    adata = mark_low_quality_cell(
        adata,
        config=config.marking_config,
        sample_key=config.sample_key,
    )

    # --- 5. Filtering ---
    adata_filtered = filter_cells(adata, config=config.filter_config, copy=True)

    # --- 6. Reporting ---
    generate_qc_report(
        adata_filtered,
        save_dir=results_path / "report",
        sample_key=config.sample_key,
        adata_before=adata,
    )

    adata_filtered.uns.setdefault("qc", {})["workflow"] = "standard"
    adata_filtered.uns["qc"]["config_used"] = config.to_dict()
    return adata_filtered


def run_advanced_qc(
    adata_in: AnnData,
    config: QCWorkflowConfig,
    overwrite: bool = False,
) -> AnnData:
    """
    Run an advanced, fully configurable single-cell RNA-seq QC workflow.

    This workflow is entirely controlled by the provided QCWorkflowConfig object,
    allowing fine-grained control over every step.

    Args:
        adata_in: Input AnnData object.
        config: A fully populated QCWorkflowConfig object.
        overwrite: If True, overwrite existing results directory.

    Returns:
        Filtered AnnData object after QC.
    """
    config.validate()
    adata, results_path = _setup_workflow(adata_in, config.results_dir, overwrite)
    log.info(f"QC workflow started, results will be saved in: {results_path}")

    # --- 1. Metrics Calculation ---
    adata = calculate_qc_metric(
        adata,
        sample_key=config.sample_key,
        reporting_config=config.metrics_reporting_config,
        calculate_cell_cycle=True,
        cell_cycle_species=config.species,
        # You can add advanced parameters such as extra_gene_sets in config.metrics_reporting_config if desired
    )

    # --- 2. Doublet Detection ---
    config.doublet_config.save_dir = str(results_path / "doublet")
    # First, analyze lineage coexpression to provide transparency and report
    adata = predict_doublets(
        adata, config=config.doublet_config, sample_key=config.sample_key
    )

    # --- 3. (Optional) Threshold Suggestion (can be toggled in config) ---
    # if getattr(config, "suggest_thresholds", False):
    #     suggested_thresholds = suggest_qc_thresholds(adata, method="mad")
    #     log.info(f"Suggested thresholds: {suggested_thresholds.to_dict()}")

    # --- 4. Marking and Filtering ---
    config.marking_config.save_dir = str(results_path / "low_quality")
    adata = mark_low_quality_cell(
        adata,
        config=config.marking_config,
        sample_key=config.sample_key,
    )

    # The filtering logic is completely controlled by filter_config
    adata_filtered = filter_cells(adata, config=config.filter_config, copy=True)

    # --- 5. Reporting ---
    generate_qc_report(
        adata_filtered,
        save_dir=results_path / "report",
        sample_key=config.sample_key,
        adata_before=adata,
    )

    adata_filtered.uns.setdefault("qc", {})["workflow"] = "advanced"
    adata_filtered.uns["qc"]["config_used"] = config.to_dict()
    log.info("Advanced QC workflow completed successfully.")
    return adata_filtered
