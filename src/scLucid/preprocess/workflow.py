"""
High-level workflow for orchestrating the scLucid preprocessing pipeline.
"""

import logging
from pathlib import Path  # NEW: Use pathlib for robust path handling
from typing import Optional

import scanpy as sc
from anndata import AnnData

from .config import PreprocessingWorkflowConfig  # CHANGED: Import main config
from .integrate import batch_correction
from .hvg import find_hvgs, select_hvg_sets
from .normalize import normalize_data
from .scale import scale_data

log = logging.getLogger(__name__)

__all__ = ["run_preprocessing"]


def run_preprocessing(
    adata: AnnData,
    config: Optional[PreprocessingWorkflowConfig] = None,
    results_dir: Optional[str] = None,
    force: bool = False,
) -> AnnData:
    """
    Run the full preprocessing workflow using a centralized configuration object.

    This orchestrated pipeline executes the following steps:
    1. Normalization
    2. Highly Variable Gene (HVG) selection
    3. Data scaling and regression
    4. Principal Component Analysis (PCA)
    5. Batch correction / Integration (optional)
    6. Neighborhood graph construction and UMAP embedding

    Args:
        adata: AnnData object after the QC stage.
        config: A PreprocessingWorkflowConfig object to control all steps.
                If None, default parameters will be used.
        results_dir: Directory to save plots and reports. If None, plots are not saved.
        force: If True, rerun steps even if their outputs already exist in the AnnData object.

    Returns:
        A new, fully preprocessed AnnData object ready for downstream analysis.
    """
    if config is None:
        log.info("No PreprocessingWorkflowConfig provided, using default parameters.")
        config = PreprocessingWorkflowConfig()

    # --- 0. Setup and Initialization ---
    log.info("=" * 50)
    log.info("=== Starting Preprocessing Workflow ===")
    log.info("=" * 50)
    
    # Work on a copy to prevent modifying the original object
    adata = adata.copy()

    # Setup results directory using pathlib
    if results_dir:
        results_path = Path(results_dir)
        results_path.mkdir(parents=True, exist_ok=True)
        log.info(f"Results will be saved to: {results_path.resolve()}")
    else:
        results_path = None

    # --- 1. Normalization ---
    log.info("Step 1: Normalization")
    # CHANGED: Pass the entire normalization config object directly
    adata = normalize_data(
        adata,
        config=config.normalization,
        input_layer=config.counts_layer,
        output_layer=config.normalized_layer,
        force=force,
    )

    # --- 2. Highly Variable Gene (HVG) Selection ---
    log.info("Step 2: Highly Variable Gene (HVG) Selection")
    # CHANGED: Pass the entire HVG config object directly
    adata = find_hvgs(
        adata, 
        config=config.hvg, 
        force=force
    )

    # Store a copy of the full data in .raw before subsetting to HVGs
    adata.raw = adata.copy()
    log.info("Stored full dataset in `adata.raw`.")

    # Subset the AnnData object to only the selected HVGs for downstream steps
    hvg_mask = adata.var[adata.uns['sclucid']['preprocess']['hvg']['output_key']]
    adata = adata[:, hvg_mask].copy()
    log.info(f"Subsetted data to {adata.n_vars} highly variable genes.")

    # --- 3. Scaling ---
    log.info("Step 3: Data Scaling and Regression")
    # CHANGED: Pass the entire scaling config object directly
    adata = scale_data(
        adata,
        config=config.scaling,
        input_layer=config.normalized_layer,
        output_layer=config.scaled_layer,
        force=force,
    )

    # --- 4. Principal Component Analysis (PCA) ---
    log.info(f"Step 4: PCA (using {config.graph.n_pcs} components)")
    # Use the scaled data layer for PCA
    sc.tl.pca(adata, n_comps=config.graph.n_pcs, use_highly_variable=False)
    
    # CHANGED: Save plot directly to the correct path
    if results_path:
        pca_plot_path = results_path / "pca_variance_ratio.png"
        sc.pl.pca_variance_ratio(adata, log=True, save=str(pca_plot_path), show=False)
        # Scanpy's save adds a prefix, so we rename for a clean filename
        default_save_path = Path(f"./figures/{pca_plot_path.name}")
        if default_save_path.exists():
            default_save_path.rename(pca_plot_path)

    # --- 5. Integration / Batch Correction ---
    downstream_rep = "X_pca"  # Default representation for downstream analysis
    if config.integration.method and config.integration.batch_key:
        log.info(f"Step 5: Batch Correction using '{config.integration.method}'")
        # CHANGED: Pass the entire integration config object
        adata = batch_correction(
            adata,
            config=config.integration,
            use_rep="X_pca",  # Integration is typically done on PCA
            plot=True,
            save_dir=str(results_path / "integration") if results_path else None,
            force=force,
        )
        downstream_rep = f"X_{config.integration.method}"
        log.info(f"Downstream analysis will use integrated embedding: '{downstream_rep}'")
    else:
        log.info("Step 5: Skipping batch correction (method or batch_key not specified).")

    # --- 6. Neighborhood Graph and UMAP ---
    log.info("Step 6: Computing neighborhood graph and UMAP embedding")
    sc.pp.neighbors(
        adata,
        n_pcs=config.graph.n_pcs,
        n_neighbors=config.graph.n_neighbors,
        use_rep=downstream_rep,
    )
    sc.tl.umap(adata)

    # --- 7. Final Visualizations ---
    if results_path:
        log.info("Generating final UMAP plots...")
        # Identify valid columns to color by
        color_vars = [
            var for var in [config.integration.batch_key, "phase"] 
            if var and var in adata.obs.columns
        ]
        if color_vars:
            # CHANGED: Save UMAP plot directly
            umap_plot_path = results_path / "final_umap.png"
            sc.pl.umap(
                adata,
                color=color_vars,
                save=str(umap_plot_path),
                show=False,
                dpi=300
            )
            default_umap_path = Path(f"./figures/{umap_plot_path.name}")
            if default_umap_path.exists():
                default_umap_path.rename(umap_plot_path)

    # --- 8. Finalization ---
    # Store the complete configuration used for this run for perfect reproducibility
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["workflow_config"] = config.to_dict()

    log.info("=" * 50)
    log.info("=== Preprocessing Workflow Complete! ===")
    log.info(f"Final AnnData object shape: {adata.n_obs} cells × {adata.n_vars} genes")
    log.info("=" * 50)
    
    return adata