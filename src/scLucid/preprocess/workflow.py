import os
import logging
from typing import Optional

import scanpy as sc
from anndata import AnnData

from .normalize import normalize_data
from .hvg import find_hvgs, select_hvg_sets, evaluate_hvg_stability, plot_hvg_metrics
from .scale import scale_data
from .integrate import batch_correction
from .config import PreprocessingConfig

log = logging.getLogger(__name__)

__all__ = ["run_preprocessing"]

def run_preprocessing(
    adata: AnnData,
    config: Optional[PreprocessingConfig] = None,
    results_dir: Optional[str] = None,
    force: bool = False,
) -> AnnData:
    """
    Run the full preprocessing workflow: normalization -> HVG -> scaling -> PCA -> integration -> neighbors/UMAP.

    Args:
        adata: AnnData object after QC.
        config: PreprocessingConfig object. If None, use default.
        results_dir: Directory to save plots/results.
        force: Rerun steps even if already done.

    Returns:
        Preprocessed AnnData object ready for downstream analysis.
    """
    if config is None:
        log.info("No PreprocessingConfig provided, using default parameters.")
        config = PreprocessingConfig()

    log.info("="*50)
    log.info("=== Starting Preprocessing Workflow ===")
    log.info("="*50)
    adata = adata.copy()

    # --- 1. Normalization ---
    log.info("Step 1: Normalization")
    os.makedirs(results_dir, exist_ok=True) if results_dir else None
    adata = normalize_data(
        adata,
        method=config.normalization.method,
        layer=config.counts_layer,
        output_layer=config.normalized_layer,
        target_sum=config.normalization.target_sum,
        exclude_highly_expressed=config.normalization.exclude_highly_expressed,
        max_fraction=config.normalization.max_fraction,
        plot_global_distribution=config.normalization.plot_global_distribution,
        save_dir=results_dir,
        force=force,
    )

    # --- 2. HVG selection ---
    log.info("Step 2: Highly Variable Gene (HVG) Selection")
    adata = find_hvgs(
        adata,
        config=config.hvg,
        force=force,
        plot=True,
        report=True,
        max_cells_plot=40000,
    )

    # --- 3. HVG set selection and scaling ---
    log.info("Step 3: HVG set selection and scaling")
    # Optionally allow intersection/union of multiple HVG masks, or just use default
    hvg_key = f"highly_variable_{config.hvg.method}_{config.hvg.flavor}" if config.hvg.method == "scanpy" else f"highly_variable_{config.hvg.method}"
    adata = select_hvg_sets(
        adata,
        hvg_keys=hvg_key,
        mode="direct",
        subset=True,
        output_key="highly_variable_final",
        show_stats=True,
    )
    adata.raw = adata  # Store for marker/DE

    adata = scale_data(
        adata,
        input_layer=config.normalized_layer,
        output_layer=config.scaled_layer,
        max_value=config.scaling.max_value,
        vars_to_regress=config.scaling.vars_to_regress,
        subset_highly_variable=False,  # Already subsetted
        plot=True,
        save_dir=results_dir,
        force=force,
    )

    # --- 4. PCA ---
    log.info("Step 4: PCA")
    sc.tl.pca(adata, n_comps=config.graph.n_pcs, use_highly_variable=False)
    if results_dir:
        sc.pl.pca_variance_ratio(adata, log=True, save=os.path.join(results_dir, "pca_variance_ratio.png"), show=False)

    # --- 5. Integration/Batch Correction ---
    use_rep_downstream = "X_pca"
    if config.integration.method and config.integration.batch_key:
        log.info(f"Step 5: Batch Correction ({config.integration.method})")
        adata = batch_correction(
            adata,
            batch_key=config.integration.batch_key,
            method=config.integration.method,
            use_rep="X_pca",
            plot=True,
            save_dir=os.path.join(results_dir, "integration") if results_dir else None,
            force=force,
            **(config.integration.harmony_params if config.integration.method == "harmony" else {}),
            **(config.integration.scvi_params if config.integration.method == "scvi" else {}),
        )
        use_rep_downstream = f"X_{config.integration.method}"
        log.info(f"Downstream analysis will use the integrated embedding: '{use_rep_downstream}'")
    else:
        log.info("Step 5: Skipping batch correction/integration.")

    # --- 6. Neighbors & UMAP ---
    log.info("Step 6: Neighbors graph and UMAP")
    sc.pp.neighbors(
        adata, n_pcs=config.graph.n_pcs, n_neighbors=config.graph.n_neighbors, use_rep=use_rep_downstream
    )
    sc.tl.umap(adata)

    # --- Save final plots ---
    if results_dir:
        color_vars = [v for v in [config.integration.batch_key, "phase"] if v and v in adata.obs.columns]
        if color_vars:
            sc.pl.umap(adata, color=color_vars, save=os.path.join(results_dir, "umap.png"), show=False, dpi=300)
        sc.pl.pca_variance_ratio(adata, log=True, save=os.path.join(results_dir, "pca_variance_ratio.png"), show=False, dpi=300)

    # --- Store config for traceability ---
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["workflow_config"] = config

    log.info("="*50)
    log.info("=== Preprocessing Workflow Complete! ===")
    log.info("="*50)
    return adata