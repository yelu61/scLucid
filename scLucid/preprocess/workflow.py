

import os
import logging
from dataclasses import dataclass, field
from typing import Optional, List, Literal

import scanpy as sc
from anndata import AnnData

from .normalize import normalize_data
from .hvg import find_highly_variable_genes
from .scale import scale_data
from .integrate import batch_correction

log = logging.getLogger(__name__)

__all__ = ["PreprocessingConfig", "run_preprocessing"]


@dataclass
class PreprocessingConfig:
    """
    Configuration for the standard preprocessing workflow.

    This class encapsulates all parameters needed to run the preprocessing
    pipeline from normalized data to a final UMAP embedding.
    """
    # Normalization
    normalization_method: str = "standard"
    
    # Highly Variable Genes (HVG)
    n_top_genes: int = 2000
    hvg_flavor: str = "seurat_v3"
    
    # Scaling & Regression
    vars_to_regress: Optional[List[str]] = field(default_factory=lambda: ["total_counts", "pct_counts_mt"])
    
    # Dimensionality Reduction
    n_pcs: int = 50
    
    # Integration / Batch Correction
    batch_key: Optional[str] = "sampleID"
    integration_method: Optional[Literal["harmony", "scanorama", "scvi"]] = "harmony"
    
    # Final Graph and Embedding
    n_neighbors: int = 15
    
    def __post_init__(self):
        """Validate configuration."""
        if self.n_top_genes <= 0:
            raise ValueError("n_top_genes must be positive.")
        if self.n_pcs <= 0:
            raise ValueError("n_pcs must be positive.")
        if self.batch_key and not self.integration_method:
            log.warning("`batch_key` is provided, but `integration_method` is None. No batch correction will be performed.")

def run_preprocessing(
    adata: AnnData,
    config: Optional[PreprocessingConfig] = None,
    results_dir: Optional[str] = None,
    force: bool = False,
) -> AnnData:
    """
    Execute a comprehensive preprocessing workflow.

    This high-level function orchestrates the entire preprocessing pipeline,
    from normalization to UMAP embedding, using a flexible configuration.

    Args:
        adata: AnnData object after QC.
        config: A PreprocessingConfig object. If None, default parameters will be used.
        results_dir: Directory to save intermediate plots and results.
        force: Whether to force re-computation of existing steps.

    Returns:
        The preprocessed AnnData object, ready for clustering and downstream analysis.
    """
    if config is None:
        log.info("No PreprocessingConfig provided, using default parameters.")
        config = PreprocessingConfig()

    log.info("\n" + "="*50)
    log.info("=== Starting Preprocessing Workflow ===")
    log.info("="*50)
    
    adata = adata.copy() # Work on a copy to avoid modifying the original object
    
    # --- 1. Normalization ---
    log.info("\n---- Step 1: Normalizing Data ----")
    adata = normalize_data(
        adata,
        method=config.normalization_method,
        layer="counts", # Assumes raw counts are in this layer
        output_layer="log1p_norm",
        force=force
    )

    # --- 2. Highly Variable Gene (HVG) Selection ---
    log.info("\n---- Step 2: Finding Highly Variable Genes ----")
    adata = find_highly_variable_genes(
        adata,
        layer="log1p_norm",
        n_top_genes=config.n_top_genes,
        flavor=config.hvg_flavor,
        batch_key=config.batch_key, # Use batch key for more robust HVG selection
        plot=True,
        save_dir=results_dir,
        force=force
    )
    
    # --- 3. Subset to HVGs and Scale Data ---
    log.info("\n---- Step 3: Subsetting to HVGs and Scaling ----")
    # Store the full normalized data in .raw before subsetting
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable].copy()
    log.info(f"Subsetted data to {adata.n_vars} highly variable genes.")
    
    # Scale the data, optionally regressing out covariates
    adata = scale_data(
        adata,
        layer=None, # Use adata.X now, which contains log1p_norm of HVGs
        output_layer="scaled", # Store scaled data in a new layer
        vars_to_regress=config.vars_to_regress,
        force=force
    )
    # The scaled data is now in adata.X after scale_data finishes, ready for PCA

    # --- 4. Dimensionality Reduction (PCA) ---
    log.info("\n---- Step 4: Principal Component Analysis (PCA) ----")
    sc.tl.pca(adata, n_comps=config.n_pcs, use_highly_variable=False) # Use all genes in the current adata
    
    if results_dir:
        sc.pl.pca_variance_ratio(adata, log=True, save="_pca_variance.png", show=False)
        # Move the saved plot to the correct directory
        if os.path.exists("./figures/pca_variance_ratio_pca_variance.png"):
             os.rename("./figures/pca_variance_ratio_pca_variance.png", os.path.join(results_dir, "pca_variance_ratio.png"))


    # --- 5. Batch Correction / Integration (Optional) ---
    use_rep_downstream = "X_pca" # Default embedding for neighbor calculation
    if config.batch_key and config.integration_method:
        log.info(f"\n---- Step 5: Batch Correction using {config.integration_method.upper()} ----")
        adata = batch_correction(
            adata,
            batch_key=config.batch_key,
            method=config.integration_method,
            use_rep="X_pca",
            plot=True,
            save_dir=os.path.join(results_dir, "integration") if results_dir else None,
            force=force
        )
        use_rep_downstream = f"X_{config.integration_method}"
        log.info(f"Downstream analysis will use the integrated embedding: '{use_rep_downstream}'")
    else:
        log.info("\n---- Step 5: Skipping Batch Correction ----")

    # --- 6. Neighborhood Graph and UMAP Embedding ---
    log.info(f"\n---- Step 6: Computing Neighborhood Graph and UMAP ----")
    log.info(f"Using '{use_rep_downstream}' for graph construction.")
    sc.pp.neighbors(adata, n_pcs=config.n_pcs, n_neighbors=config.n_neighbors, use_rep=use_rep_downstream)
    sc.tl.umap(adata)
    
    # Final visualization
    if results_dir:
        save_path = os.path.join(results_dir, "final_preprocessing_umaps.png")
        color_vars = [v for v in [config.batch_key, 'phase'] if v in adata.obs.columns]
        if color_vars:
            sc.pl.umap(adata, color=color_vars, save="_umaps.png", show=False)
            if os.path.exists("./figures/umap_umaps.png"):
                os.rename("./figures/umap_umaps.png", save_path)
    
    # --- Store final configuration ---
    adata.uns.setdefault('scrnatk', {})
    adata.uns['scrnatk'].setdefault('preprocess', {})
    adata.uns['scrnatk']['preprocess']['workflow_config'] = config.__dict__

    log.info("\n" + "="*50)
    log.info("=== Preprocessing Workflow Complete! ===")
    log.info("="*50)
    
    return adata