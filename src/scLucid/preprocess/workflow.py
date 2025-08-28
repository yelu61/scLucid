import logging
from dataclasses import asdict
from pathlib import Path
from typing import Optional

import scanpy as sc
from anndata import AnnData

from .config import WorkflowConfig
from .hvg import find_hvgs, select_hvg_sets
from .integrate import batch_correction
from .normalize import normalize_data
from .scale import regress_out, scale_data

log = logging.getLogger(__name__)

__all__ = ["run_preprocessing"]


def run_preprocessing(
    adata: AnnData,
    config: Optional[WorkflowConfig] = None,
    results_dir: Optional[str] = None,
    force: bool = False,
) -> AnnData:
    """
    Run the final, fully corrected preprocessing workflow.
    """
    if config is None:
        config = WorkflowConfig()

    log.info("=" * 50)
    log.info("=== Starting Preprocessing Workflow ===")
    log.info("=" * 50)

    adata = adata.copy()
    results_path = Path(results_dir) if results_dir else None
    if results_path:
        results_path.mkdir(parents=True, exist_ok=True)

    # --- 1. Normalization ---
    log.info("Step 1: Normalization")
    adata = normalize_data(
        adata,
        config=config.normalization,
        force=force,
        save_dir=str(results_path / "normalization") if results_path else None,
    )

    # --- 2. Set .raw with normalized data BEFORE regression ---
    log.info(f"Step 2: Storing data from '{config.normalized_layer}' in .raw")
    raw_adata = AnnData(
        X=adata.layers[config.normalized_layer].copy(),
        var=adata.var.copy(),
        obs=adata.obs.copy(),
    )
    adata.raw = raw_adata

    # --- 3. Regression (Optional) ---
    hvg_input_layer = config.normalized_layer
    if config.scaling.vars_to_regress:
        log.info("Step 3: Regressing out covariates")
        adata = regress_out(
            adata,
            config=config.scaling,
            input_layer=config.normalized_layer,
            output_layer=config.regressed_layer,
        )
        hvg_input_layer = config.regressed_layer
    else:
        log.info("Step 3: Skipping regression.")

    # --- 4. HVG Selection ---
    log.info(f"Step 4: HVG Selection (from layer '{hvg_input_layer}')")
    adata = find_hvgs(
        adata,
        config=config.hvg,
        input_layer=hvg_input_layer,
        force=force,
        save_dir=str(results_path / "hvg") if results_path else None,
    )

    # --- 5. Subset to HVGs ---
    log.info("Step 5: Subsetting data to final HVG set")
    hvg_key = adata.uns["sclucid"]["preprocess"]["hvg"]["output_key"]
    adata = select_hvg_sets(
        adata,
        hvg_keys=[hvg_key],
        mode="direct",
        subset=True,
        keep_raw=False,  # .raw is already correctly set
    )

    # --- 6. Scaling ---
    log.info("Step 6: Scaling data (on HVGs)")
    adata = scale_data(adata, config=config.scaling)

    # --- 7. PCA ---
    log.info(f"Step 7: PCA (using {config.graph.n_pcs} components)")
    sc.tl.pca(adata, n_comps=config.graph.n_pcs)
    if results_path:
        try:
            sc.pl.pca_variance_ratio(
                adata, log=True, save="_variance_ratio.png", show=False
            )
            Path("./figures/pca_variance_ratio.png").rename(
                results_path / "pca_variance_ratio.png"
            )
        except Exception:
            log.warning("Could not save PCA variance plot.")

    # --- 8. Integration/Batch Correction ---
    use_rep_downstream = "X_pca"
    if config.integration.method and config.integration.batch_key:
        log.info(f"Step 8: Batch Correction ({config.integration.method})")
        adata = batch_correction(
            adata,
            config=config.integration,
            save_dir=str(results_path / "integration") if results_path else None,
        )
        use_rep_downstream = adata.uns["sclucid"]["preprocess"]["integration"][
            "output_key"
        ]
    else:
        log.info("Step 8: Skipping batch correction.")

    # --- 9. Neighbors & UMAP ---
    log.info("Step 9: Neighbors graph and UMAP")
    sc.pp.neighbors(
        adata,
        n_pcs=config.graph.n_pcs,
        n_neighbors=config.graph.n_neighbors,
        use_rep=use_rep_downstream,
    )
    sc.tl.umap(adata)
    if results_path:
        try:
            color_vars = [
                v
                for v in [config.integration.batch_key, "phase"]
                if v and v in adata.obs.columns
            ]
            if color_vars:
                sc.pl.umap(
                    adata, color=color_vars, save="_final.png", show=False, dpi=300
                )
                Path("./figures/umap_final.png").rename(results_path / "final_umap.png")
        except Exception:
            log.warning("Could not save final UMAP plot.")

    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})[
        "workflow_config"
    ] = asdict(config)

    log.info("=" * 50)
    log.info("=== Preprocessing Workflow Complete! ===")
    log.info("=" * 50)
    return adata
