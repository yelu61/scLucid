"""Manual preprocessing step-by-step example.

This teaching example shows how to call low-level scLucid preprocessing APIs
while still producing a workflow-compatible review contract. For a shorter
stage-level template, prefer ``examples/02_simple_api/qc_preprocess_review.py``.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

import scLucid as scl

DATA_PATH = Path("results/qc_filtered.h5ad")  # output from qc_step_by_step.py or qc_preprocess_review.py
OUTPUT_DIR = Path("results/preprocess_step_by_step")
OUTPUT_PATH = OUTPUT_DIR / "preprocessed.h5ad"


def ensure_counts_layer(adata):
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    return adata


def finalize_manual_preprocess_review(adata, config, steps_executed, save_dir: Path):
    """Write the same preprocess review envelope used by workflow examples."""

    base_summary = {
        "manual_workflow": True,
        "manual_workflow_note": (
            "This example uses low-level preprocessing APIs and finalizes the "
            "standard preprocess review contract explicitly."
        ),
    }
    enriched = scl.pp.enrich_preprocessing_review_summary(
        base_summary,
        adata=adata,
        config=config,
        successful_steps=steps_executed,
        tissue_type="unknown",
        keep_intermediate_layers=True,
    )
    review_summary = scl.ut.normalize_review_summary(
        enriched,
        module="preprocess",
        workflow_name="manual_step_by_step",
        adata=adata,
        steps_executed=steps_executed,
        config=config.to_dict(),
        warnings=enriched.get("preprocess_readiness", {}).get("review_reasons", []),
    )
    scl.ut.validate_review_summary_schema(review_summary, module="preprocess", raise_on_error=True)
    scl.pp.validate_preprocessing_review_summary(review_summary, raise_on_error=True)

    scl.ut.save_result(adata, "preprocess", "workflow_config", config.to_dict())
    scl.ut.save_result(adata, "preprocess", "steps_executed", steps_executed)
    scl.ut.save_result(adata, "preprocess", "review_summary", review_summary)
    scl.ut.export_review_summary(
        review_summary,
        save_dir=save_dir,
        module="preprocess",
        title="Manual Preprocessing Review Summary",
        adata=adata,
    )
    return review_summary


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(DATA_PATH)
    adata = ensure_counts_layer(adata)

    config = scl.pp.PreprocessingWorkflowConfig(
        normalization=scl.pp.NormalizationConfig(
            method="standard",
            input_layer="counts",
            output_layer="normalized",
            target_sum=1e4,
            save_dir=str(OUTPUT_DIR / "normalization"),
        ),
        hvg=scl.pp.HVGConfig(
            method="scanpy",
            n_top_genes=2000,
            flavor="auto",
            exclude_gene_types=["mitochondrial", "ribosomal"],
        ),
        scaling=scl.pp.ScalingConfig(max_value=10),
        integration=scl.pp.IntegrationConfig(method=None, batch_key=None),
        graph=scl.pp.GraphConfig(n_pcs=50, n_neighbors=15),
        run_gene_filtering=True,
        run_regression=False,
        run_integration=False,
        save_dir=str(OUTPUT_DIR),
        n_jobs=1,
    )

    steps_executed = []

    print("Filtering low-detection genes...")
    initial_genes = int(adata.n_vars)
    sc.pp.filter_genes(adata, min_cells=config.min_cells_per_gene)
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["gene_filtering"] = {
        "source": "counts" if "counts" in adata.layers else "X",
        "min_cells_per_gene": int(config.min_cells_per_gene),
        "initial_genes": initial_genes,
        "final_genes": int(adata.n_vars),
        "removed_genes": int(initial_genes - adata.n_vars),
        "manual_step": True,
    }
    steps_executed.append("gene_filtering")

    print("Normalizing data...")
    adata = scl.pp.normalize_data(adata, config=config.normalization, force=True)
    steps_executed.append("normalization")

    print("Storing full-gene normalized expression in adata.raw...")
    adata.raw = adata.copy()
    steps_executed.append("set_raw")

    print("Finding highly variable genes...")
    adata = scl.pp.find_hvgs(
        adata,
        config=config.hvg,
        input_layer="normalized",
        force=True,
        save_dir=str(OUTPUT_DIR / "hvg"),
    )
    steps_executed.append("hvg_selection")

    print("Subsetting to selected HVGs...")
    hvg_keys = [k for k in adata.var.columns if k.startswith("highly_variable") and not k.endswith(("_means", "_dispersions", "_dispersions_norm"))]
    if not hvg_keys:
        hvg_keys = ["highly_variable"]
    adata = scl.pp.select_hvg_sets(
        adata,
        hvg_keys=hvg_keys,
        mode="direct",
        subset=True,
        keep_raw=True,
        output_key="highly_variable_selected",
        plot_venn=False,
        save_dir=str(OUTPUT_DIR / "hvg"),
    )
    steps_executed.append("subset_hvg")

    print("Scaling data...")
    if "normalized" in adata.layers:
        adata.X = adata.layers["normalized"].copy()
    adata = scl.pp.scale_data(adata, config=config.scaling, output_layer="scaled")
    steps_executed.append("scaling")

    print("Running PCA...")
    n_pcs = min(config.graph.n_pcs, adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=n_pcs)
    steps_executed.append("pca")

    print("Skipping batch correction in this light manual path.")

    print("Computing neighbors and UMAP...")
    sc.pp.neighbors(
        adata,
        n_pcs=min(config.graph.n_pcs, adata.obsm["X_pca"].shape[1]),
        n_neighbors=config.graph.n_neighbors,
    )
    sc.tl.umap(adata)
    steps_executed.append("neighbors_umap")

    review_summary = finalize_manual_preprocess_review(
        adata,
        config=config,
        steps_executed=steps_executed,
        save_dir=OUTPUT_DIR,
    )
    compact = scl.pp.summarize_preprocess_review_summary(review_summary)

    adata.write_h5ad(OUTPUT_PATH, compression="gzip")
    print("Preprocessing complete")
    print(f"Final shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"Raw shape: {adata.raw.shape if adata.raw is not None else 'missing'}")
    print(f"Readiness: {compact['readiness_status']} ({compact['readiness_score']})")
    print(f"Saved to: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
