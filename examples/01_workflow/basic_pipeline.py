"""Workflow layer example: run the supported scLucid baseline pipeline.

This is the recommended first-pass script for a new dataset. It uses the
highest-level workflow API and relies on scLucid to record QC, preprocessing,
analysis, contracts, and review summaries under ``adata.uns["sclucid"]``.
"""

from __future__ import annotations

from pathlib import Path

import scanpy as sc

import scLucid as scl

DATA_PATH = Path("data/pbmc3k.h5ad")
OUTPUT_DIR = Path("results/examples/workflow_basic")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(DATA_PATH)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    if "sampleID" not in adata.obs.columns:
        adata.obs["sampleID"] = "pbmc3k"

    qc_config = scl.qc.QCWorkflowConfig(
        save_dir=str(OUTPUT_DIR / "qc"),
        species="human",
        tissue_type="normal_tissue",
    )
    preprocess_config = scl.pp.PreprocessingWorkflowConfig(
        save_dir=str(OUTPUT_DIR / "preprocess"),
    )
    analysis_config = scl.al.AnalysisWorkflowConfig(
        save_dir=str(OUTPUT_DIR / "analysis"),
    )

    adata = scl.run_pipeline(
        adata,
        stages=["qc", "preprocess", "analysis"],
        dataset_type="pbmc_or_blood",
        species="human",
        qc_config=qc_config,
        preprocess_config=preprocess_config,
        analysis_config=analysis_config,
        show_progress=True,
    )

    adata.write_h5ad(OUTPUT_DIR / "pbmc3k_workflow_result.h5ad")

    qc_summary = adata.uns["sclucid"]["qc"]["review_summary"]
    preprocess_summary = adata.uns["sclucid"]["preprocess"]["review_summary"]
    analysis_summary = adata.uns["sclucid"]["analysis"]["review_summary"]

    print("Workflow complete")
    print(f"Final shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"QC steps: {qc_summary['steps_executed']}")
    print(f"Preprocessing steps: {preprocess_summary['steps_executed']}")
    print(f"Analysis steps: {analysis_summary['steps_executed']}")
    analysis_compact = scl.al.summarize_analysis_review_summary(analysis_summary)
    print(f"Analysis readiness: {analysis_compact['readiness_status']}")
    print(f"Post-hoc QC review required: {analysis_compact['posthoc_qc_review_required']}")
    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    if "hierarchical_annotation_plan" in annotation_ns:
        print("Hierarchical annotation plan recorded under analysis.annotation")
    if "subset_annotation_refinement_reconciliation" in annotation_ns:
        print("Subset annotation reconciliation is available for manual review")
    if analysis_compact.get("malignancy_enabled"):
        print(
            "Suspect/malignant fraction: "
            f"{analysis_compact.get('suspect_or_malignant_fraction')}"
        )
    print(f"Saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
