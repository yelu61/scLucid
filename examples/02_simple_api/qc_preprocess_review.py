"""Simple API example: make QC and preprocessing decisions reviewable.

This script focuses on the two modules that should become scLucid's benchmark
modules. It uses stage-level APIs instead of the unified pipeline so the user
can inspect each stage's review summary before moving on.
"""

from pathlib import Path

import scanpy as sc

import scLucid as scl
from scLucid.preprocess import WorkflowConfig
from scLucid.qc import QCWorkflowConfig


DATA_PATH = Path("data/pbmc3k.h5ad")
OUTPUT_DIR = Path("results/examples/qc_preprocess_review")


def print_review_header(name: str, summary: dict) -> None:
    """Print the stable review-summary fields shared by scLucid modules."""
    print(f"\n{name}")
    print("-" * len(name))
    print(f"Module: {summary.get('module')}")
    print(f"Workflow: {summary.get('workflow_name')}")
    print(f"Steps: {summary.get('steps_executed')}")
    print(f"Warnings: {len(summary.get('warnings', []))}")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(DATA_PATH)
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()

    if "sampleID" not in adata.obs.columns:
        adata.obs["sampleID"] = "pbmc3k"

    qc_config = QCWorkflowConfig(
        save_dir=str(OUTPUT_DIR / "qc"),
        species="human",
        tissue_type="normal_tissue",
        use_recommendations=True,
        threshold_mode="hierarchical",
        use_parallel=False,
        n_jobs=1,
    )

    adata = scl.qc.run_standard_qc(
        adata,
        config=qc_config,
        tissue_type="normal_tissue",
        show_progress=True,
    )
    qc_summary = adata.uns["sclucid"]["qc"]["review_summary"]
    print_review_header("QC review summary", qc_summary)
    qc_validation = scl.qc.validate_qc_module_completeness(adata)
    qc_compact = scl.qc.summarize_qc_review_summary(qc_summary)
    print(f"QC module valid: {qc_validation['valid']}")
    print(f"QC maturity: {qc_compact['maturity_status']}")
    print(f"QC readiness: {qc_compact['readiness_status']} ({qc_compact['readiness_score']})")
    print(f"QC retained cells: {qc_compact['final_cells']}")

    preprocess_config = WorkflowConfig.quick(
        n_top_genes=1000,
        run_regression=False,
        run_integration=False,
        save_dir=str(OUTPUT_DIR / "preprocess"),
        n_jobs=1,
    )
    preprocess_config.graph.n_pcs = 30
    preprocess_config.graph.n_neighbors = 15

    adata = scl.pp.run_preprocessing(
        adata,
        config=preprocess_config,
        tissue_type="normal_tissue",
        show_progress=True,
    )
    preprocess_summary = adata.uns["sclucid"]["preprocess"]["review_summary"]
    print_review_header("Preprocessing review summary", preprocess_summary)
    pp_validation = scl.pp.validate_preprocess_module_completeness(adata)
    pp_compact = scl.pp.summarize_preprocess_review_summary(preprocess_summary)
    print(f"Preprocess module valid: {pp_validation['valid']}")
    print(f"Preprocess maturity: {pp_compact['maturity_status']}")
    print(
        "Preprocess readiness: "
        f"{pp_compact['readiness_status']} ({pp_compact['readiness_score']})"
    )
    print(f"QC input available: {pp_compact['qc_input_available']}")
    print(f"HVG selected: {pp_compact['n_hvg_selected']}")
    print(f"PCA components: {pp_compact['actual_n_pcs']}")
    print(f"UMAP computed: {pp_compact['umap_computed']}")
    print(f"Step status counts: {pp_compact['step_status_counts']}")
    print(f"Review-required steps: {pp_compact['review_required_steps']}")

    adata.write(OUTPUT_DIR / "pbmc3k_qc_preprocess_result.h5ad")
    print(f"\nSaved reviewed QC/preprocessing result to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
