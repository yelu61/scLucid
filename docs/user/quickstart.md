# Quick Start

This page shows the recommended minimal path for using scLucid on a new
dataset. For longer runnable scripts, see `examples`. For full narrative
analyses, see `notebooks`.

## Recommended Learning Order

1.  Describe the project with `scLucid.ProjectContext`
2.  Inspect `scLucid.plan_analysis()` before running
3.  Run the supported QC -> preprocessing -> analysis path
4.  Use `scLucid.review_run()` to resolve `BLOCKED` and `REVIEW` items
5.  Drop down to a stage workflow only for decisions that need adjustment
6.  Export reviewer-facing summaries before making biological claims

## Minimal End-To-End Example

``` python
import scanpy as sc
import scLucid as scl

adata = sc.read_h5ad("data/pbmc3k.h5ad")
adata.layers["counts"] = adata.X.copy()

context = scl.ProjectContext(
    dataset_type="pbmc_or_blood",
    species="human",
    sample_key="sample",
    study_objective="broad cell atlas",
)
plan = scl.plan_analysis(adata, context=context)
print(plan.to_frame())

adata = scl.run_pipeline(
    adata,
    plan=plan,
    preprocess_save_dir="results/preprocess",
    show_progress=True,
)

review = scl.review_run(adata)
print(review.to_frame())
print(review.show_next_actions())

adata.write("results/final_annotated.h5ad")
```

Replace `"sample"` with the actual biological-sample column in
`adata.obs`. For treatment or response projects, also provide
`condition_key`, `experimental_unit_key`, and `paired_key` when applicable.

## What This Path Gives You

- QC trace and `qc_reviewer_table` under `adata.uns["sclucid"]["qc"]`
- QC review sidecars when `save_dir` is set
- standard preprocessing outputs such as normalized layers, HVG
  metadata, PCA, and neighbors/UMAP
- preprocessing layer contract and reviewer table describing
  `counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph`
- clustering labels in `adata.obs`
- annotation evidence, analysis output contract, inference policy, and
  analysis reviewer table in `adata.obs` and `adata.uns`
- one cross-stage `run_review` with status, rationale, next action, and
  rerun scope under `adata.uns["sclucid"]["run_review"]`

## How To Read The Outcome

- `BLOCKED`: a structural prerequisite is missing; do not rely on downstream
  interpretation until it is fixed
- `REVIEW`: the workflow may continue as a first pass, but a biological or
  statistical decision still needs confirmation
- `READY`: the recorded contract contains no unresolved blocker for that stage
- `NOT_RUN`: the stage was absent; this is not automatically an error

`READY` means ready for the declared handoff, not proof that a biological claim
is true. See [Reviewing Results](reviewing_results.md) for the full action
contract.

## Light Default, Optional Enhancements

The recommended path is intentionally light by default. It avoids
mandatory R dependencies and avoids aggressive correction steps unless
the data provide a reason to use them.

Default preprocessing:

- filter genes detected in too few cells after QC and before
  normalization
- normalize raw counts with library-size normalization and `log1p`
- select HVGs with `flavor="auto"`, which resolves to dependency-light
  `seurat` on log-normalized inputs
- scale, run PCA, build neighbors, and compute UMAP
- skip regression and batch correction unless explicitly enabled

Optional enhancements:

- `normalization.method="scran"` keeps Scanpy's
  `scanpy.external.pp.scran_normalize` path for users who already have a
  working R/scran environment
- `hvg.flavor="seurat_v3"` can be used on raw-count inputs when
  `scikit-misc` is installed
- `run_integration=True` with Harmony/scVI/scANVI/BBKNN/ComBat should be
  used only after inspecting batch effects and over-correction risk

## Choosing Between Default And Intelligent Preprocessing

Use `PreprocessingWorkflowConfig.default()` when:

- you want the canonical light-dependency package path
- your dataset is standard scRNA-seq with familiar batch structure
- you value stability, signal preservation, and simplicity over
  parameter search

Use `scl.recommendation.run_intelligent_preprocessing()` when:

- you want data-driven parameter suggestions
- you want a reviewer-facing summary before applying recommendations
- you need help choosing HVG / PCA / neighbors / integration settings

## Related Repository Entry Points

- `examples/01_workflow/basic_pipeline.py`: shortest maintained
  workflow-layer script
- `examples/02_simple_api/qc_step_by_step.py`: composable QC inspection
  path
- `examples/02_simple_api/preprocess_step_by_step.py`: composable
  preprocessing path
- `examples/03_advanced_notebooks/`: full notebook analyses with richer
  outputs
- `scripts/run_pbmc_golden_path.py`: real-data acceptance script with
  manifest output

## When To Use Stage-Specific Functions

Use `scl.run_qc()`, `scl.pp.run_preprocessing()`, and
`scl.analysis.run_standard_analysis()` when rerunning a complete stage.
`run_standard_qc()` is the compatibility/step-control QC path; low-level
functions such as `cluster_cells()` and `run_annotation()` belong to the
simple API layer. The unified plan -> run -> review path remains the
recommended first screen.

See [Usage Layers](usage_layers.md) for the full product-layer model and
[QC And Preprocess Maturity](qc_preprocess_maturity.md) for the
QC/preprocessing hardening standard.
