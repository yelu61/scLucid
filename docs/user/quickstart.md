# Quick Start

This page shows the recommended minimal path for using scLucid on a new
dataset. For longer runnable scripts, see `examples`. For full narrative
analyses, see `notebooks`.

## Recommended Learning Order

1.  Describe the project with `scLucid.ProjectContext`
2.  Review QC without mutating the input
3.  Explicitly apply the reviewed QC policy
4.  Review preprocessing for a declared consumer
5.  Explicitly apply the reviewed preprocessing policy
6.  Do not proceed to biological claims while a blocker remains

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
    condition_key="condition",
    input_provenance="filtered_counts",
)

qc_review = scl.recommend_qc_policy(adata, context=context)
print(qc_review.status, qc_review.reason, qc_review.next_action)
qc_result = scl.apply_qc_policy(adata, qc_review.policy)

pp_review = scl.recommend_preprocess_policy(
    qc_result.adata,
    context,
    consumer="exploration",
)
print(pp_review.status, pp_review.reason, pp_review.next_action)
pp_result = scl.apply_preprocess_policy(qc_result.adata, pp_review.policy)

pp_result.adata.write("results/qc_preprocess_reviewed.h5ad")
```

Replace `"sample"` with the actual biological-sample column in
`adata.obs`. For treatment or response projects, also provide
`condition_key`, `experimental_unit_key`, and `paired_key` when applicable.

## What This Path Gives You

- a compact `DecisionCard` before either stage mutates data
- an immutable, fingerprinted `QCPolicy` and `PreprocessPolicy`
- separate counts, full-gene interpretation, discovery, and optional
  integration spaces
- executed `RunEvidence` with status, claim boundary, and limitations

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

## Current Scientific Maturity

QC and Preprocess remain `REVIEW`, not `CORE`. A contract-complete run does not
establish superiority. Analysis and Tumor feature development remains frozen
until the locked blinded/held-out validation and all three real-project UX gates
pass.

## Related Repository Entry Points

- `examples/02_simple_api/qc_preprocess_review.py`: canonical four-action path
- `examples/01_workflow/basic_pipeline.py`: compatibility pipeline
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
simple API layer. These functions are compatibility or teaching paths; the
four-action recommend/apply workflow is the recommended first screen.

See [Usage Layers](usage_layers.md) for the full product-layer model and
[QC And Preprocess Maturity](qc_preprocess_maturity.md) for the
QC/preprocessing hardening standard.
