# Best Practices

This guide defines the recommended division of labor between the major
scLucid entrypoints and the repository artifacts around them.

## Recommended Pipeline Policy

The maintained QC -\> Preprocess -\> Analysis path is light by default
and optionally extensible:

- QC uses reviewer-first iterative filtering, Python-native metrics,
  adaptive recommendations, conservative multi-evidence decisions,
  canonical contamination/stress/doublet fields, and optional
  Scrublet/heuristic/external doublet evidence.
- Preprocessing preserves counts in `adata.layers["counts"]`, filters
  low-detection genes, runs standard log-normalization, uses
  dependency-light HVG selection, then PCA/neighbors/UMAP. Regression
  and batch correction are explicit opt-ins.
- Analysis consumes the unambiguous layers and embeddings from
  preprocessing for clustering review, marker/evidence tables,
  annotation consensus, DE, and proportion summaries. Analysis also
  records post-hoc QC review evidence for doublet-heavy,
  high-mitochondrial, or stress-high clusters without deleting cells
  automatically.

Optional enhancements such as `scanpy.external.pp.scran_normalize`,
`seurat_v3` HVGs, Harmony, scVI/scANVI, BBKNN, SOLO, or DoubletDetection
are available when their dependencies and biological rationale are
present.

Project-level ambient RNA correction is diagnostic-only by default, but
external CellBender/SoupX/DecontX-style outputs can be registered into
the canonical QC schema and counts-layer contract when available. Custom
rpy2 execution branches are not part of the recommended default path.

## Project Decision Contract

For a new dataset, start with:

```python
context = scl.ProjectContext(
    dataset_type="tumor_tissue",
    sample_key="sample",
    condition_key="condition",
    experimental_unit_key="patient",
    paired_key="patient",
    study_objective="paired treatment response",
)
plan = scl.plan_analysis(adata, context=context)
adata = scl.run_pipeline(adata, plan=plan)
review = scl.review_run(adata)
```

The plan is a structural and study-design check, not an automatic claim that
its parameters are biologically optimal. Resolve missing metadata and inspect
all `REVIEW` items before promoting automated annotations, differential
signals, proportion changes, or malignant labels to final results.

Use the stage entrypoints consistently:

- full core workflow: `scl.run_pipeline()`
- reviewer-first QC stage: `scl.run_qc()` or `scl.qc.run_qc()`
- preprocessing stage: `scl.pp.run_preprocessing()`
- analysis stage: `scl.analysis.run_standard_analysis()`
- tumor interpretation stage: `scl.tumor.run_tumor_analysis()`
- cross-stage decision surface: `scl.review_run()`

## Recommended Pipeline Shape

scLucid's maintained path is a light-dependency default with optional
enhancements. The default should be easy to install, reproducible on
ordinary Python single-cell environments, and conservative about
removing biological signal. Heavier methods remain available when the
project has enough evidence to justify them.

Light default:

- QC: metric calculation, adaptive/hierarchical threshold suggestions,
  reviewer-first `qc_decision` filtering, canonical ambient/doublet/cell
  probability fields, and Scrublet/heuristic doublet support
- Preprocessing: count preservation, library-size normalization + log1p,
  low-detection gene filtering, dependency-light HVG selection, scaling,
  PCA, neighbors, and UMAP
- Analysis: clustering review, marker/evidence tables, annotation
  consensus, differential expression, and proportion summaries

Optional enhancements:

- `normalization.method="scran"` through
  `scanpy.external.pp.scran_normalize` when an R/scran environment is
  already available
- `hvg.flavor="seurat_v3"` for raw-count HVG selection when
  `scikit-misc` is installed
- Harmony, Scanorama, scVI/scANVI, BBKNN, or ComBat when batch effects
  are documented and biological groups are not fully confounded with
  batch
- heavier doublet methods such as SOLO or DoubletDetection when their
  optional dependencies are installed and the dataset size justifies
  them

Review evidence, not automatic deletion:

- ambient RNA and empty-droplet diagnostics are heuristic review
  prompts; external CellBender/SoupX/DecontX/EmptyDrops-style evidence
  should be registered into `ambient_evidence_summary`,
  `ambient_layer_contract`, `cell_probability`, and
  `ambient_fraction` when available
- cell-cycle regression diagnostics report associations and group
  imbalance; regression should be enabled only when the
  biology/technical tradeoff is explicit
- curated marker/pathway/tumor genes belong in a labeled sensitivity
  analysis; they are not forced into the default unsupervised discovery set
- integration diagnostics, post-hoc QC cluster review, and doublet
  evidence should guide manual review before removing cells or
  collapsing biology

Removed from the recommended path:

- project-level ambient RNA correction as an automatic default
  preprocess stage
- custom rpy2/Bioconductor execution paths outside Scanpy's optional
  scran bridge

## QC

Recommended default:

- use `recommend_qc_policy(adata, context)` as the read-only first screen
- inspect the `DecisionCard`, candidate disagreement, affected samples/cells,
  and missing evidence before calling `apply_qc_policy()`
- use `run_qc()` and `run_standard_qc()` only for compatibility, explicit
  step control, resume, or legacy-threshold workflows
- keep filtering conservative by requiring multiple independent
  low-quality criteria before removing cells
- review the stored QC trace and summary outputs before finalizing
  thresholds
- use the `DecisionCard` as the first screen; detailed reviewer tables remain
  annex evidence for thresholds, sources, affected cells, and biological risk
- inspect `ambient_evidence_summary`, `doublet_evidence_summary`,
  `post_annotation_qc_review`, and `qc_benchmark_scorecard` before
  claiming benchmark-grade QC readiness

When to override defaults:

- use pooled thresholds only when samples are intentionally treated as
  one shared population
- set explicit thresholds when project constraints require reproducible
  fixed cutoffs
- treat tumor-aware behavior as a cautionary layer, not a substitute for
  human review

## Preprocessing

Recommended default:

- use `recommend_preprocess_policy(adata, context, consumer=...)` and inspect
  its representation and integration decisions before explicit application
- use `apply_preprocess_policy()` for the canonical light-dependency baseline
- keep `run_preprocessing()` and the older recommendation engine as
  compatibility or sensitivity tools, not as the product decision surface
- do not regress out `total_counts` or `pct_counts_mt` by default; enable
  regression only when diagnostics show a technical covariate dominates
  the biological signal
- do not run batch correction by default; enable integration only after
  inspecting batch mixing, sample structure, and the risk of
  over-correction

When to compare preprocessing candidates:

- when batch correction is uncertain
- when HVG / PCA / neighbors settings need reviewable justification
- when you want the exported preprocessing review summary before
  applying a workflow

Recommended stage handoff:

- preserve raw counts in `adata.layers["counts"]`
- filter genes detected in too few cells at the start of preprocessing
  (default `min_cells_per_gene=3`), after cell-level QC and before
  normalization/HVG selection
- store log-normalized full-gene expression in
  `adata.layers["normalized_full"]` and `adata.raw`
- mark unsupervised discovery features without subsetting the persistent
  interpretation space; scaling is a temporary PCA intermediate
- expect the review and RunEvidence to document `counts`, `normalized_full`,
  `discovery_rep`, and optional `integrated_rep` separately
- use discovery features for PCA/neighbors, but use `normalized_full` for
  marker/program interpretation and `counts` for count models
- use integrated embeddings for visualization/clustering only when
  justified; keep unintegrated normalized expression for marker and DE
  interpretation

Analysis defaults:

- treat cell-level marker tests as discovery/exploratory evidence
- use sample-level pseudobulk for condition DE when biological
  replicates exist
- use CLR sample-level tests for cell-type proportions; scLucid closes
  count, proportion, percentage, and sub-composition inputs before CLR
  transformation
- when using `linear_model_logcpm` for pseudobulk DE, the default
  covariance estimator is HC3 robust standard error; set
  `robust_cov_type="nonrobust"` only for compatibility with ordinary OLS
  behavior

## Stage-Wise Uns Compaction

`adata.uns["sclucid"]` is useful during review, but it can grow after
QC, preprocessing, and analysis because it may contain benchmark
summaries, diagnostic tables, differential-expression tables,
integration evidence, and other intermediate objects. After exporting
sidecar reports, compact the object before handing it to the next stage
or writing a portable `.h5ad`:

``` python
# Inspect first.
scl.compact_sclucid_uns(adata, dry_run=True)

# Keep review summaries, configs, artifacts, context, warnings, and contract
# metadata; remove heavier intermediate result payloads.
scl.compact_sclucid_uns(adata)
```

Use `keep_keys=[...]` when a specific top-level module payload is still
needed inside `uns`. For example, keep `"de"` during active marker
exploration, but drop it after exporting DE tables to CSV/Parquet
sidecars.

## Examples Vs Docs Vs Notebooks

The repository should be read in this order:

- **docs**: stable explanation of the recommended path
- **examples**: short runnable scripts for each supported usage pattern
- **notebooks**: complete analyses with intermediate decisions and
  richer plots

Do not use examples as the authoritative source for package policy. The
authoritative recommended path should live in docs, and examples should
implement that policy rather than redefine it.

## Suggested Usage Patterns

Use `examples/` for:

- quick copy-and-adapt scripts
- testing one workflow in isolation
- minimal reporting examples

Use `notebooks/` for:

- full project walkthroughs
- real-data or publication-style analyses
- exploration with rich intermediate figures

Use `docs/` for:

- installation and onboarding
- recommended defaults
- method selection guidance
- public API reference

## Reproducibility

- store raw counts in `adata.layers["counts"]` before preprocessing
- prefer config objects over scattered keyword arguments
- keep workflow outputs under stable `save_dir` locations
- review sidecar outputs before treating automated decisions as final

## Workflow Hardening

When improving scLucid, prefer vertical workflow slices over broad
module-by-module expansion. Start with a small reproducible dataset, run
the full supported path, inspect the review summaries, and polish the
module boundary that fails first.

Recommended validation tiers:

- `data/pbmc3k.h5ad` for the fast normal-tissue baseline
- `data/lin2020.pdac.h5ad` for the first tumor golden path
- `data/schlesinger2020.pdac.h5ad` for tumor generalization
- active project data for final product acceptance and biological
  plausibility

See `workflow_hardening` for the detailed execution plan.

## QC And Preprocessing First

For productization, QC and preprocessing should be treated as the first
two benchmark modules. They are the foundation for annotation,
differential expression, tumor analysis, and visualization.

Recommended order:

- stabilize the shared workflow contract
- harden QC on PBMC and tumor data
- harden preprocessing on the QC output
- then move to analysis and annotation polish

See `qc_preprocess_maturity` for the module-level maturity checklist.
