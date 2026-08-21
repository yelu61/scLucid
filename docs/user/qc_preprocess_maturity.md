# QC And Preprocessing Maturity Plan

> **Current release boundary (2026-08-20):** QC and Preprocess are `REVIEW`,
> not `CORE`. A review bundle may be contract-complete while scientific
> superiority remains unvalidated. Analysis and Tumor feature development is
> frozen until `validation/qc_preprocess/acceptance_contract.json` and the three
> real-project usability gates pass.

The P0 dataset portfolio is an additional fail-closed gate. Its accessions,
download/license boundaries, metadata requirements, and preregistered metrics
are maintained in [Validation Dataset Portfolio](validation_dataset_portfolio.md).
Availability is not validation: every required dataset-specific endpoint must
have an intact `RunEvidence` artifact before QC or Preprocess can become `CORE`.

This page started as a QC/preprocess maturity plan and now also records
nearby module maturity context. For the compact current module map and
documentation source-of-truth hierarchy, see `module_features_and_plan`.

QC and preprocessing should become the first benchmark modules for
scLucid. They define the user's trust in the rest of the workflow: if
filtering, normalization, feature selection, and graph construction are
not explainable, the later annotation and tumor modules cannot be
trusted.

## Current Maturity Assessment

The table below is a snapshot of where each module stands today. It
should be updated as modules move through the maturity gates described
later in this document.

| Area | Current Level | What Works Now | Main Gaps |
|----|----|----|----|
| QC | **REVIEW; scientific core withheld** | Read-only `DecisionCard`, sample-level fail-closed gate, four candidate families, separate evidence heads, fingerprinted policy execution | Independent blinded labels, real miQC/SampleQC adapters, passing grouped-bootstrap superiority gate, all three real-project UX runs |
| Preprocessing | **REVIEW; scientific core withheld** | Four-space representation contract, simple baseline executor, controlled mixology held-out regret gate, integration Pareto review | External tumor-project method ranking, leave-one-project-out validation, all three real-project UX runs |
| Analysis | **Feature development frozen** | Existing compatibility workflows remain callable | Do not expand until QC and Preprocess pass locked gates |
| Marker Resources | Strong architectural direction | Unified `Manager`, human/mouse registry resources, tissue/tumor marker views, artifact/program/tumor routing, curation SOP | Source provenance at scale, mouse tissue/tumor parity, atlas-derived marker review |
| Tumor Module | **Feature development frozen** | Existing CNV, malignancy, TME, therapy, and heterogeneity code remains callable for compatibility | Resume only after QC and Preprocess pass locked gates |
| Plotting | Useful foundation | Publication-style themes and domain plots | Top-journal figure templates, richer multi-panel reports, visual regression checks |
| Tools / Evidence Modules | Expanding tumor support | Python-facing wrappers, bulk deconvolution, bulk/spatial clean-room utilities, R parity scaffolds | Selective method validation, dependency isolation, bulk/spatial tumor use cases |
| Documentation / Examples | Good skeleton | Three usage layers, advanced notebooks, golden-path scripts | Keep docs synchronized with maturity contracts and real-data acceptance results |

## Target State

The target is not simply "more functions". A mature QC or preprocessing
module should satisfy five standards:

- stable public entrypoints for workflow, simple API, and advanced usage
- explicit AnnData input and output contracts
- review summaries that explain recommendations, applied parameters, and
  user overrides
- real-data golden-path validation on PBMC and tumor datasets
- clear docs and examples that match the maintained code path

The maintained code path is a light-dependency default with optional
enhancements. Defaults should not require R, scVI-tools, Scanorama,
DoubletDetection, or other heavy extras. Those methods remain supported
as explicit opt-in extensions when the dataset and environment justify
them.

## QC As The First Benchmark Module

QC should answer four user questions:

1.  What thresholds were recommended?
2.  Why were they recommended?
3.  What was actually applied?
4.  What biological or technical risk remains after filtering?
5.  What exactly is safe to hand off to preprocessing, and which cells
    must remain in review or sensitivity records?

Required QC outputs:

- `adata.obs` metrics such as `n_genes_by_counts`, `total_counts`, and
  `pct_counts_mt` when available
- low-quality and doublet flags
- `adata.uns["sclucid"]["qc"]["workflow_config"]`
- `adata.uns["sclucid"]["qc"]["review_summary"]`
- `review_summary["policy_flow"]` describing profile -\> threshold
  proposal -\> biological-risk scoring -\> policy choice -\> reviewer
  table -\> optional apply
- `review_summary["doublet_evidence_summary"]` with prediction rates,
  score ranges, doublet risk metadata, and external-evidence notes when
  present
- `review_summary["ambient_evidence_summary"]` with ambient risk,
  correction status, counts-layer contract, `cell_probability`, and
  `ambient_fraction` availability
- `review_summary["post_annotation_qc_review"]` with retained
  stress/ambient/high-MT/doublet signals stratified by annotation and
  sample when labels are available
- `review_summary["qc_benchmark_scorecard"]` summarizing threshold,
  doublet, ambient, retention, and post-annotation sensitivity evidence
- `review_summary["qc_handoff_readiness"]` declaring the recommended
  preprocessing counts layer, retained decision columns, review/sensitivity
  cell fractions, downstream safe-to-continue flags, and blockers
- optional report sidecars under the configured `save_dir`

QC reporting boundary:

- `scLucid.qc.generate_qc_report` is implemented in the reporting layer
  and exported from `scLucid.qc.reporting`.
- Filtering code stays focused on cell filtering. Threshold recommendation
  lives in `scLucid.qc.policy.thresholds`, and report generation lives in
  `scLucid.qc.reporting`.

QC entries that can usually be compacted after review sidecars are
exported:

- detailed benchmark payloads duplicated inside `review_summary`
- doublet diagnostic parameter/evidence payloads when exported tables
  exist
- raw recommendation objects when applied/original configs and review
  summaries are retained

QC hardening tasks:

- keep `recommend_qc_policy` as the canonical read-only first screen and
  `apply_qc_policy` as the explicit execution entrypoint
- keep `run_qc`, `run_standard_qc`, and `recommend_intelligent_qc` callable
  for compatibility, step control, or sensitivity analysis, without presenting
  them as the scientific decision surface
- keep ambient RNA correction diagnostic-only by default, but support
  explicit external CellBender/SoupX/DecontX-style result registration
  into canonical obs columns and the ambient layer contract
- make user overrides explicit in the review summary
- make the threshold `decision_table` reviewer-readable with
  recommended, applied, source, confidence, evidence, review-required
  status, affected cells, biological guardrails, strategy rank,
  recommended-policy status, and risk note
- keep doublet benchmark evidence connected to normal QC reports through
  `doublet_evidence_summary`, rather than requiring users to inspect raw
  doublet parameter payloads
- keep `qc_handoff_readiness` as the QC-to-preprocess contract so real
  projects can distinguish removed cells, review-required cells,
  sensitivity-only cells, and tumor-fragile states before preprocessing
- treat the former QC evidence package as exploratory work under
  `validation_outputs/work/qc_evidence_package/`; only explicitly bound
  endpoint evidence under `validation_outputs/current/` can affect maturity
- test tumor-aware behavior on PDAC data where high mitochondrial
  content may be a warning rather than an automatic removal criterion
- test edge cases: missing mitochondrial genes, single-sample data,
  small cell counts, sparse matrices, and absent `sampleID`
- maintain real-data evidence runners under `validation/qc/` for
  threshold comparison, tumor biological fidelity, doublet evidence, and
  ambient/empty droplet contracts

Current QC evidence status:

- Auditability and execution contracts are implemented, but contract
  completeness is not scientific validation.
- Lin 2020 is a development failure-control dataset, not a high-quality lockbox:
  `GSM4679533` has a catastrophic joint pattern and must be blocked.
- Existing tumor marker/program retention results are proxy evidence only and
  cannot establish correct cell-level QC labels or superiority.
- Kang 2018 supports a limited doublet calibration question, mainly for
  genotype-detectable multiplets; doublets remain a separate review head.
- Public mixology supplies cross-protocol cell-line identity truth for
  preprocessing fidelity. Its released objects were already sample-QC filtered,
  so they are not used as low-quality-cell truth.
- Ambient RNA remains `NOT_EVALUABLE` without appropriate source evidence and
  is never converted into an automatic cell-deletion pass.

## Preprocessing As The Second Benchmark Module

Preprocessing should answer four user questions:

1.  Which expression layer was used at each step?
2.  Why were HVG, PCA, neighbors, and batch-correction parameters
    chosen?
3.  Was biological signal protected from over-correction?
4.  Is the output ready for clustering and annotation?

Required preprocessing outputs:

- authoritative `adata.layers["counts"]`
- full-gene `adata.layers["normalized_full"]` and `adata.raw`
- `adata.var["discovery_feature"]` or the recorded discovery-feature key
- `adata.obsm["X_pca"]`
- neighbors graph and `adata.obsm["X_umap"]` when graph steps run
- `adata.uns["sclucid"]["preprocess"]["workflow_config"]`
- `adata.uns["sclucid"]["preprocess"]["review_summary"]`

The preprocessing review summary should include both a compact
`layer_transition_summary` and a row-wise `layer_transition_table` so a
reviewer can inspect each step's input layer, output slot, `adata.X`
semantics, `adata.raw` semantics, review-required status, and risk note.

Preprocessing entries that can usually be compacted after the handoff to
analysis:

- normalization/HVG/scaling/integration diagnostic dictionaries once
  their summaries are in `review_summary`
- large integration evaluation payloads and temporary
  intelligent-preprocess recommendation details after sidecar export
- redundant layer statistics when counts, normalized layer, PCA,
  neighbors, and review summary are retained

Analysis entries that can usually be compacted before long-term storage:

- full DE, enrichment, proportion, and annotation-review tables after
  exporting them to sidecar files
- intermediate clustering-resolution grids after the chosen cluster
  labels and review summary are retained

Preprocessing hardening tasks:

- keep `recommend_preprocess_policy` and `apply_preprocess_policy` as the
  canonical review/apply entrypoints; retain `run_preprocessing` for compatibility
- make layer transitions explicit in the review summary
- make HVG selection evidence inspectable
- keep regression and batch correction opt-in by default, and document
  when to enable them
- document when to skip regression and why persistent HVG subsetting is not
  part of the four-space policy path
- keep `scanpy.external.pp.scran_normalize` as the only supported
  optional scran path; avoid custom rpy2/Bioconductor execution branches
- warn when tumor data are batch-corrected in a way that may remove
  malignant, clone, patient, or microenvironment signal
- test small datasets where PCA components and neighbors must be clipped
  safely
- maintain real-data evidence runners under
  `validation/preprocess/` for layer contracts, HVG
  marker/program preservation, batch-correction diagnostics, and graph
  handoff stability

## Recommended Implementation Order

1.  Freeze the minimal shared contract. Keep `adata.uns["sclucid"]`,
    review summary envelopes, and canonical keys stable before
    refactoring internals.
2.  Stabilize QC. Run PBMC and PDAC slices, then fix every unclear
    threshold, warning, or report field discovered by those runs.
3.  Stabilize preprocessing. Use the QC output as input, then harden
    normalization, HVG, PCA, neighbors, and optional integration until
    the handoff to analysis is predictable.
4.  Only after the locked QC and preprocessing gates pass, unfreeze Analysis
    and Tumor feature development.
5.  Build the combined QC/preprocess evidence package from the
    stabilized review summaries, then run it inside active real projects
    as an acceptance loop. Every mismatch between reviewer output and
    analyst expectation should become a contract, default, or document
    improvement.

## Definition Of Done For A Module

A module can be treated as contract-complete when all of the following are
true. It becomes scientific `CORE` only after the additional locked acceptance
contract and real-project usability gates pass:

- the workflow entrypoint passes lightweight tests
- at least one real-data golden path exercises the module
- public examples cover workflow and simple API use
- docs explain default behavior, override behavior, and review artifacts
- review summaries contain enough evidence for another analyst to audit
  the run
- output can be serialized to `.h5ad` after compacting heavy artifacts
- known limitations are documented rather than hidden

## Suggested Test Tiers

- smoke tests: imports, public API availability, config construction
- unit tests: threshold logic, layer selection, HVG behavior, edge cases
- lightweight integration: synthetic or small AnnData through the full
  module
- golden path: PBMC and PDAC subsets with saved artifacts
- project acceptance: active real-world datasets with biological
  plausibility review
