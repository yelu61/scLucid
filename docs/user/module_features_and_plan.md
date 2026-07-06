# Module Features And Stage Plan

This page is the current bridge between the stable root `README.md` and
the more detailed design, roadmap, and developer-audit documents under
`docs/`. It should describe the current shape of scLucid without turning
the root README into a frequently changing implementation log.

## Source Of Truth Layers

Use these layers in order when documents disagree:

1.  Code, tests, and executed notebooks define what is implemented.
2.  `docs/api/*.md` and maintained workflow guides define the
    public user-facing contract.
3.  This page and `docs/CURRENT_IMPLEMENTATION_AND_DOCS_POLICY.md`
    summarize the current module design and documentation policy.
4.  `docs/roadmap/` contains phase plans and submission-oriented
    execution goals.
5.  `docs/dev/` and `docs/archive/` preserve audit history and
    provenance; they are useful context, but not the final current
    contract.

## Core Workflow Modules

| Module | Distinctive role | Maintained entrypoints | Review surface | Current focus |
|----|----|----|----|----|
| QC | Reviewer-first filtering with sample-aware, tumor-aware, ambient/stress, and benchmark scorecard guardrails | `run_qc`, `run_iterative_qc`, `run_standard_qc`, `run_qc_threshold_decision` | `qc_reviewer_table`, `qc_decision`, `qc_reason`, `qc_confidence`, `ambient_evidence_summary`, `doublet_evidence_summary`, `post_annotation_qc_review`, `qc_benchmark_scorecard` | Real raw-matrix ambient evidence and larger homotypic/solid-tissue doublet validation |
| Preprocess | Explicit layer contract from counts to graph | `run_preprocessing`, `run_embedding_pipeline` | `normalization_decision_policy`, `preprocess_layer_contract`, `preprocess_decision_summary`, `preprocess_reviewer_table` | Larger multi-sample validation and clearer batch/HVG evidence |
| Analysis | Conservative inference boundary for clustering, markers, annotation, and DE | `run_standard_analysis`, `run_pseudobulk_de` | `analysis_inference_policy`, `analysis_output_contract`, `analysis_decision_summary`, `analysis_reviewer_table` | Pseudobulk DE review summary and real-data annotation acceptance |
| Annotation | Multi-evidence consensus instead of single black-box labels | `run_annotation`, `run_annotation_evidence`, `build_annotation_consensus` | annotation review tables and consensus summaries | Clearer lineage, subtype, state, and uncertainty examples |
| Tumor | Tumor-context interpretation over QC/preprocess/analysis evidence | `run_tumor_analysis`, `run_malignancy_interpretation` | malignancy, CNV, TME, therapy, and tumor-stage review summaries | Tighter consumption of analysis contracts and tumor case-study validation |
| Marker resources | Source-aware marker knowledge infrastructure | marker resource managers and registry loaders | quality summaries, curation status, and nomenclature rules | Provenance completion and mouse/tumor parity |
| Tools / Evidence modules | Supporting bulk, spatial, proportion, and external evidence | `scLucid.tools.bulk`, `scLucid.tools.spatial`, proportion and DE helpers | method-specific evidence tables | Dependency isolation and selective method validation |

## Workflow Contract Spine

The current core workflow should remain readable as a chain of
contracts:

`QC decision -> preprocess layer contract -> analysis inference policy -> tumor/annotation evidence`

The practical layer contract for preprocessing is:

`counts -> normalized -> raw -> HVG -> scaled -> PCA -> graph`

Each stage should write a compact reviewer-facing summary that explains:

- what was recommended
- what was actually applied
- why the decision was made
- confidence and review-required status
- affected cells, genes, clusters, or comparisons where applicable
- remaining biological or technical risk

## Stage Plan

Stage 0: Documentation hygiene  
Keep root README stable; make MkDocs pages under `docs/` the user-facing current
contract; mark roadmap and developer audit documents as plan/provenance
when they are not the current implementation truth.

Stage 1: Core contract stabilization  
Keep QC, preprocess, and analysis review summaries aligned around stable
schema fields. New outputs should extend reviewer tables instead of
creating isolated `adata.uns` payloads.

Stage 2: Validation hardening  
Expand PBMC and tumor real-data validation for QC thresholds, doublets,
ambient contamination, layer contracts, HVG/PCA/batch handling,
clustering, annotation, and pseudobulk DE.

Stage 3: Tumor interpretation integration  
Make tumor modules consume stable analysis outputs directly, preserve
uncertainty, and validate malignancy/TME/program calls in public tumor
case studies.

Stage 4: Evidence bridges and release assets  
Connect bulk, spatial, clinical, marker resources, and manuscript
figures as support evidence around the main tumor single-cell workflow.

## What To Update When Code Changes

- Public entrypoint or output schema changes: update `docs/api/*.md` and the
  relevant workflow guide.
- Module positioning or maturity changes: update this page and
  `docs/CURRENT_IMPLEMENTATION_AND_DOCS_POLICY.md`.
- Submission or validation phase changes: update `docs/roadmap/`.
- API uncertainty, deprecation, or compatibility decisions: update
  `docs/dev/API_TRIAGE.md` or the relevant developer audit.
- Historical drafts should move to `docs/archive/` once they are
  superseded and no longer useful as active planning documents.
