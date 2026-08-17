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
| Decision layer | Project design and cross-stage action contract | `ProjectContext`, `plan_analysis`, `review_run` | `analysis_plan`, `run_review`, READY/REVIEW/BLOCKED action table | Real-project acceptance and applied-override feedback |
| QC | Reviewer-first filtering with sample-aware, tumor-aware, ambient/stress, and benchmark scorecard guardrails | `run_qc`, `run_iterative_qc`, `run_standard_qc`, `run_qc_threshold_decision` | `qc_reviewer_table`, `qc_decision`, `qc_reason`, `qc_confidence`, `ambient_evidence_summary`, `doublet_evidence_summary`, `post_annotation_qc_review`, `qc_benchmark_scorecard` | Real raw-matrix ambient evidence and larger homotypic/solid-tissue doublet validation |
| Preprocess | Explicit layer contract from counts to graph | `run_preprocessing`, `run_embedding_pipeline` | `normalization_decision_policy`, `preprocess_layer_contract`, `preprocess_decision_summary`, `preprocess_reviewer_table` | Larger multi-sample validation and clearer batch/HVG evidence |
| Analysis | Conservative inference boundary for clustering, markers, annotation, and DE | `run_standard_analysis`, `run_pseudobulk_de` | `analysis_inference_policy`, `analysis_output_contract`, `analysis_decision_summary`, `analysis_reviewer_table` | Pseudobulk DE review summary and real-data annotation acceptance |
| Annotation | Multi-evidence consensus instead of single black-box labels | `run_annotation`, `run_annotation_evidence`, `build_annotation_consensus` | annotation review tables and consensus summaries | Clearer lineage, subtype, state, and uncertainty examples |
| Tumor | Tumor-context interpretation over QC/preprocess/analysis evidence | `run_tumor_analysis`, `run_malignancy_interpretation` | malignancy, CNV, TME, therapy, and tumor-stage review summaries | Tighter consumption of analysis contracts and tumor case-study validation |
| Marker resources | Source-aware marker knowledge infrastructure | marker resource managers and registry loaders | quality summaries, curation status, and nomenclature rules | Provenance completion and mouse/tumor parity |
| Tools / Evidence modules | Supporting bulk, spatial, proportion, and external evidence | `scLucid.tools.bulk`, `scLucid.tools.spatial`, proportion and DE helpers | method-specific evidence tables | Dependency isolation and selective method validation |

## Differentiation Spine

scLucid's differentiator is not a larger checklist of single-cell methods. The
project should keep turning routine analysis steps into context-aware scientific
decisions:

| Layer | Current contract | Near-term deepening | Long-term direction |
|----|----|----|----|
| Biological context | Dataset context, tumor-aware QC guardrails, marker/tumor resources | Make tissue, species, cancer type, sample structure, and analysis goal more visible in reviewer tables | Analysis-intent ontology for rare-cell discovery, TME interpretation, trajectory, therapy response, and other project goals |
| QC interpretation | Conservative reviewer-first filtering, ambient/doublet/stress summaries, post-annotation QC review | Report biological impact of competing QC policies, especially tumor purity and immune/TME composition changes | Quality-state attribution and sensitivity analysis for cells that could be low-quality, stressed, hypoxic, doublet-like, or biologically meaningful |
| Inference semantics | `claim_level`, `inference_level`, review summaries, layer contracts | Harmonize evidence, claim, and inference terms across QC, preprocess, analysis, annotation, tumor, bulk, and spatial outputs | Unified evidence ontology that downstream reports and agent interfaces can consume without overclaiming |
| Tumor interpretation | Malignancy, CNV, TME, therapy, and program-oriented scaffolds | Make tumor modules consume stable analysis contracts and write tumor-stage review summaries | Sample-level ecosystem and ecotype-style interpretation with stability, limitations, and support evidence |
| Product surface | `ProjectContext`, analysis plans, unified run review, HTML/report artifacts, reviewer tables, notebooks, examples | Real-project decision cards, applied overrides, and Methods-ready summaries | Interactive review and natural-language interfaces over mature evidence bundles |

Planned items in the right columns should remain roadmap claims until they are
backed by code, tests, executed notebooks, or validation outputs.

## Workflow Contract Spine

The current core workflow should remain readable as a chain of
contracts:

`ProjectContext -> analysis plan -> QC decision -> preprocess layer contract -> analysis inference policy -> tumor/annotation evidence -> run review`

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
creating isolated `adata.uns` payloads. Use `review_run()` as the common
product-facing read model rather than duplicating module evidence.

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
