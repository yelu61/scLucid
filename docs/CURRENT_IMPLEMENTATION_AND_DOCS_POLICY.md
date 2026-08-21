# Current Implementation And Documentation Policy

**Status**: current documentation governance  
**Purpose**: keep scLucid documentation from mixing current implementation,
active plans, developer audits, and historical design notes at the same level.

## Source-Of-Truth Layers

Use this order when documents appear to disagree:

1. **Code and tests**
   - `src/scLucid/**`
   - `tests/**`
   - These define what is actually implemented.

2. **User-facing API and best-practice docs**
   - `docs/api/*.md`
   - `docs/user/best_practices.md`
   - `docs/user/quickstart.md`
   - These define the maintained user contract.

3. **Current implementation snapshot**
   - this file
   - `docs/README.md`
   - These explain how to read the documentation tree.

4. **Strategic and phase plans**
   - `docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md`
   - `docs/roadmap/*.md`
   - These describe intended direction and acceptance gates, not guaranteed
     implementation.

5. **Developer audits**
   - `docs/dev/*.md`
   - These are useful engineering notes but can become stale after
     implementation rounds. Treat unchecked gap items as prompts to verify
     against code, not as current truth.

6. **Archive**
   - `docs/archive/**`
   - Historical provenance only.

The root `README.md` should stay stable: identity, core differentiators,
quickstart, and documentation links. It should not track module-by-module
implementation status.

## Current Module Feature Map

| Module | Core Differentiator | Current Maintained Surface | Next Focus |
|--------|---------------------|----------------------------|------------|
| Decision layer | Explicit project design and prioritized cross-stage actions | `ProjectContext`, `plan_analysis`, `review_run`, analysis plan and run-review audit panels | AK112/LPJ-like applied-override acceptance feedback |
| QC | Reviewer-first filtering with tumor-aware guardrails | `run_qc`, QC review summary, decision/reviewer tables, doublet and contamination evidence; `run_standard_qc` for compatibility/step control | More real-data ambient/doublet calibration and report polish |
| Preprocess | Explicit layer and inference handoff | `run_preprocessing`, normalization policy, layer contract, preprocess reviewer table | More real-data layer/HVG/batch benchmarks and docs examples |
| Analysis | Conservative interpretation boundary | `run_standard_analysis`, replicate-aware proportion/pseudobulk design audits, fail-closed metadata and contrast contracts, analysis inference policy, output contract, decision summary, reviewer table | Dedicated pseudobulk reviewer summary, richer annotation evidence, and broader real-data acceptance |
| Annotation | Multi-evidence consensus rather than black-box labeling | marker/reference/LLM evidence merge, consensus labels, review table | clearer examples for lineage/subtype/state workflows |
| Tumor | Tumor-context interpretation on top of stable analysis outputs | malignancy/CNV/program/TME helpers and tumor workflow scaffold | consume analysis contracts more tightly and add tumor-stage review summaries |
| Tools / Bulk / Spatial | Support evidence rather than core product sprawl | selected bulk, spatial, deconvolution, R/Python parity helpers | dependency isolation and evidence-contract integration |
| Marker Resources | Source-aware biological knowledge | marker registry, curation contract, quality summaries | provenance completion and mouse/tumor parity |

## Stage Plan

### Stage 0: Documentation Hygiene

- Keep `README.md` stable and high-level.
- Keep implementation details in `docs/api/*.md` and
  `docs/user/best_practices.md`.
- Add a banner to stale `docs/dev` audits when they are superseded by current
  code.
- Prefer updating an existing roadmap or this policy over creating another
  parallel plan.

### Stage 1: Core Contract Stabilization

- Keep QC, preprocess, and analysis review summaries serializable to `.h5ad`.
- Keep reviewer tables as the notebook/API/report-facing surface.
- Keep `adata.uns["sclucid"]` as the shared audit namespace.
- Keep `review_run()` as a derived product read model; module review summaries
  remain the source evidence.

### Stage 2: Module-Specific Hardening

- QC: strengthen ambient correction evidence, doublet calibration, and
  tumor-aware validation.
- Preprocess: benchmark normalization/HVG/layer handoff across PBMC, tumor,
  low-RNA, and multi-sample settings.
- Analysis: build a dedicated pseudobulk DE reviewer summary on the frozen
  sample/condition/experimental-unit contract and validate annotation evidence
  on additional real datasets.

### Stage 3: Tumor Interpretation And Evidence Bridges

- Use stable analysis outputs as inputs to tumor interpretation.
- Keep malignancy, CNV, TME, therapy, bulk, and spatial outputs labelled with
  confidence, evidence source, and inference level.

## Cleanup Rules

- Do not delete historical notes unless they are duplicated generated clutter.
  Move or label them as archive instead.
- Do not let `docs/dev` gap audits serve as user-facing truth after a gap is
  implemented; update the audit status or link to the current API docs.
- Do not add new phase documents when an existing phase can be updated.
- Every new core-module feature should update:
  - the relevant API `.md`
  - `docs/user/best_practices.md` if behavior changes user guidance
  - this policy if it changes the module feature map or stage plan
  - tests that validate the stored review contract

## Immediate Cleanup Queue

- Review `docs/dev/QC_MODULE_BENCHMARK_GAP_AUDIT.md` and mark implemented
  items as superseded by current QC review/decision code.
- Review `docs/dev/CORE_WORKFLOW_COVERAGE_MATRIX.md` after each public API
  reorganization.
- Keep `docs/user/qc_preprocess_maturity.md` as the benchmark-maturity page,
  but avoid duplicating every implementation detail already covered by API docs.
- Keep `docs/validation/qc_preprocess_evidence_pilot.md` only as a historical
  proxy-evidence log; `validation/evidence_run_index.json` and
  `validation_outputs/current/` define the active evidence surface.
