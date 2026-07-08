# QC / Preprocess API Cleanup Audit

> Developer note: this audit records cleanup decisions from a specific
> hardening round. It should guide compatibility maintenance, but current
> user-facing contracts live in `docs/api/*.md`,
> `docs/user/module_features_and_plan.md`, and code/tests.

This audit records low-value wrappers and boundary risks found while building
real-data validation benchmarks. The goal is to keep scLucid's public surface
small without breaking older notebooks abruptly.

## Current Decision

Do not delete compatibility aliases immediately. Keep them importable, hide them
from `__all__`, and prefer canonical entrypoints in docs, tests, and workflow
code.

## Canonical Entrypoints

| Area | Canonical API | Notes |
|---|---|---|
| QC workflow | `scLucid.qc.run_qc` | Recommended default reviewer-first QC workflow. |
| QC advanced workflow | `scLucid.qc.run_iterative_qc` | Explicit reviewer-first advanced entrypoint. |
| QC compatibility workflow | `scLucid.qc.run_standard_qc` | Compatibility, step control, and resume entrypoint. |
| QC threshold recommendation | `scLucid.qc.recommend_qc_thresholds` | Structured recommendation bundle with candidate evidence. |
| QC threshold decision | `scLucid.qc.run_qc_threshold_decision` | Canonical recommendation-to-marking chain. |
| QC reporting | `scLucid.qc.generate_qc_report` from `qc.reporting` | Reporting is no longer exported from filtering modules. |
| Preprocess workflow | `scLucid.preprocess.run_preprocessing` | Main maintained preprocessing workflow. |
| Preprocess iterative workflow | `scLucid.preprocess.run_iterative_preprocessing` | Reviewer-first real-project workflow for HVG audit, diagnostic embedding, integration decision, and final graph construction. |
| Preprocess embedding | `scLucid.preprocess.run_embedding_pipeline` | Replaces `run_embedding_workflow`. |
| Preprocess normalization | `scLucid.preprocess.normalize_data` | High-level workflow entry; low-level adaptive helpers live in `scLucid.preprocess.adaptive_normalize`. |
| Preprocess recommendation | `scLucid.recommendation.preprocess` / `scLucid.recommendation.run_intelligent_preprocessing` | Preprocess-specific parameter recommendation now belongs to the recommendation layer, not `scl.pp`. |

## Legacy Aliases Removed From Top-Level API

| Symbol | Current status | Recommended action |
|---|---|---|
| `scLucid.preprocess.run_embedding_workflow` | Removed | Use `scLucid.preprocess.run_embedding_pipeline`. |
| `scLucid.preprocess.adaptive_normalize` | Removed from `scl.pp`; still available from `scLucid.preprocess.adaptive_normalize` | Use `normalize_data()` or `run_iterative_preprocessing()` for workflows; import low-level helpers explicitly when needed. |
| `scLucid.preprocess.quality_aware_normalize` | Removed from `scl.pp`; still available from `scLucid.preprocess.adaptive_normalize` | Import explicitly from the low-level module in advanced notebooks. |
| `apply_gene_biotype_strategy`, `get_gene_biotype_cache_dir`, `list_gene_biotype_resources` | Removed from `scl.pp`; still available from `scLucid.preprocess.gene_biotype` | Keep as explicit low-level utilities. |
| `scLucid.preprocess.intelligent` and `scl.pp.run_intelligent_preprocessing` | Moved to `scLucid.recommendation.preprocess` and `scl.recommendation.run_intelligent_preprocessing` | Keep preprocess execution APIs separate from recommendation/advisor APIs. |

## Boundary Fixes Already Applied

| Issue | Resolution |
|---|---|
| QC report generation lived beside threshold helpers | Canonical implementation moved to `qc.reporting`; filtering no longer exports report helpers. |
| QC decision table lacked reviewer impact fields | Added `affected_cells`, `affected_fraction`, `review_required`, and `risk_note` in `qc.trace`. |
| QC benchmark evidence was scattered across threshold, tumor, doublet, and ambient tables | Added `validation/qc/build_qc_evidence_package.py` to generate unified QC source data and `qc_claim_scorecard.tsv`. |
| Doublet benchmark evidence stayed outside normal reports | Added `doublet_evidence_summary.benchmark_decision` and surfaced recommended default mode, primary method, candidate `algorithm_weight`, and review status in QC reports. |
| Preprocess layer contract was nested and hard to scan | Added `layer_transition_table` with row-wise layer/slot/`.X`/`.raw` semantics. |
| HDF5 sanitizer converted review table rows into dicts | Added `layer_transition_table` to review sequence restoration keys. |
| `run_advanced_qc` duplicated `run_standard_qc` without independent semantics | Removed the wrapper and the top-level alias; use `run_qc`, `run_iterative_qc`, or `run_standard_qc` according to workflow depth. |
| Preprocess compatibility aliases blurred public API boundaries | Removed `run_embedding_workflow` and top-level adaptive/gene-biotype aliases; tests now assert they are absent from `scl.pp`. |
| Preprocess intelligent recommendation lived under the execution namespace | Moved the package to `recommendation.preprocess`; `scl.pp` no longer exports recommender classes or `run_intelligent_preprocessing`. |

## Cleanup Rules Going Forward

- New report/rendering functions should live in `qc.reporting` or the relevant
  reporting module, not in filtering/threshold modules.
- New threshold/filter helpers should not write reports directly.
- New preprocessing validation should write review-summary fields rather than
  ad hoc `uns` payloads when the evidence is part of the public contract.
- New benchmark scripts should either write panel-specific source tables or be
  consumed by the QC/preprocess package builders; avoid creating one-off
  evidence files with no claim-level scorecard entry.
- Legacy threshold entrypoints should be removed from package-level exports once
  a canonical replacement is available; avoid keeping multiple public ways to
  perform the same QC threshold step.

## Next Cleanup Candidates

| Candidate | Why not delete now | Next step |
|---|---|---|
| Low-level adaptive normalization helpers | Useful as algorithm-level APIs. | Keep in `scLucid.preprocess.adaptive_normalize`; do not re-export from `scl.pp`. |
