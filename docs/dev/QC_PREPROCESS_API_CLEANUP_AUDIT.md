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
| QC workflow | `scLucid.qc.run_standard_qc` | Main maintained QC workflow. |
| QC threshold decision | `scLucid.qc.run_qc_threshold_decision` | Replaces `run_qc_decision_workflow`. |
| QC reporting | `scLucid.qc.generate_qc_report` from `qc.reporting` | `qc.filtering.generate_qc_report` is compatibility only. |
| Preprocess workflow | `scLucid.preprocess.run_preprocessing` | Main maintained preprocessing workflow. |
| Preprocess embedding | `scLucid.preprocess.run_embedding_pipeline` | Replaces `run_embedding_workflow`. |
| Preprocess normalization | `scLucid.preprocess.normalize_data` | High-level workflow entry; low-level adaptive helpers remain hidden. |

## Compatibility Aliases To Keep Hidden

| Symbol | Current status | Recommended action |
|---|---|---|
| `scLucid.qc.run_advanced_qc` | Importable compatibility wrapper, omitted from `qc.__all__` | Keep for one deprecation cycle; docs should use `run_standard_qc`. |
| `scLucid.qc.run_qc_decision_workflow` | Importable compatibility wrapper, omitted from `qc.__all__` | Keep for one deprecation cycle; docs should use `run_qc_threshold_decision`. |
| `scLucid.qc.filtering.generate_qc_report` | Thin compatibility wrapper to `qc.reporting.generate_qc_report` | Keep importable, but do not call from workflow internals. |
| `scLucid.preprocess.run_embedding_workflow` | Importable compatibility alias, omitted from `preprocess.__all__` | Keep for one deprecation cycle; docs should use `run_embedding_pipeline`. |
| `scLucid.preprocess.adaptive_normalize` | Importable hidden low-level helper | Keep as low-level algorithm API; do not present as canonical workflow. |
| `scLucid.preprocess.quality_aware_normalize` | Importable hidden low-level helper | Keep as low-level algorithm API; docs should clarify it is not the canonical workflow. |
| `apply_gene_biotype_strategy`, `get_gene_biotype_cache_dir`, `list_gene_biotype_resources` | Importable hidden utilities | Keep hidden unless a documented workflow needs them. |

## Boundary Fixes Already Applied

| Issue | Resolution |
|---|---|
| QC report generation lived in `filtering/suggestions.py` | Canonical implementation moved to `qc.reporting`; filtering wrapper retained only for compatibility. |
| QC decision table lacked reviewer impact fields | Added `affected_cells`, `affected_fraction`, `review_required`, and `risk_note` in `qc.trace`. |
| QC benchmark evidence was scattered across threshold, tumor, doublet, and ambient tables | Added `validation/qc/build_figure2_qc_evidence_package.py` to generate unified Figure 2 source data and `qc_claim_scorecard.tsv`. |
| Doublet benchmark evidence stayed outside normal reports | Added `doublet_evidence_summary.benchmark_decision` and surfaced recommended default mode, primary method, candidate `algorithm_weight`, and review status in QC reports. |
| Preprocess layer contract was nested and hard to scan | Added `layer_transition_table` with row-wise layer/slot/`.X`/`.raw` semantics. |
| HDF5 sanitizer converted review table rows into dicts | Added `layer_transition_table` to review sequence restoration keys. |

## Cleanup Rules Going Forward

- New report/rendering functions should live in `qc.reporting` or the relevant
  reporting module, not in filtering/threshold modules.
- New threshold/filter helpers should not write reports directly.
- New preprocessing validation should write review-summary fields rather than
  ad hoc `uns` payloads when the evidence is part of the public contract.
- New benchmark scripts should either write panel-specific source tables or be
  consumed by the Figure 2/Figure 3 package builders; avoid creating one-off
  evidence files with no claim-level scorecard entry.
- Compatibility aliases should be covered by tests that assert they are
  importable but absent from `__all__`.
- A compatibility alias can be deleted only after docs, notebooks, and tests no
  longer reference it and a deprecation note has existed for at least one
  release cycle.

## Next Cleanup Candidates

| Candidate | Why not delete now | Next step |
|---|---|---|
| `run_advanced_qc` | Older notebooks may still import it. | Keep warning; verify notebooks no longer use it before removal. |
| `run_qc_decision_workflow` | Older threshold-decision examples may use it. | Keep warning; migrate examples to `run_qc_threshold_decision`. |
| `run_embedding_workflow` | Older preprocess notebooks may use it. | Keep warning; migrate docs/examples to `run_embedding_pipeline`. |
| Hidden adaptive normalization helpers | Useful as algorithm-level APIs. | Document as low-level opt-in helpers, not workflow entrypoints. |
