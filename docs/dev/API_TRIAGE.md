# API Triage Staging Area

> Developer note: this is a triage log for uncertain, deprecated, or
> private-but-exposed symbols. It is not the current public API contract. Verify
> active user-facing behavior against code, tests, and `docs/api/*.md`.

This document tracks symbols whose public/internal status is uncertain. Each
entry follows the template below. Decisions are made by maintainers and
recorded here; resolved entries are moved to the **Resolved** section.

**Template for new entries:**

| Field | Value |
|-------|-------|
| **Symbol** | `name` |
| **Location** | `module.path` |
| **Current Status** | public / internal / deprecated / uncertain |
| **Concern** | Why is this borderline? |
| **Proposed Action** | promote / demote / keep / remove |
| **Decision Deadline** | YYYY-MM-DD |
| **Decision** | (filled when resolved) |
| **Decision Date** | (filled when resolved) |

---

## P0.5 Surface Classification Snapshot

**Audited:** 2026-08-14
**Evidence:** current `__all__` declarations, `scripts/audit_public_api.py`, API
docs, import/contract tests, and deprecation warnings. The generated inventory
contains 1,099 exported symbols; its syntactic "stable" label means exported
and unflagged, not that every optional scientific method is release-stable.

| Surface | Stable contract | Experimental / optional | Deprecated / compatibility | Internal or unresolved |
|---|---|---|---|---|
| Top level (`scLucid`) | Context/plan/run-review objects, I/O helpers, `run_pipeline`, and canonical QC/preprocess/analysis entrypoints | Optional submodule aliases when their dependencies are unavailable | None removed in P0.5 | Dynamic placeholder-to-import pattern remains T009 |
| `qc` | `run_qc`, `run_standard_qc`, QC configs, decisions, traces, review summaries | External correction backends and optional R integrations | Compatibility config aliases | Backend probes and workflow implementation helpers |
| `preprocess` | `run_preprocessing`, `run_iterative_preprocessing`, layer contracts, graph/integration diagnostics | Specialized adaptive/quality-aware normalization helpers | `results_dir`, renamed configuration fields, and method aliases | Workflow implementation and backend helpers |
| `analysis` | Clustering, annotation, scoring, proportion, pseudobulk-first and review contracts | Optional proportion backends and support-evidence integrations | Analysis-level malignancy wrapper and `analysis.bulk` shim | Private workflow helpers; malignancy belongs to `tumor` |
| `tumor` | Tumor workflow, malignancy/CNV/TME review surfaces | Therapy/evolution utilities and external-data-dependent helpers | Deprecated aliases recorded by the generated inventory | `AnalysisStep` adapters (T002) and TCGA placeholder (T011) |
| `recommendation` | Public schema and result contract | Recommendation engine remains evidence-expanding; treat outputs as advice requiring review | No removal in P0.5 | Engine implementation modules |
| `tools` | Importable optional support-evidence namespaces | Bulk, spatial, pyMonocle3, deconvolution and communication methods remain optional support evidence | Bulk legacy aliases and `analysis.bulk` compatibility path | Backend-specific implementation modules; unsupported limma branch (T012) |

P0.5 removes no callable API. Symbols that might have external callers remain
importable and are triaged for a versioned deprecation cycle instead of being
deleted from source.

---

## Active Triage Items

### T001: `run_malignancy_interpretation` in `scLucid.analysis`

| Field | Value |
|-------|-------|
| **Symbol** | `run_malignancy_interpretation` |
| **Location** | `scLucid.analysis` |
| **Current Status** | deprecated (public, emits `FutureWarning`) |
| **Concern** | Re-exported from `scLucid.tumor.malignancy` with a deprecation wrapper. Exists for backward compatibility since v0.1. |
| **Proposed Action** | remove in v0.3.0 |
| **Decision Deadline** | 2026-09-01 |
| **Decision** | — |
| **Decision Date** | — |

### T002: Tumor `AnalysisStep` adapters in `scLucid.tumor.steps`

| Field | Value |
|-------|-------|
| **Symbol** | `CNVInferenceStep`, `MalignancyScoringStep`, `MalignancyInterpretationStep`, `TMEDeconvolutionStep`, `TherapyPredictionStep` |
| **Location** | `scLucid.tumor.steps` |
| **Current Status** | public |
| **Concern** | Internal adapters to make tumor functions work with the `AnalysisStep` interface. They are in `__all__` but are not documented. Users may confuse them with user-facing classes. |
| **Proposed Action** | demote to internal (remove from `__all__`, keep importable from submodule) or document clearly |
| **Decision Deadline** | 2026-08-01 |
| **Decision** | — |
| **Decision Date** | — |

### T003: Trace/contract symbols in public `__all__`

| Field | Value |
|-------|-------|
| **Symbol** | `QC_TRACE_SCHEMA_VERSION`, `ANALYSIS_REQUIRED_REVIEW_SECTIONS`, `PREPROCESS_STABLE_ENTRYPOINTS`, `TUMOR_TRACE_SCHEMA_VERSION`, etc. |
| **Location** | `scLucid.qc`, `scLucid.analysis`, `scLucid.preprocess`, `scLucid.tumor` |
| **Current Status** | public |
| **Concern** | Schema constants and validation helpers for the review-summary system. Useful for advanced users but clutter the public API. ~24 symbols across subpackages. |
| **Proposed Action** | keep public but consider grouping under a `trace` sub-namespace (e.g., `scLucid.qc.trace.QC_TRACE_SCHEMA_VERSION`) |
| **Decision Deadline** | 2026-08-15 |
| **Decision** | — |
| **Decision Date** | — |

### T004: `_get_cancer_markers` in `scLucid.utils`

| Field | Value |
|-------|-------|
| **Symbol** | `_get_cancer_markers` |
| **Location** | `scLucid.utils` |
| **Current Status** | private-but-exposed (starts with `_`, in `__all__`) |
| **Concern** | Only symbol in `utils.__all__` with `_` prefix. Likely slipped in by mistake. |
| **Proposed Action** | demote (remove from `__all__`, keep importable from `scLucid.utils.manager`) |
| **Decision Deadline** | 2026-07-01 |
| **Decision** | — |
| **Decision Date** | — |

### T005: `analysis.bulk` shim module

| Field | Value |
|-------|-------|
| **Symbol** | `diagnose_bulk_data_quality`, `BulkDiagnosticsConfig`, `run_bulk_de`, etc. |
| **Location** | `scLucid.analysis.bulk` |
| **Current Status** | public via `_export(..., optional=True)` fallback |
| **Concern** | Canonical implementation moved to `scLucid.tools.bulk`. The `analysis.bulk` submodule is a backward-compatible shim. |
| **Proposed Action** | deprecate and remove in v0.3.0; direct users to `scLucid.tools.bulk` |
| **Decision Deadline** | 2026-09-01 |
| **Decision** | — |
| **Decision Date** | — |

### T006: `tools.bulk` legacy aliases

| Field | Value |
|-------|-------|
| **Symbol** | `differential_abundance`, `run_deconvolution` |
| **Location** | `scLucid.tools.bulk` |
| **Current Status** | public (imported from `_legacy.py`) |
| **Concern** | Aliases to canonical functions (`run_bulk_abundance_test`, `deconvolve_bulk`). Kept for backward compatibility. |
| **Proposed Action** | deprecate in v0.2.0, remove in v0.3.0 |
| **Decision Deadline** | 2026-07-15 |
| **Decision** | — |
| **Decision Date** | — |

### T007: `plotting.main` module

| Field | Value |
|-------|-------|
| **Symbol** | All functions in `plotting.main` |
| **Location** | `scLucid.plotting.main` |
| **Current Status** | internal (module docstring says DEPRECATED) |
| **Concern** | Module still exists but is not in `plotting.__all__`. It may contain stale code. |
| **Proposed Action** | remove module in v0.3.0 |
| **Decision Deadline** | 2026-08-01 |
| **Decision** | — |
| **Decision Date** | — |

### T008: `recommendation` subpackage visibility

| Field | Value |
|-------|-------|
| **Symbol** | `recommend_analysis_parameters`, `RecommendationEngine` |
| **Location** | `scLucid.recommendation` |
| **Current Status** | public |
| **Concern** | Small subpackage (6 symbols). Maturity relative to core workflow is uncertain. Currently aliased as `rc` at the top level. |
| **Proposed Action** | keep public; monitor usage |
| **Decision Deadline** | 2026-09-01 |
| **Decision** | — |
| **Decision Date** | — |

### T009: Top-level workflow placeholders

| Field | Value |
|-------|-------|
| **Symbol** | `run_standard_qc`, `run_preprocessing`, `run_annotation`, `run_standard_analysis`, `run_tumor_analysis`, etc. |
| **Location** | `scLucid` (top-level) |
| **Current Status** | public (set to `None`, then reassigned via `_import_optional`) |
| **Concern** | If optional dependencies fail, these are `None` at runtime. Users get confusing errors when calling them. |
| **Proposed Action** | keep current pattern but add `__getattr__` lazy loading with helpful error messages in a future refactor |
| **Decision Deadline** | 2026-10-01 |
| **Decision** | — |
| **Decision Date** | — |

### T010: `utils` mega-subpackage

| Field | Value |
|-------|-------|
| **Symbol** | ~120 symbols in `scLucid.utils` |
| **Location** | `scLucid.utils` |
| **Current Status** | public |
| **Concern** | `utils` is the largest subpackage by symbol count. It mixes contracts, validation, storage, profiling, marker management, I/O, and workflow utilities. |
| **Proposed Action** | keep for now; consider exposing thematic sub-namespaces (e.g., `utils.contracts`, `utils.storage`) as public subpackages in v0.3.0 |
| **Decision Deadline** | 2026-12-01 |
| **Decision** | — |
| **Decision Date** | — |

### T011: `query_tcga_data` placeholder

| Field | Value |
|-------|-------|
| **Symbol** | `query_tcga_data` |
| **Location** | `scLucid.tumor.utils` |
| **Current Status** | experimental public placeholder |
| **Concern** | Returns a one-row explanatory DataFrame without querying TCGA. Treating it as evidence would violate the source/provenance contract. |
| **Proposed Action** | keep behavior unchanged in P0.5; decide whether to replace it with a connector-backed implementation or deprecate the export |
| **Decision Deadline** | maintainer decision before v0.2.0 |
| **Decision** | — |
| **Decision Date** | — |

### T012: bulk `limma` method option

| Field | Value |
|-------|-------|
| **Symbol** | `BulkDEConfig.method="limma"` |
| **Location** | `scLucid.tools.bulk` |
| **Current Status** | explicit unsupported option (`NotImplementedError`) |
| **Concern** | The config surface accepts a method that is not implemented; older design notes can be read as implying parity. |
| **Proposed Action** | keep the explicit error in P0.5; either implement with validation evidence in a later support-evidence milestone or remove through a documented breaking-change cycle |
| **Decision Deadline** | maintainer decision before v0.2.0 |
| **Decision** | — |
| **Decision Date** | — |

### T013: CNV reference-signature comparison placeholder

| Field | Value |
|-------|-------|
| **Symbol** | `CNVSigExtractor.compare_to_reference` |
| **Location** | `scLucid.tumor.cnv` |
| **Current Status** | experimental public method with placeholder comparison values |
| **Concern** | Returns `best_match="unknown"` and `correlation=0.0` for every reference. It must not be presented as a completed biological comparison. |
| **Proposed Action** | keep unchanged in P0.5 because removal changes public behavior; either implement a validated comparison or deprecate the method with a versioned migration |
| **Decision Deadline** | maintainer decision before tumor validation claims expand |
| **Decision** | — |
| **Decision Date** | — |

---

## Resolved Triage Items

No resolved items yet. When a triage item is closed, move it here with the
resolution date and commit reference.
