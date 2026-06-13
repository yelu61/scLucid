# API Triage Staging Area

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
| **Symbol** | `run_standard_qc`, `run_advanced_qc`, `run_preprocessing`, `run_annotation`, `run_standard_analysis`, `run_tumor_analysis`, etc. |
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

---

## Resolved Triage Items

*(None yet — populate as decisions are made.)*
