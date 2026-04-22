# FIX Summary — scLucid Review-Fix Iteration

## Overview

This document summarizes the modifications from the review-fix iteration driven by `REVIEW.md` and `PLAN.md`. The focus is on making QC recommendations trustworthy, hierarchical thresholds biologically valid, tests intention-revealing, and documentation factually consistent with the repository.

---

## 1. Historical Context: Tumor Workflow Phase A+B (not restored)

Phase A+B tumor workflow changes were partially implemented earlier but the full test suite (`tests/tumor/`) was never restored to disk. The current repository **does not** contain:

- `tests/tumor/` directory or its tests
- `run_tumor_analysis` / `run_tumor_analysis_expert` exports from `src/scLucid/tumor/__init__.py`

The repository **does** contain `src/scLucid/tumor/config.py` and `src/scLucid/tumor/workflow.py`, but the public exports and test coverage are still incomplete. Because the public surface and tests are missing, verification commands referencing a restored `tests/tumor/` suite will still fail. The tumor-related code that **does** exist and is functional is limited to:

- `src/scLucid/tumor/` submodules (`cnv`, `malignancy`, `microenvironment`, etc.) imported as domain libraries
- Tumor-aware QC behavior in `src/scLucid/qc/workflow.py` (`tissue_type` hint, MT filtering relaxation, flag storage)

---

## 2. QC Module Fixes

### 2.1 Recommendation Application Semantics

**Problem:** QC recommendations were applied unconditionally, overwriting caller-specified config values. `model_copy()` was not deep, so the caller's nested config objects could be mutated as a side effect. The trace could not reliably distinguish recommendation-applied values from genuine user overrides.

**Fixes in `src/scLucid/qc/workflow.py`:**
- `_apply_qc_recommendations()` now returns `(applied_config, original_config_snapshot)` and uses `model_copy(deep=True)` to guarantee isolation from the caller's object.
- Recommendations are applied **only** to fields that were **not explicitly set** by the user, checked via `model_fields_set` at each nested level (`QCWorkflowConfig` → `MarkingConfig` → `QCThresholds`).
- Explicit user-provided thresholds (e.g., `min_genes=999`) are preserved even when `use_recommendations=True`.
- `_run_qc_workflow()` now returns the modified `applied_config` so that downstream steps (`filter_cells`, trace storage) actually use the post-recommendation, post-tumor-adjustment configuration.
- The unified trace gained a new `original_config` key, and `user_overrides` now diffs the **original user config** against the recommendation (not the applied config against itself).

**New tests:**
- `tests/qc/test_qc_recommendation_executable.py::test_caller_config_not_mutated`
- `tests/qc/test_qc_recommendation_executable.py::test_explicit_user_thresholds_survive_recommendations`

### 2.2 Domain-Aware Bounds for Hierarchical Thresholds

**Problem:** `AdaptiveThresholdCalculator._suggest_adaptive_thresholds()` computed `lower` / `upper` as `adjusted_mean ± z_score * adjusted_std` with no clipping, producing biologically invalid values such as negative count thresholds and MT percentages above 100.

**Fix in `src/scLucid/qc/filtering.py`:**
- After computing thresholds, metric-aware clipping is applied:
  - Count-like metrics (`n_genes_by_counts`, `total_counts`): `lower` clipped to `>= 0.0`
  - Percentage metrics (`pct_counts_mt`, `pct_counts_hb`, etc.): `lower` clipped to `>= 0.0`, `upper` clipped to `<= 100.0`

**Validation evidence:**
- `validation_outputs/qc_lin2020/execution_summary.json` no longer contains negative `n_genes_by_counts.lower` or `pct_counts_mt.upper > 100`.
- Suspect samples (`GSM4679533`, `GSM4679535`) still surface with elevated MT bounds, but now within valid domains.

**New test:**
- `tests/qc/test_qc_threshold_mode.py::test_hierarchical_thresholds_are_clipped_to_valid_ranges`

### 2.3 Correct PBMC Threshold-Mode Test Intent

**Problem:** `test_hierarchical_and_pooled_yield_different_filtering_counts` asserted a hard inequality between hierarchical and pooled final cell counts on `pbmc3k.h5ad`. On a near-homogeneous baseline dataset, this assertion checked incidental divergence rather than intended behavior.

**Fix in `tests/qc/test_qc_threshold_mode.py`:**
- Replaced the hard inequality with a stability assertion: both modes must retain `> 50%` of the original cells. This encodes the real product requirement that hierarchical mode should not produce pathological removal on well-behaved data.

---

## 3. Preprocess Module Fixes

### 3.1 Pydantic v2.11 Compatibility

**Problem:** `config.model_fields` instance access is deprecated in Pydantic v2.11 and would raise errors when `warnings.filterwarnings('error')` is active.

**Fix in `src/scLucid/preprocess/normalize.py`:**
- Changed `config.model_fields.keys()` → `type(config).model_fields.keys()`
- Changed `active_config.model_fields` → `type(active_config).model_fields`

### 3.2 Adaptive Normalize Public API

**Fix in `src/scLucid/preprocess/__init__.py`:**
- Added missing exports for `AdaptiveNormalizationConfig`, `adaptive_normalize`, `estimate_cell_size_factors`, `quality_aware_normalize`

### 3.3 Core Preprocess Test Coverage

**New test files:**
- `tests/preprocess/test_normalize.py` — 23 tests covering standard/CLR/Pearson normalization, layer handling, sparse input, validation, metadata storage, config immutability
- `tests/preprocess/test_scale.py` — 19 tests covering z-score/robust/minmax scaling (dense + sparse), regress_out, inline regression, metadata, helper functions
- `tests/preprocess/test_integrate.py` — 9 tests covering batch_correction (harmony/scanorama/combat mocks, no-method early return), evaluate_integration
- `tests/preprocess/test_neighbors.py` — 5 tests for optimize_neighbors_pcs grid search, metadata storage, missing PCA error
- `tests/preprocess/test_backend.py` — 12 tests for backend abstraction (get/set, ScanpyBackend, RapidsBackend skip, custom subclass)

---

## 4. Cross-Module Polish (Preprocess, Analysis)

### 4.1 Unified kwargs Handling with Deep Copy

**Problem:** `normalize.py`, `scale.py`, `integrate.py`, `neighbors.py`, `analysis/workflow.py`, `analysis/clustering.py`, and `analysis/differential_expression/de_core.py` each implemented their own ad-hoc kwargs override logic. Most used shallow `model_copy()` with `hasattr` + `setattr` loops. None used `deep=True`, leaving caller configs vulnerable to nested mutation.

**Fix:**
- Added `apply_config_overrides(config, *, ignored_keys, **kwargs)` to `src/scLucid/base_config.py`.
  - Filters kwargs against `type(config).model_fields`
  - Warns on unknown parameters
  - Silently skips keys in `ignored_keys` (e.g. `force`)
  - Returns `config.model_copy(update=..., deep=True)`
- Replaced all ad-hoc override blocks with a single call to `apply_config_overrides`.

**Files changed:**
- `src/scLucid/preprocess/normalize.py`
- `src/scLucid/preprocess/scale.py`
- `src/scLucid/preprocess/integrate.py`
- `src/scLucid/preprocess/neighbors.py`
- `src/scLucid/analysis/workflow.py`
- `src/scLucid/analysis/clustering.py`
- `src/scLucid/analysis/differential_expression/de_core.py`
- `src/scLucid/preprocess/config.py` (re-exports from `base_config.py`)

**New tests:**
- `tests/preprocess/test_config_helpers.py` — 5 tests for deep copy, unknown param warnings, ignored keys, valid overrides

### 4.2 Shared Input Validation Utility

**Problem:** `normalize.py` had a private `_validate_normalization_input()` that checked empty matrix, NaN/Inf, and negative values. `scale.py` and `integrate.py` had no equivalent validation, meaning invalid matrices could propagate downstream with opaque errors.

**Fix:**
- Created `src/scLucid/preprocess/utils.py` with `validate_matrix_input(data, name, *, allow_negative=True)`.
- Replaced `_validate_normalization_input()` body with a call to `validate_matrix_input(..., allow_negative=False)` for count-based methods.
- Added `validate_matrix_input(adata.X, name="adata.X", allow_negative=True)` at the start of `scale_data()`.

### 4.3 Factory Methods for AnalysisWorkflowConfig

**Problem:** `AnalysisWorkflowConfig` lacked the convenient factory methods that `PreprocessingWorkflowConfig` already had (`from_simple_dict`, `quick`).

**Fix in `src/scLucid/analysis/config.py`:**
- Added `from_simple_dict(cls, simple_config)` — extracts `clustering_*`, `annotation_*`, `de_*` prefixed keys into nested configs
- Added `quick(cls, clustering_method="leiden", resolution=1.0, run_annotation=True, **kwargs)` — one-liner for standard analyses

### 4.4 Standardized Error Messages

**Fix:** Unified error message templates across preprocess and analysis modules to match a cross-module pattern:
- Preprocess failures: `"[preprocess] <Operation> failed: {e}. Check input data format..."`
- Analysis failures: `"[analysis] Workflow failed at step '{current_step}': {e}"`
- Unknown method: `"[<module>] Unknown <thing> method '<method>'. Expected one of: ..."`

### 4.5 Bug Fixes Discovered During Testing

**Fix in `src/scLucid/analysis/differential_expression/de_core.py`:**
- Missing imports: `_to_frac`, `_safe_filename` from `de_utils`; `plot_volcano` from `de_plots`. These caused `NameError` in `compare_groups` and `compare_conditions`.

**Fix in `src/scLucid/analysis/config.py`:**
- `annotation` field changed to `Optional[AnnotationConfig]` to support `quick(run_annotation=False)`.
- `CompareConditionsConfig` added missing `n_top_genes: int = Field(default=50, ge=1)` field, referenced by `compare_conditions`.

**Fix in `src/scLucid/analysis/workflow.py`:**
- Hard-coded `adata.obs['leiden']` replaced with dynamic `cluster_key` lookup in `run_custom_analysis` and `compare_clustering_resolutions`.
- `find_resolution` return value unpacking fixed (`eval_df, recommended_res` instead of dict access).
- Annotation config propagation: explicitly checks `isinstance(config.annotation, AnnotationConfig)` and passes `config=config.annotation` to `run_annotation()`. Dict-based inputs are converted to `AnnotationConfig(**config.annotation)` before passing.
- Removed invalid `save_path=` kwarg from `characterize_clusters()` call.

---

## 5. Annotation Config Propagation in Analysis Workflow

**Problem:** `src/scLucid/analysis/workflow.py` passed annotation configuration via `**config.annotation if isinstance(config.annotation, dict) else {}`, which silently dropped `AnnotationConfig` objects and fell back to `run_annotation()` defaults.

**Fixes in `src/scLucid/analysis/workflow.py`:**
- Annotation step now explicitly checks `isinstance(config.annotation, AnnotationConfig)` and passes `config=config.annotation` to `run_annotation()`.
- Dict-based inputs are converted to `AnnotationConfig(**config.annotation)` before passing.

**New tests:**
- `tests/analysis/test_workflow.py::test_run_standard_analysis_passes_annotation_config`
- `tests/analysis/test_workflow.py::test_run_standard_analysis_passes_dict_annotation_config`

---

## 6. Modified Files (Cumulative)

| File | Change |
|------|--------|
| `src/scLucid/base_config.py` | Added shared `apply_config_overrides` helper for cross-module config immutability |
| `src/scLucid/qc/workflow.py` | Deep-copy recommendation application; preserve user overrides; return `applied_config` from `_run_qc_workflow`; store `original_config` in trace; added `_build_qc_review_summary()` and `_export_qc_review_summary()` |
| `src/scLucid/qc/filtering.py` | Clip hierarchical thresholds to biologically valid domains (`>=0` for counts, `[0,100]` for percentages) |
| `src/scLucid/analysis/workflow.py` | Config immutability via `apply_config_overrides`; pass `AnnotationConfig` correctly; standardized error messages; dynamic `cluster_key` lookup |
| `src/scLucid/analysis/clustering.py` | Unified kwargs via `apply_config_overrides`; standardized error messages |
| `src/scLucid/analysis/differential_expression/de_core.py` | Unified kwargs via `apply_config_overrides`; added missing imports (`_to_frac`, `_safe_filename`, `plot_volcano`) |
| `src/scLucid/analysis/config.py` | Added `from_simple_dict` and `quick` factory methods; `annotation` made Optional; added `n_top_genes` to `CompareConditionsConfig` |
| `src/scLucid/preprocess/normalize.py` | Pydantic v2.11 fix; unified kwargs via `apply_config_overrides`; use shared `validate_matrix_input`; standardized error messages |
| `src/scLucid/preprocess/scale.py` | Unified kwargs via `apply_config_overrides`; added `validate_matrix_input`; standardized error messages |
| `src/scLucid/preprocess/integrate.py` | Unified kwargs via `apply_config_overrides`; standardized error messages |
| `src/scLucid/preprocess/neighbors.py` | Unified kwargs via `apply_config_overrides` |
| `src/scLucid/preprocess/config.py` | Re-exports `apply_config_overrides`; added `default()` factory method for canonical preprocessing path |
| `src/scLucid/preprocess/utils.py` | New file: `validate_matrix_input` shared validation utility |
| `src/scLucid/preprocess/__init__.py` | Added adaptive normalize exports |
| `src/scLucid/preprocess/intelligent/data_classes.py` | Added `to_review_summary()` method on `PreprocessingStrategy` |
| `src/scLucid/preprocess/intelligent/recommender.py` | Added `_export_preprocess_review_summary()`; stores `intelligent_review_summary` in trace; exports JSON/Markdown sidecars |
| `src/scLucid/qc/workflow.py` | Added `_build_qc_review_summary()` and `_export_qc_review_summary()`; review summary stored in trace and exported when `save_dir` is set |
| `tests/qc/test_qc_trace.py` | Extended with review-summary contract, filtering consistency, tumor-aware flag, and disk-export tests |
| `tests/preprocess/test_preprocess_semantics.py` | Added tests for `default()` factory, `to_review_summary()`, review summary storage, and disk export |
| `examples/preprocessing.py` | Added commented usage patterns for `default()` and intelligent preprocessing |
| `examples/README.md` | Updated preprocessing example description to mention `default()` and intelligent review summary |
| `tests/preprocess/test_normalize.py` | New file: 23 normalization tests |
| `tests/preprocess/test_scale.py` | New file: 19 scaling tests |
| `tests/preprocess/test_integrate.py` | New file: 9 integration tests |
| `tests/preprocess/test_neighbors.py` | New file: 5 neighbors tests |
| `tests/preprocess/test_backend.py` | New file: 12 backend tests |
| `tests/preprocess/test_config_helpers.py` | New file: 5 config helper tests |
| `tests/qc/test_qc_trace.py` | Added `original_config` to required trace keys |
| `tests/qc/test_qc_threshold_mode.py` | Stability-oriented assertion for PBMC; added clipping test |
| `tests/qc/test_qc_recommendation_executable.py` | Added caller-config immutability and user-override survival tests |
| `tests/analysis/test_workflow.py` | Extended: error recovery, custom analysis, resolution comparison, kwargs override, AnnotationConfig propagation |
| `tests/analysis/test_config.py` | New file: 15 tests for `AnalysisWorkflowConfig` factory methods and validation |
| `tests/analysis/test_de_core.py` | New file: 19 tests for `find_markers`, `filter_markers`, `compare_groups`, `compare_conditions` |

---

## 7. Real-Dataset Validation Results (Refreshed)

| Dataset | Mode | Before | After | Notes |
|---------|------|--------|-------|-------|
| pbmc3k | hierarchical | 2,700 | 2,523 | Per-sample thresholds computed; all values within valid domain |
| schlesinger2020.pdac | hierarchical | 6,499 | 5,892 | Tumor-aware: `outlier_mt` excluded from filtering; only `outlier_min_genes` applied |
| lin2020.pdac | pooled | 9,621 | 8,613 | Global thresholds; tumor-aware MT exclusion |
| lin2020.pdac | hierarchical | 9,621 | 9,593 | Sample-specific thresholds preserve more cells; suspect samples (`GSM4679533`, `GSM4679535`) now have clipped, valid bounds |

**Key clipping evidence (lin2020):**
- `GSM4679533` `pct_counts_mt.upper`: previously `127.6` → now `100.0`
- `GSM4679533` `n_genes_by_counts.lower`: previously `-1094.9` → now `0.0`
- `GSM4679535` `pct_counts_mt.lower`: previously `-23.89` → now `0.0`

---

## 8. Verification Commands

```bash
# Full QC test suite
pytest tests/qc/ -v

# Full preprocess test suite
pytest tests/preprocess/ -v

# Full analysis test suite
pytest tests/analysis/ -v

# Real-dataset validations
python scripts/validate_qc_pbmc3k.py
python scripts/validate_qc_schlesinger2020.py
python scripts/validate_qc_lin2020.py
```

---

## 9. Design Notes

### How user overrides interact with recommendations
- If the caller explicitly sets a threshold field (detected via `model_fields_set`), the recommendation **never overwrites it**.
- If the caller leaves a field at its default (default factory), the recommendation **may fill it**.
- The trace stores four artifacts:
  1. `original_config` — exactly what the user passed
  2. `recommendation` — what the engine suggested
  3. `applied_config` — what was actually executed
  4. `user_overrides` — fields where `original_config != recommendation`

### Hierarchical threshold clipping policy
- **Counts** (`n_genes_by_counts`, `total_counts`): lower bound floor at `0.0`.
- **Percentages** (`pct_counts_mt`, `pct_counts_hb`, etc.): lower bound floor at `0.0`, upper bound ceiling at `100.0`.
- Clipping is applied **after** shrinkage computation but **before** the thresholds are stored in the trace, so the trace reflects the values that were actually used for filtering.

---

## 10. Deferred / Out of Scope

- Restoration of the full `tests/tumor/` test suite and tumor workflow public surface
- CNV-augmented QC logic
- Preprocessing or annotation workflow redesign beyond config propagation
- Publication-ready report styling

---

## 11. Preprocess Module Fixes (Review-Fix Iteration Continuation)

### 11.1 Default Main Path Factory Method

**Problem:** `PreprocessingWorkflowConfig.quick()` disables regression and integration by default (`run_regression=False`, `run_integration=False`), which does not represent the canonical standard preprocessing path. Users had no single, well-documented entry point for the "true default" pipeline.

**Fix in `src/scLucid/preprocess/config.py`:**
- Added `PreprocessingWorkflowConfig.default(**kwargs)` class method.
- Returns a config with all major steps enabled:
  `run_regression=True`, `run_scaling=True`, `run_pca=True`, `run_neighbors=True`, `run_integration=True`.
- Docstring explicitly lists the 9-step canonical pipeline.

**New test:**
- `tests/preprocess/test_preprocess_semantics.py::test_workflow_config_default_sets_all_run_flags_true`

### 11.2 Intelligent Preprocessing Review Summary

**Problem:** `run_intelligent_preprocessing()` stored the raw `PreprocessingStrategy.to_dict()` in `adata.uns`, but there was no consolidated, reviewer-facing digest. The raw strategy dict is verbose and mixes evidence with recommendations, making manual review difficult.

**Fix in `src/scLucid/preprocess/intelligent/data_classes.py`:**
- Added `PreprocessingStrategy.to_review_summary()` which returns a structured dict with:
  - `data_profile` (cells, genes, sparsity, quality score, strategy type)
  - `hvg` / `pca` / `neighbors` / `resolution` recommendations (value + confidence + CI)
  - `batch_correction` decision (needs_correction, severity, recommended_method, alternatives)
  - `overall_confidence`, `concerns`, `recommendations`

**Fix in `src/scLucid/preprocess/intelligent/recommender.py`:**
- Added `_export_preprocess_review_summary()` helper that writes:
  - `preprocess_review_summary.json` (machine-readable)
  - `preprocess_review_summary.md` (human-readable Markdown with tables)
- `run_intelligent_preprocessing()` now:
  1. Calls `strategy.to_review_summary()`
  2. Stores it in `adata.uns["sclucid"]["preprocess"]["intelligent_review_summary"]`
  3. Exports JSON/Markdown sidecars when `save_dir` is set

**New tests:**
- `tests/preprocess/test_preprocess_semantics.py::test_strategy_to_review_summary_structure`
- `tests/preprocess/test_preprocess_semantics.py::test_run_intelligent_preprocessing_stores_review_summary`
- `tests/preprocess/test_preprocess_semantics.py::test_run_intelligent_preprocessing_exports_review_summary_to_disk`

### 11.3 Example and Documentation Updates

**Fix in `examples/preprocessing.py`:**
- Added three commented usage patterns at the top:
  - Option A: `PreprocessingWorkflowConfig.default()` one-liner
  - Option B: `run_intelligent_preprocessing()` with review summary
  - Option C: Manual step-by-step (existing code)

**Fix in `examples/README.md`:**
- Updated `preprocessing.py` description to mention `default()` and intelligent preprocessing with review summary export.

### 11.4 Verification Commands

See Section 8 above for consolidated verification commands.

