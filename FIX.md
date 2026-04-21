# FIX Summary — scLucid QC Module Review-Fix Iteration

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

## 2. QC Module Fixes (Current Iteration)

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

### 2.4 Annotation Config Propagation in Analysis Workflow

**Problem:** `src/scLucid/analysis/workflow.py` passed annotation configuration via `**config.annotation if isinstance(config.annotation, dict) else {}`, which silently dropped `AnnotationConfig` objects and fell back to `run_annotation()` defaults.

**Fixes in `src/scLucid/analysis/workflow.py`:**
- Annotation step now explicitly checks `isinstance(config.annotation, AnnotationConfig)` and passes `config=config.annotation` to `run_annotation()`.
- Dict-based inputs are converted to `AnnotationConfig(**config.annotation)` before passing.
- Removed invalid `save_path=` kwarg from `characterize_clusters()` call (pre-existing bug).

**New tests:**
- `tests/analysis/test_workflow.py::test_run_standard_analysis_passes_annotation_config`
- `tests/analysis/test_workflow.py::test_run_standard_analysis_passes_dict_annotation_config`

---

## 3. Modified Files

| File | Change |
|------|--------|
| `src/scLucid/qc/workflow.py` | Deep-copy recommendation application; preserve user overrides; return `applied_config` from `_run_qc_workflow`; store `original_config` in trace; wire applied config through `filter_cells` and reporting |
| `src/scLucid/qc/filtering.py` | Clip hierarchical thresholds to biologically valid domains (`>=0` for counts, `[0,100]` for percentages) |
| `src/scLucid/analysis/workflow.py` | Pass `AnnotationConfig` objects correctly to `run_annotation`; remove invalid `save_path` kwarg from `characterize_clusters` |
| `tests/qc/test_qc_trace.py` | Added `original_config` to required trace keys |
| `tests/qc/test_qc_threshold_mode.py` | Stability-oriented assertion for PBMC; added clipping test |
| `tests/qc/test_qc_recommendation_executable.py` | Added caller-config immutability and user-override survival tests |
| `tests/analysis/test_workflow.py` | New file: verifies `AnnotationConfig` propagation |

---

## 4. Real-Dataset Validation Results (Refreshed)

| Dataset | Mode | Before | After | Notes |
|---------|------|--------|-------|-------|
| pbmc3k | hierarchical | 2,700 | 2,523 | Per-sample thresholds computed; all values within valid domain |
| schlesinger2020.pdac | hierarchical | 6,499 | 5,898 | Tumor-aware: `outlier_mt` excluded from filtering; only `outlier_min_genes` applied |
| lin2020.pdac | pooled | 9,621 | 8,606 | Global thresholds; tumor-aware MT exclusion |
| lin2020.pdac | hierarchical | 9,621 | 9,593 | Sample-specific thresholds preserve more cells; suspect samples (`GSM4679533`, `GSM4679535`) now have clipped, valid bounds |

**Key clipping evidence (lin2020):**
- `GSM4679533` `pct_counts_mt.upper`: previously `127.6` → now `100.0`
- `GSM4679533` `n_genes_by_counts.lower`: previously `-1094.9` → now `0.0`
- `GSM4679535` `pct_counts_mt.lower`: previously `-23.89` → now `0.0`

---

## 5. Verification Commands

```bash
# Full QC test suite (50 tests)
pytest tests/qc/ -v

# Analysis workflow config propagation tests
pytest tests/analysis/test_workflow.py -v

# Real-dataset validations
python scripts/validate_qc_pbmc3k.py
python scripts/validate_qc_schlesinger2020.py
python scripts/validate_qc_lin2020.py
```

---

## 6. Design Notes

### How user overrides interact with recommendations
- If the caller explicitly sets a threshold field (detected via `model_fields_set`), the recommendation **never overwrites it**.
- If the caller leaves a field at its default (default factory), the recommendation **may fill it**.
- The trace stores three artifacts:
  1. `original_config` — exactly what the user passed
  2. `recommendation` — what the engine suggested
  3. `applied_config` — what was actually executed
  4. `user_overrides` — fields where `original_config != recommendation`

### Hierarchical threshold clipping policy
- **Counts** (`n_genes_by_counts`, `total_counts`): lower bound floor at `0.0`.
- **Percentages** (`pct_counts_mt`, `pct_counts_hb`, etc.): lower bound floor at `0.0`, upper bound ceiling at `100.0`.
- Clipping is applied **after** shrinkage computation but **before** the thresholds are stored in the trace, so the trace reflects the values that were actually used for filtering.

---

## 7. Deferred / Out of Scope

- Restoration of the full `tests/tumor/` test suite and tumor workflow public surface
- CNV-augmented QC logic
- Preprocessing or annotation workflow redesign beyond config propagation
- Publication-ready report styling
