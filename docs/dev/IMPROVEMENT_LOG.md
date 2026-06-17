# scLucid Improvement Log

This document captures concrete problems, gaps, and enhancement ideas that emerge
from real-world use of scLucid in analysis projects. It is the bridge between
**project-driven discovery** and **scLucid source-code polishing**.

## How to Use This Log

1. **During analysis** (in your project directory):
   - Do not modify scLucid source code.
   - Record each issue using the template below.
   - Attach only minimal, structured context (not full AnnData objects or long
     notebooks).

2. **Sync to this log** (in scLucid directory):
   - Periodically move high-priority items into this file.
   - Group related issues into themes.

3. **Source-code polishing** (in scLucid directory):
   - Pick 1–3 related items per session.
   - Use Plan mode for behavior changes; direct edits for obvious bugs.
   - Add or update tests.
   - Commit with a reference to the log entry.

4. **Verify in project** (back in project directory):
   - Install the updated scLucid (`pip install -e /path/to/scLucid`).
   - Run the analysis cell or script that exposed the issue.
   - Update the log entry status.

## Relationship to Other Docs

- `API_TRIAGE.md` — tracks uncertain *public API symbols* (promote/demote/remove).
- `CORE_WORKFLOW_COVERAGE_MATRIX.md` — tracks *workflow-stage coverage*.
- `IMPROVEMENT_LOG.md` (this file) — tracks *behavioral improvements, bugs, and
  usability issues* discovered during real analyses.

## Entry Template

```markdown
### Ixxx: Short title

| Field | Value |
|-------|-------|
| **Module** | qc / preprocess / analysis / tumor / tools / utils / plotting |
| **Function(s)** | `run_standard_qc`, `suggest_qc_thresholds`, ... |
| **Discovered in** | project name / dataset / date |
| **Severity** | blocker / high / medium / low |
| **Status** | reported / triaged / in_progress / resolved / wontfix |
| **Phenomenon** | What went wrong or felt awkward? |
| **Expected behavior** | What should happen instead? |
| **Minimal repro** | 5–15 lines of code that trigger the issue |
| **Key output / error** | One error message, one number, or one small table |
| **Proposed fix** | Optional rough direction |
| **Resolution commit** | Filled when fixed |
| **Verification** | Filled when verified in a real project |
```

## Active Items

*(All active items from this batch have been moved to **Resolved Items**.)*

## Resolved Items

### I001: Add `is_raw_count_matrix` utility for raw-count semantics guard

| Field | Value |
|-------|-------|
| **Module** | utils / preprocess |
| **Function(s)** | `is_raw_count_matrix` |
| **Status** | resolved |
| **Resolution summary** | Implemented in `src/scLucid/utils/validation.py` and re-exported from `utils`. `_looks_like_counts` and `_matrix_looks_like_counts` now share the same canonical diagnostics. |
| **Resolution commit** | TBD |
| **Verification** | `tests/utils/test_validation.py`; smoke + affected module tests pass. |

### I002: Add `build_metadata_dicts` helper for multi-sample loading

| Field | Value |
|-------|-------|
| **Module** | utils |
| **Function(s)** | `build_metadata_dicts` |
| **Status** | resolved |
| **Resolution summary** | Implemented in `src/scLucid/utils/helpers.py`. Converts `{sample: value}` group/batch dicts into `{column: {sample: value}}` for `read_10x` / `load_10x_data`. |
| **Resolution commit** | TBD |
| **Verification** | `tests/utils/test_helpers.py`; smoke + affected module tests pass. |

### I003: Add `audit_filtering` helper for QC retention audit

| Field | Value |
|-------|-------|
| **Module** | qc |
| **Function(s)** | `audit_filtering` |
| **Status** | resolved |
| **Resolution summary** | Implemented in `src/scLucid/qc/filtering/core.py`. Compares cell counts before/after filtering by sample and optional group. |
| **Resolution commit** | TBD |
| **Verification** | `tests/qc/test_filtering.py`; smoke + affected module tests pass. |

### I004: Add `audit_doublets` helper for post-filter doublet check

| Field | Value |
|-------|-------|
| **Module** | qc / doublet |
| **Function(s)** | `audit_doublets` |
| **Status** | resolved |
| **Resolution summary** | Implemented in `src/scLucid/qc/doublet/core.py`. Summarizes remaining doublet predictions and score distributions after filtering. |
| **Resolution commit** | TBD |
| **Verification** | `tests/qc/test_doublet.py`; smoke + affected module tests pass. |

### I005: Detect batch-biology confounding before integration

| Field | Value |
|-------|-------|
| **Module** | preprocess |
| **Function(s)** | `detect_integration_confounding`, `diagnose_integration_risk` |
| **Status** | resolved |
| **Resolution summary** | Implemented in `src/scLucid/preprocess/integrate.py`. Detects one-to-one confounding and produces a structured risk assessment. `IntegrationConfig.auto_decide` wires the check into `batch_correction`. |
| **Resolution commit** | TBD |
| **Verification** | `tests/preprocess/test_integrate.py`; smoke + affected module tests pass. |

### I006: Add `resolve_qc_thresholds` for merging threshold sources

| Field | Value |
|-------|-------|
| **Module** | qc |
| **Function(s)** | `resolve_qc_thresholds` |
| **Status** | resolved |
| **Resolution summary** | Implemented in `src/scLucid/qc/filtering/suggestions.py`. Merges intelligent/MAD/manual thresholds with configurable policy. Wired into `run_standard_qc` threshold application; user-explicit thresholds are authoritative. |
| **Resolution commit** | TBD |
| **Verification** | `tests/qc/test_filtering.py`, `tests/qc/test_qc_recommendation_executable.py`; smoke + affected module tests pass. |

### I007: Add `decide_integration` for auto integration decision

| Field | Value |
|-------|-------|
| **Module** | preprocess |
| **Function(s)** | `decide_integration` |
| **Status** | resolved |
| **Resolution summary** | Implemented in `src/scLucid/preprocess/integrate.py`. Returns `(run, warnings, risk_dict)` for `run_integration="auto"` or explicit bool. |
| **Resolution commit** | TBD |
| **Verification** | `tests/preprocess/test_integrate.py`; smoke + affected module tests pass. |

### I008: Add `set_raw` option to `normalize_data`

| Field | Value |
|-------|-------|
| **Module** | preprocess |
| **Function(s)** | `normalize_data(..., set_raw=True)` |
| **Status** | resolved |
| **Resolution summary** | `NormalizationConfig` gained `set_raw`; `normalize_data` optionally sets `adata.raw` automatically after normalization. |
| **Resolution commit** | TBD |
| **Verification** | `tests/preprocess/test_normalize.py`; smoke + affected module tests pass. |

## Theme Backlog

Use this section to group related issues that should be tackled together in a
single polishing session.

### QC robustness

- *no entries yet*

### Preprocess defaults

- *no entries yet*

### Analysis / annotation usability

- *no entries yet*

### Tumor interpretation edge cases

- *no entries yet*

### Documentation / discoverability

- *no entries yet*

## Polishing Session Checklist

When you sit down to polish scLucid:

- [ ] Review active items and pick 1–3 related ones.
- [ ] Confirm each item has a minimal repro or clear expected behavior.
- [ ] Implement the fix.
- [ ] Add/update unit or integration tests.
- [ ] Update relevant docs (docstrings, this log, coverage matrix if needed).
- [ ] Run `python scripts/audit_public_api.py --write` if public API changed.
- [ ] Run `pytest` for affected modules.
- [ ] Commit with `Refs Ixxx` in the message.
- [ ] Verify in the project that reported the issue.
- [ ] Move the item to **Resolved Items** with resolution commit and verification.
