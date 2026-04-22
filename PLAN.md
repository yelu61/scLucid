# scLucid PLAN.md

## Review-Fix Iteration Plan

### 1. Purpose

This iteration is a focused repair cycle driven by `REVIEW.md`.

The goal is not to add new capability. The goal is to make the recently added QC and tumor-related changes trustworthy, internally consistent, and accurately documented.

This cycle is complete only if:
- recommendation application no longer corrupts user intent,
- hierarchical thresholds are biologically valid,
- tests reflect intended product behavior,
- annotation config propagation is actually fixed,
- `FIX.md` matches the real repository state.

---

### 2. Scope Of This Iteration

This cycle is limited to five repair targets:

1. QC recommendation application semantics
2. Hierarchical threshold validity
3. QC threshold-mode test correctness
4. Analysis annotation config propagation
5. `FIX.md` accuracy

This cycle does **not** include:
- redesign of QC architecture,
- new tumor features,
- new reporting surfaces,
- new benchmark claims,
- broader workflow expansion.

---

### 3. Required Changes

#### 3.1 Fix QC Recommendation Application Semantics

Current problem:
- QC recommendations are applied unconditionally,
- caller-specified config can be overwritten,
- trace fields such as `user_overrides` are no longer reliable,
- shallow-copy behavior may mutate the original config object.

Required implementation changes:
- recommendation application must operate on an isolated applied-config object, not mutate the caller’s config in place,
- explicit user-provided values must take precedence over recommended values,
- only unset or default-valued fields may be auto-filled from recommendation,
- the trace must distinguish:
  - original user config,
  - recommendation,
  - final applied config,
  - actual user overrides relative to recommendation.

Acceptance criteria:
- the caller’s config object remains unchanged after `run_standard_qc()` / `run_advanced_qc()`,
- user-specified thresholds are preserved when recommendations are enabled,
- `user_overrides` records genuine user-vs-recommendation divergence rather than recommendation-induced self-diff.

#### 3.2 Add Domain-Aware Bounds To Hierarchical Thresholds

Current problem:
- hierarchical threshold generation can emit negative lower bounds,
- percentage thresholds can exceed 100,
- validation artifacts already show biologically invalid results.

Required implementation changes:
- clip count-like lower bounds to `>= 0`,
- clip percentage-like bounds to `[0, 100]`,
- apply metric-aware bounds consistently for:
  - `n_genes_by_counts`
  - `total_counts`
  - `pct_counts_mt`
  - any future percentage metric routed through the same helper,
- preserve raw pre-clipping values only if needed for debugging, but never use them for execution.

Acceptance criteria:
- no executed QC threshold in trace is biologically invalid,
- `validation_outputs/qc_lin2020/execution_summary.json` no longer contains negative count thresholds or MT thresholds above 100,
- suspect samples in `lin2020.pdac.h5ad` still surface as poor-quality candidates after clipping.

#### 3.3 Correct QC Threshold-Mode Test Intent

Current problem:
- the PBMC threshold-mode test currently asserts that hierarchical and pooled results must differ,
- this conflicts with the intended product behavior for a relatively homogeneous baseline dataset.

Required implementation changes:
- replace the hard inequality assertion with a correctness-oriented expectation,
- for `pbmc3k.h5ad`, test should validate:
  - hierarchical mode produces per-sample thresholds,
  - pooled mode produces no per-sample thresholds,
  - final results are stable and not wildly divergent,
  - hierarchical mode does not invent pathological divergence on near-homogeneous samples.

Acceptance criteria:
- the revised test encodes intended behavior, not incidental implementation noise,
- PBMC baseline validation still passes after hierarchical threshold clipping.

#### 3.4 Actually Fix Annotation Config Propagation In Analysis Workflow

Current problem:
- `FIX.md` claims annotation config passing was fixed,
- current `analysis/workflow.py` still drops `AnnotationConfig` objects and falls back to defaults.

Required implementation changes:
- `run_standard_analysis()` must pass annotation configuration explicitly to `run_annotation()`,
- both dict-based and Pydantic-config-based inputs must work,
- the configured annotation method, cluster key, and output key must actually be honored during workflow execution.

Acceptance criteria:
- analysis workflow no longer silently ignores `AnnotationConfig`,
- at least one test covers the object-config path end to end,
- the implementation behavior matches the claim in `FIX.md`.

#### 3.5 Bring `FIX.md` Back In Sync With Reality

Current problem:
- `FIX.md` still contains stale tumor-workflow claims and dead verification commands,
- it overstates some fixes and references files/tests not present in the current tree.

Required implementation changes:
- remove or rewrite any claim that is not true in the current repo,
- remove references to missing paths such as `tests/tumor/` if they are not restored,
- correct statements about exported tumor workflow surface if `src/scLucid/tumor/__init__.py` still does not export it,
- keep the QC iteration section, but revise any result summary if threshold clipping changes outputs.

Acceptance criteria:
- every verification command in `FIX.md` must run or be explicitly marked historical,
- every claimed file/test/export must exist in the current repository,
- `FIX.md` should be readable as a factual execution summary, not a mixed historical log.

---

### 4. Required Test And Validation Work

The implementation must update and rerun the following:

#### 4.1 QC Tests
- `tests/qc/test_qc_trace.py`
- `tests/qc/test_qc_threshold_mode.py`
- `tests/qc/test_qc_recommendation_executable.py`
- `tests/qc/test_qc_tumor_aware.py`

Additional required coverage:
- a test proving caller config is not mutated by recommendation application,
- a test proving explicit user thresholds survive when recommendations are enabled,
- a test proving clipped hierarchical thresholds stay within valid ranges.

#### 4.2 Analysis Workflow Test
- add or update a test that exercises `run_standard_analysis()` with a real `AnnotationConfig` object and verifies the config is honored.

#### 4.3 Real-Data Validation Refresh
- rerun:
  - `scripts/validate_qc_pbmc3k.py`
  - `scripts/validate_qc_schlesinger2020.py`
  - `scripts/validate_qc_lin2020.py`
- refresh `validation_outputs/` to reflect the repaired threshold behavior.

Minimum review evidence required:
- new execution summary for `lin2020.pdac.h5ad`,
- proof that suspect samples `GSM4679533` and `GSM4679535` still surface in QC output,
- no invalid threshold values in the stored QC trace.

---

### 5. Deliverables For Review

When this cycle is implemented, review input must include:
- updated `FIX.md`,
- updated `PLAN.md` only if the scope changes again,
- changed files list,
- test commands and outputs,
- refreshed `validation_outputs/`,
- a short note explaining how user overrides now interact with recommendations,
- a short note explaining hierarchical threshold clipping policy.

---

### 6. Completion Standard

This cycle is complete only when all of the following are true:
- QC recommendation application preserves user intent and does not mutate caller config,
- hierarchical thresholds are clipped to valid domains,
- PBMC threshold-mode tests reflect intended behavior,
- `run_standard_analysis()` truly propagates `AnnotationConfig`,
- `FIX.md` is factually consistent with the repository,
- refreshed tests and validation artifacts support the updated claims.

---

## QC And Preprocess Polishing Plan

### Summary

`qc` and `preprocess` should now be treated as basically feature-complete foundations rather than active feature-expansion fronts.

The next cycle for both modules is a refinement cycle focused on:
- trust,
- API clarity,
- reviewability,
- validation,
- and documentation coherence.

Net-new feature expansion should shift primarily to `analysis`. Work in `qc` and `preprocess` should be accepted when it improves reliability, interpretability, or default-path usability.

### QC Module Assessment

Current strengths:
- sample-aware workflow behavior with fallback handling for missing sample keys,
- tumor-aware QC logic that favors cautionary flagging over blind filtering,
- multi-evidence doublet detection combining algorithmic and marker-based heuristic evidence,
- recommendation and reporting surfaces for threshold review,
- strong config-driven execution and stored traceability.

Current weaknesses:
- overlapping control surfaces make the default user path heavier than necessary,
- heuristic quality depends strongly on marker quality and tissue context,
- intelligent recommendation logic still needs stronger validation to build trust,
- override semantics must remain easy to explain and inspect,
- reviewer-facing outputs exist but are not yet fully standardized.

### QC Polishing Directions

- Simplify the default path so most users can rely on `run_standard_qc()` with a shorter, safer configuration surface.
- Preserve expert controls, but separate advanced knobs from default knobs more clearly in docs and examples.
- Standardize QC review artifacts so every run can emit:
  - recommendation summary,
  - applied-threshold summary,
  - user-override summary,
  - sample-level threshold table,
  - tumor-aware caution note when applicable.
- Strengthen trust in intelligent QC by validating representative homogeneous, heterogeneous, and tumor-like datasets.
- Tighten and document recommendation-vs-user-override semantics consistently across code, trace outputs, and examples.
- Do not add new QC algorithms unless validation evidence shows a real gap in current behavior.

The next QC cycle is about trust, clarity, and reviewability, not feature breadth.

### Preprocess Module Assessment

Current strengths:
- a complete multi-step preprocessing workflow,
- flexible step control, resume/error recovery, and memory-management options,
- multiple integration methods with traceability,
- intelligent preprocessing recommendation support,
- backend abstraction and advanced utility support.

Current weaknesses:
- the capability surface is broad and can feel fragmented,
- the relationship between the main workflow, intelligent recommendation, and experimental normalization paths is not yet clear enough,
- many preprocessing and integration options exist without a strongly opinionated primary path,
- benchmark and documentation support for method choice is thinner than the implementation surface,
- advanced options can overwhelm standard users.

### Preprocess Polishing Directions

- Establish one clear default preprocessing path and explicitly label other routes as advanced or experimental.
- Reduce ambiguity around normalization, HVG, integration, and neighbor-selection defaults.
- Standardize intelligent preprocessing outputs so users can inspect:
  - recommended config summary,
  - rationale trace,
  - method-selection note,
  - integration decision summary.
- Unify naming and documentation across workflow config, intelligent strategy objects, and public API examples.
- Validate the default path and common integration modes on small, medium, and multi-batch scenarios.
- Do not add more preprocessing methods unless they replace or clearly outperform an existing path.

The next preprocess cycle is about path consolidation and decision support, not option expansion.

## Current Non-Goals For QC And Preprocess

- No broad expansion of new QC heuristics.
- No addition of more batch-correction methods without benchmark justification.
- No major architecture rewrite unless required by simplification work.
- No expansion of experimental branches into default workflow without validation.
- No effort to maximize option count at the expense of clearer defaults.

## Priority Position In The Roadmap

- `analysis` is now the main expansion target.
- `qc` and `preprocess` should receive targeted polishing work in parallel or between analysis milestones.
- Work on these modules should be accepted when it improves reliability, interpretability, documentation quality, or public API coherence.
