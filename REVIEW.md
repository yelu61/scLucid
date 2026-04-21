# REVIEW.md

## Review Scope

This review compares the claims in `FIX.md` against the current repository state, with emphasis on:
- the original tumor workflow Phase A+B summary,
- the later QC module iteration summary,
- whether the claimed fixes are actually implemented and verifiable.

---

## Findings

### 1. High — QC recommendations overwrite user intent and make `user_overrides` effectively untrustworthy

The QC workflow now applies recommendations unconditionally at the start of `_run_qc_workflow()`, without preserving any explicit user-specified thresholds. See [src/scLucid/qc/workflow.py](/Users/luye/Scripts/scLucid/src/scLucid/qc/workflow.py:447). This directly conflicts with the stated goal of capturing user overrides and auditing recommendation-vs-actual behavior.

The problem is compounded in `_apply_qc_recommendations()`, which uses `config.model_copy()` and then mutates nested config objects such as `marking_config.thresholds` and `doublet_config` in place; because `model_copy()` is shallow by default in Pydantic v2, this can mutate the caller’s original config object rather than creating an isolated applied-config snapshot. See [src/scLucid/qc/workflow.py](/Users/luye/Scripts/scLucid/src/scLucid/qc/workflow.py:303). The result is:
- explicit user thresholds are silently replaced when `use_recommendations=True`,
- the trace cannot reliably distinguish recommendation-applied values from caller-specified overrides,
- the config object passed by the caller may be mutated as a side effect.

This is a correctness issue, not just a trace-quality issue.

### 2. High — The new hierarchical QC policy can emit biologically invalid thresholds, and the validation output already shows it

The hierarchical branch in `AdaptiveThresholdCalculator._suggest_adaptive_thresholds()` computes `lower` and `upper` as `adjusted_mean ± z_score * adjusted_std` with no clipping or domain-aware bounds. See [src/scLucid/qc/filtering.py](/Users/luye/Scripts/scLucid/src/scLucid/qc/filtering.py:179). For count-like metrics this can produce negative lower bounds, and for percentages it can produce values above 100.

The shipped validation artifact already demonstrates this failure mode. In [validation_outputs/qc_lin2020/execution_summary.json](/Users/luye/Scripts/scLucid/validation_outputs/qc_lin2020/execution_summary.json), the hierarchical QC output includes examples such as:
- `GSM4679533` with `n_genes_by_counts.lower = -1094.9`
- `GSM4679533` with `pct_counts_mt.upper = 127.6`
- `GSM4679535` with `pct_counts_mt.lower = -23.89`

These thresholds are not biologically meaningful and make the “hierarchical policy preserves far more cells” result hard to trust. The implementation needs metric-aware clipping and plausibility constraints before the validation summary can be treated as evidence of success.

### 3. Medium — The new PBMC threshold-mode test locks in the wrong product behavior

The QC-specialized `FIX.md` and `PLAN.md` both position `pbmc3k` as a near-homogeneous multi-sample baseline where hierarchical and pooled policies should behave similarly. However, the test in [tests/qc/test_qc_threshold_mode.py](/Users/luye/Scripts/scLucid/tests/qc/test_qc_threshold_mode.py:45) asserts that hierarchical and pooled modes **must** produce different final cell counts.

That assertion checks incidental divergence rather than intended behavior. On a baseline dataset with similar sample distributions, a correct implementation could reasonably produce identical or near-identical filtering results. This test therefore encourages unstable threshold behavior rather than validating the product requirement.

### 4. Medium — The earlier tumor workflow summary in `FIX.md` still overstates what was actually fixed

`FIX.md` claims that `analysis/workflow.py` was fixed to “pass config correctly instead of unpacking a Pydantic object with `**`”, but the current code still does exactly that fallback. In the annotation step, [src/scLucid/analysis/workflow.py](/Users/luye/Scripts/scLucid/src/scLucid/analysis/workflow.py:203) calls:

```python
run_annotation(adata, **config.annotation if isinstance(config.annotation, dict) else {})
```

When `config.annotation` is a Pydantic `AnnotationConfig`, the workflow passes no config at all and silently falls back to `run_annotation()` defaults. That means the tumor workflow still cannot reliably carry annotation recommendations or user-provided annotation config through the standard analysis stage.

This remains a functional gap in the earlier Phase A+B implementation summary.

### 5. Medium — The earlier tumor verification section in `FIX.md` is not reproducible from the current repo state

The original tumor portion of `FIX.md` lists `tests/tumor/` and recommends `pytest tests/tumor/ -v`, but `tests/tumor` does not currently exist in the repository. Running that command now fails immediately with “file or directory not found”. That makes the original verification section stale at best.

`FIX.md` also says [src/scLucid/tumor/__init__.py](/Users/luye/Scripts/scLucid/src/scLucid/tumor/__init__.py:1) exports the tumor workflow surface and configs, but the module currently exports domain helpers only; it does not export `run_tumor_analysis`, `run_tumor_analysis_expert`, `TumorWorkflowConfig`, or `TumorAnalysisConfig`. So the summary still overstates the public-surface cleanup that actually landed.

---

## Verified Positives

The QC-focused iteration did land real, testable improvements:
- `QCWorkflowConfig` now exposes `threshold_mode`, `use_recommendations`, and `tissue_type`.
- `DoubletConfig` now has an executable `score_threshold` field.
- The QC trace contract under `adata.uns["sclucid"]["qc"]` is implemented and the new trace-focused tests pass.
- The four new QC test files named in `FIX.md` are present and currently passing.
- Real-data validation artifacts exist for:
  - `pbmc3k.h5ad`
  - `schlesinger2020.pdac.h5ad`
  - `lin2020.pdac.h5ad`

These are meaningful improvements. The main remaining issues are around policy correctness, trace semantics, and stale claims in `FIX.md`.

---

## Verification Performed

I directly verified the following during review:

- Read and compared `FIX.md` against current code and file layout.
- Inspected:
  - [src/scLucid/qc/config.py](/Users/luye/Scripts/scLucid/src/scLucid/qc/config.py)
  - [src/scLucid/qc/workflow.py](/Users/luye/Scripts/scLucid/src/scLucid/qc/workflow.py)
  - [src/scLucid/qc/filtering.py](/Users/luye/Scripts/scLucid/src/scLucid/qc/filtering.py)
  - [src/scLucid/qc/doublet.py](/Users/luye/Scripts/scLucid/src/scLucid/qc/doublet.py)
  - [src/scLucid/analysis/workflow.py](/Users/luye/Scripts/scLucid/src/scLucid/analysis/workflow.py)
  - [src/scLucid/tumor/__init__.py](/Users/luye/Scripts/scLucid/src/scLucid/tumor/__init__.py)
- Checked validation artifacts under [validation_outputs](/Users/luye/Scripts/scLucid/validation_outputs).
- Ran:

```bash
pytest tests/qc/test_qc_trace.py \
       tests/qc/test_qc_threshold_mode.py \
       tests/qc/test_qc_recommendation_executable.py \
       tests/qc/test_qc_tumor_aware.py -q
```

Result:
- `11 passed` in about `40s`

- Ran:

```bash
pytest tests/tumor -q
```

Result:
- immediate failure because `tests/tumor` does not exist

---

## Recommended Next Fixes

Priority order:

1. Make QC recommendation application preserve explicit user overrides and avoid mutating caller configs.
2. Add domain-aware clipping to hierarchical thresholds before treating validation outputs as evidence.
3. Replace the PBMC threshold-mode assertion with a bounded-divergence or stability-oriented assertion.
4. Correct `analysis/workflow.py` so `AnnotationConfig` objects are actually passed through.
5. Update `FIX.md` to match the current repository state, especially around tumor exports and test commands.

