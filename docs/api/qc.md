# QC API

This page is generated at build time with `mkdocstrings`. The Markdown source is
intentionally small; the rendered documentation below is pulled from the current
Python objects and docstrings.

The stable user entrypoints are the workflow functions and configuration
objects. Lower-level modules are included so reviewers can inspect the decision
engine, contamination scoring, doublet calling, and filtering contracts.

## Recommended Workflow Contract

For most analyses, use `scLucid.qc.run_standard_qc()` or the package-level
workflow wrapper instead of calling threshold internals directly. The maintained
QC path is:

1. calculate QC metrics and optional cell-cycle evidence;
2. recommend threshold policy from count/percentage metrics;
3. resolve threshold decisions and mark evidence columns;
4. build `qc_decision`, `qc_remove`, `qc_reason`, and `qc_confidence`;
5. call `filter_cells()` only after final evidence exists;
6. review `qc_review_summary.json`, `qc_review_summary.md`, and
   `qc_benchmark.md` when `save_dir` is set.

The benchmark and review summaries are evidence for human review, not proof that
the filtered object is biologically correct. Treat `status="fail"` as a stop
signal, `status="review_required"` as a manual-review signal, and `status="pass"`
as supporting evidence to archive with downstream preprocessing records.

Important review locations:

- `adata.uns["sclucid"]["qc"]["review_summary"]["data"]`
- `review_summary["qc_reviewer_table"]`
- `review_summary["decision_table"]`
- `review_summary["benchmark_summary"]["assessment"]`
- `review_summary["evidence_bundle"]`

## Workflow

::: scLucid.qc.workflow
    options:
      show_root_heading: true

## Artifact Contract

::: scLucid.qc.artifacts
    options:
      show_root_heading: true

## Configuration

::: scLucid.qc.config
    options:
      show_root_heading: true

## Policy And Decisions

::: scLucid.qc.policy.thresholds
    options:
      show_root_heading: true

::: scLucid.qc.policy.decisions
    options:
      show_root_heading: true

::: scLucid.qc.policy.benchmark
    options:
      show_root_heading: true

## Metrics And Filtering

::: scLucid.qc.metrics
    options:
      show_root_heading: true

::: scLucid.qc.filtering
    options:
      show_root_heading: true

## Contamination

::: scLucid.qc.ambient
    options:
      show_root_heading: true
