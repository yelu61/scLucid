# Planning And Run Review API

The decision API is the product-facing layer above module-specific configs and
review summaries.

Recommended sequence:

```python
context = scl.ProjectContext(...)
plan = scl.plan_analysis(adata, context=context)
adata = scl.run_pipeline(adata, plan=plan)
review = scl.review_run(adata)
```

`ProjectContext` is the public alias of `AnalysisContext`; existing code using
`AnalysisContext` remains compatible.

## Project Context

::: scLucid.utils.context.AnalysisContext

## Analysis Plan

::: scLucid.decision.AnalysisPlan

## Decision Card

::: scLucid.decision.DecisionCard

## Run Review

::: scLucid.decision.RunReview

## Plan Analysis

::: scLucid.decision.plan_analysis

## Review Run

::: scLucid.decision.review_run
