# Parameter Profiles

Profiles reduce the number of choices shown at the start of a project. They are
context labels for conservative planning; they do not silently auto-apply a
scientifically optimal parameter set.

```python
plan = scl.plan_analysis(adata, context=context, profile="auto")
print(plan.profile)
```

## Current Profiles

| Profile | Intended context | Conservative starting policy |
|---|---|---|
| `baseline` | Standard non-tumor dataset | Preserve counts, light preprocessing, no integration by default |
| `tumor_conservative` | Single tumor dataset | Tumor-aware QC review, broad annotation first, separate malignancy review |
| `multi_sample_tumor` | Multi-sample/patient tumor cohort | Per-sample QC, batch-condition audit, sample-aware inference, multi-evidence malignancy |
| `treatment_response` | Paired or grouped intervention study | Explicit condition/experimental unit/pairing; sample-level DE and composition |
| `cell_line` | Cell-line experiment | Replicate identity remains required; tumor-tissue assumptions are disabled |

With `profile="auto"`, scLucid selects a label from the declared
`ProjectContext`. The selected profile is recorded in the analysis plan and can
be overridden explicitly.

## Parameters Safe To Expose Early

Most users should initially see only:

- project profile;
- sample, condition, batch, experimental-unit, and paired keys;
- whether counts semantics are confirmed;
- whether integration should be evaluated;
- the intended annotation depth;
- the formal comparison/contrast.

Exact QC thresholds, HVG flavor/count, PC count, neighbor count, clustering
resolution, annotation evidence methods, and tumor scoring parameters belong in
stage review cards or expert configs.

## Decisions That Must Not Be Hidden In A Profile

- filtering a tumor/stress state because it resembles low-quality cells;
- integrating when batch and condition are confounded;
- selecting a final clustering resolution;
- treating automated annotation as final;
- choosing the statistical contrast and experimental unit;
- calling malignant cells from a single evidence source.

Profiles reduce setup friction. They do not replace the plan -> run -> review ->
override -> rerun evidence loop.
