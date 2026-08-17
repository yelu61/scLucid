# Reviewing Results

`review_run()` converts module-specific review summaries into one product-facing
decision table. It answers four questions:

1. Can this stage be handed downstream?
2. What still requires human review?
3. Which evidence should be inspected?
4. From which stage should the project be rerun after a change?

## Basic Use

```python
review = scl.review_run(adata)

print(review.overall_status)
print(review.to_frame())
print(review.show_next_actions())
```

The serializable result is also stored at:

```python
adata.uns["sclucid"]["run_review"]
```

and appears near the top of `scl.export_audit_report()`.

## Status Contract

| Status | Meaning | User action |
|---|---|---|
| `BLOCKED` | A structural prerequisite or required output is missing | Fix it before relying on downstream interpretation |
| `REVIEW` | The first pass can exist, but a biological/statistical choice is unresolved | Inspect evidence, accept or override, then rerun the affected stage |
| `READY` | The recorded stage contract has no unresolved blocker | Continue to the next declared handoff |
| `NOT_RUN` | No review summary exists for that stage | Run it only if required by the project objective |

`READY` is a workflow handoff claim, not biological truth. A ready annotation
stage may still contain exploratory labels; a ready DE stage is only formal if
its result contract records a sample-aware estimand and adequate replication.

## Decision Table Fields

| Field | Interpretation |
|---|---|
| `stage` | QC, preprocess, analysis, or tumor |
| `status` | BLOCKED, REVIEW, or READY |
| `decision` | Parameter, handoff, or review action |
| `recommended` | System recommendation or required target state |
| `applied` | Value actually used or observed |
| `reason` | Why the row exists |
| `evidence` | Source fields, plots, or tables to inspect |
| `next_action` | Concrete action to take |
| `rerun_scope` | Earliest stage that must be rerun after a change |

## Stage Review Checklist

### QC

Inspect raw-count semantics, per-sample retention, mitochondrial/count/gene
thresholds, doublet evidence, ambient evidence, and tumor-aware sensitivity
cells. A threshold recommendation is not enough: compare recommended versus
applied values and verify which samples or biological states are affected.

### Preprocess

Confirm the counts-to-normalized layer contract, HVGs, PCs, graph stability,
and any integration decision. If integration is used, compare integrated and
unintegrated representations and retain unintegrated normalized expression for
marker and DE interpretation.

### Analysis

Review clustering resolution, broad-lineage annotation, fine-label conflicts,
post-hoc cluster QC, and the inference contract. Cell-level marker tests are
exploratory. Condition DE and composition claims require sample-aware methods
and the correct experimental unit.

### Tumor

Review the claim boundary as well as readiness. Malignant labels require CNV
plus transcriptional/sample-aware evidence. TME, therapy, and state results
remain exploratory when their upstream annotation or malignancy calls are under
review.

## Accepting Or Overriding A Decision

For every material override, retain:

- the recommended value;
- the applied value;
- the evidence inspected;
- the biological or technical rationale;
- the affected stage and rerun scope;
- any limitation that remains on downstream claims.

After rerunning a stage, call `review_run()` again. Do not manually edit the
unified table as if it were source evidence; it is a read model derived from the
module review summaries.

## Highest-Risk Failure Modes

- reading every warning but resolving none of the recorded actions;
- treating cells as independent replicates;
- integrating a batch that is confounded with treatment or patient;
- promoting automated cell labels without multi-marker review;
- inferring malignancy from epithelial markers alone;
- changing a QC/preprocessing choice without rerunning dependent stages.

The smallest defensible next step is the first `BLOCKED` action, otherwise the
first `REVIEW` action returned by `show_next_actions()`.
