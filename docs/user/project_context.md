# Project Context

`ProjectContext` is the minimum biological and study-design contract supplied
before scLucid recommends or runs a workflow. Its purpose is to prevent a
technically valid pipeline from silently using the wrong sample, batch,
condition, or inference unit.

## Minimal Context

```python
import scLucid as scl

context = scl.ProjectContext(
    dataset_type="tumor_tissue",
    species="human",
    tissue="lung",
    cancer_type="NSCLC",
    sample_key="sample_id",
    batch_key="library_batch",
    condition_key="treatment_group",
    experimental_unit_key="patient_id",
    paired_key="patient_id",
    study_objective="paired treatment response",
)

plan = scl.plan_analysis(adata, context=context)
print(plan.to_frame())
```

Every referenced key must be an actual `adata.obs` column. Do not use a file
name, Cell Ranger library ID, or sequencing run as `sample_key` unless it truly
identifies the biological sample.

## Field Semantics

| Field | Meaning | Typical value |
|---|---|---|
| `dataset_type` | Biological dataset class | `tumor_tissue`, `pbmc_or_blood`, `cell_line` |
| `sample_key` | Biological sample represented by cells | `sample_id` |
| `batch_key` | Technical processing batch | `library_batch` |
| `condition_key` | Group or condition to compare | `treatment_group` |
| `experimental_unit_key` | Independent unit used for inference | `patient_id`, `donor_id`, `mouse_id` |
| `paired_key` | Unit linking repeated/paired observations | `patient_id` |
| `study_objective` | Plain-language analysis goal | `paired treatment response` |
| `cell_type_key` | Existing annotation column, when available | `cell_type_final` |

`sample_key` and `experimental_unit_key` are often, but not always, identical.
For a pre/post patient study, `sample_key` may identify each biopsy while
`experimental_unit_key` and `paired_key` identify the patient.

## Treatment-Response Project Pattern

AK112/LPJ-like projects should explicitly separate the biopsy, patient,
condition, and technical batch:

```python
context = scl.ProjectContext(
    dataset_type="tumor_tissue",
    sample_key="biopsy_id",
    condition_key="timepoint_or_response",
    batch_key="sequencing_batch",
    experimental_unit_key="patient_id",
    paired_key="patient_id",
    study_objective="treatment response and tumor ecosystem change",
)
```

Before accepting the plan, inspect:

```python
adata.obs[[
    "biopsy_id",
    "patient_id",
    "timepoint_or_response",
    "sequencing_batch",
]].drop_duplicates().sort_values(["patient_id", "biopsy_id"])
```

The table should make the number of independent patients, missing pairs, and
batch-condition confounding immediately visible.

## What scLucid Infers Conservatively

When fields are omitted, scLucid may detect common column names such as
`sample`, `patient`, `donor`, `condition`, and `batch`. Detected values are
recorded as assumptions and must still be confirmed.

- repeated conditions within the same experimental unit are flagged as paired
- multiple batches trigger integration review, not automatic integration
- a comparative objective without biological sample metadata is `BLOCKED`
- fewer than two experimental units in a condition is not treated as formal
  replicated inference
- tumor context adds a separate multi-evidence malignancy review boundary

## Required Handoff Metadata

Before formal DE or composition testing, confirm at minimum:

1. raw counts are preserved in `adata.layers["counts"]`;
2. `sample_key` identifies biological samples;
3. `experimental_unit_key` matches the independent replicate;
4. `condition_key` encodes the intended contrast;
5. `paired_key` is explicit for repeated or paired observations;
6. batch is not silently substituted for sample or condition;
7. patient/sample identifiers remain unchanged through subsetting and merging.

The smallest defensible next step is to correct this metadata table before
tuning QC, integration, clustering, or DE parameters.

`run_pipeline(plan=plan)` rejects a `BLOCKED` plan by default. For an explicitly
exploratory or prerequisite-only run, `allow_blocked_plan=True` can bypass that
guard, but the blocker remains recorded and downstream claims stay limited.
