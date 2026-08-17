# Cell Proportion Analysis Method Guide

This guide explains how scLucid frames cell-proportion analysis. The central
rule is simple: formal condition inference should happen at the biological
sample level whenever replicates are available.

!!! warning
    Cell proportions are compositional data. Raw-proportion `t-test` or
    `wilcoxon` paths are legacy exploratory tools. For publication-grade
    inference, prefer sample-level CLR, covariate-aware pseudobulk models, or
    a validated compositional backend.

## Method Selection

| Method | Best use | Inference level | Status |
|---|---|---|---|
| Pseudobulk / CLR | Cell-type abundance testing with biological replicates | sample-level | recommended default |
| Pseudobulk / logCPM model | Gene-level condition DE with batch, patient, or paired design terms | sample-level | recommended for condition DE |
| scCODA | Small-sample Bayesian compositional analysis | sample-level compositional model | optional backend |
| Milo | Neighborhood-level abundance discovery | cell-neighborhood exploratory | planned interface |

## Automatic Recommendation

```python
from scLucid.analysis import analyze_celltype_proportion

result = analyze_celltype_proportion(
    adata,
    sample_col="sample_id",
    condition_col="condition",
)
```

To inspect the recommendation before running:

```python
from scLucid.analysis import recommend_method

method = recommend_method(
    adata,
    sample_col="sample_id",
    condition_col="condition",
)

print(method.value)
```

## Pseudobulk Cell-Proportion Testing

Use this path when each condition has biological samples and the question is
whether a cell type changes in abundance.

```python
from scLucid.analysis import ProportionConfig, analyze_celltype_proportion

config = ProportionConfig(
    test_method="clr-t-test",
    plot_types=["bar", "box", "volcano"],
    out_dir="./results",
)

prop_df, stat_df = analyze_celltype_proportion(
    adata,
    method="pseudobulk",
    sample_col="sample_id",
    condition_col="condition",
    cell_type_col="cell_type",
    config=config,
)
```

Review result fields such as `inference_level`,
`compositional_data_warning`, and adjusted p-values before making biological
claims.

## Pseudobulk Condition DE

When the question is gene-level condition DE, use `run_pseudobulk_de`. This
keeps formal inference at the sample level and avoids treating cells as
independent biological replicates.

```python
from scLucid.analysis import PseudobulkDEConfig, run_pseudobulk_de

config = PseudobulkDEConfig(
    sample_col="donor_condition_sample",
    condition_key="condition",
    experimental_unit_col="patient_id",
    groupby="cell_type_final",
    group_names=["T cell", "B cell"],
    contrasts=[("control", "treated")],
    method="linear_model_logcpm",
    design_covariates=["batch"],
    block_col="patient_id",
    min_cells_per_sample=10,
)

de_df = run_pseudobulk_de(adata, config)
```

Formal results should report `inference_level == "sample_level"` and
`valid_for_publication_inference == True`. If only one biological sample exists
per group, scLucid should return descriptive effect-size results rather than
formal p-values.

`sample_col` identifies the aggregate pseudobulk row; it is not automatically
the biological replicate. For paired or repeated designs, set both
`experimental_unit_col` and `block_col`. scLucid rejects unobserved contrast
levels, repeated units without a block, unresolved technical replicates, and
covariates that the selected backend ignores or cannot identify.

## scCODA

Use scCODA only when the optional backend is installed and you have checked the
reference cell type, model assumptions, and sampling settings.

```python
from scLucid.analysis import analyze_celltype_proportion

adata_result = analyze_celltype_proportion(
    adata,
    method="sccoda",
    sample_col="sample_id",
    condition_col="condition",
    cell_type_col="cell_type",
    reference_cell_type="T_cells",
    reference_level="control",
    n_samples=25000,
    out_dir="./results",
)
```

## Interpretation Boundary

- Use pseudobulk or compositional sample-level methods for formal condition
  claims.
- Treat cell-level comparisons as exploratory marker discovery unless a method
  explicitly models sample-level replication.
- Record method choice, sample counts, design covariates, confidence, and
  review-required status in the analysis review summary.
