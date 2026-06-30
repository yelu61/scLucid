Cell Proportion Analysis Method Selection Guide
================================================

This guide explains how to choose and use the three cell-proportion analysis
methods available in scLucid:

- **Pseudo-bulk / CLR**: aggregate to biological-sample level; prefer
  CLR-transformed compositional data for testing.
- **Pseudo-bulk / covariate-aware logCPM**: use a sample-level linear model
  when batch, patient, or paired designs are present.
- **scCODA**: optional Bayesian compositional-data-analysis backend.
- **Milo**: cell-level neighborhood-based analysis; still a planned interface.

Quick Start
-----------

.. warning::

   Cell proportions are compositional data. Raw-proportion ``t-test`` /
   ``wilcoxon`` tests are kept only as legacy exploratory paths in scLucid.
   For formal inference, prefer sample-level CLR, DESeq2-style count models,
   or covariate-aware sample-level models.

Automatic Method Recommendation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest way to choose a method is to let the system recommend one:

```python
from scLucid.analysis import analyze_celltype_proportion

# Recommend and analyze automatically.
result = analyze_celltype_proportion(
    adata,
    sample_col="sample_id",
    condition_col="condition",
)

# The system prints the recommendation and its rationale, e.g.:
# INFO: recommended method: sccoda
# INFO: reason: n_samples=3 < 5, batch_effect=True
```

Inspect the Recommendation Without Running
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

```python
from scLucid.analysis import recommend_method

method = recommend_method(
    adata,
    sample_col="sample",
    condition_col="condition",
)

print(f"Recommended method: {method.value}")
# Output: 'pseudobulk' or 'sccoda' (Milo is still a planned interface)
```

Compare All Methods
^^^^^^^^^^^^^^^^^^^

```python
from scLucid.analysis import compare_methods

comparison = compare_methods(adata)
print(comparison[["method", "overall_score", "recommendation"]])

# Example output:
#          method  overall_score  recommendation
#      pseudobulk           0.80  strongly recommended
#          sccoda           0.65  usable
#            milo           0.45  not recommended
```

---

Method Details
--------------

1. Pseudo-bulk Methods
^^^^^^^^^^^^^^^^^^^^^^

**Principle**: aggregate to biological-sample level. By default, use a
compositional-aware CLR test; when batch or patient design factors exist, use
``linear_model_logcpm`` for sample-level covariate modeling.

``composition_transform(method="clr")`` first checks whether input is already a
closed composition (rows sum to ~1), percentages (~100), raw counts, or a 0–1
subcomposition that is not closed. Except for already-closed compositions, all
non-negative inputs are row-closed per sample before CLR. Negative or non-finite
values raise an error instead of silently producing invalid compositional
results.

**Strengths**:

- Mature and widely accepted in the literature.
- High statistical power from sample-level aggregation.
- Avoids treating individual cells as independent biological replicates.
- Supports CLR, paired/batch-aware paths, and FDR.
- Easy to interpret.

**Weaknesses**:

- Ignores cell-to-cell heterogeneity.
- Loses single-cell resolution.
- Raw-proportion tests are not suitable for formal compositional inference.

**When to use**:

- Biological replicates are available per group.
- ``sample``, ``condition``, and optional ``batch``/``patient`` metadata exist.
- Cell-type annotations are complete.
- The focus is cell-type-level abundance change.

**Example**:

```python
from scLucid.analysis import analyze_celltype_proportion, ProportionConfig

config = ProportionConfig(
    test_method="clr-t-test",  # recommended: sample-level CLR test
    plot_types=["bar", "box", "volcano"],
    out_dir="./results",
)

prop_df, stat_df = analyze_celltype_proportion(
    adata,
    method="pseudobulk",
    config=config,
)

print(stat_df[stat_df["padj"] < 0.05])  # significant cell types

# Result tables include audit fields such as inference_level and
# compositional_data_warning.
```

Pseudo-bulk Conditional DE: Batch/Patient-Aware Sample-Level Models
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the goal is gene-level condition DE rather than cell-proportion change, use
``run_pseudobulk_de``. If batch or patient blocking is present, prefer the pure
Python ``linear_model_logcpm`` path:

```python
from scLucid.analysis import PseudobulkDEConfig, run_pseudobulk_de

config = PseudobulkDEConfig(
    sample_col="sampleID",
    condition_key="condition",
    groupby="cell_type_final",
    group_names=["T cell", "B cell"],
    contrasts=[("control", "treated")],
    method="linear_model_logcpm",
    design_covariates=["batch"],
    block_col="patient_id",  # optional; added as a categorical covariate
    min_cells_per_sample=10,
)

de_df = run_pseudobulk_de(adata, config)

# Default robust_cov_type="HC3" uses heteroskedasticity-robust standard errors.
# To reproduce ordinary OLS standard errors, set robust_cov_type="nonrobust".

# Check formal results for:
# - inference_level == "sample_level"
# - valid_for_publication_inference == True
# - design_formula / design_covariates / block_col / covariance_type
```

When only one biological sample is available per group, scLucid returns
``descriptive_single_sample`` effect-size-only results without formal p-values.
Forcing a cell-level fallback is flagged as ``exploratory_cell_level``.

---

2. scCODA Methods
^^^^^^^^^^^^^^^^^

**Principle**: Bayesian compositional model designed for compositional abundance
changes. In scLucid this path is an optional backend; confirm dependencies,
reference cell type, and sampling parameters before using it.

**Strengths**:

- Models batch effects.
- Works for small samples (N < 5).
- Provides credible intervals.
- Friendly to multi-condition comparisons.

**Weaknesses**:

- MCMC sampling is slower.
- Bayesian model tuning is more complex.
- Newer method with less literature uptake.

**When to use**:

- Fewer than 5 samples per group.
- Batch effects are present.
- Bayesian credible intervals are desired.
- Multi-condition comparisons are needed.

**Example**:

```python
from scLucid.analysis import analyze_celltype_proportion

adata_result = analyze_celltype_proportion(
    adata,
    method="sccoda",
    reference_cell_type="T_cells",
    reference_level="control",
    n_samples=25000,  # MCMC samples
    out_dir="./results",
)

sccoda_results = adata_result.uns["sclucid"]["sccoda"]
print(sccoda_results["final_results"])
```

**scCODA-specific functions**:

```python
from scLucid.tools import (
    run_sccoda,
    summarize_sccoda,
    plot_sccoda_proportion_with_significance,
)

adata = run_sccoda(
    adata,
    cell_type_col="cell_type",
    sample_col="sample_id",
    condition_col="condition",
)

summary = summarize_sccoda(adata)
plot_sccoda_proportion_with_significance(
    adata,
    condition="condition",
    save_path="./sccoda_plot.pdf",
)
```

---

3. Milo Method
^^^^^^^^^^^^^^

**Status**: **Not yet implemented** (planned). Do not rely on this interface for
production analysis.

**Principle**: define neighborhoods in UMAP/PCA space and test for changes in
neighborhood cell composition.

**Strengths**:

- Preserves single-cell resolution.
- Detects subpopulation-level changes.
- Does not require predefined cell types.
- Visualizes spatial patterns.

**Weaknesses**:

- Computationally expensive.
- Requires tuning (neighborhood size).
- Results are harder to interpret.

**When to use** (future):

- Subpopulation changes are of interest.
- Cell-type annotations may be incomplete.
- Spatial distribution patterns matter.
- Discovery of new subpopulations is a goal.

**Planned example**:

```python
# Not yet implemented
from scLucid.analysis import analyze_celltype_proportion

adata_result = analyze_celltype_proportion(
    adata,
    method="milo",
    n_neighbors=30,
    n_pcs=30,
    alpha=0.1,
    out_dir="./results",
)
```

---

Method Comparison
-----------------

| Feature | Pseudo-bulk | scCODA | Milo |
|---------|-------------|--------|------|
| Sample requirement | biological replicates | small samples OK | N ≥ 3 / group |
| Batch effects | ✅ explicit ``linear_model_logcpm`` | ✅ modeled | ⚠️ partial |
| Single-cell resolution | ❌ no | ❌ no | ✅ yes |
| Speed | ⚡⚡⚡ fast | ⚡⚡ medium | ⚡ slow |
| Statistical power | ✅ high | ⚠️ medium | ⚠️ medium |
| Interpretability | ✅ simple | ⚠️ complex | ⚠️ complex |
| Literature acceptance | ✅ high | ⚠️ medium | ⚠️ medium |
| Maturity | ✅ mature | ⚠️ newer | ⚠️ newer |

---

Workflow Examples
-----------------

Workflow 1: Automated Analysis (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Best when you are unsure which method to use.

```python
from scLucid.analysis import analyze_celltype_proportion, recommend_method

# Step 1: inspect recommendation
method = recommend_method(adata)
print(f"Recommended method: {method.value}")

# Step 2: run analysis using the recommended method
result = analyze_celltype_proportion(adata)

# Step 3: extract results
if isinstance(result, tuple):
    prop_df, stat_df = result
    sig_celltypes = stat_df[stat_df["padj"] < 0.05]
else:  # scCODA returns AnnData
    sccoda_results = result.uns["sclucid"]["sccoda"]
```

Workflow 2: Cross-Method Validation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Best when you need to verify result consistency.

```python
from scLucid.analysis import analyze_all_methods

results = analyze_all_methods(
    adata,
    methods=["pseudobulk", "sccoda"],
    out_dir="./comparison",
    compare=True,
)

# Comparison report is saved to ./comparison/method_comparison.csv
pb_prop, pb_stat = results["pseudobulk"]
adata_sccoda = results["sccoda"]

import matplotlib.pyplot as plt
plt.scatter(pb_stat["pval"], sccoda_stat["pval"])
plt.xlabel("Pseudo-bulk p-value")
plt.ylabel("scCODA p-value")
plt.savefig("./comparison/pval_correlation.pdf")
```

Workflow 3: Large-Sample Standard Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Best for large cohorts (N ≥ 10 / group) with clear batch/paired information.

```python
from scLucid.analysis import analyze_celltype_proportion, ProportionConfig

config = ProportionConfig(
    test_method="clr-t-test",
    plot_types=["bar", "box", "heatmap", "volcano"],
)

prop_df, stat_df = analyze_celltype_proportion(
    adata,
    method="pseudobulk",
    config=config,
)
```

Workflow 4: Small-Sample Batch Correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Best for small cohorts (N < 5 / group) with batch effects.

```python
from scLucid.analysis import analyze_celltype_proportion

adata_result = analyze_celltype_proportion(
    adata,
    method="sccoda",
    reference_cell_type="T_cells",
    reference_level="control",
    n_samples=25000,
    n_burnin=5000,
)
```

---

FAQ
---

Q1: Which method should I choose?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the automatic recommendation:

```python
from scLucid.analysis import recommend_method, compare_methods

method = recommend_method(adata)
comparison = compare_methods(adata)
print(comparison[["method", "overall_score", "recommendation"]])
```

Q2: Can I use multiple methods at once?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes. Use ``analyze_all_methods()``:

```python
from scLucid.analysis import analyze_all_methods

results = analyze_all_methods(
    adata,
    methods=["pseudobulk", "sccoda"],
    out_dir="./comparison",
    compare=True,
)
```

Q3: What if pseudo-bulk and scCODA disagree?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is normal because the methods detect different types of changes:

- **Pseudo-bulk / CLR**: sample-level compositional change.
- **scCODA**: compositional change with compositional constraints.

Recommendations:

1. Prioritize changes significant in both methods (higher confidence).
2. Use biological knowledge to judge which result is more plausible.
3. Consider sample size and batch effects.

Q4: When will Milo be implemented?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Milo is on the development roadmap. Current recommendations:

- **Short term**: use pseudo-bulk + manual subpopulation analysis.
- **Medium term**: use clustering + Milo (requires separate implementation).
- **Long term**: integrated into the unified scLucid interface.

You can open a GitHub issue to raise the priority of Milo support.

---

References
----------

1. **DESeq2**: Love, Huber, and Anders. *Genome Biology* 2014.
2. **scCODA**: Büttner et al. *Nature Methods* 2021.
3. **Milo**: Dann et al. *Nature Methods* 2022.
4. **scCODA tutorial**: https://github.com/theislab/scCODA
5. **Milo tutorial**: https://github.com/MarioniLab/milo

---

Feedback
--------

For questions or suggestions:

- Open a GitHub issue.
- See the scLucid documentation.
- Contact the development team.
