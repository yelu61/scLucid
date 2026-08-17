# scLucid Advanced Notebooks

This directory contains **complete advanced notebook workflows**. These are
project-style templates, not short examples. They are the reference for users
who need manual checkpoints while preserving the same scLucid contracts used by
workflow scripts.

The intended split in this repository is:

- `docs/`: stable explanations, API reference, recommended defaults, and layer policy
- `examples/01_workflow/`: one-call or stage-level workflow scripts
- `examples/02_simple_api/`: composable stage-level scripts and teaching references
- `examples/03_advanced_notebooks/`: full analysis narratives with richer intermediate review

For a publication-oriented package, notebooks should be used to show:

- real-data or project-style analyses
- full QC → preprocess → analysis flows
- reviewer-facing outputs and interpretation checkpoints
- longer result narratives that would be too verbose for `docs/` or short examples
- module maturity checks, compact audit summaries, and step-level evidence when a notebook implements a benchmark-grade module path

## Recommended Product-Facing Sequence

Use this split sequence for new project templates:

1. `Step1A-QC_Audit.ipynb`
2. `Step1B-Preprocessing_Audit.ipynb`
3. `Step2-Annotation_and_Malignancy.ipynb`
4. `Step3-Standard_Downstream.ipynb`
5. `Step4-Signature_and_Target_Analysis.ipynb`

The split Step1A/Step1B sequence is preferred over a single combined QC +
preprocessing notebook because it creates a clean handoff between quality-control
decisions and preprocessing decisions.

## Current Notebook Set

| Notebook | Current role |
|---|---|
| `Step1A-QC_Audit.ipynb` | Canonical advanced QC audit notebook. It starts from raw combined data and writes `data/processed/Step1-sce_cleaned.h5ad`. |
| `Step1B-Preprocessing_Audit.ipynb` | Canonical advanced preprocessing audit notebook. It starts from `Step1-sce_cleaned.h5ad` and writes `data/processed/Step2-sce_preprocessed.h5ad`. |
| `Step2-Annotation_and_Malignancy.ipynb` | Evidence-first analysis acceptance shell. It calls `scripts/run_analysis_acceptance.py` and writes `data/processed/Step3-sce_annotated.h5ad`. |
| `Step3-Standard_Downstream.ipynb` | Standard downstream composition, proportion, differential expression, and enrichment analyses. |
| `Step4-Signature_and_Target_Analysis.ipynb` | Project-specific signatures, focused cell-state analysis, and target-oriented exports. |

## Manual Review Contract

Advanced notebooks may intentionally bypass one-call workflow functions so users
can inspect thresholds, plots, handoff state, and biological assumptions. That is
valid scLucid usage only if the notebook still writes the same module-level
review contract as the workflow layer:

```python
adata.uns["sclucid"][module]["workflow_config"]
adata.uns["sclucid"][module]["steps_executed"]
adata.uns["sclucid"][module]["review_summary"]
```

Use `scLucid.utils.finalize_manual_review_summary()` when a notebook executes a
stage manually. It normalizes, validates, stores, and optionally exports the
same review contract as the workflow layer. Existing notebook-local finalizers
are migration references and should not be copied into new projects. After the
manual stages are finalized, call `scLucid.review_run()` for one cross-stage
action table.

## Step Boundaries

Step1A and Step1B are stable audit/handoff notebooks. They calculate, document,
and review QC/preprocessing decisions using light-dependency package defaults.
They should not absorb project-specific ambient RNA correction,
CellBender/SoupX/DecontX execution, ScDblFinder, or other R bridge workflows.
Those tools can be used upstream or in project-specific expert appendices when
there is a clear rationale.

Tumor purity, malignancy evidence, CNV interpretation, stress-state biology, and
doublet-heavy cluster interpretation belong in Step2 and the analysis/tumor
modules, where they can be reviewed with annotation and biological context.
Step2 should treat these signals as review evidence first; automatic deletion
belongs only after project-specific manual confirmation.

## Recommended Maintenance Rules

- keep notebooks narrative and result-oriented
- keep package policy and recommended defaults in `docs/`, not only in notebooks
- keep short runnable scripts in `examples/01_workflow/` and `examples/02_simple_api/`
- use real or representative datasets, and make expected inputs explicit near the top of each notebook
- when a notebook bypasses one-call workflow functions, still write the same `adata.uns["sclucid"]` review contracts used by the package workflow layer
- keep Step2 synchronized with `scripts/run_analysis_acceptance.py`; the notebook should inspect acceptance artifacts rather than reimplement the workflow
- mark legacy notebooks as references when they are retained for comparison
- do not let notebook-local helpers drift into de facto package APIs without moving them into `src/scLucid/`
