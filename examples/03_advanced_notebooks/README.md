# scLucid Advanced Notebooks

This directory contains **complete advanced notebook workflows**. These are
project-style templates, not short examples.

The intended split in this repository is:

- `docs/`: stable explanations, API reference, recommended defaults
- `examples/01_workflow/`: one-call workflow scripts
- `examples/02_simple_api/`: composable stage-level scripts
- `examples/03_advanced_notebooks/`: full analysis narratives with richer intermediate results

For a publication-oriented package, notebooks should be used to show:

- real-data or project-style analyses
- full QC → preprocess → analysis flows
- reviewer-facing outputs and interpretation checkpoints
- longer result narratives that would be too verbose for `docs/` or `examples/`
- module maturity checks, compact audit summaries, and step-level evidence when
  a notebook implements a benchmark-grade module path

## Current Notebook Set

- `Step1A-QC_Audit.ipynb` - benchmark-grade QC audit notebook. It starts from
  raw combined data and writes `data/processed/Step1-sce_cleaned.h5ad`.
- `Step1B-Preprocessing_Audit.ipynb` - benchmark-grade preprocessing audit
  notebook. It starts from `Step1-sce_cleaned.h5ad` and writes
  `data/processed/Step2-sce_preprocessed.h5ad`.
- `Step2-Annotation_and_Malignancy.ipynb` - evidence-first analysis acceptance
  shell for clustering review, marker/CellTypist/LLM annotation evidence,
  consensus labels, optional malignancy interpretation, and reviewable artifacts.
  It starts from `Step2-sce_preprocessed.h5ad`, calls
  `scripts/run_analysis_acceptance.py`, and writes
  `data/processed/Step3-sce_annotated.h5ad`.
- `Step3-Standard_Downstream.ipynb` - standard downstream composition,
  proportion, differential expression, and enrichment analyses. It starts from
  `Step3-sce_annotated.h5ad`.
- `Step4-Signature_and_Target_Analysis.ipynb` - project-specific signatures,
  focused cell-state analysis, and target-oriented exports. It starts from
  `Step3-sce_annotated.h5ad`.
- `Step1-QC_and_Preprocessing.ipynb` - legacy unsplit QC + preprocessing
  reference retained for comparison.
- `Step2-Celltype_annotation.ipynb` - legacy unsplit project notebook retained
  for comparison.
- `04_advanced_topics.ipynb`
- `04_differential_expression.ipynb`
- `05_trajectory_inference.ipynb`
- `06_advanced_tools.ipynb`

## Recommended Run Order

1. `Step1A-QC_Audit.ipynb`
2. `Step1B-Preprocessing_Audit.ipynb`
3. `Step2-Annotation_and_Malignancy.ipynb`
4. `Step3-Standard_Downstream.ipynb`
5. `Step4-Signature_and_Target_Analysis.ipynb`

## Recommended Maintenance Rules

- keep notebooks narrative and result-oriented
- keep package policy and recommended defaults in `docs/`, not only in notebooks
- keep `examples/` short and script-like, and reserve notebooks for deeper walkthroughs
- use real or representative datasets, and make the expected inputs explicit near the top of each notebook
- when a notebook bypasses one-call workflow functions, still write the same
  `adata.uns["sclucid"]` review contracts used by the package workflow layer
- keep Step2 synchronized with `scripts/run_analysis_acceptance.py`; the notebook
  should inspect acceptance artifacts rather than reimplement the workflow
