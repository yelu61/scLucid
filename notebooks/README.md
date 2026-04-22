# scLucid Notebooks

This directory is reserved for **complete notebook workflows**.

The intended split in this repository is:

- `docs/`: stable explanations, API reference, recommended defaults
- `examples/`: short runnable scripts
- `notebooks/`: full analysis narratives with richer intermediate results

For a publication-oriented package, notebooks should be used to show:

- real-data or project-style analyses
- full QC → preprocess → analysis flows
- reviewer-facing outputs and interpretation checkpoints
- longer result narratives that would be too verbose for `docs/` or `examples/`

## Current Notebook Set

- `01_quality_control.ipynb`
- `02_preprocessing.ipynb`
- `03_clustering_annotation.ipynb`
- `04_advanced_topics.ipynb`
- `04_differential_expression.ipynb`
- `05_trajectory_inference.ipynb`
- `06_advanced_tools.ipynb`

## Recommended Maintenance Rules

- keep notebooks narrative and result-oriented
- keep package policy and recommended defaults in `docs/`, not only in notebooks
- keep `examples/` short and script-like, and reserve notebooks for deeper walkthroughs
- use real or representative datasets, and make the expected inputs explicit near the top of each notebook
