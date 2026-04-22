# scLucid Examples

This directory contains short runnable scripts.

Use these files when you want a template to adapt quickly.
Use `docs/` when you want the official recommended path and API guidance.
Use `notebooks/` when you want a full analysis walkthrough with richer outputs.

## Current Scripts

`quickstart.py`
- Minimal end-to-end path: QC → preprocessing → clustering → annotation.
- Best first script to copy for a new project.

`preprocessing.py`
- Shows three preprocessing modes:
  - standard default path
  - intelligent preprocessing with review summary
  - manual step-by-step path

`intelligent_qc_example.py`
- Focused demonstration of intelligent QC recommendations.
- Useful when reviewing threshold recommendation behavior.

`intelligent_preprocessing_example.py`
- Focused demonstration of intelligent preprocessing recommendations and summaries.

`annotation_report.py`
- Shows how to export reviewer-facing annotation reports and sidecars.

`curated_annotation_workflow.py`
- Shows how to replace long notebook-specific annotation code with reusable APIs.

## Scope Rules

Examples should stay:

- short
- runnable
- scenario-based
- aligned with the package's documented defaults

Examples should not become the place where package policy is defined.
If the recommended workflow changes, update `docs/source/quickstart.rst` and
`docs/source/best_practices.rst` first, then keep examples consistent with them.
