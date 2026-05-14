# scLucid Examples

This directory contains runnable scripts organized by **analysis layer**.

scLucid is designed around three user-facing layers. Choose the one that matches your needs:

## 01_workflow/ — One-Line Analysis

For **beginners and standard projects**. Load data, configure, run.

| Script | What it shows |
|--------|---------------|
| `basic_pipeline.py` | Minimal end-to-end: QC → preprocess → cluster → annotate in one call |
| `prepare_data.py` | Loading 10x data and attaching metadata |
| `plugin_development.py` | How to extend scLucid with custom plugins |

**When to use**: You have a standard dataset and want results fast.
**API**: `scl.run_pipeline()`, `scl.run_standard_qc()`, `scl.run_preprocessing()`

## 02_simple_api/ — Composable Steps

For **analysts who need control**. Inspect, tweak, or replace individual stages.

| Script | What it shows |
|--------|---------------|
| `preprocess_step_by_step.py` | Manual normalization → HVG → scaling → PCA → integration |
| `qc_preprocess_review.py` | Stage-level QC + preprocessing with review-summary inspection |
| `intelligent_qc.py` | Data-driven threshold recommendations with confidence intervals |
| `intelligent_preprocess.py` | Smart parameter selection with review summaries |
| `annotation_workflow.py` | Customizable cell-type annotation pipeline |
| `annotation_report.py` | Export reviewer-facing reports |
| `qc_evaluation.py` | Evaluate QC decisions with benchmarks |

**When to use**: You want to understand what each step does and adjust parameters.
**API**: `scl.qc.calculate_qc_metric()`, `scl.pp.normalize_data()`, `scl.pl.plot_embedding()`, etc.

## 03_advanced_notebooks/ — Full Transparency

For **real exploratory projects** where every decision must be auditable.

Use the split advanced sequence when presenting a real project-style analysis:

| Notebook | What it shows |
|----------|---------------|
| `Step1A-QC_Audit.ipynb` | QC benchmark path, threshold evidence, module maturity, and `Step1-sce_cleaned.h5ad` |
| `Step1B-Preprocessing_Audit.ipynb` | QC handoff, layer audit, preprocessing parameter/layer evidence, and `Step2-sce_preprocessed.h5ad` |
| `Step2-Annotation_and_Malignancy.ipynb` | Clustering, annotation, malignancy review, CNV-aware interpretation, and `Step3-sce_annotated.h5ad` |
| `Step3-Standard_Downstream.ipynb` | Composition, proportion, differential expression, and enrichment |
| `Step4-Signature_and_Target_Analysis.ipynb` | Project-specific signatures, focused cell states, and target-oriented exports |

The legacy unsplit `Step1-QC_and_Preprocessing.ipynb` and
`Step2-Celltype_annotation.ipynb` are retained as project references, but the
split sequence is the recommended product-facing demonstration.

**When to use**: You are doing research where every threshold, diagnostic, and override must remain visible and reviewable.
**Format**: Jupyter notebooks with step-by-step parameter blocks, decision-support tools, and audit trails.

## 04_publication_figures/ — Journal-Ready Figures

For **manuscript figure preparation**. Self-contained scripts that produce one publication-quality PDF each, with TrueType-embedded fonts (`pdf.fonttype=42`) so every label can be edited in Illustrator before submission.

| Script | Figure type |
|--------|-------------|
| `01_umap_annotation.py` | UMAP scatter colored by cell type — Nature-themed |
| `02_marker_heatmap.py` | Per-cell-type marker expression heatmap (z-score) |
| `03_volcano_de.py` | Differential expression volcano with top-hit labels |
| `04_cnv_heatmap.py` | Chromosome-ordered CNV profile by cell group |

**When to use**: You have analyzed data and need camera-ready figures. Copy a script, swap in your AnnData, tweak palette / sizes / labels.
**API**: `scl.pl.plot_embedding`, `scl.pl.plot_marker_heatmap`, `scl.pl.plot_volcano`, plus direct matplotlib for custom layouts.

## Scope Rules

Examples should stay:

- short
- runnable
- scenario-based
- aligned with the package's documented defaults

Examples should not become the place where package policy is defined.
If the recommended workflow changes, update `docs/source/quickstart.rst` and
`docs/source/best_practices.rst` first, then keep examples consistent with them.
