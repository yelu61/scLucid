# scLucid Examples

This directory contains runnable scripts and project-style notebooks organized by
**analysis layer**. The examples are meant to demonstrate how the same scLucid
contracts can serve users with different levels of control.

The layers should not become separate products. Workflow scripts, simple API
scripts, and advanced notebooks must all preserve the same core handoff:

```python
adata.uns["sclucid"][module]["workflow_config"]
adata.uns["sclucid"][module]["steps_executed"]
adata.uns["sclucid"][module]["review_summary"]
```

## Recommended Path Through the Examples

| User need | Start here | Why |
|---|---|---|
| Recommended first screen | `02_simple_api/qc_preprocess_review.py` | Uses the four-action read-only review → explicit apply contract. |
| Compatibility pipeline | `01_workflow/basic_pipeline.py` | Preserves the previous QC → preprocessing → analysis entrypoint during migration. |
| Fully auditable project notebook | `03_advanced_notebooks/Step1A-QC_Audit.ipynb` then `Step1B-Preprocessing_Audit.ipynb` | Shows manual audit checkpoints while still writing scLucid review contracts. |
| Evidence-first annotation | `03_advanced_notebooks/Step2-Annotation_and_Malignancy.ipynb` | Delegates the analysis acceptance path to `scripts/run_analysis_acceptance.py`. |

## 00_data_io/ — Input Preparation

For users who need to turn raw files into a scLucid-ready AnnData before running
any workflow.

| Script | Current role |
|---|---|
| `prepare_data.py` | Data I/O reference for 10x inputs and existing `.h5ad` objects. It is not the canonical workflow example and should not define package policy. |

## 01_workflow/ — One-Call Baselines

For **beginners and standard projects**. Load data, configure, run.

| Script | Current role |
|---|---|
| `basic_pipeline.py` | Compatibility minimal end-to-end example: QC → preprocess → cluster → annotate in one call. |
| `plugin_development.py` | Extension reference. It should not imply arbitrary registered plugins are executed by `run_custom_analysis()` unless the package implements that registry bridge. |

**When to use**: You have a standard dataset and want a reproducible first pass.
**Primary APIs**: `scl.run_pipeline()`, `scl.qc.run_standard_qc()`, `scl.pp.run_preprocessing()`, `scl.al.run_standard_analysis()`.

Workflow examples use the maintained light-dependency QC/preprocess path. R
bridges, ScDblFinder wrappers, CellBender/SoupX/DecontX, and project-level
ambient RNA correction are not part of these defaults.

## 02_simple_api/ — Composable, Reviewable Steps

For **analysts who need control**. Inspect, tweak, or replace individual stages.

| Script | Current role |
|---|---|
| `qc_preprocess_review.py` | Canonical DecisionCard/QCPolicy/PreprocessPolicy/RunEvidence example. |
| `qc_step_by_step.py` | Manual QC teaching example. It should be paired with a manual review finalizer before being used as a project template. |
| `preprocess_step_by_step.py` | Manual preprocessing teaching example. It preserves full-gene `counts`/`normalized_full`, limits scaling to the temporary discovery matrix, and writes a review contract. |
| `intelligent_qc.py` | Read-only DecisionCard context-sensitivity example; candidate impacts are not a scalar quality score. |
| `intelligent_preprocess.py` | Read-only PreprocessPolicy example with optional explicit application and RunEvidence. |
| `annotation_workflow.py` | Curated annotation recipe for marker/enrichment evidence, manual mapping, module scoring, and composition plots. |
| `annotation_report.py` | Reviewer-facing annotation report export. |
| `qc_evaluation.py` | QC decision evaluation and benchmark-style reporting. |

**When to use**: You want to understand what each step does and adjust parameters.
**Primary policy APIs**: `scl.recommend_qc_policy()`, `scl.apply_qc_policy()`, `scl.recommend_preprocess_policy()`, `scl.apply_preprocess_policy()`.

The low-level teaching scripts additionally show maintained composable calls
such as `calculate_qc_metric()`, `filter_cells()`, `normalize_data()`, and
`find_hvgs()`, but they do not define product defaults.

Manual simple-API scripts must not stop at producing modified AnnData. If a
manual path is promoted to a project template, it should finalize the same review
contract as the workflow layer. Until package-level finalizers exist, use the
advanced notebooks as the reference for how to normalize and export manual
review summaries.

Optional enhancements such as scran, Harmony, scVI/scANVI, BBKNN, SOLO, or
DoubletDetection should be enabled only when their dependencies and biological
rationale are explicit.

## 03_advanced_notebooks/ — Full Transparency

For **real exploratory projects** where every decision must be auditable.

Use the split advanced sequence when presenting a real project-style analysis:

| Notebook | What it shows |
|---|---|
| `Step1A-QC_Audit.ipynb` | QC benchmark path, threshold evidence, module maturity, and `Step1-sce_cleaned.h5ad`. |
| `Step1B-Preprocessing_Audit.ipynb` | QC handoff, layer audit, preprocessing parameter/layer evidence, and `Step2-sce_preprocessed.h5ad`. |
| `Step2-Annotation_and_Malignancy.ipynb` | Evidence-first analysis acceptance via `scripts/run_analysis_acceptance.py`: clustering review, annotation evidence, consensus labels, post-hoc QC cluster review, optional malignancy interpretation, and `Step3-sce_annotated.h5ad`. |
| `Step3-Standard_Downstream.ipynb` | Sample-level composition/proportion, pseudobulk differential expression, covariate-aware downstream inference, and enrichment. |
| `Step4-Signature_and_Target_Analysis.ipynb` | Project-specific signatures, focused cell states, and target-oriented exports. |

The legacy unsplit notebooks are retained as references, not as the
product-facing path:

- `Step1-QC_and_Preprocessing.ipynb`
- `Step2-Celltype_annotation.ipynb`

QC and preprocessing are `REVIEW` layers, not scientifically locked `CORE`
modules. Analysis/Tumor feature development is frozen until their acceptance
gates pass. Doublet-heavy,
high-mitochondrial, stress-high, and low tumor-purity signals should first be
surfaced as review evidence, then acted on after project-specific manual
confirmation.

Downstream examples must keep inference levels explicit:

- cell-level marker discovery is for annotation and exploratory biology, not
  formal condition DE
- publication-grade condition DE should use `run_pseudobulk_de()` with
  biological samples as replicates
- when batch, donor, patient, or paired structure exists, use
  `method="linear_model_logcpm"` with `design_covariates` / `block_col`
- raw proportion `t-test` / `wilcoxon` outputs are legacy exploratory summaries;
  prefer CLR sample-level proportion tests or a compositional backend
- ambient RNA and empty-droplet checks in examples are diagnostics only unless a
  project explicitly registers an external correction result

**When to use**: You are doing research where every threshold, diagnostic, and
override must remain visible and reviewable.

## 04_publication_figures/ — Journal-Ready Figures

For **manuscript figure preparation**. Self-contained scripts that produce one
publication-quality PDF each, with TrueType-embedded fonts (`pdf.fonttype=42`) so
every label can be edited in Illustrator before submission.

| Script | Figure type |
|---|---|
| `01_umap_annotation.py` | UMAP scatter colored by cell type — Nature-themed. |
| `02_marker_heatmap.py` | Per-cell-type marker expression heatmap. |
| `03_volcano_de.py` | Differential expression volcano with top-hit labels. |
| `04_cnv_heatmap.py` | Chromosome-ordered CNV profile by cell group. |

## Scope Rules

Examples should stay:

- short, runnable, and scenario-based
- aligned with package defaults documented in `docs/`
- explicit about required inputs near the top of the file
- clear about whether they are canonical, teaching references, or legacy references

Examples should not become the only place where package policy is defined. If a
recommended workflow changes, update `docs/source/usage_layers.rst` first, then
keep examples consistent with it.

## Maintenance Checklist

Before treating an example as product-facing, verify that it:

- runs against an installed scLucid package without ad hoc `sys.path` edits
- has at least one lightweight smoke test or golden-path test
- writes or inspects the expected scLucid review contract
- does not imply unsupported plugin or registry behavior
- avoids project-specific heavy tools unless explicitly framed as optional
