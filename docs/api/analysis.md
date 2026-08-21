# Analysis API

This page is generated at build time with `mkdocstrings`. The Markdown source is
intentionally small; the rendered documentation below is pulled from the current
Python objects and docstrings.

The stable user entrypoint is `run_standard_analysis`. Lower-level modules are
included for clustering review, annotation, differential expression, proportion
analysis, and analysis review contracts.

## Workflow

::: scLucid.analysis.workflow
    options:
      show_root_heading: true

## Configuration

::: scLucid.analysis.config
    options:
      show_root_heading: true

## Clustering

::: scLucid.analysis.clustering
    options:
      show_root_heading: true

## Annotation

::: scLucid.analysis.annotation
    options:
      show_root_heading: true

## Differential expression and enrichment

Use these public APIs in notebooks instead of redefining large helper functions:

- `PseudobulkDEConfig` + `run_pseudobulk_de` for sample-level pseudobulk DEG.
- `EnrichmentConfig` + `run_enrichment` for ORA/GSEA-style pathway review.
- `visualize_markers`, `plot_grouped_marker_dotplot`, and
  `plot_categorized_gene_heatmap` for marker/program visualization.

::: scLucid.analysis.differential_expression
    options:
      show_root_heading: true

::: scLucid.analysis.enrichment
    options:
      show_root_heading: true

## Composition analysis

Use `ProportionConfig` + `analyze_celltype_proportion` for composition testing.
Use `summarize_composition_shift`, `plot_composition_shift_bubble`, and
`plot_composition_shift_effect` for LVA-style effect display. These plotting
helpers visualize existing sample-level results; they do not perform a new
statistical test.

::: scLucid.analysis.proportion
    options:
      show_root_heading: true
