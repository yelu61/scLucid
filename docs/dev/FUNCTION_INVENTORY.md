# scLucid Public API Inventory

**Regenerate:** `python scripts/audit_public_api.py --write`

This document lists every public symbol in the scLucid API. The auto-generated
section below is maintained by `scripts/audit_public_api.py`. You may add notes
above the `AUTO-GENERATED` markers; they will be preserved across regenerations.

## Legend

| Tag | Meaning |
|-----|---------|
| `[A]` | Alias — points to another symbol or module |
| `[C]` | Config class |
| `[W]` | Workflow orchestrator |
| `[T]` | Trace / contract / schema constant |
| `[D]` | Deprecated |
| `[P]` | Private-but-exposed (starts with `_` but in `__all__`) |
| `[O]` | Optional / depends on extra dependencies |
| `[?]` | Uncertain / could not be traced cleanly |

## Maintainer Notes

- To add a symbol to the public API, add it to `__all__` in the subpackage `__init__.py`
  or register it via `_export()`.
- To deprecate a symbol, emit a `FutureWarning`/`DeprecationWarning`; the audit script
  will detect it automatically.
- To remove a symbol from the public API, remove it from `__all__` and run the audit
  script with `--write`.

<!-- AUTO-GENERATED INVENTORY START -->
<!-- Generated: 2026-06-17 18:54:52 by scripts/audit_public_api.py -->
<!-- Total public symbols: 978 -->

## scLucid

### Stable APIs

#### Workflow Orchestrator

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_annotation` | Workflow Orchestrator | `src/scLucid/__init__.py` | [W] dynamically resolved from submodule workflow |
| `run_custom_analysis` | Workflow Orchestrator | `src/scLucid/__init__.py` | [W] dynamically resolved from submodule workflow |
| `run_pipeline` | Workflow Orchestrator | `src/scLucid/__init__.py` | [W] |
| `run_preprocessing` | Workflow Orchestrator | `src/scLucid/__init__.py` | [W] dynamically resolved from submodule workflow |
| `run_standard_analysis` | Workflow Orchestrator | `src/scLucid/__init__.py` | [W] dynamically resolved from submodule workflow |
| `run_standard_qc` | Workflow Orchestrator | `src/scLucid/__init__.py` | [W] dynamically resolved from submodule workflow |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `AnalysisContext` | Function | `src/scLucid/__init__.py` |  |
| `AnalysisStep` | Function | `src/scLucid/__init__.py` |  |
| `AnalysisStepFactory` | Function | `src/scLucid/__init__.py` |  |
| `CellAnnotator` | Function | `src/scLucid/__init__.py` |  |
| `characterize_clusters` | Function | `src/scLucid/__init__.py` | dynamically resolved from submodule workflow |
| `compact_sclucid_uns` | Function | `src/scLucid/__init__.py` |  |
| `DatasetProfile` | Function | `src/scLucid/__init__.py` |  |
| `export_audit_report` | Function | `src/scLucid/__init__.py` |  |
| `get_config` | Function | `src/scLucid/__init__.py` |  |
| `infer_analysis_context` | Function | `src/scLucid/__init__.py` |  |
| `infer_dataset_profile` | Function | `src/scLucid/__init__.py` |  |
| `is_interactive_mode` | Function | `src/scLucid/__init__.py` |  |
| `normalize_dataset_type` | Function | `src/scLucid/__init__.py` |  |
| `PlottingBackend` | Function | `src/scLucid/__init__.py` |  |
| `ProportionAnalysisMethod` | Function | `src/scLucid/__init__.py` |  |
| `QCFilter` | Function | `src/scLucid/__init__.py` |  |
| `read_10x` | Function | `src/scLucid/__init__.py` |  |
| `read_h5ad` | Function | `src/scLucid/__init__.py` |  |
| `recommend_analysis_parameters` | Function | `src/scLucid/__init__.py` | dynamically resolved from submodule workflow |
| `reset_config` | Function | `src/scLucid/__init__.py` |  |
| `ScoringMethod` | Function | `src/scLucid/__init__.py` |  |
| `set_config` | Function | `src/scLucid/__init__.py` |  |

#### Alias

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `al` | Alias | `src/scLucid/__init__.py` | [A] from `analysis` |
| `analysis` | Alias | `src/scLucid/__init__.py` | [A] module analysis |
| `pl` | Alias | `src/scLucid/__init__.py` | [A] from `plotting` |
| `plotting` | Alias | `src/scLucid/__init__.py` | [A] module plotting |
| `pp` | Alias | `src/scLucid/__init__.py` | [A] from `preprocess` |
| `preprocess` | Alias | `src/scLucid/__init__.py` | [A] module preprocess |
| `qc` | Alias | `src/scLucid/__init__.py` | [A] module qc |
| `rc` | Alias | `src/scLucid/__init__.py` | [A] from `recommendation` |
| `recommendation` | Alias | `src/scLucid/__init__.py` | [A] module recommendation |
| `reset_figure_params` | Alias | `src/scLucid/__init__.py` | [A] from `_settings_unavailable` |
| `set_figure_params` | Alias | `src/scLucid/__init__.py` | [A] from `_settings_unavailable` |
| `set_interactive_mode` | Alias | `src/scLucid/__init__.py` | [A] from `_settings_unavailable` |
| `setup_logging` | Alias | `src/scLucid/__init__.py` | [A] from `_settings_unavailable` |
| `tl` | Alias | `src/scLucid/__init__.py` | [A] from `tools` |
| `tools` | Alias | `src/scLucid/__init__.py` | [A] module tools |
| `ut` | Alias | `src/scLucid/__init__.py` | [A] from `utils` |
| `utils` | Alias | `src/scLucid/__init__.py` | [A] module utils |

#### Constant

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `FONT_CELL` | Constant | `src/scLucid/__init__.py` |  |
| `FONT_NATURE` | Constant | `src/scLucid/__init__.py` |  |
| `FONT_TRADITIONAL` | Constant | `src/scLucid/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

#### Deprecated

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_advanced_qc` | Deprecated | `src/scLucid/__init__.py` | [D] dynamically resolved from submodule workflow |

**Summary:** 49 symbols (48 stable, 1 flagged). workflow=6, config=0, class=0, function=22, alias=17, constant=3, trace=0, deprecated=1, uncertain=0, private_but_exposed=0.

## scLucid.analysis

### Stable APIs

#### Workflow Orchestrator

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_annotation` | Workflow Orchestrator | `src/scLucid/analysis/__init__.py` | [W] from `annotation` |
| `run_custom_analysis` | Workflow Orchestrator | `src/scLucid/analysis/__init__.py` | [W] from `workflow` |
| `run_standard_analysis` | Workflow Orchestrator | `src/scLucid/analysis/__init__.py` | [W] from `workflow` |

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `AnalysisWorkflowConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `AnnotationConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `BulkAbundanceConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `bulk`; optional |
| `BulkClinicalAssociationConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `bulk`; optional |
| `BulkDEConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `bulk`; optional |
| `BulkDeconvolutionConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `bulk`; optional |
| `BulkDiagnosticsConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `bulk`; optional |
| `BulkNormalizationConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `bulk`; optional |
| `BulkTraitAssociationConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `bulk`; optional |
| `ClusteringConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `CompareConditionsConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `CompareGroupsConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `DifferentialConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `EnrichmentConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `FilterMarkersConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `MergeClustersConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `MethodSelectionConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `proportion`; optional |
| `ProportionConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] [O] from `config` |
| `PseudobulkDEConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `ResolutionSearchConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |
| `ScoringConfig` | Config Class | `src/scLucid/analysis/__init__.py` | [C] from `config` |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `AnalysisWorkflowError` | Function | `src/scLucid/analysis/__init__.py` | from `workflow` |
| `analyze_all_methods` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `analyze_celltype_proportion` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `annotate_clusters` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `apply_annotation_mapping` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `apply_final_annotation` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `apply_subset_annotation_reconciliation` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `batch_compare_scores` | Function | `src/scLucid/analysis/__init__.py` | from `scoring` |
| `batch_plot_delta_heatmap` | Function | `src/scLucid/analysis/__init__.py` | from `scoring` |
| `build_annotation_consensus` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `build_hierarchical_annotation_plan` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `build_llm_annotation_bundle` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `build_subset_annotation_reconciliation` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `calculate_signature_matrix` | Function | `src/scLucid/analysis/__init__.py` | from `scoring` |
| `celltype_proportion_analysis` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `cluster_cells` | Function | `src/scLucid/analysis/__init__.py` | from `clustering` |
| `compare_clustering_resolutions` | Function | `src/scLucid/analysis/__init__.py` | from `workflow` |
| `compare_methods` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `compute_celltype_proportion` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `correlate_abundance_with_clinical` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `deconvolve_bulk` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `deduplicate_var_names` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `diagnose_bulk_data_quality` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `estimate_size_factors_median_ratio` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `evaluate_annotation` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `export_analysis_data` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `filter_bulk_genes` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `filter_marker_table_for_annotation` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `flag_suspect_clusters` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `merge_annotation_evidence` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `merge_clusters` | Function | `src/scLucid/analysis/__init__.py` | from `clustering` |
| `normalize_bulk_counts` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `PartialAnalysisResult` | Function | `src/scLucid/analysis/__init__.py` | from `workflow` |
| `pb_analysis` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_batch_effect` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_box_summary` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_cell_counts` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_celltype_alluvial` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_celltype_correlation` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_celltype_variability` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_composition` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_composition_pca` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_composition_transform_heatmap` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_delta_heatmap` | Function | `src/scLucid/analysis/__init__.py` | from `scoring` |
| `plot_diff_stats` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_effect_size_volcano` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_grouped_celltype_counts` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_grouped_proportion_bar` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_individual_boxplots` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_paired_proportion_shifts` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_proportion_bar` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_proportion_heatmap` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_proportion_shifts` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_proportion_timeseries` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_proportion_with_ci` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `plot_score_violin_with_stats` | Function | `src/scLucid/analysis/__init__.py` | from `scoring` |
| `plot_signature_heatmap` | Function | `src/scLucid/analysis/__init__.py` | from `scoring` |
| `ProportionMethod` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `recommend_method` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `remap_labels` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `run_annotation_evidence` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `run_bulk_abundance_test` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `run_bulk_de` | Function | `src/scLucid/analysis/__init__.py` | [O] from `bulk`; optional |
| `run_celltypist` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `run_clustering_review` | Function | `src/scLucid/analysis/__init__.py` | from `clustering` |
| `run_lineage_state_annotation` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `run_marker_annotation_evidence` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `run_module_scoring_workflow` | Function | `src/scLucid/analysis/__init__.py` | from `scoring` |
| `run_program_annotation_evidence` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `run_sccoda` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `run_statistical_test` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |
| `run_subset_annotation_refinement` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `score_by_gene_sets` | Function | `src/scLucid/analysis/__init__.py` | from `scoring` |
| `score_cell_types` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `standardize_cluster_marker_table` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `summarize_annotation_evidence` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `transfer_labels` | Function | `src/scLucid/analysis/__init__.py` | from `annotation` |
| `transform_composition` | Function | `src/scLucid/analysis/__init__.py` | [O] from `proportion`; optional |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `ANALYSIS_REQUIRED_REVIEW_SECTIONS` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] from `trace` |
| `build_annotation_review_table` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] from `annotation` |
| `build_bulk_review_summary` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] [O] from `bulk`; optional |
| `build_posthoc_qc_review_summary` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] from `trace` |
| `enrich_analysis_review_summary` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] from `trace` |
| `get_analysis_module_contract` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] from `trace` |
| `summarize_analysis_review_summary` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] from `trace` |
| `validate_analysis_module_completeness` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] from `trace` |
| `validate_analysis_review_summary` | Trace / Contract | `src/scLucid/analysis/__init__.py` | [T] from `trace` |

### Deprecated / Uncertain / Private-but-Exposed

#### Deprecated

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_malignancy_interpretation` | Deprecated | `src/scLucid/analysis/__init__.py` | [D] |

**Summary:** 112 symbols (111 stable, 1 flagged). workflow=3, config=21, class=0, function=78, alias=0, constant=0, trace=9, deprecated=1, uncertain=0, private_but_exposed=0.

## scLucid.analysis.annotation

### Stable APIs

#### Workflow Orchestrator

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_annotation` | Workflow Orchestrator | `src/scLucid/analysis/annotation/__init__.py` | [W] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `annotate_clusters` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `apply_annotation_mapping` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `apply_final_annotation` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `apply_subset_annotation_reconciliation` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `build_annotation_consensus` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `build_hierarchical_annotation_plan` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `build_llm_annotation_bundle` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `build_subset_annotation_reconciliation` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `evaluate_annotation` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `evaluate_annotation_benchmark` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `filter_marker_table_for_annotation` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `flag_suspect_clusters` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `merge_annotation_evidence` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `remap_labels` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `run_annotation_evidence` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `run_celltypist` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `run_lineage_state_annotation` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `run_marker_annotation_evidence` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `run_program_annotation_evidence` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `run_subset_annotation_refinement` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `score_cell_types` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `standardize_cluster_marker_table` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `summarize_annotation_evidence` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |
| `transfer_labels` | Function | `src/scLucid/analysis/annotation/__init__.py` |  |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `ANALYSIS_REVIEW_SUMMARY_SCHEMA` | Trace / Contract | `src/scLucid/analysis/annotation/__init__.py` | [T] |
| `ANNOTATION_REVIEW_SCHEMA` | Trace / Contract | `src/scLucid/analysis/annotation/__init__.py` | [T] |
| `build_annotation_review_table` | Trace / Contract | `src/scLucid/analysis/annotation/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 28 symbols (28 stable, 0 flagged). workflow=1, config=0, class=0, function=24, alias=0, constant=0, trace=3, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.analysis.bulk

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `BulkAbundanceConfig` | Config Class | `src/scLucid/analysis/bulk/__init__.py` | [C] |
| `BulkClinicalAssociationConfig` | Config Class | `src/scLucid/analysis/bulk/__init__.py` | [C] |
| `BulkDEConfig` | Config Class | `src/scLucid/analysis/bulk/__init__.py` | [C] |
| `BulkDeconvolutionConfig` | Config Class | `src/scLucid/analysis/bulk/__init__.py` | [C] |
| `BulkDiagnosticsConfig` | Config Class | `src/scLucid/analysis/bulk/__init__.py` | [C] |
| `BulkNormalizationConfig` | Config Class | `src/scLucid/analysis/bulk/__init__.py` | [C] |
| `BulkTraitAssociationConfig` | Config Class | `src/scLucid/analysis/bulk/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `correlate_abundance_with_clinical` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `deconvolve_bulk` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `deduplicate_var_names` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `diagnose_bulk_data_quality` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `differential_abundance` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `estimate_size_factors_median_ratio` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `filter_bulk_genes` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `normalize_bulk_counts` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `run_bulk_abundance_test` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `run_bulk_de` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |
| `run_deconvolution` | Function | `src/scLucid/analysis/bulk/__init__.py` |  |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `build_bulk_review_summary` | Trace / Contract | `src/scLucid/analysis/bulk/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 19 symbols (19 stable, 0 flagged). workflow=0, config=7, class=0, function=11, alias=0, constant=0, trace=1, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.analysis.differential_expression

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `batch_celltype_deg_enrichment` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `characterize_clusters` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `compare_conditions` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `compare_groups` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `export_enrichment_results` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `filter_markers` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `find_markers` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `get_conserved_markers` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `load_results` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `plot_multi_cluster_deg` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `plot_volcano` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `ResultManager` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `run_enrichment` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `run_pseudobulk_de` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `save_results` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `summarize_markers_and_enrichment` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |
| `visualize_markers` | Function | `src/scLucid/analysis/differential_expression/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 17 symbols (17 stable, 0 flagged). workflow=0, config=0, class=0, function=17, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.analysis.proportion

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `MethodSelectionConfig` | Config Class | `src/scLucid/analysis/proportion/__init__.py` | [C] |
| `ProportionConfig` | Config Class | `src/scLucid/analysis/proportion/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `analyze_all_methods` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `analyze_celltype_proportion` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `celltype_proportion_analysis` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `compare_methods` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `composition_transform` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `compute_celltype_proportion` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `export_analysis_data` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_batch_effect` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_box_summary` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_cell_counts` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_celltype_alluvial` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_celltype_correlation` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_celltype_variability` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_composition` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_composition_pca` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_composition_transform_heatmap` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_diff_stats` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_effect_size_volcano` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_grouped_celltype_counts` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_grouped_proportion_bar` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_individual_boxplots` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_paired_proportion_shifts` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_proportion_bar` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_proportion_heatmap` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_proportion_shifts` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_proportion_timeseries` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `plot_proportion_with_ci` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `ProportionMethod` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `recommend_method` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `run_sccoda` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `run_statistical_test` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `summarize_sccoda` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |
| `transform_composition` | Function | `src/scLucid/analysis/proportion/__init__.py` |  |

#### Alias

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `pb_analysis` | Alias | `src/scLucid/analysis/proportion/__init__.py` | [A] from `celltype_proportion_analysis` |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 36 symbols (36 stable, 0 flagged). workflow=0, config=2, class=0, function=33, alias=1, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.plotting

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `apply_theme` | Function | `src/scLucid/plotting/__init__.py` |  |
| `build_color_palette` | Function | `src/scLucid/plotting/__init__.py` |  |
| `build_obs_palette` | Function | `src/scLucid/plotting/__init__.py` |  |
| `export_annotation_report` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_annotation_evidence_panel` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_coexpression` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_differential_abundance` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_dotplot` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_embedding` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_faceted_embedding` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_faceted_feature` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_feature_correlation` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_marker_expression` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_marker_heatmap` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_ridge` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_split_violin_with_stats` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_stacked_violin` | Function | `src/scLucid/plotting/__init__.py` |  |
| `plot_volcano` | Function | `src/scLucid/plotting/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 18 symbols (18 stable, 0 flagged). workflow=0, config=0, class=0, function=18, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.preprocess

### Stable APIs

#### Workflow Orchestrator

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_preprocessing` | Workflow Orchestrator | `src/scLucid/preprocess/__init__.py` | [W] |

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `AdaptiveNormalizationConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `GraphConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `HVGConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `IntegrationConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `IntelligentPreprocessConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `NeighborsConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `NormalizationConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `PreprocessingWorkflowConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `ScalingConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |
| `WorkflowConfig` | Config Class | `src/scLucid/preprocess/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `adaptive_normalize` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `annotate_gene_biotypes` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `apply_gene_biotype_strategy` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `batch_correction` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `BatchCorrectionRecommendation` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `build_qc_input_context` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `build_step_evidence_summary` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `DataProfile` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `decide_integration` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `detect_integration_confounding` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `diagnose_cell_cycle_regression` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `diagnose_integration_risk` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `estimate_cell_size_factors` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `evaluate_hvg_stability` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `evaluate_integration` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `filter_genes_by_biotype` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `find_hvgs` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `get_backend` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `get_biotype_statistics` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `get_gene_biotype_cache_dir` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `HVGRecommendation` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `IntelligentPreprocessRecommender` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `list_available_backends` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `list_gene_biotype_resources` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `load_gene_biotypes` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `NeighborsRecommendation` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `normalize_data` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `optimize_neighbors_pcs` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `PartialWorkflowResult` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `PCARecommendation` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `plot_hvg_metrics` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `plot_normalization_effect` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `plot_scaling_effect` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `PreprocessingBackend` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `PreprocessingStrategy` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `quality_aware_normalize` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `RapidsBackend` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `recommend_biotype_strategy` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `recommend_intelligent_preprocessing` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `regress_out` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `ResolutionRecommendation` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `run_intelligent_preprocessing` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `scale_data` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `ScanpyBackend` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `select_and_audit_hvgs` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `select_hvg_sets` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `set_backend` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `suggest_hvg_choice` | Function | `src/scLucid/preprocess/__init__.py` |  |
| `WorkflowError` | Function | `src/scLucid/preprocess/__init__.py` |  |

#### Constant

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `WORKFLOW_STEPS` | Constant | `src/scLucid/preprocess/__init__.py` |  |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `build_preprocess_module_maturity_assessment` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `enrich_preprocessing_review_summary` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `get_preprocess_module_contract` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `PREPROCESS_MODULE_MATURITY_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `PREPROCESS_REQUIRED_REVIEW_SECTIONS` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `PREPROCESS_STABLE_ENTRYPOINTS` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `PREPROCESS_TRACE_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `summarize_preprocess_review_summary` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `validate_preprocess_module_completeness` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |
| `validate_preprocessing_review_summary` | Trace / Contract | `src/scLucid/preprocess/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

#### Deprecated

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_embedding_workflow` | Deprecated | `src/scLucid/preprocess/__init__.py` | [D] |

**Summary:** 72 symbols (71 stable, 1 flagged). workflow=1, config=10, class=0, function=49, alias=0, constant=1, trace=10, deprecated=1, uncertain=0, private_but_exposed=0.

## scLucid.preprocess.intelligent

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `IntelligentPreprocessConfig` | Config Class | `src/scLucid/preprocess/intelligent/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `BatchCorrectionRecommendation` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `DataProfile` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `HVGRecommendation` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `IntelligentPreprocessRecommender` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `NeighborsRecommendation` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `PCARecommendation` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `PreprocessingStrategy` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `recommend_intelligent_preprocessing` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `ResolutionRecommendation` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |
| `run_intelligent_preprocessing` | Function | `src/scLucid/preprocess/intelligent/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 11 symbols (11 stable, 0 flagged). workflow=0, config=1, class=0, function=10, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.qc

### Stable APIs

#### Workflow Orchestrator

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_standard_qc` | Workflow Orchestrator | `src/scLucid/qc/__init__.py` | [W] from `workflow` |

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `DoubletConfig` | Config Class | `src/scLucid/qc/__init__.py` | [C] from `config` |
| `FilterConfig` | Config Class | `src/scLucid/qc/__init__.py` | [C] from `config` |
| `MarkerConfig` | Config Class | `src/scLucid/qc/__init__.py` | [C] from `config` |
| `MarkingConfig` | Config Class | `src/scLucid/qc/__init__.py` | [C] from `config` |
| `MetricsReportingConfig` | Config Class | `src/scLucid/qc/__init__.py` | [C] from `config` |
| `QCWorkflowConfig` | Config Class | `src/scLucid/qc/__init__.py` | [C] from `config` |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `AdaptiveThresholdLearner` | Function | `src/scLucid/qc/__init__.py` | [O] from `adaptive_threshold`; optional |
| `audit_doublets` | Function | `src/scLucid/qc/__init__.py` | from `doublet` |
| `audit_filtering` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `build_qc_benchmark_assessment` | Function | `src/scLucid/qc/__init__.py` | from `benchmark` |
| `build_qc_decision_table` | Function | `src/scLucid/qc/__init__.py` | from `trace` |
| `calculate_qc_metric` | Function | `src/scLucid/qc/__init__.py` | from `metrics` |
| `compute_marker_fidelity` | Function | `src/scLucid/qc/__init__.py` | from `benchmark` |
| `compute_retention_metrics` | Function | `src/scLucid/qc/__init__.py` | from `benchmark` |
| `create_custom_marker_dict` | Function | `src/scLucid/qc/__init__.py` | from `doublet` |
| `diagnose_ambient_rna` | Function | `src/scLucid/qc/__init__.py` | from `ambient` |
| `diagnose_empty_droplets` | Function | `src/scLucid/qc/__init__.py` | from `ambient` |
| `EnhancedQCReport` | Function | `src/scLucid/qc/__init__.py` | [O] from `reporting`; optional |
| `evaluate_qc_benchmark` | Function | `src/scLucid/qc/__init__.py` | from `benchmark` |
| `export_qc_benchmark_report` | Function | `src/scLucid/qc/__init__.py` | from `benchmark` |
| `filter_cells` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `generate_doublet_rates` | Function | `src/scLucid/qc/__init__.py` | from `doublet` |
| `generate_qc_html_report` | Function | `src/scLucid/qc/__init__.py` | [O] from `reporting`; optional |
| `generate_qc_report` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `identify_outliers` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `infer_qc_benchmark_profile` | Function | `src/scLucid/qc/__init__.py` | from `benchmark` |
| `IntelligentQCRecommender` | Function | `src/scLucid/qc/__init__.py` | [O] from `intelligent_qc`; optional |
| `InteractiveReportGenerator` | Function | `src/scLucid/qc/__init__.py` | [O] from `reporting`; optional |
| `mark_low_quality_cell` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `mark_low_quality_cells_adaptive` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `MultiMetricAdaptiveLearner` | Function | `src/scLucid/qc/__init__.py` | [O] from `adaptive_threshold`; optional |
| `predict_doublets` | Function | `src/scLucid/qc/__init__.py` | from `doublet` |
| `predict_doublets_with_profiling` | Function | `src/scLucid/qc/__init__.py` | from `doublet` |
| `QCRecommendation` | Function | `src/scLucid/qc/__init__.py` | [O] from `intelligent_qc`; optional |
| `QCThresholds` | Function | `src/scLucid/qc/__init__.py` | from `config` |
| `recommend_intelligent_qc` | Function | `src/scLucid/qc/__init__.py` | [O] from `intelligent_qc`; optional |
| `record_ambient_correction_status` | Function | `src/scLucid/qc/__init__.py` | from `ambient` |
| `register_external_ambient_result` | Function | `src/scLucid/qc/__init__.py` | from `ambient` |
| `render_qc_benchmark_compact_markdown` | Function | `src/scLucid/qc/__init__.py` | from `benchmark` |
| `render_qc_benchmark_markdown` | Function | `src/scLucid/qc/__init__.py` | from `benchmark` |
| `resolve_qc_thresholds` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `run_qc_threshold_decision` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `score_cell_cycle` | Function | `src/scLucid/qc/__init__.py` | from `cycle` |
| `StrategyType` | Function | `src/scLucid/qc/__init__.py` | [O] from `intelligent_qc`; optional |
| `suggest_qc_thresholds` | Function | `src/scLucid/qc/__init__.py` | from `filtering` |
| `ThresholdRecommendation` | Function | `src/scLucid/qc/__init__.py` | [O] from `intelligent_qc`; optional |

#### Constant

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `BENCHMARK_PROFILES` | Constant | `src/scLucid/qc/__init__.py` | from `benchmark` |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `build_qc_module_maturity_assessment` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `enrich_qc_review_summary` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `get_qc_module_contract` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `QC_BENCHMARK_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `benchmark` |
| `QC_MODULE_MATURITY_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `QC_REQUIRED_OBS_METRICS` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `QC_REQUIRED_REVIEW_SECTIONS` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `QC_STABLE_ENTRYPOINTS` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `QC_TRACE_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `summarize_qc_review_summary` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `validate_qc_module_completeness` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |
| `validate_qc_review_summary` | Trace / Contract | `src/scLucid/qc/__init__.py` | [T] from `trace` |

### Deprecated / Uncertain / Private-but-Exposed

#### Deprecated

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_advanced_qc` | Deprecated | `src/scLucid/qc/__init__.py` | [D] from `workflow` |
| `run_qc_decision_workflow` | Deprecated | `src/scLucid/qc/__init__.py` | [D] from `filtering` |

**Summary:** 62 symbols (60 stable, 2 flagged). workflow=1, config=6, class=0, function=40, alias=0, constant=1, trace=12, deprecated=2, uncertain=0, private_but_exposed=0.

## scLucid.recommendation

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `RecommendationConfig` | Config Class | `src/scLucid/recommendation/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `ParameterRecommendation` | Function | `src/scLucid/recommendation/__init__.py` |  |
| `recommend_analysis_parameters` | Function | `src/scLucid/recommendation/__init__.py` |  |
| `RecommendationEngine` | Function | `src/scLucid/recommendation/__init__.py` |  |
| `RecommendationSection` | Function | `src/scLucid/recommendation/__init__.py` |  |
| `WorkflowRecommendations` | Function | `src/scLucid/recommendation/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 6 symbols (6 stable, 0 flagged). workflow=0, config=1, class=0, function=5, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.recommendation.engine

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `recommend_analysis_parameters` | Function | `src/scLucid/recommendation/engine/__init__.py` |  |
| `RecommendationEngine` | Function | `src/scLucid/recommendation/engine/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 2 symbols (2 stable, 0 flagged). workflow=0, config=0, class=0, function=2, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tools

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `BulkAbundanceConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `bulk` |
| `BulkClinicalAssociationConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `bulk` |
| `BulkDEConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `bulk` |
| `BulkDeconvolutionConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `bulk` |
| `BulkDiagnosticsConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `bulk` |
| `BulkNormalizationConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `bulk` |
| `BulkTraitAssociationConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `bulk` |
| `PrismConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `pyBayesPrism` |
| `SpatialAutocorrConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `spatial` |
| `SpatialDiagnosticsConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `spatial` |
| `SpatialNeighborsConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `spatial` |
| `SpatialWindowConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `spatial` |
| `SVGConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `spatial` |
| `TissueZonesConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `spatial` |
| `VisiumIOConfig` | Config Class | `src/scLucid/tools/__init__.py` | [C] from `spatial` |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `analyze_scenic_results` | Function | `src/scLucid/tools/__init__.py` | from `pySCENIC` |
| `BayesPrism` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `BayesPrismEmbedding` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `BayesPrismReference` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `build_spatial_neighbors` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `CellChat` | Function | `src/scLucid/tools/__init__.py` | from `pyCellChat` |
| `CellChatDB` | Function | `src/scLucid/tools/__init__.py` | from `pyCellChat` |
| `CellDataSet` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `cleanup_genes` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `cluster_cells` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `compute_correlation` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `compute_moran_i` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `compute_rmse` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `compute_spatial_autocorr` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `correlate_abundance_with_clinical` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `create_cds_from_scanpy` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `create_cellchat_from_scanpy` | Function | `src/scLucid/tools/__init__.py` | from `pyCellChat` |
| `create_pseudo_bulk` | Function | `src/scLucid/tools/__init__.py` | from `pyDWLS` |
| `crop_visium` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `CrossValidator` | Function | `src/scLucid/tools/__init__.py` | from `pyDWLS` |
| `DampenedWLS` | Function | `src/scLucid/tools/__init__.py` | from `pyDWLS` |
| `deconvolve_bulk` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `deduplicate_var_names` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `diagnose_bulk_data_quality` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `diagnose_spatial_data_quality` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `differential_abundance` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `estimate_size_factors_median_ratio` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `export_scenic_report` | Function | `src/scLucid/tools/__init__.py` | from `pySCENIC` |
| `export_spatial_report` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `export_to_scanpy` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `filter_bulk_genes` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `filter_genes` | Function | `src/scLucid/tools/__init__.py` | from `pyDWLS` |
| `find_spatially_variable_genes` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `find_tissue_zones` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `get_default_database` | Function | `src/scLucid/tools/__init__.py` | from `pyCellChat` |
| `GibbsSampler` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `graph_test` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `learn_graph` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `MarkerSelector` | Function | `src/scLucid/tools/__init__.py` | from `pyDWLS` |
| `new_cell_data_set` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `normalize_bulk_counts` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `normalize_data` | Function | `src/scLucid/tools/__init__.py` | from `pyDWLS` |
| `order_cells` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `plot_cells` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `plot_correlation` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `plot_fraction` | Function | `src/scLucid/tools/__init__.py` | from `pyBayesPrism` |
| `plot_heatmap` | Function | `src/scLucid/tools/__init__.py` | from `pyCellChat` |
| `plot_spatial` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `preprocess_cds` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `read_visium_10x` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `reduce_dimension` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |
| `rotate_visium` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `run_bulk_abundance_test` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `run_bulk_de` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `run_cellphonedb` | Function | `src/scLucid/tools/__init__.py` | from `cellphonedb` |
| `run_cellphonedb_batch` | Function | `src/scLucid/tools/__init__.py` | from `cellphonedb` |
| `run_cellphonedb_by_group` | Function | `src/scLucid/tools/__init__.py` | from `cellphonedb` |
| `run_deconvolution` | Function | `src/scLucid/tools/__init__.py` | from `bulk` |
| `run_scenic` | Function | `src/scLucid/tools/__init__.py` | from `pySCENIC` |
| `run_scenic_batch` | Function | `src/scLucid/tools/__init__.py` | from `pySCENIC` |
| `run_scenic_by_group` | Function | `src/scLucid/tools/__init__.py` | from `pySCENIC` |
| `run_spatial_analysis` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `run_spatial_batch` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `SignatureBuilder` | Function | `src/scLucid/tools/__init__.py` | from `pyDWLS` |
| `solve_nnls` | Function | `src/scLucid/tools/__init__.py` | from `pyDWLS` |
| `subset_spatial_window` | Function | `src/scLucid/tools/__init__.py` | from `spatial` |
| `summarize_cellphonedb` | Function | `src/scLucid/tools/__init__.py` | from `cellphonedb` |
| `top_markers` | Function | `src/scLucid/tools/__init__.py` | from `pyMonocle3` |

#### Constant

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `DWLS` | Constant | `src/scLucid/tools/__init__.py` | from `pyDWLS` |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `build_bulk_review_summary` | Trace / Contract | `src/scLucid/tools/__init__.py` | [T] from `bulk` |
| `build_spatial_review_summary` | Trace / Contract | `src/scLucid/tools/__init__.py` | [T] from `spatial` |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 86 symbols (86 stable, 0 flagged). workflow=0, config=15, class=0, function=68, alias=0, constant=1, trace=2, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tools.bulk

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `BulkAbundanceConfig` | Config Class | `src/scLucid/tools/bulk/__init__.py` | [C] |
| `BulkClinicalAssociationConfig` | Config Class | `src/scLucid/tools/bulk/__init__.py` | [C] |
| `BulkDEConfig` | Config Class | `src/scLucid/tools/bulk/__init__.py` | [C] |
| `BulkDeconvolutionConfig` | Config Class | `src/scLucid/tools/bulk/__init__.py` | [C] |
| `BulkDiagnosticsConfig` | Config Class | `src/scLucid/tools/bulk/__init__.py` | [C] |
| `BulkNormalizationConfig` | Config Class | `src/scLucid/tools/bulk/__init__.py` | [C] |
| `BulkTraitAssociationConfig` | Config Class | `src/scLucid/tools/bulk/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `associate_tme_with_response` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `bulk_immune_landscape` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `correlate_abundance_with_clinical` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `deconvolve_bulk` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `deconvolve_tumor_tme` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `deduplicate_var_names` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `diagnose_bulk_data_quality` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `differential_abundance` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `estimate_size_factors_median_ratio` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `estimate_tumor_purity_from_bulk` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `filter_bulk_genes` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `normalize_bulk_counts` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `run_bulk_abundance_test` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `run_bulk_de` | Function | `src/scLucid/tools/bulk/__init__.py` |  |
| `run_deconvolution` | Function | `src/scLucid/tools/bulk/__init__.py` |  |

#### Alias

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `correlate_abundance_with_clinical_alias` | Alias | `src/scLucid/tools/bulk/__init__.py` | [A] from `clinical.correlate_abundance_with_clinical` |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `build_bulk_review_summary` | Trace / Contract | `src/scLucid/tools/bulk/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 24 symbols (24 stable, 0 flagged). workflow=0, config=7, class=0, function=15, alias=1, constant=0, trace=1, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tools.pyBayesPrism

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `DeconvolutionConfig` | Config Class | `src/scLucid/tools/pyBayesPrism/__init__.py` | [C] |
| `PrismConfig` | Config Class | `src/scLucid/tools/pyBayesPrism/__init__.py` | [C] |
| `ReferenceConfig` | Config Class | `src/scLucid/tools/pyBayesPrism/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `batch_correct` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `BayesPrism` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `BayesPrismEmbedding` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `BayesPrismReference` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `cleanup_genes` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `compute_correlation` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `compute_rmse` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `find_outlier_genes` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `GibbsSampler` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `normalize_expression` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `plot_correlation` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `plot_cv` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `plot_fraction` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `plot_gene_programs` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `plot_program_usage` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `plot_stacked_bar` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `plot_validation_scatter` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |
| `subsample_cells` | Function | `src/scLucid/tools/pyBayesPrism/__init__.py` |  |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `validate_inputs` | Trace / Contract | `src/scLucid/tools/pyBayesPrism/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 22 symbols (22 stable, 0 flagged). workflow=0, config=3, class=0, function=18, alias=0, constant=0, trace=1, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tools.pyCellChat

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `CellChatConfig` | Config Class | `src/scLucid/tools/pyCellChat/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `CellChat` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `CellChatDB` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `compare_cellchat_objects` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `compute_centrality` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `compute_network_similarity` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `create_cellchat_from_scanpy` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `get_default_database` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `identify_conserved_pathways` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `identify_differential_pathways` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `identify_roles` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `merge_cellchat_objects` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `plot_bubble` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `plot_chord_diagram` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `plot_circle_network` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `plot_contribution` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `plot_heatmap` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |
| `plot_signaling_gene_expression` | Function | `src/scLucid/tools/pyCellChat/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 18 symbols (18 stable, 0 flagged). workflow=0, config=1, class=0, function=17, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tools.pyDWLS

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `align_data` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |
| `create_pseudo_bulk` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |
| `CrossValidator` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |
| `DampenedWLS` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |
| `filter_genes` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |
| `MarkerSelector` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |
| `normalize_data` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |
| `SignatureBuilder` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |
| `solve_nnls` | Function | `src/scLucid/tools/pyDWLS/__init__.py` |  |

#### Constant

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `DWLS` | Constant | `src/scLucid/tools/pyDWLS/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 10 symbols (10 stable, 0 flagged). workflow=0, config=0, class=0, function=9, alias=0, constant=1, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tools.pyMonocle3

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `aggregate_gene_expression` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `align_cds` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `calculate_gene_modules` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `CellDataSet` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `choose_graph_segments` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `cluster_cells` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `compare_genes` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `convert_to_dense` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `convert_to_sparse` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `create_cds_from_scanpy` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `detect_genes` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `detect_sparse_type` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `estimate_memory_usage` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `estimate_size_factors` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `export_to_scanpy` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `find_cluster_markers` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `graph_test` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `group_cells` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `learn_graph` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `merge_datasets` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `new_cell_data_set` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `normalize_expression` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `order_cells` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `partition_cells` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `plot_cells` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `plot_genes_by_group` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `plot_pseudotime_heatmap` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `plot_trajectory` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `preprocess_cds` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `pseudotime_de` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `reduce_dimension` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `run_pca` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `run_umap` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `select_highly_variable_genes` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `subsample_cells` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |
| `top_markers` | Function | `src/scLucid/tools/pyMonocle3/__init__.py` |  |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `validate_cds` | Trace / Contract | `src/scLucid/tools/pyMonocle3/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 37 symbols (37 stable, 0 flagged). workflow=0, config=0, class=0, function=36, alias=0, constant=0, trace=1, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tools.spatial

### Stable APIs

#### Config Class

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `SpatialAutocorrConfig` | Config Class | `src/scLucid/tools/spatial/__init__.py` | [C] |
| `SpatialDiagnosticsConfig` | Config Class | `src/scLucid/tools/spatial/__init__.py` | [C] |
| `SpatialNeighborsConfig` | Config Class | `src/scLucid/tools/spatial/__init__.py` | [C] |
| `SpatialWindowConfig` | Config Class | `src/scLucid/tools/spatial/__init__.py` | [C] |
| `SVGConfig` | Config Class | `src/scLucid/tools/spatial/__init__.py` | [C] |
| `TissueZonesConfig` | Config Class | `src/scLucid/tools/spatial/__init__.py` | [C] |
| `VisiumIOConfig` | Config Class | `src/scLucid/tools/spatial/__init__.py` | [C] |

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `analyze_spatial_niches` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `build_spatial_neighbors` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `compute_immune_infiltration_score` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `compute_moran_i` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `compute_spatial_autocorr` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `crop_visium` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `diagnose_spatial_data_quality` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `export_spatial_report` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `find_spatially_variable_genes` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `find_tissue_zones` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `find_tumor_stroma_boundary` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `infer_spatial_platform` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `plot_spatial` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `read_visium_10x` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `rotate_visium` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `run_spatial_analysis` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `run_spatial_batch` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `spatial_ici_response_signature` | Function | `src/scLucid/tools/spatial/__init__.py` |  |
| `subset_spatial_window` | Function | `src/scLucid/tools/spatial/__init__.py` |  |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `build_spatial_review_summary` | Trace / Contract | `src/scLucid/tools/spatial/__init__.py` | [T] |
| `validate_spatial_coords` | Trace / Contract | `src/scLucid/tools/spatial/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 28 symbols (28 stable, 0 flagged). workflow=0, config=7, class=0, function=19, alias=0, constant=0, trace=2, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tumor

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `align_progression_trajectories` | Function | `src/scLucid/tumor/__init__.py` |  |
| `analyze_cell_interactions` | Function | `src/scLucid/tumor/__init__.py` |  |
| `analyze_dissemination` | Function | `src/scLucid/tumor/__init__.py` |  |
| `analyze_ecosystem_composition` | Function | `src/scLucid/tumor/__init__.py` |  |
| `analyze_immune_infiltration` | Function | `src/scLucid/tumor/__init__.py` |  |
| `analyze_regional_heterogeneity` | Function | `src/scLucid/tumor/__init__.py` |  |
| `analyze_treatment_response_trajectory` | Function | `src/scLucid/tumor/__init__.py` |  |
| `analyze_tumor_progression` | Function | `src/scLucid/tumor/__init__.py` |  |
| `assign_cnv_signature` | Function | `src/scLucid/tumor/__init__.py` |  |
| `build_phylogenetic_tree` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_clonal_diversity` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_cnv_score` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_diversity_indices` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_proliferation_index` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_regional_expression_differences` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_signature_scores` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_stemness_score` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_transcriptional_diversity` | Function | `src/scLucid/tumor/__init__.py` |  |
| `calculate_tumor_microenvironment_score` | Function | `src/scLucid/tumor/__init__.py` |  |
| `classify_malignant_cells` | Function | `src/scLucid/tumor/__init__.py` |  |
| `CloneAnalyzer` | Function | `src/scLucid/tumor/__init__.py` |  |
| `CNVInferenceStep` | Function | `src/scLucid/tumor/__init__.py` |  |
| `CNVSigExtractor` | Function | `src/scLucid/tumor/__init__.py` |  |
| `compare_ecosystems` | Function | `src/scLucid/tumor/__init__.py` |  |
| `compare_primary_vs_metastasis` | Function | `src/scLucid/tumor/__init__.py` |  |
| `compare_resistance_between_groups` | Function | `src/scLucid/tumor/__init__.py` |  |
| `compare_stemness_between_groups` | Function | `src/scLucid/tumor/__init__.py` |  |
| `deconvolve_tme` | Function | `src/scLucid/tumor/__init__.py` |  |
| `detect_clonal_sweep` | Function | `src/scLucid/tumor/__init__.py` |  |
| `discover_therapeutic_targets` | Function | `src/scLucid/tumor/__init__.py` |  |
| `estimate_intratumoral_heterogeneity` | Function | `src/scLucid/tumor/__init__.py` |  |
| `estimate_metastatic_potential` | Function | `src/scLucid/tumor/__init__.py` |  |
| `estimate_stromal_content` | Function | `src/scLucid/tumor/__init__.py` |  |
| `evaluate_biomarker` | Function | `src/scLucid/tumor/__init__.py` |  |
| `extract_cnv_signatures` | Function | `src/scLucid/tumor/__init__.py` |  |
| `find_dominant_interactions` | Function | `src/scLucid/tumor/__init__.py` |  |
| `find_tumor` | Function | `src/scLucid/tumor/__init__.py` |  |
| `find_tumor_cells` | Function | `src/scLucid/tumor/__init__.py` |  |
| `get_drug_targets` | Function | `src/scLucid/tumor/__init__.py` |  |
| `get_immune_markers` | Function | `src/scLucid/tumor/__init__.py` |  |
| `get_stromal_markers` | Function | `src/scLucid/tumor/__init__.py` |  |
| `get_tumor_markers` | Function | `src/scLucid/tumor/__init__.py` |  |
| `identify_cancer_stem_cells` | Function | `src/scLucid/tumor/__init__.py` |  |
| `identify_clones` | Function | `src/scLucid/tumor/__init__.py` |  |
| `identify_resistance_mechanisms` | Function | `src/scLucid/tumor/__init__.py` |  |
| `identify_spatial_patterns` | Function | `src/scLucid/tumor/__init__.py` |  |
| `identify_transition_states` | Function | `src/scLucid/tumor/__init__.py` |  |
| `infer_clonal_phylogeny` | Function | `src/scLucid/tumor/__init__.py` |  |
| `infer_cnv` | Function | `src/scLucid/tumor/__init__.py` |  |
| `load_hallmark_signatures` | Function | `src/scLucid/tumor/__init__.py` |  |
| `MalignancyInterpretationStep` | Function | `src/scLucid/tumor/__init__.py` |  |
| `MalignancyScoringStep` | Function | `src/scLucid/tumor/__init__.py` |  |
| `predict_metastasis_risk` | Function | `src/scLucid/tumor/__init__.py` |  |
| `predict_therapy_response` | Function | `src/scLucid/tumor/__init__.py` |  |
| `prioritize_druggable_genes` | Function | `src/scLucid/tumor/__init__.py` |  |
| `query_cancer_gene_census` | Function | `src/scLucid/tumor/__init__.py` |  |
| `root_tree` | Function | `src/scLucid/tumor/__init__.py` |  |
| `run_cnv_analysis` | Function | `src/scLucid/tumor/__init__.py` |  |
| `score_drug_resistance` | Function | `src/scLucid/tumor/__init__.py` |  |
| `score_immune_interactions` | Function | `src/scLucid/tumor/__init__.py` |  |
| `score_malignancy` | Function | `src/scLucid/tumor/__init__.py` |  |
| `score_malignancy_potential` | Function | `src/scLucid/tumor/__init__.py` |  |
| `stratify_patients` | Function | `src/scLucid/tumor/__init__.py` |  |
| `suggest_targeted_therapies` | Function | `src/scLucid/tumor/__init__.py` |  |
| `TherapyPredictionStep` | Function | `src/scLucid/tumor/__init__.py` |  |
| `TMEDeconvolutionStep` | Function | `src/scLucid/tumor/__init__.py` |  |
| `track_temporal_dynamics` | Function | `src/scLucid/tumor/__init__.py` |  |

#### Alias

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `identify_clones_from_cnv` | Alias | `src/scLucid/tumor/__init__.py` | [A] from `cnv.clone_analysis.identify_clones` |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `enrich_tumor_review_summary` | Trace / Contract | `src/scLucid/tumor/__init__.py` | [T] |
| `get_tumor_module_contract` | Trace / Contract | `src/scLucid/tumor/__init__.py` | [T] |
| `TUMOR_MODULE_MATURITY_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/tumor/__init__.py` | [T] |
| `TUMOR_REVIEW_SUMMARY_REQUIRED_KEYS` | Trace / Contract | `src/scLucid/tumor/__init__.py` | [T] |
| `TUMOR_TRACE_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/tumor/__init__.py` | [T] |
| `validate_tumor_module_completeness` | Trace / Contract | `src/scLucid/tumor/__init__.py` | [T] |
| `validate_tumor_review_summary` | Trace / Contract | `src/scLucid/tumor/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

#### Deprecated

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_malignancy_interpretation` | Deprecated | `src/scLucid/tumor/__init__.py` | [D] |

**Summary:** 76 symbols (75 stable, 1 flagged). workflow=0, config=0, class=0, function=67, alias=1, constant=0, trace=7, deprecated=1, uncertain=0, private_but_exposed=0.

## scLucid.tumor.cnv

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `assign_cnv_signature` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `calculate_clonal_diversity` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `calculate_cnv_score` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `CloneAnalyzer` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `CNVAnalyzer` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `CNVSigExtractor` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `extract_cnv_signatures` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `find_tumor` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `find_tumor_cells` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `identify_clones` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `infer_clonal_phylogeny` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `infer_cnv` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `plot_aneuploid_proportion` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `plot_cnv_distribution` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `plot_cnv_heatmap` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `plot_per_chromosome_scores` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |
| `run_cnv_analysis` | Function | `src/scLucid/tumor/cnv/__init__.py` |  |

#### Alias

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `identify_clones_from_cnv` | Alias | `src/scLucid/tumor/cnv/__init__.py` | [A] from `clone_analysis.identify_clones` |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 18 symbols (18 stable, 0 flagged). workflow=0, config=0, class=0, function=17, alias=1, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tumor.evolution

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `analyze_dissemination` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |
| `analyze_tumor_progression` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |
| `build_phylogenetic_tree` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |
| `identify_transition_states` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |
| `MetastasisTracker` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |
| `PhylogenyBuilder` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |
| `predict_metastasis_risk` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |
| `ProgressionAnalyzer` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |
| `root_tree` | Function | `src/scLucid/tumor/evolution/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 9 symbols (9 stable, 0 flagged). workflow=0, config=0, class=0, function=9, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tumor.heterogeneity

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `analyze_regional_heterogeneity` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |
| `analyze_treatment_response_trajectory` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |
| `calculate_diversity_indices` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |
| `DiversityAnalyzer` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |
| `estimate_intratumoral_heterogeneity` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |
| `identify_spatial_patterns` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |
| `RegionalAnalyzer` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |
| `TemporalAnalyzer` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |
| `track_temporal_dynamics` | Function | `src/scLucid/tumor/heterogeneity/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 9 symbols (9 stable, 0 flagged). workflow=0, config=0, class=0, function=9, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tumor.malignancy

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `calculate_proliferation_index` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `calculate_stemness_score` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `classify_malignant_cells` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `compare_stemness_between_groups` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `estimate_metastatic_potential` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `identify_cancer_stem_cells` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `MalignancyClassifier` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `MalignancyScorer` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `score_malignancy` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `score_malignancy_potential` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |
| `StemnessAnalyzer` | Function | `src/scLucid/tumor/malignancy/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

#### Deprecated

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `run_malignancy_interpretation` | Deprecated | `src/scLucid/tumor/malignancy/__init__.py` | [D] |

**Summary:** 12 symbols (11 stable, 1 flagged). workflow=0, config=0, class=0, function=11, alias=0, constant=0, trace=0, deprecated=1, uncertain=0, private_but_exposed=0.

## scLucid.tumor.microenvironment

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `analyze_cell_interactions` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `analyze_ecosystem_composition` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `analyze_immune_infiltration` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `calculate_tumor_microenvironment_score` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `compare_ecosystems` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `deconvolve_tme` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `EcosystemAnalyzer` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `estimate_stromal_content` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `find_dominant_interactions` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `InteractionAnalyzer` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `score_immune_interactions` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |
| `TMEProfiler` | Function | `src/scLucid/tumor/microenvironment/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 12 symbols (12 stable, 0 flagged). workflow=0, config=0, class=0, function=12, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tumor.therapy

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `discover_therapeutic_targets` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |
| `identify_resistance_mechanisms` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |
| `predict_therapy_response` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |
| `prioritize_druggable_genes` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |
| `ResistanceAnalyzer` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |
| `ResponsePredictor` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |
| `score_drug_resistance` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |
| `stratify_patients` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |
| `TargetDiscovery` | Function | `src/scLucid/tumor/therapy/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 9 symbols (9 stable, 0 flagged). workflow=0, config=0, class=0, function=9, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.tumor.utils

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `calculate_signature_scores` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `get_drug_targets` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `get_emt_markers` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `get_immune_markers` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `get_proliferation_markers` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `get_stromal_markers` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `get_tumor_markers` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `HallmarkCalculator` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `load_hallmark_signatures` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `query_cancer_gene_census` | Function | `src/scLucid/tumor/utils/__init__.py` |  |
| `query_tcga_data` | Function | `src/scLucid/tumor/utils/__init__.py` |  |

### Deprecated / Uncertain / Private-but-Exposed

*No flagged symbols.*

**Summary:** 11 symbols (11 stable, 0 flagged). workflow=0, config=0, class=0, function=11, alias=0, constant=0, trace=0, deprecated=0, uncertain=0, private_but_exposed=0.

## scLucid.utils

### Stable APIs

#### Function

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `AnalysisContext` | Function | `src/scLucid/utils/__init__.py` |  |
| `api_layer_contract_to_dict` | Function | `src/scLucid/utils/__init__.py` |  |
| `APILayerContract` | Function | `src/scLucid/utils/__init__.py` |  |
| `AssayKeys` | Function | `src/scLucid/utils/__init__.py` |  |
| `assert_analysis_ready` | Function | `src/scLucid/utils/__init__.py` |  |
| `assert_preprocessing_ready` | Function | `src/scLucid/utils/__init__.py` |  |
| `assert_qc_ready` | Function | `src/scLucid/utils/__init__.py` |  |
| `assert_trusted_resources` | Function | `src/scLucid/utils/__init__.py` |  |
| `audit_curation_index` | Function | `src/scLucid/utils/__init__.py` |  |
| `audit_geneset_resources` | Function | `src/scLucid/utils/__init__.py` |  |
| `audit_marker_entry_quality` | Function | `src/scLucid/utils/__init__.py` |  |
| `audit_marker_resources` | Function | `src/scLucid/utils/__init__.py` |  |
| `audit_resource_manifest` | Function | `src/scLucid/utils/__init__.py` |  |
| `BaseWorkflow` | Function | `src/scLucid/utils/__init__.py` |  |
| `BenchmarkRunner` | Function | `src/scLucid/utils/__init__.py` |  |
| `build_config_lineage` | Function | `src/scLucid/utils/__init__.py` |  |
| `build_metadata_dicts` | Function | `src/scLucid/utils/__init__.py` |  |
| `build_qc_preprocess_validation` | Function | `src/scLucid/utils/__init__.py` |  |
| `build_resource_trust_report` | Function | `src/scLucid/utils/__init__.py` |  |
| `canonicalize_marker_label` | Function | `src/scLucid/utils/__init__.py` |  |
| `CellType` | Function | `src/scLucid/utils/__init__.py` |  |
| `check_layer_consistency` | Function | `src/scLucid/utils/__init__.py` |  |
| `classify_literature_resource_utility` | Function | `src/scLucid/utils/__init__.py` |  |
| `clear_sclucid_results` | Function | `src/scLucid/utils/__init__.py` |  |
| `clear_storage` | Function | `src/scLucid/utils/__init__.py` |  |
| `compact_sclucid_uns` | Function | `src/scLucid/utils/__init__.py` |  |
| `ContractError` | Function | `src/scLucid/utils/__init__.py` |  |
| `ContractValidationResult` | Function | `src/scLucid/utils/__init__.py` |  |
| `DatasetProfile` | Function | `src/scLucid/utils/__init__.py` |  |
| `DatasetType` | Function | `src/scLucid/utils/__init__.py` |  |
| `DecisionRecord` | Function | `src/scLucid/utils/__init__.py` |  |
| `effective_n_jobs` | Function | `src/scLucid/utils/__init__.py` |  |
| `ensure_sclucid_namespace` | Function | `src/scLucid/utils/__init__.py` |  |
| `estimate_adata_memory` | Function | `src/scLucid/utils/__init__.py` |  |
| `EvidenceBundle` | Function | `src/scLucid/utils/__init__.py` |  |
| `EvidenceItem` | Function | `src/scLucid/utils/__init__.py` |  |
| `EvidenceLevel` | Function | `src/scLucid/utils/__init__.py` |  |
| `export_audit_report` | Function | `src/scLucid/utils/__init__.py` |  |
| `export_review_summary` | Function | `src/scLucid/utils/__init__.py` |  |
| `filter_by_species` | Function | `src/scLucid/utils/__init__.py` |  |
| `filter_by_tissue_type` | Function | `src/scLucid/utils/__init__.py` |  |
| `filter_marker_dict` | Function | `src/scLucid/utils/__init__.py` |  |
| `finalize_manual_review_summary` | Function | `src/scLucid/utils/__init__.py` |  |
| `flatten_marker_dict` | Function | `src/scLucid/utils/__init__.py` |  |
| `format_contract_error` | Function | `src/scLucid/utils/__init__.py` |  |
| `get_api_layer_spec` | Function | `src/scLucid/utils/__init__.py` |  |
| `get_dataset_info` | Function | `src/scLucid/utils/__init__.py` |  |
| `get_gene_display_aliases` | Function | `src/scLucid/utils/__init__.py` |  |
| `get_marker_aliases` | Function | `src/scLucid/utils/__init__.py` |  |
| `get_marker_manager` | Function | `src/scLucid/utils/__init__.py` |  |
| `get_memory_usage` | Function | `src/scLucid/utils/__init__.py` |  |
| `get_progress_bar` | Function | `src/scLucid/utils/__init__.py` |  |
| `get_storage` | Function | `src/scLucid/utils/__init__.py` |  |
| `has_result` | Function | `src/scLucid/utils/__init__.py` |  |
| `infer_analysis_context` | Function | `src/scLucid/utils/__init__.py` |  |
| `infer_anndata_semantics` | Function | `src/scLucid/utils/__init__.py` |  |
| `infer_dataset_profile` | Function | `src/scLucid/utils/__init__.py` |  |
| `is_ci_environment` | Function | `src/scLucid/utils/__init__.py` |  |
| `is_multi_sample_hint` | Function | `src/scLucid/utils/__init__.py` |  |
| `is_raw_count_matrix` | Function | `src/scLucid/utils/__init__.py` |  |
| `LayerKeys` | Function | `src/scLucid/utils/__init__.py` |  |
| `LayerSemanticKeys` | Function | `src/scLucid/utils/__init__.py` |  |
| `list_results` | Function | `src/scLucid/utils/__init__.py` |  |
| `list_sclucid_modules` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_10x_data` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_all_datasets` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_config` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_gene_set_manager` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_gene_sets` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_luad` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_marker_aliases` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_marker_curation_literature_index` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_melanoma` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_pbmc3k` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_reference_index` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_resource_manifest` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_result` | Function | `src/scLucid/utils/__init__.py` |  |
| `load_workflow_result` | Function | `src/scLucid/utils/__init__.py` |  |
| `Manager` | Function | `src/scLucid/utils/__init__.py` |  |
| `memory_tracker` | Function | `src/scLucid/utils/__init__.py` |  |
| `merge_obs_metadata` | Function | `src/scLucid/utils/__init__.py` |  |
| `merge_partial_results` | Function | `src/scLucid/utils/__init__.py` |  |
| `migrate_legacy_storage` | Function | `src/scLucid/utils/__init__.py` |  |
| `ModalityContractResult` | Function | `src/scLucid/utils/__init__.py` |  |
| `ModalityKeys` | Function | `src/scLucid/utils/__init__.py` |  |
| `model_to_dict` | Function | `src/scLucid/utils/__init__.py` |  |
| `module_namespace` | Function | `src/scLucid/utils/__init__.py` |  |
| `Modules` | Function | `src/scLucid/utils/__init__.py` |  |
| `normalize_dataset_type` | Function | `src/scLucid/utils/__init__.py` |  |
| `normalize_review_summary` | Function | `src/scLucid/utils/__init__.py` |  |
| `ObsKeys` | Function | `src/scLucid/utils/__init__.py` |  |
| `ObsmKeys` | Function | `src/scLucid/utils/__init__.py` |  |
| `PartialResultManager` | Function | `src/scLucid/utils/__init__.py` |  |
| `PerformanceProfiler` | Function | `src/scLucid/utils/__init__.py` |  |
| `PerformanceStats` | Function | `src/scLucid/utils/__init__.py` |  |
| `print_dataset_summary` | Function | `src/scLucid/utils/__init__.py` |  |
| `print_sample_crosstab` | Function | `src/scLucid/utils/__init__.py` |  |
| `profile_function` | Function | `src/scLucid/utils/__init__.py` |  |
| `profile_performance` | Function | `src/scLucid/utils/__init__.py` |  |
| `progress_decorator` | Function | `src/scLucid/utils/__init__.py` |  |
| `read_10x` | Function | `src/scLucid/utils/__init__.py` |  |
| `read_h5ad` | Function | `src/scLucid/utils/__init__.py` |  |
| `record_artifact` | Function | `src/scLucid/utils/__init__.py` |  |
| `record_config_lineage` | Function | `src/scLucid/utils/__init__.py` |  |
| `record_contract_result` | Function | `src/scLucid/utils/__init__.py` |  |
| `record_error` | Function | `src/scLucid/utils/__init__.py` |  |
| `RecoveryError` | Function | `src/scLucid/utils/__init__.py` |  |
| `register_anndata_semantics` | Function | `src/scLucid/utils/__init__.py` |  |
| `ReviewAction` | Function | `src/scLucid/utils/__init__.py` |  |
| `rollup_step_status` | Function | `src/scLucid/utils/__init__.py` |  |
| `run_joblib_or_sequential` | Function | `src/scLucid/utils/__init__.py` |  |
| `sanitize_for_hdf5` | Function | `src/scLucid/utils/__init__.py` |  |
| `save_result` | Function | `src/scLucid/utils/__init__.py` |  |
| `save_workflow_result` | Function | `src/scLucid/utils/__init__.py` |  |
| `setup_runtime_environment` | Function | `src/scLucid/utils/__init__.py` |  |
| `stage_contract_to_dict` | Function | `src/scLucid/utils/__init__.py` |  |
| `StageContract` | Function | `src/scLucid/utils/__init__.py` |  |
| `step_results_from_storage` | Function | `src/scLucid/utils/__init__.py` |  |
| `step_results_to_storage` | Function | `src/scLucid/utils/__init__.py` |  |
| `StepError` | Function | `src/scLucid/utils/__init__.py` |  |
| `StepResult` | Function | `src/scLucid/utils/__init__.py` |  |
| `StepStatus` | Function | `src/scLucid/utils/__init__.py` |  |
| `subset_adata` | Function | `src/scLucid/utils/__init__.py` |  |
| `subset_from_annotations` | Function | `src/scLucid/utils/__init__.py` |  |
| `summarize_step_results` | Function | `src/scLucid/utils/__init__.py` |  |
| `UnsKeys` | Function | `src/scLucid/utils/__init__.py` |  |
| `use_layer_as_X` | Function | `src/scLucid/utils/__init__.py` |  |
| `validation_table_to_dataframe` | Function | `src/scLucid/utils/__init__.py` |  |
| `ValidationError` | Function | `src/scLucid/utils/__init__.py` |  |
| `VarKeys` | Function | `src/scLucid/utils/__init__.py` |  |
| `with_error_recovery` | Function | `src/scLucid/utils/__init__.py` |  |
| `WorkflowCheckpoint` | Function | `src/scLucid/utils/__init__.py` |  |
| `WorkflowError` | Function | `src/scLucid/utils/__init__.py` |  |
| `WorkflowStepIterator` | Function | `src/scLucid/utils/__init__.py` |  |
| `write_h5ad_safe` | Function | `src/scLucid/utils/__init__.py` |  |
| `write_validation_outputs` | Function | `src/scLucid/utils/__init__.py` |  |

#### Constant

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `API_LAYER_ORDER` | Constant | `src/scLucid/utils/__init__.py` |  |
| `COMPARATIVE_READINESS_LABEL` | Constant | `src/scLucid/utils/__init__.py` |  |
| `KNOWN_SPECIES` | Constant | `src/scLucid/utils/__init__.py` |  |
| `MARKER_FORMATS` | Constant | `src/scLucid/utils/__init__.py` |  |
| `REVIEW_SUMMARY_RECOMMENDED_KEYS` | Constant | `src/scLucid/utils/__init__.py` |  |
| `SCHEMA_VERSION` | Constant | `src/scLucid/utils/__init__.py` |  |
| `SCLUCID_ROOT` | Constant | `src/scLucid/utils/__init__.py` |  |
| `STAGE_ORDER` | Constant | `src/scLucid/utils/__init__.py` |  |
| `STORAGE_ROOT` | Constant | `src/scLucid/utils/__init__.py` |  |
| `VALID_MODULES` | Constant | `src/scLucid/utils/__init__.py` |  |
| `VALIDATION_SCOPE` | Constant | `src/scLucid/utils/__init__.py` |  |

#### Trace / Contract

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `API_LAYER_CONTRACTS` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `EVIDENCE_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `get_contract_spec` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `get_minimal_workflow_contract` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `get_stage_contract` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `MINIMAL_WORKFLOW_CONTRACT` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `REVIEW_SUMMARY_REQUIRED_KEYS` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `STAGE_CONTRACTS` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `validate_adata` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `validate_all_stage_contracts` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `validate_analysis_results` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `validate_config` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `validate_modality_contract` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `validate_review_summary_schema` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `validate_stage_contract` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `validate_workflow_contract` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |
| `VALIDATION_SCAFFOLD_SCHEMA_VERSION` | Trace / Contract | `src/scLucid/utils/__init__.py` | [T] |

### Deprecated / Uncertain / Private-but-Exposed

#### Private-but-Exposed

| Symbol | Kind | Source | Notes |
|--------|------|--------|-------|
| `_get_cancer_markers` | Private-but-Exposed | `src/scLucid/utils/__init__.py` | [P] |

**Summary:** 165 symbols (164 stable, 1 flagged). workflow=0, config=0, class=0, function=136, alias=0, constant=11, trace=17, deprecated=0, uncertain=0, private_but_exposed=1.

## Global Summary

| Kind | Count |
|------|-------|
| Workflow Orchestrator | 12 |
| Config Class | 81 |
| Class | 0 |
| Function | 772 |
| Alias | 21 |
| Constant | 18 |
| Trace / Contract | 66 |
| Deprecated | 7 |
| Uncertain | 0 |
| Private-but-Exposed | 1 |
| **Total** | **978** |

<!-- AUTO-GENERATED INVENTORY END -->
