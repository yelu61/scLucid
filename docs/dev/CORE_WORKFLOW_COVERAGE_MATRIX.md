# Core Workflow Coverage Matrix

This document maps the six core workflow stages of scLucid to their entry
functions, supporting APIs, configuration classes, review/trace functions, and
known gaps. It is manually curated; update it when adding new workflow stages
or reorganizing modules.

**Related:** [Function Inventory](FUNCTION_INVENTORY.md)

---

## Stage Overview

| Stage | Entry Function | Config Class | Status |
|-------|---------------|--------------|--------|
| QC | [`run_standard_qc`](FUNCTION_INVENTORY.md#scLucidqc) | `QCWorkflowConfig` | Stable |
| Preprocess | [`run_preprocessing`](FUNCTION_INVENTORY.md#scLucidpreprocess) | `PreprocessingWorkflowConfig` | Stable |
| Analysis | [`run_standard_analysis`](FUNCTION_INVENTORY.md#scLucidanalysis) | `AnalysisWorkflowConfig` | Stable |
| Annotation | [`run_annotation`](FUNCTION_INVENTORY.md#scLucidanalysisannotation) | `AnnotationConfig` | Stable |
| Tumor Interpretation | [`run_tumor_analysis`](FUNCTION_INVENTORY.md#scLucidtumor) | `TumorAnalysisConfig` | Evolving |
| Ecosystem Modeling | (none yet) | — | Emerging |

---

## QC Stage

### Entry Functions
- [`run_standard_qc`](FUNCTION_INVENTORY.md#scLucidqc) — Standard QC workflow
- [`run_qc_threshold_decision`](FUNCTION_INVENTORY.md#scLucidqc) — Reusable threshold decision / marking / optional filtering layer for interactive QC
- [`run_advanced_qc`](FUNCTION_INVENTORY.md#scLucidqc) — Advanced QC with adaptive thresholds *(deprecated alias for `run_standard_qc` with a full config)*

### Key Supporting Functions

| Function | Role | Status |
|----------|------|--------|
| `calculate_qc_metric` | Core metric computation | Stable |
| `suggest_qc_thresholds` | Distribution-based threshold recommendation | Stable |
| `resolve_qc_thresholds` | Merge intelligent / MAD / manual thresholds | Stable |
| `predict_doublets` | Doublet detection | Stable |
| `audit_doublets` | Post-filter doublet audit | Stable |
| `diagnose_ambient_rna` | Ambient RNA diagnosis | Stable |
| `mark_low_quality_cell` | Low-quality marking (single-threshold) | Stable |
| `mark_low_quality_cells_adaptive` | Adaptive low-quality marking | Stable |
| `filter_cells` | Final filtering | Stable |
| `audit_filtering` | Retention audit before/after filtering | Stable |
| `generate_qc_report` | QC report generation | Stable |

### Config Classes
- `QCWorkflowConfig`
- `QCThresholds`
- `DoubletConfig`
- `FilterConfig`
- `MarkerConfig`
- `MetricsReportingConfig`

### Review / Trace Functions
- `build_qc_decision_table`
- `enrich_qc_review_summary`
- `validate_qc_review_summary`
- `summarize_qc_review_summary`
- `validate_qc_module_completeness`

### Gaps / TODO
- [x] No dedicated `run_qc_review` workflow orchestrator (review is currently manual). *Addressed by `run_qc_threshold_decision` for threshold/marking decisions; full review summary is automated.*
- [ ] Benchmark utilities (`build_qc_benchmark_assessment`) are public but not integrated into `run_standard_qc`.
- [x] Intelligent QC (`recommend_intelligent_qc`) is optional and may be promoted to a core workflow step. *Addressed by `resolve_qc_thresholds` being wired into threshold application when no explicit thresholds are provided.*

---

## Preprocessing Stage

### Entry Function
- [`run_preprocessing`](FUNCTION_INVENTORY.md#scLucidpreprocess)

### Key Supporting Functions

| Function | Role | Status |
|----------|------|--------|
| `normalize_data` | Normalization | Stable |
| `find_hvgs` | HVG selection | Stable |
| `select_hvg_sets` | HVG set operations (intersection/union/difference) | Stable |
| `select_and_audit_hvgs` | HVG selection + stability + audit | Stable |
| `scale_data` | Scaling | Stable |
| `batch_correction` | Batch correction | Stable |
| `detect_integration_confounding` | Batch-biology confounding check | Stable |
| `diagnose_integration_risk` | Pre-integration risk assessment | Stable |
| `decide_integration` | Auto/skip integration decision | Stable |
| `evaluate_integration` | Post-integration quality metrics | Stable |
| `regress_out` | Covariate regression | Stable |
| `optimize_neighbors_pcs` | Neighbor/PCA optimization | Stable |
| `run_embedding_pipeline` | Optimized graph + named UMAP generation | Stable |
| `is_raw_count_matrix` | Raw-count semantic guard | Stable |
| `build_metadata_dicts` | Build `metadata_dicts` for multi-sample loaders | Stable |

### Config Classes
- `PreprocessingWorkflowConfig`
- `NormalizationConfig`
- `AdaptiveNormalizationConfig`
- `HVGConfig`
- `IntegrationConfig`
- `ScalingConfig`
- `NeighborsConfig`
- `GraphConfig`
- `IntelligentPreprocessConfig`
- `PreprocessingStrategy`

### Review / Trace Functions
- `enrich_preprocessing_review_summary`
- `validate_preprocessing_review_summary`
- `summarize_preprocess_review_summary`
- `validate_preprocess_module_completeness`

### Gaps / TODO
- [ ] `run_intelligent_preprocessing` is public but not wired into `run_preprocessing`.
- [ ] Backend abstraction (`PreprocessingBackend`, `ScanpyBackend`, `RapidsBackend`) is public but under-documented.
- [x] Gene biotype helpers are stable but optional in many workflows. *Stable helpers (`annotate_gene_biotypes`, `filter_genes_by_biotype`) exist; not yet default in `run_preprocessing`.*
- [x] Integration risk diagnostics are public but not wired into `run_preprocessing`. *Addressed by `IntegrationConfig.auto_decide` and `.evaluate`, which `batch_correction` honors.*

---

## Analysis Stage

### Entry Functions
- [`run_standard_analysis`](FUNCTION_INVENTORY.md#scLucidanalysis)
- [`run_custom_analysis`](FUNCTION_INVENTORY.md#scLucidanalysis)

### Key Supporting Functions

| Function | Role | Status |
|----------|------|--------|
| `cluster_cells` | Clustering | Stable |
| `find_markers` | Marker discovery | Stable |
| `annotate_clusters` | Cell type annotation | Stable |
| `score_cell_types` | Annotation scoring | Stable |
| `compare_groups` | Group comparison | Stable |
| `run_enrichment` | Functional enrichment | Stable |
| `score_by_gene_sets` | Gene-set scoring | Stable |

### Config Classes
- `AnalysisWorkflowConfig`
- `ClusteringConfig`
- `AnnotationConfig`
- `DifferentialConfig`
- `ScoringConfig`
- `PseudobulkDEConfig`
- `EnrichmentConfig`

### Review / Trace Functions
- `build_posthoc_qc_review_summary`
- `validate_analysis_review_summary`
- `summarize_analysis_review_summary`
- `validate_analysis_module_completeness`

### Gaps / TODO
- [ ] `run_malignancy_interpretation` is deprecated in `scLucid.analysis`; canonical location is `scLucid.tumor`.
- [ ] Proportion analysis (`analyze_celltype_proportion`) is large but not in the default workflow.
- [ ] Bulk shim (`analysis.bulk`) is legacy; canonical location is `tools.bulk`.

---

## Annotation Stage

### Entry Function
- [`run_annotation`](FUNCTION_INVENTORY.md#scLucidanalysisannotation)

### Key Supporting Functions

| Function | Role | Status |
|----------|------|--------|
| `annotate_clusters` | Core annotation | Stable |
| `score_cell_types` | Scoring support | Stable |
| `run_celltypist` | CellTypist integration | Stable |
| `transfer_labels` | Label transfer | Stable |
| `evaluate_annotation` | Quality evaluation | Stable |
| `build_annotation_consensus` | Consensus building | Stable |
| `flag_suspect_clusters` | Quality flagging | Stable |

### Config Classes
- `AnnotationConfig`

### Review / Trace Functions
- `summarize_annotation_evidence`
- `build_annotation_review_table`
- `build_annotation_consensus`

### Gaps / TODO
- [ ] LLM annotation bundle (`build_llm_annotation_bundle`) is experimental.
- [ ] Subset annotation refinement (`run_subset_annotation_refinement`) is powerful but complex.
- [ ] Lineage-state annotation (`run_lineage_state_annotation`) needs clearer integration examples.

---

## Tumor Interpretation Stage

### Entry Functions
- [`run_tumor_analysis`](FUNCTION_INVENTORY.md#scLucidtumor) — Unified tumor workflow
- [`run_malignancy_interpretation`](FUNCTION_INVENTORY.md#scLucidtumor) — Malignancy interpretation

### Key Supporting Functions

| Function | Role | Status |
|----------|------|--------|
| `score_malignancy` | Malignancy scoring | Stable |
| `classify_malignant_cells` | Classification | Stable |
| `infer_cnv` | CNV inference | Stable |
| `deconvolve_tme` | TME deconvolution | Stable |
| `predict_therapy_response` | Therapy prediction | Evolving |
| `analyze_cell_interactions` | Cell-cell interaction | Evolving |
| `calculate_diversity_indices` | Heterogeneity indices | Evolving |

### Config Classes
- `TumorAnalysisConfig`
- `TumorWorkflowConfig`

### Review / Trace Functions
- `enrich_tumor_review_summary`
- `validate_tumor_review_summary`
- `validate_tumor_module_completeness`

### Gaps / TODO
- [ ] `AnalysisStep` adapters in `tumor.steps` are public but their role is unclear to users (see [API_TRIAGE.md](API_TRIAGE.md#t002-tumor-analysisstep-adapters)).
- [ ] Therapy prediction is marked `exploratory` evidence level.
- [ ] Evolution tracking (`build_phylogenetic_tree`) is not integrated into `run_tumor_analysis`.
- [ ] No unified tumor report generator beyond the review summary helpers.

---

## Ecosystem Modeling Stage

### Entry Function
- None yet — this stage is emerging.

### Key Supporting Functions

| Function | Role | Status |
|----------|------|--------|
| `analyze_ecosystem_composition` | Ecosystem profiling | Emerging |
| `compare_ecosystems` | Cross-sample comparison | Emerging |
| `score_immune_interactions` | Interaction scoring | Emerging |
| `analyze_cell_interactions` | Cell-cell interaction | Emerging |
| `find_dominant_interactions` | Dominant interaction detection | Emerging |

### Config Classes
- None dedicated yet.

### Review / Trace Functions
- None dedicated yet.

### Gaps / TODO
- [ ] No unified workflow orchestrator for ecosystem modeling.
- [ ] CellPhoneDB integration (`tools.cellphonedb`) is external-tool dependent.
- [ ] Spatial ecosystem analysis (`tools.spatial`) is separate from tumor TME.
- [ ] Ecotype / microenvironment archetype APIs are not yet defined.

---

## Cross-Cutting Concerns

| Concern | Current Coverage | Notes |
|---------|------------------|-------|
| Audit / evidence bundle | Strong | `adata.uns["sclucid"]` + review summary + HTML report |
| Config serialization | Strong | Pydantic-based config classes |
| Optional dependency handling | Strong | `_export()` and `_import_optional()` patterns |
| Public API discoverability | Improving | Function inventory now automated |
| Core workflow integration | Partial | Some advanced helpers are public but not wired into default workflows |
