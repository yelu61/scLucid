# Analysis Module Specification

## Purpose

Provide tools for single-cell data analysis, including cell clustering, annotation, differential expression, enrichment analysis, and cell type proportion statistics.

## Requirements

### Requirement: Cell Clustering
The system SHALL perform cell clustering using multiple algorithms.

#### Scenario: Leiden clustering
- **GIVEN** AnnData with neighborhood graph and `ClusteringConfig(method="leiden", resolution=1.0)`
- **WHEN** calling `cluster_cells(adata, config)`
- **THEN** system:
  - Runs Leiden clustering
  - Stores cluster labels in `adata.obs[config.key_added]`
  - Returns clustering results

#### Scenario: Louvain clustering
- **GIVEN** `ClusteringConfig(method="louvain", resolution=0.8)`
- **WHEN** calling `cluster_cells(adata, config)`
- **THEN** system runs Louvain clustering with specified resolution

#### Scenario: K-means clustering
- **GIVEN** `ClusteringConfig(method="kmeans", n_clusters=10)`
- **WHEN** calling `cluster_cells(adata, config)`
- **THEN** system:
  - Runs k-means on specified embedding
  - Requires `n_clusters` parameter

#### Scenario: HDBSCAN clustering
- **GIVEN** `ClusteringConfig(method="hdbscan")`
- **WHEN** calling `cluster_cells(adata, config)`
- **THEN** system runs density-based clustering

#### Scenario: Invalid clustering method
- **GIVEN** `ClusteringConfig(method="invalid")`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (method must be one of: leiden, louvain, kmeans, hdbscan)

#### Scenario: Invalid resolution raises error
- **GIVEN** `ClusteringConfig(method="leiden", resolution=-1.0)`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (resolution must be > 0)

#### Scenario: K-means requires n_clusters
- **GIVEN** `ClusteringConfig(method="kmeans", n_clusters=None)`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (n_clusters required for kmeans)

### Requirement: Resolution Optimization
The system SHALL optimize clustering resolution using multiple metrics.

#### Scenario: Find optimal resolution
- **GIVEN** AnnData with clusters and `ResolutionSearchConfig`
- **WHEN** calling `find_resolution(adata, config)`
- **THEN** system:
  - Tests resolutions in specified range
  - Computes metrics: silhouette, marker abundance, stability
  - Suggests optimal resolution

#### Scenario: Configure resolution search
- **GIVEN** `ResolutionSearchConfig(resolution_range=(0.2, 2.0, 10))`
- **WHEN** running resolution search
- **THEN** system tests 10 resolutions between 0.2 and 2.0

### Requirement: Cluster Merging
The system SHALL merge similar clusters based on marker gene overlap.

#### Scenario: Merge clusters
- **GIVEN** AnnData with clusters and `MergeClustersConfig`
- **WHEN** calling `merge_clusters(adata, config)`
- **THEN** system:
  - Identifies clusters with high marker overlap
  - Merges similar clusters
  - Updates cluster labels

### Requirement: Cell Type Annotation
The system SHALL annotate cell types using multiple methods.

#### Scenario: Manual marker-based annotation
- **GIVEN** a marker database and `AnnotationConfig`
- **WHEN** calling `annotate_clusters(adata, marker_manager, config)`
- **THEN** system:
  - Scores clusters for each cell type
  - Assigns labels based on scores
  - Stores labels in `adata.obs`

#### Scenario: CellTypist annotation
- **GIVEN** `AnnotationConfig(method="celltypist", model="immune_human")`
- **WHEN** calling `run_celltypist(adata, config)`
- **THEN** system:
  - Runs CellTypist prediction
  - Returns predicted labels and confidence scores

#### Scenario: Ensemble annotation
- **GIVEN** `AnnotationConfig(method="ensemble", methods=["markers", "celltypist"])`
- **WHEN** calling `run_annotation(adata, config)`
- **THEN** system:
  - Combines predictions from multiple methods
  - Resolves conflicts using confidence scores

#### Scenario: Invalid annotation method
- **GIVEN** `AnnotationConfig(method="invalid")`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (method must be valid)

#### Scenario: Invalid confidence range
- **GIVEN** `AnnotationConfig(min_confidence=1.5)`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (min_confidence must be <= 1.0)

### Requirement: Cell Type Scoring
The system SHALL score cells for cell type signatures.

#### Scenario: Score cells with markers
- **GIVEN** AnnData and marker gene list
- **WHEN** calling `score_cell_types(adata, markers, "T cells")`
- **THEN** system:
  - Calculates average expression of marker genes
  - Stores scores in `adata.obs['score_T_cells']`

#### Scenario: Configure scoring
- **GIVEN** `ScoringConfig(use_raw=True, score_method="mean")`
- **WHEN** scoring cells
- **THEN** system uses `.raw` data and mean scoring method

### Requirement: Differential Expression Analysis
The system SHALL perform differential expression analysis between groups.

#### Scenario: Find cluster markers
- **GIVEN** AnnData with clusters and `DifferentialConfig`
- **WHEN** calling `find_markers(adata, groupby="cluster", config=config)`
- **THEN** system:
  - Performs DE analysis (Wilcoxon by default)
  - Returns marker genes with p-values and fold changes
  - Stores results in `adata.uns['sclucid']['analysis']['de']`

#### Scenario: Compare two groups
- **GIVEN** `CompareGroupsConfig(group1="cluster_0", group2="cluster_1")`
- **WHEN** calling `compare_groups(adata, config)`
- **THEN** system performs DE between specified groups

#### Scenario: Compare conditions
- **GIVEN** samples from multiple conditions
- **WHEN** calling `compare_conditions(adata, condition_col="condition")`
- **THEN** system performs condition-level DE analysis

#### Scenario: Filter markers
- **GIVEN** DE results and `FilterMarkersConfig`
- **WHEN** calling `filter_markers(results, config)`
- **THEN** system:
  - Filters by p-value threshold (`max_padj`)
  - Filters by log2FC threshold (`min_log2fc`)
  - Filters by expression percentage (`min_pct`)
  - Returns filtered markers

#### Scenario: Find conserved markers
- **GIVEN** multiple batches and `ConservedMarkersConfig`
- **WHEN** calling `get_conserved_markers(adata, config)`
- **THEN** system:
  - Finds markers consistent across batches
  - Returns conserved marker genes

#### Scenario: Invalid DE method raises error
- **GIVEN** `DifferentialConfig(method="invalid")`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (method must be valid)

### Requirement: Gene Set Enrichment Analysis
The system SHALL perform enrichment analysis on gene sets.

#### Scenario: GSEA enrichment
- **GIVEN** gene list and `EnrichmentConfig(method="gsea")`
- **WHEN** calling `run_enrichment(gene_list, config)`
- **THEN** system:
  - Runs GSEA with specified permutations
  - Returns enriched pathways with NES and p-values

#### Scenario: Over-representation analysis
- **GIVEN** `EnrichmentConfig(method="ora", database="go_bp")`
- **WHEN** running enrichment
- **THEN** system performs ORA using GO Biological Process database

#### Scenario: Configure enrichment parameters
- **GIVEN** `EnrichmentConfig(max_padj=0.05, min_size=15, max_size=500)`
- **WHEN** running enrichment
- **THEN** system:
  - Filters results by adjusted p-value
  - Filters pathways by size range

#### Scenario: Invalid permutation count
- **GIVEN** `EnrichmentConfig(method="gsea", n_perm=10)`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (n_perm must be >= 100)

### Requirement: Cell Type Proportion Analysis
The system SHALL analyze cell type proportions across samples.

#### Scenario: Calculate proportions
- **GIVEN** AnnData with samples and cell types
- **WHEN** calling function with `ProportionConfig`
- **THEN** system:
  - Calculates cell type proportions per sample
  - Returns proportion table
  - Optionally performs statistical tests

### Requirement: Marker Gene Summarization
The system SHALL summarize marker genes and enrichment results.

#### Scenario: Summarize markers
- **GIVEN** DE and enrichment results
- **WHEN** calling `summarize_markers_and_enrichment(adata)`
- **THEN** system:
  - Combines marker genes and pathway enrichment
  - Returns comprehensive summary table
  - Generates visualizations if configured

### Requirement: Cluster Characterization
The system SHALL provide comprehensive cluster characterization.

#### Scenario: Characterize clusters
- **GIVEN** AnnData with clusters
- **WHEN** calling `characterize_clusters(adata, marker_manager)`
- **THEN** system:
  - Finds cluster markers
  - Performs enrichment analysis
  - Annotates cell types
  - Returns characterization summary

### Requirement: Configuration Serialization
All analysis configs MUST support serialization.

#### Scenario: Serialize clustering config
- **GIVEN** a `ClusteringConfig` instance
- **WHEN** calling `config.to_dict()`
- **THEN** returns dictionary with all clustering parameters

#### Scenario: Serialize DE config
- **GIVEN** a `DifferentialConfig` instance
- **WHEN** storing to `adata.uns`
- **THEN** config is serialized via `to_dict()`

### Requirement: Analysis Workflow
The system SHALL provide integrated analysis workflow.

#### Scenario: Run complete analysis
- **GIVEN** preprocessed AnnData
- **WHEN** using `AnalysisWorkflowConfig`
- **THEN** system can execute:
  1. Clustering
  2. Resolution optimization (optional)
  3. Marker finding
  4. Cell type annotation
  5. Differential expression
  6. Enrichment analysis
  7. Proportion analysis
