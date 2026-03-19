# Preprocessing Module Specification

## Purpose

Provide complete preprocessing workflow for single-cell RNA-seq data, including normalization, HVG selection, scaling, dimensionality reduction, batch correction, and UMAP embedding.

## Requirements

### Requirement: Data Normalization
The system SHALL normalize raw UMI counts using library size normalization and log1p transformation.

#### Scenario: Standard normalization
- **GIVEN** an AnnData object with raw counts in `.X`
- **WHEN** calling `normalize_data(adata, config)`
- **THEN** system:
  - Normalizes to `target_sum` (default: 1e4)
  - Applies log1p transformation
  - Stores result in `adata.layers[config.normalized_layer]` (default: "normalized")
  - Preserves raw counts in `adata.layers["counts"]`

#### Scenario: Configure normalization parameters
- **GIVEN** a `NormalizationConfig` with `target_sum=1e5, max_fraction=0.05`
- **WHEN** running normalization
- **THEN** system:
  - Normalizes to 1e5 counts per cell
  - Clips genes expressing in >5% of cells

#### Scenario: Invalid target_sum raises error
- **GIVEN** a `NormalizationConfig` with `target_sum=-1`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (target_sum must be > 0)

#### Scenario: Reserved layer name protection
- **GIVEN** a `NormalizationConfig` with `output_layer="counts"`
- **WHEN** instantiating the config
- **THEN** system raises `ValueError` ("counts" is reserved for raw counts)

### Requirement: Highly Variable Gene Selection
The system SHALL identify highly variable genes using dispersion-based methods.

#### Scenario: Find HVGs using Seurat method
- **GIVEN** normalized data and `HVGConfig` with `method="seurat"`
- **WHEN** calling `find_hvgs(adata, config)`
- **THEN** system:
  - Calculates dispersion for each gene
  - Selects top `n_top_genes` (default: 2000)
  - Stores HVGs in `adata.var['highly_variable']`
  - Saves results with key `output_key` in `adata.uns['sclucid']['preprocess']['hvg']`

#### Scenario: Find HVGs with batch-specific mode
- **GIVEN** data with multiple batches and `batch_key="batch"`
- **WHEN** using `method="seurat_batch"`
- **THEN** system calculates batch-specific dispersions

#### Scenario: Configure HVG parameters
- **GIVEN** an `HVGConfig` with `n_top_genes=3000, min_mean=0.0125, max_mean=3`
- **WHEN** running HVG selection
- **THEN** system:
  - Selects 3000 genes
  - Filters genes by mean expression range
  - Uses specified `span` for loess smoothing

#### Scenario: Invalid n_top_genes raises error
- **GIVEN** an `HVGConfig` with `n_top_genes=-1`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (n_top_genes must be > 0)

#### Scenario: Invalid span value raises error
- **GIVEN** an `HVGConfig` with `span=2.0`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (span must be <= 1.0)

### Requirement: Subset to HVGs
The system SHALL subset AnnData to highly variable genes.

#### Scenario: Subset to HVGs
- **GIVEN** an AnnData with HVGs identified
- **WHEN** calling `select_hvg_sets(adata, hvg_keys=['hvg'], subset=True)`
- **THEN** system:
  - Filters `adata.var` to HVGs only
  - Updates `adata.X` and layers accordingly
  - Preserves `.raw` if `keep_raw=True`

#### Scenario: Keep raw data in .raw
- **GIVEN** an AnnData before subsetting
- **WHEN** calling `select_hvg_sets` with `keep_raw=True`
- **THEN** system stores full data in `adata.raw` before subsetting

### Requirement: Data Scaling
The system SHALL scale data to zero mean and unit variance.

#### Scenario: Scale data
- **GIVEN** normalized data and `ScalingConfig`
- **WHEN** calling `scale_data(adata, config)`
- **THEN** system:
  - Z-score scales each gene
  - Clips values to `max_value` (default: 10)
  - Stores result in `adata.layers[config.scaled_layer]`

#### Scenario: Regress out covariates
- **GIVEN** a `ScalingConfig` with `vars_to_regress=['pct_counts_mt', 'phase']`
- **WHEN** calling `regress_out(adata, config)`
- **THEN** system:
  - Regresses out specified covariates
  - Stores result in `config.regressed_layer`
  - Preserves normalized data in original layer

#### Scenario: Scale after regression
- **GIVEN** regressed data
- **WHEN** calling `scale_data`
- **THEN** system scales regressed values

### Requirement: Dimensionality Reduction (PCA)
The system SHALL perform Principal Component Analysis.

#### Scenario: Run PCA
- **GIVEN** scaled data and `GraphConfig` with `n_pcs=50`
- **WHEN** calling `sc.tl.pca(adata, n_comps=n_pcs)`
- **THEN** system:
  - Computes PCA
  - Stores PCs in `adata.obsm['X_pca']`
  - Stores variance in `adata.uns['pca']`
  - Stores loadings in `adata.varm['PCs']`

#### Scenario: Configure PCA parameters
- **GIVEN** `GraphConfig` with `n_pcs=100, use_highly_variable=True`
- **WHEN** running PCA
- **THEN** system uses only HVGs for PCA and computes 100 components

### Requirement: Batch Correction
The system SHALL integrate data from multiple batches using various methods.

#### Scenario: Harmony integration
- **GIVEN** an AnnData with PCA and `IntegrationConfig` with `method="harmony", batch_key="batch"`
- **WHEN** calling `batch_correction(adata, config)`
- **THEN** system:
  - Runs Harmony integration
  - Stores corrected embeddings in `adata.obsm[output_key]`
  - Returns corrected AnnData

#### Scenario: Scanorama integration
- **GIVEN** multiple batches and `method="scanorama"`
- **WHEN** running batch correction
- **THEN** system:
  - Corrects PCs using Scanorama
  - Returns integrated AnnData

#### Scenario: BBKNN integration
- **GIVEN** `method="bbknn"`
- **WHEN** running batch correction
- **THEN** system:
  - Constructs batch-balanced neighborhood graph
  - Does not modify `adata.obsm`

#### Scenario: Skip integration if not needed
- **GIVEN** `IntegrationConfig` with `method=None`
- **WHEN** running preprocessing workflow
- **THEN** system skips batch correction step

### Requirement: Neighborhood Graph Construction
The system SHALL construct k-nearest neighbor graph for clustering.

#### Scenario: Build neighborhood graph
- **GIVEN** PCA (or integrated) embeddings
- **WHEN** calling `sc.pp.neighbors(adata, n_pcs=50, n_neighbors=15)`
- **THEN** system:
  - Constructs KNN graph
  - Stores in `adata.obsp['connectivities']` and `adata.obsp['distances']`
  - Stores parameters in `adata.uns['neighbors']`

#### Scenario: Use integrated embeddings
- **GIVEN** batch-corrected data with `use_rep="X_harmony"`
- **WHEN** building neighbors
- **THEN** system uses Harmony embeddings instead of PCA

### Requirement: UMAP Embedding
The system SHALL compute 2D UMAP embedding for visualization.

#### Scenario: Compute UMAP
- **GIVEN** neighborhood graph
- **WHEN** calling `sc.tl.umap(adata)`
- **THEN** system:
  - Computes UMAP coordinates
  - Stores in `adata.obsm['X_umap']`

### Requirement: Preprocessing Workflow
The system SHALL provide end-to-end preprocessing workflow.

#### Scenario: Run complete workflow
- **GIVEN** AnnData with raw counts
- **WHEN** calling `run_preprocessing(adata, config, results_dir)`
- **THEN** system executes in order:
  1. Normalization
  2. Store normalized data in `.raw`
  3. Regression (if vars_to_regress specified)
  4. HVG selection
  5. Subset to HVGs
  6. Scaling
  7. PCA
  8. Batch correction (if method specified)
  9. Neighbors graph
  10. UMAP

#### Scenario: Workflow with checkpointing
- **GIVEN** `run_preprocessing_v2` with checkpointing enabled
- **WHEN** workflow is interrupted
- **THEN** system:
  - Saves intermediate results
  - Can resume from last checkpoint
  - Stores config in `adata.uns['sclucid']['preprocess']['workflow_config']`

### Requirement: Configuration Serialization
All preprocessing configs MUST support serialization.

#### Scenario: Serialize workflow config
- **GIVEN** a `WorkflowConfig` instance
- **WHEN** storing to `adata.uns`
- **THEN** system calls `config.to_dict()` for serialization

#### Scenario: Save/load sub-configs
- **GIVEN** a `NormalizationConfig` instance
- **WHEN** calling `to_json_file()` and `from_json_file()`
- **THEN** config is saved and loaded with validation

### Requirement: Layer Naming Convention
The system SHALL follow consistent layer naming conventions.

#### Scenario: Standard layer names
- **GIVEN** preprocessing workflow
- **WHEN** storing intermediate results
- **THEN** layers are named:
  - `counts`: Raw UMI counts (preserved)
  - `normalized`: Log1p normalized data
  - `regressed`: Regression-corrected data (if applicable)
  - `scaled`: Z-score scaled data

#### Scenario: Custom layer names
- **GIVEN** config with custom `normalized_layer="norm"`
- **WHEN** running normalization
- **THEN** system uses custom layer name "norm"
