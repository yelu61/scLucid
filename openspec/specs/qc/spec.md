# Quality Control (QC) Module Specification

## Purpose

Provide comprehensive quality control workflows for single-cell RNA-seq data, including metric calculation, doublet detection, cell filtering, and cell cycle scoring.

## Requirements

### Requirement: QC Metrics Calculation
The system SHALL calculate standard QC metrics for single-cell data.

#### Scenario: Calculate basic metrics
- **GIVEN** an AnnData object with raw counts
- **WHEN** calling `calculate_qc_metric(adata, config)`
- **THEN** system adds QC metrics to `adata.obs`:
  - `n_genes_by_counts`: Number of genes detected per cell
  - `total_counts`: Total UMI counts per cell
  - `pct_counts_mt`: Percentage of mitochondrial counts
  - `pct_counts_ribo`: Percentage of ribosomal counts
  - Other standard metrics

#### Scenario: Configure metric reporting
- **GIVEN** a `MetricsReportingConfig` instance
- **WHEN** setting `plot_violin=False` and `export_stats=True`
- **THEN** system skips violin plots and exports CSV statistics

### Requirement: QC Filtering Thresholds
The system SHALL validate and apply QC filtering thresholds using Pydantic.

#### Scenario: Create valid thresholds
- **GIVEN** a `QCThresholds` config with `min_genes=200, pc_mt=20.0`
- **WHEN** instantiating the config
- **THEN** validation succeeds and config is created

#### Scenario: Invalid threshold raises error
- **GIVEN** a `QCThresholds` config with `min_genes=-1`
- **WHEN** instantiating the config
- **THEN** system raises `pydantic.ValidationError` (field must be >= 0)

#### Scenario: Thresholds must be consistent
- **GIVEN** a `QCThresholds` config with `min_genes=500, max_genes=200`
- **WHEN** instantiating the config
- **THEN** system raises `ValueError` (min_genes cannot be > max_genes)

#### Scenario: Percentage thresholds bounded
- **GIVEN** a `QCThresholds` config with `pc_mt=150.0`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (percentage must be 0-100)

### Requirement: Doublet Detection
The system SHALL detect doublets using multiple algorithms.

#### Scenario: Scrublet doublet detection
- **GIVEN** a `DoubletConfig` with `method="scrublet"`
- **WHEN** calling `predict_doublets(adata, config)`
- **THEN** system:
  - Runs Scrublet algorithm
  - Adds `doublet_score` and `predicted_doublet` to `adata.obs`
  - Returns filtered AnnData

#### Scenario: Doublet detection with custom parameters
- **GIVEN** a `DoubletConfig` with `expected_doublet_rate=0.1`
- **WHEN** running doublet detection
- **THEN** system uses specified expected rate for algorithm

#### Scenario: Invalid doublet method
- **GIVEN** a `DoubletConfig` with `method="invalid_method"`
- **WHEN** instantiating the config
- **THEN** system raises `ValidationError` (method must be one of: scrublet, solo, doubletdetection)

### Requirement: Cell Filtering
The system SHALL filter cells based on QC thresholds.

#### Scenario: Filter by thresholds
- **GIVEN** an AnnData with QC metrics and `QCThresholds` config
- **WHEN** calling `filter_cells(adata, thresholds)`
- **THEN** system:
  - Identifies cells passing all thresholds
  - Returns filtered AnnData
  - Stores filtering results in `adata.uns['qc']`

#### Scenario: Adaptive threshold learning
- **GIVEN** an AnnData with QC metrics
- **WHEN** using `AdaptiveThresholdLearner`
- **THEN** system:
  - Fits Gaussian Mixture Model to metric distribution
  - Suggests optimal thresholds
  - Returns threshold values

### Requirement: Cell Cycle Scoring
The system SHALL score cells for cell cycle phases.

#### Scenario: Score cell cycle
- **GIVEN** an AnnData with normalized counts
- **WHEN** calling `score_cell_cycle(adata, species='human')`
- **THEN** system:
  - Loads species-specific cell cycle genes
  - Calculates S-phase and G2M-phase scores
  - Adds scores to `adata.obs`
  - Assigns `phase` (G1, S, G2M) based on scores

### Requirement: QC Workflow
The system SHALL provide end-to-end QC workflows.

#### Scenario: Run standard QC workflow
- **GIVEN** an AnnData with raw counts
- **WHEN** calling `run_standard_qc(adata, config)`
- **THEN** system executes:
  1. Calculate QC metrics
  2. Score cell cycle
  3. Detect doublets
  4. Suggest filtering thresholds
  5. Filter cells
  6. Generate QC report

#### Scenario: Configure workflow with custom thresholds
- **GIVEN** a `QCWorkflowConfig` with custom thresholds
- **WHEN** running workflow
- **THEN** system uses custom thresholds instead of defaults

### Requirement: Configuration Serialization
All QC configs MUST support serialization to dict/JSON.

#### Scenario: Serialize thresholds
- **GIVEN** a `QCThresholds` instance
- **WHEN** calling `thresholds.to_dict()`
- **THEN** returns Python dictionary with all fields

#### Scenario: Serialize workflow config
- **GIVEN** a `QCWorkflowConfig` instance
- **WHEN** storing to `adata.uns['qc']['config_used']`
- **THEN** config is serialized via `config.to_dict()`

#### Scenario: Save config to file
- **GIVEN** a `QCThresholds` instance
- **WHEN** calling `thresholds.to_json_file('qc_config.json')`
- **THEN** config is saved as JSON file

#### Scenario: Load config from file
- **GIVEN** a saved config JSON file
- **WHEN** calling `QCThresholds.from_json_file('qc_config.json')`
- **THEN** config is loaded and validated

### Requirement: Cache Support
The system SHALL support caching of QC results.

#### Scenario: Enable caching
- **GIVEN** a `CacheConfig` with `enabled=True, cache_dir="~/.sclucid/qc_cache"`
- **WHEN** running QC with caching
- **THEN** system:
  - Generates MD5-based cache keys
  - Stores results in cache directory
  - Reuses cached results on subsequent runs

#### Scenario: Cache expiration
- **GIVEN** cached results older than `max_cache_age_days`
- **WHEN** running QC with same parameters
- **THEN** system recomputes and updates cache
