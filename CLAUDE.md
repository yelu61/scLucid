
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

scLucid is a Python toolkit for single-cell RNA-seq analysis emphasizing **evidence-driven analysis** with traceable, configurable, and biologically interpretable steps.

**Key Architectural Principles:**
- **Configuration-driven**: All operations use **Pydantic** config objects for validation and reproducibility
- **Layer-aware**: Preserves raw counts (`counts` layer) and maintains intermediate results in layers (`normalized`, `scaled`)
- **Hierarchical marker system**: Uses `Manager` class from `utils/manager.py` for biology-aware annotation
- **Plugin architecture**: Abstract base classes in `base_interfaces.py` enable extensibility without core modifications

## Development Commands

### Installation
```bash
pip install -e ".[all]"      # Full development install
pip install -e ".[dev]"      # Core + dev dependencies only
```

### Testing
```bash
pytest                        # Run all tests
pytest -m unit               # Unit tests only
pytest -m "not slow"         # Skip slow tests
pytest tests/qc/test_qc_workflow.py                    # Specific file
pytest tests/qc/test_intelligent_qc.py::test_recommend_intelligent_qc  # Specific test
pytest --cov=src/scLucid --cov-report=html             # Coverage report
```

### Code Quality
```bash
black src/ tests/            # Format code
black --check src/ tests/    # Check formatting
pre-commit run --all-files   # Run pre-commit hooks
```

### Build
```bash
python -m build              # Build distribution
```

### OpenSpec CLI (Spec Management)
```bash
openspec list                  # List active changes
openspec list --specs          # List specifications
openspec show [item]           # Display change or spec
openspec validate [item]       # Validate changes or specs
openspec archive <change-id>   # Mark change complete (use --yes for automation)
openspec spec list --long      # List all specs with details
```

## Core Architecture

### Configuration System (Pydantic)

All configs inherit from `SclucidBaseConfig` (`base_config.py`):
- **Automatic validation**: Type checking and constraint enforcement
- **Serialization**: Built-in `to_dict()`, `to_json()`, `from_json()` methods
- **Path auto-creation**: `save_dir` is created automatically if specified

**Base classes:**
- `SclucidBaseConfig`: Common fields (`save_dir`, `verbose`, `plot`)
- `WorkflowConfigBase`: Extends base with `n_jobs`, `random_state`

**Module configs:**
| Module | Config File | Key Classes |
|--------|-------------|-------------|
| QC | `qc/config.py` | `QCWorkflowConfig`, `DoubletConfig`, `QCThresholds`, `AdaptiveThresholdConfig` |
| Preprocess | `preprocess/config.py` | `PreprocessingWorkflowConfig`, `NormalizationConfig`, `HVGConfig`, `IntegrationConfig` |
| Analysis | `analysis/config.py` | `ClusteringConfig`, `AnnotationConfig` |

### Abstract Base Classes (`base_interfaces.py`)

The plugin architecture is built on abstract base classes:

```python
from scLucid.base_interfaces import AnalysisStep, AnalysisStepFactory

class MyCustomQC(AnalysisStep):
    def validate_input(self, adata): ...
    def run(self, adata, **kwargs): ...
    def get_summary(self): ...

AnalysisStepFactory.register('my_qc', MyCustomQC)
qc = AnalysisStepFactory.create('my_qc')
```

**Available ABCs:**
- `AnalysisStep`: Base for all analysis steps (QC, preprocessing, clustering, etc.)
- `QCFilter`: QC filtering operations with `calculate_metric()`, `get_threshold()`, `apply_filter()`
- `CellAnnotator`: Cell type annotation with `predict()`, `get_confidence()`
- `ScoringMethod`: Functional scoring with `score()`, `normalize()`
- `PlottingBackend`: Visualization backends
- `ProportionMethod`: Cell proportion analysis

### Layer Convention

Follow this layer naming in AnnData objects:
- `counts`: Raw UMI counts (backup from `.X`)
- `normalized`: Log1p normalized data (target_sum=1e4)
- `scaled`: Z-score scaled data
- `.raw`: Set to normalized data before regression

### Marker Manager System

The `Manager` class (`utils/manager.py`) loads hierarchical markers from `.toml` files in `resources/`:

```python
from scLucid.utils.manager import get_marker_manager

mgr = get_marker_manager(species="human", tissue="Lung", case_sensitive=True)
markers = mgr.query("markers", "T cells")
```

## High-Level Workflows

**QC Workflow** (`qc/workflow.py`):
- `run_standard_qc(adata, config)`: Basic QC with defaults
- `run_advanced_qc(adata, config)`: Fully configurable QC

**Preprocessing Workflow** (`preprocess/workflow.py`):
- `run_preprocessing(adata, config)`: Full pipeline from counts to integrated UMAP

**Analysis Workflow** (`analysis/workflow.py`):
- Note: This is a stub. Use `run_annotation()` and `characterize_clusters()` directly from `analysis/annotation.py`

## Key Implementation Details

### Doublet Detection (`qc/doublet.py`)
Dual approach combining:
1. Algorithmic: Scrublet/SOLO/DoubletDetection
2. Heuristic: Lineage co-expression analysis using marker genes

### Intelligent QC (`qc/intelligent_qc.py`)
Data-driven QC recommendations using statistical modeling:
```python
from scLucid.qc import recommend_intelligent_qc
rec = recommend_intelligent_qc(adata, tissue_type='lung_tumor')
```

### Adaptive Thresholds (`qc/adaptive_threshold.py`)
Automatically determines data-driven QC thresholds using statistical methods (MAD, IQR, percentile). Use when dataset characteristics vary significantly.

### Submodule Organization
Some modules have sub-packages for better organization:
- `analysis/differential_expression/`: Differential expression and enrichment analysis
- `analysis/proportion/`: Cell proportion analysis (sccoda, compositional analysis)
- `preprocess/intelligent/`: Intelligent preprocessing methods
- `tools/pyCellChat/`, `tools/pyDWLS/`, `tools/pyMonocle3/`: Python-native tool wrappers (R-free)

### Test Fixtures
Available in `tests/fixtures/data_loader.py`:
- `pbmc`, `pbmc_small`: 10x Genomics PBMC datasets
- `nsclc`: Lung cancer (GEO GSE119911)
- `mouse_melanoma`, `mouse_lung`: Mouse cancer datasets

## Documentation Structure

**Dual-track system:**
- `docs/notebooks/`: Interactive learning (Jupyter notebooks)
- `examples/`: Production-ready scripts with Pydantic configs
- `docs/source/`: Sphinx-generated API reference

**OpenSpec:**
- `openspec/specs/`: Current requirements (source of truth)
- `openspec/changes/`: Active proposals
- `openspec/changes/archive/`: Completed changes
- `openspec/project.md`: Project conventions and context

**OpenSpec workflow:**
1. Review existing specs before making changes: `openspec list --specs`
2. Create proposals for new features/breaking changes in `openspec/changes/`
3. Validate proposals: `openspec validate <change-id> --strict`
4. Archive after deployment: `openspec archive <change-id>`

## Resource Files

Built-in marker databases in `resources/`:
- `marker_base_human.toml`: Core cell type markers
- `marker_base_mouse.toml`: Mouse equivalents
- `marker_tissue_specific_human.toml`: Tissue-context markers
- `marker_cell_state_human.toml`: Activation states

Set `metadata.doublet_lineage = true` for lineages used in doublet detection.

## Key Files for Tasks

| Task | Key Files |
|------|-----------|
| Add QC metric | `qc/metrics.py`, `qc/config.py` |
| Modify doublet detection | `qc/doublet.py`, `qc/config.py:DoubletConfig` |
| Adaptive thresholds | `qc/adaptive_threshold.py` |
| Change normalization | `preprocess/normalize.py`, `preprocess/config.py:NormalizationConfig` |
| Add integration method | `preprocess/integrate.py`, `preprocess/config.py:IntegrationConfig` |
| Add annotation method | `analysis/annotation.py`, `analysis/config.py:AnnotationConfig` |
| Differential expression | `analysis/differential_expression/` |
| Cell proportion analysis | `analysis/proportion/` |
| Custom plots | `plotting/main.py`, `plotting/advanced_plots.py` |
| Add markers | Edit `.toml` files in `resources/` |
| Web API endpoints | `web/api/routes/*.py` |
| Plugin development | `base_interfaces.py`, `docs/PLUGIN_DEVELOPMENT_GUIDE.md` |
| Naming conventions | `docs/NAMING_CONVENTIONS.md` |
| Intelligent QC | `qc/intelligent_qc.py`, `qc/strategy_decision_tree.py` |

## Naming Conventions

**Functions:** Use descriptive prefixes
- `run_*`: Workflows (`run_preprocessing()`, `run_annotation()`)
- `calculate_*`: Metrics (`calculate_qc_metrics()`)
- `find_*`: Identification (`find_markers()`, `find_hvgs()`)
- `get_*`: Data retrieval (`get_marker_manager()`)
- `plot_*`: Visualization (`plot_embedding()`)
- `score_*`: Scoring (`score_cell_types()`)

**Classes:**
- Configs: `*Config` (e.g., `QCWorkflowConfig`)
- Managers: `*Manager` (e.g., `MarkerManager`)
- Private members: Leading underscore `_`

**Variables:**
- Constants: `UPPER_CASE`
- Booleans: `is_*`, `has_*`, `can_*` prefixes

## Module Organization

**Analysis submodules:**
- `analysis/clustering.py`: Clustering algorithms (Leiden, HDBSCAN)
- `analysis/annotation.py`: Cell type annotation (CellTypist, scoring, enrichment)
- `analysis/differential_expression/`: DE analysis and gene set enrichment
- `analysis/proportion/`: Cell proportion analysis (scCODA)

**QC submodule structure:**
- `qc/metrics.py`: QC metric calculation
- `qc/filtering.py`: Cell/gene filtering
- `qc/doublet.py`: Doublet detection (algorithmic + heuristic)
- `qc/intelligent_qc.py`: Data-driven QC recommendations
- `qc/adaptive_threshold.py`: Adaptive threshold determination
- `qc/strategy_decision_tree.py`: QC strategy selection logic
- `qc/workflow.py`: QC workflow orchestration

**Preprocessing submodules:**
- `preprocess/normalize.py`: Normalization methods
- `preprocess/hvg.py`: Highly variable gene selection
- `preprocess/scale.py`: Data scaling
- `preprocess/integrate.py`: Batch correction/integration
- `preprocess/neighbors.py`: Neighbor graph construction
- `preprocess/intelligent/`: Intelligent preprocessing methods

**Tools (Python-native, R-free):**
- `tools/pyCellChat/`: Cell-cell communication analysis
- `tools/pyDWLS/`: Bulk deconvolution
- `tools/pyMonocle3/`: Trajectory inference
- `tools/pySCENIC.py`: Gene regulatory network inference
- `tools/sccoda.py`: Compositional analysis
- `tools/spatial.py`: Spatial transcriptomics analysis

**Plotting modules:**
- `plotting/main.py`: High-level plotting interface
- `plotting/embedding_plots.py`: UMAP/tSNE visualizations
- `plotting/feature_plots.py`: Gene/metric feature plots
- `plotting/marker_plots.py`: Marker gene visualization
- `plotting/advanced_plots.py`: Complex multi-panel figures

## Global Settings

```python
from scLucid import set_figure_params, setup_logging

setup_logging(level="INFO", file_path="analysis.log")
set_figure_params(dpi=100, color_theme="dark")  # FONT_NATURE, FONT_CELL, FONT_TRADITIONAL
```

## Common Patterns

### Working with Configs
All configs support automatic validation, serialization, and can be saved/loaded:

```python
from scLucid.qc.config import QCWorkflowConfig

# Create config
config = QCWorkflowConfig(
    species="human",
    save_dir="./results",
    doublet_detection=True
)

# Serialize
config_dict = config.to_dict()
config_json = config.to_json()
config.to_json_file("qc_config.json")

# Load from file
config2 = QCWorkflowConfig.from_json_file("qc_config.json")
```

### Running Workflows with Checkpoints
Workflows support checkpointing for resumable analysis:

```python
from scLucid.qc import run_advanced_qc
from scLucid.qc.config import QCWorkflowConfig

config = QCWorkflowConfig(
    save_dir="./results",
    checkpoint_file="qc_checkpoint.pkl"  # Enable checkpointing
)
adata = run_advanced_qc(adata, config=config)
# Can resume from checkpoint if interrupted
```

### Module Aliases for Interactive Use
When working interactively, use the module aliases:

```python
import scLucid as scl

# Use aliases for brevity
adata = scl.pp.normalize(adata)
adata = scl.pp.find_hvgs(adata)
adata = scl.al.cluster_cells(adata)
```

### Accessing Marker Genes
```python
from scLucid.utils.manager import get_marker_manager

mgr = get_marker_manager(species="human", tissue="Lung")
# Query markers by cell type
markers = mgr.query("markers", "T cells")
# Get state-specific markers
state_markers = mgr.query("cell_state", "activated")
```
