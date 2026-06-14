# scLucid Naming Conventions

This document defines the unified naming conventions for the scLucid project to
ensure consistency and readability across the codebase.

## Function Naming Conventions

### Public API Functions

All public functions should use descriptive names following these prefixes:

| Prefix | Purpose | Examples |
|--------|---------|----------|
| `run_*` | Run a complete workflow or analysis pipeline | `run_preprocessing()`, `run_annotation()` |
| `calculate_*` | Calculate metrics or scores | `calculate_qc_metrics()`, `calculate_signature_matrix()` |
| `find_*` | Find or identify features | `find_markers()`, `find_hvgs()` |
| `get_*` | Retrieve data or configuration | `get_marker_manager()`, `get_summary()` |
| `plot_*` | Produce visualizations | `plot_embedding()`, `plot_volcano()` |
| `score_*` | Compute scores | `score_cell_types()`, `score_by_gene_sets()` |
| `compare_*` | Comparative analyses | `compare_groups()`, `compare_conditions()` |
| `filter_*` | Filtering operations | `filter_cells()`, `filter_markers()` |
| `predict_*` | Predictive operations | `predict_doublets()` |
| `suggest_*` | Suggest parameters | `suggest_qc_thresholds()`, `suggest_hvg_choice()` |

### Private Functions

Private functions (for internal use only) should start with a single underscore:

```python
def _validate_input(adata):
    """Internal validation helper."""
    pass


def _calculate_metric(values):
    """Internal calculation helper."""
    pass
```

## Class Naming Conventions

### Configuration Classes

All configuration classes should end with `Config`:

```python
class QCWorkflowConfig(SclucidBaseConfig):
    """QC workflow configuration."""
    pass


class ClusteringConfig(SclucidBaseConfig):
    """Clustering configuration."""
    pass
```

### Manager Classes

Manager classes should end with `Manager`:

```python
class ResourceManager:
    """Resource manager."""
    pass


class CacheManager:
    """Cache manager."""
    pass
```

### Analyzer / Predictor Classes

Analyzers, predictors, and similar classes should use descriptive names:

```python
class CellAnnotator:
    """Cell type annotator."""
    pass


class DoubletPredictor:
    """Doublet predictor."""
    pass
```

### Abstract Base Classes

Abstract base classes should use descriptive names without a special prefix:

```python
class AnalysisStep(ABC):
    """Abstract base class for analysis steps."""
    pass


class QCFilter(ABC):
    """Abstract base class for QC filters."""
    pass
```

## Variable Naming Conventions

### Constants

Constants should use `UPPER_CASE_WITH_UNDERSCORES`:

```python
DEFAULT_RESOLUTION = 0.8
MAX_MARKERS = 100
```

### Regular Variables

Regular variables should use `lower_case_with_underscores`:

```python
n_cells = 1000
cluster_labels = adata.obs["leiden"]
```

### Boolean Variables

Boolean variables should start with `is_`, `has_`, `can_`, etc.:

```python
is_normalized = True
has_batch_effect = False
can_parallelize = True
```

## Module and Package Naming Conventions

### Module Files

Module files should use lowercase letters and underscores:

```
qc/filtering.py
preprocess/normalize.py
analysis/clustering.py
```

### Subpackages

Subpackages should use lowercase names and avoid version suffixes such as `_v2`:

```
analysis/proportion/               # ✓ good
analysis/differential_expression/  # ✓ good
qc/workflow_v2.py                  # ✗ avoid
```

## Best Practices

### 1. Start Function Names with a Verb

Function names should begin with a verb that clearly expresses the action:

```python
# ✓ good
find_markers(adata)
filter_cells(adata)
plot_volcano(results)

# ✗ avoid
markers(adata)
cells_filter(adata)
volcano_plot(results)
```

### 2. Keep Names Concise but Descriptive

Function names should be short yet informative:

```python
# ✓ good
run_standard_qc(adata)

# ✗ too long
run_standard_quality_control_workflow(adata)

# ✗ too short
qc(adata)
```

### 3. Avoid Abbreviations

Avoid ambiguous abbreviations:

```python
# ✓ good
calculate_n_genes(adata)
plot_umap(adata)

# ✗ avoid abbreviations
calc_n_genes(adata)
plt_umap(adata)
```

### 4. Use Consistent Terminology

Use the same terms consistently throughout the codebase:

| Concept | Preferred Term | Avoid |
|---------|---------------|-------|
| Cell type | `cell_type` | `celltype`, `CellType` |
| Sample ID | `sample_id` | `sample`, `SampleID`, `sampleID` |
| Marker genes | `markers` | `marker_genes`, `marker` |
| Doublet | `doublet` | `doublets` (plural) |

## Example

### Code that follows the conventions

```python
from scLucid.base_interfaces import AnalysisStep
from scLucid.base_config import SclucidBaseConfig


class MyAnalysisConfig(SclucidBaseConfig):
    """Custom analysis configuration."""
    threshold: float = 0.5
    max_iter: int = 100


class MyAnalyzer(AnalysisStep):
    """Custom analyzer."""

    def __init__(self, config: MyAnalysisConfig):
        self.config = config

    def validate_input(self, adata):
        """Validate input."""
        if adata.n_obs == 0:
            raise ValueError("No cells in data")
        return True

    def run(self, adata, **kwargs):
        """Run the analysis."""
        # analysis logic
        return adata

    def get_summary(self):
        """Return result summary."""
        return {"config": self.config.model_dump()}


def run_my_analysis(adata, config=None):
    """Run custom analysis."""
    if config is None:
        config = MyAnalysisConfig()

    analyzer = MyAnalyzer(config)
    analyzer.validate_input(adata)
    return analyzer.run(adata)
```

## Checklist

Before submitting code, confirm:

- [ ] Public functions use the standard prefixes (`run_`, `calculate_`, `find_`, `get_`, `plot_`, etc.)
- [ ] Private functions start with `_`
- [ ] Class names use `PascalCase`
- [ ] Configuration classes end with `Config`
- [ ] Manager classes end with `Manager`
- [ ] Constants use `UPPER_CASE_WITH_UNDERSCORES`
- [ ] Boolean variables start with `is_`, `has_`, or `can_`
- [ ] Abbreviations are avoided unless universally clear
- [ ] Terminology is consistent (`cell_type`, `sample_id`, `markers`)
- [ ] Module files use lowercase letters and underscores
