# scLucid Plugin Development Guide

> This guide describes extension points for custom scLucid analysis steps. It is
> a developer-facing document, not the primary project roadmap. For current
> product positioning and implementation priorities, see
> `docs/SCLUCID_STRATEGIC_IMPLEMENTATION_PLAN.md`.

## What Is Plugin Development?

### Traditional Approach (problematic)

```python
# Adding a custom QC method requires editing core code.
def my_custom_qc(adata, threshold):
    # ... your logic ...
    pass

# You would need to modify scLucid/qc/filtering.py.
# A bug could break the whole package, and code review becomes harder.
```

### Plugin Approach (preferred)

```python
# Create a standalone plugin file without modifying core code.
from scLucid import AnalysisStep


class MyCustomQC(AnalysisStep):
    def validate_input(self, adata):
        return True

    def run(self, adata, **kwargs):
        # ... your logic ...
        return adata


# Register with the factory.
from scLucid import AnalysisStepFactory

AnalysisStepFactory.register("my_qc", MyCustomQC)

# Use it.
qc = AnalysisStepFactory.create("my_qc")
qc.run(adata)
```

---

## Core Concept: Abstract Base Classes (ABC)

An **abstract base class** is a contract that defines:

- Which methods must be implemented
- The expected method signatures
- How the implementation interacts with the rest of the system

```
Abstract base class = blueprint
            │
    ┌───────┼───────┐
    │       │       │
  Impl 1  Impl 2  Impl 3
  (MyQC) (YourQC) (CustomQC)
    │       │       │
    └───────┴───────┘
            │
    All implementations follow
    the same blueprint
```

---

## scLucid Abstract Base Class System

```
scLucid architecture
│
├─ AnalysisStep          → base class for all analysis steps
│   ├─ QCFilter          → QC filters
│   ├─ CellAnnotator     → cell annotators
│   ├─ ScoringMethod     → scoring methods
│   ├─ PlottingBackend   → plotting backends
│   └─ ProportionMethod  → proportion-analysis methods
│
└─ AnalysisStepFactory   → factory pattern
    ├─ register()        → register a plugin
    ├─ create()          → create an instance
    └─ list_steps()      → list all registered plugins
```

---

## Practical Scenarios

### Scenario 1: Create a Custom QC Filter

```python
from scLucid import AnalysisStep, SclucidBaseConfig
from anndata import AnnData
from pydantic import Field


# Step 1: define a config class.
class StrictQCConfig(SclucidBaseConfig):
    min_genes: int = Field(default=500, ge=0)
    max_mt: float = Field(default=10.0, ge=0, le=100)


# Step 2: implement the plugin class.
class StrictQCFilter(AnalysisStep):
    def __init__(self, config: StrictQCConfig):
        self.config = config

    def validate_input(self, adata: AnnData) -> bool:
        """Check that input is valid."""
        if adata.n_obs == 0:
            raise ValueError("No cells in data")
        return True

    def run(self, adata: AnnData, **kwargs) -> AnnData:
        """Run the QC filter."""
        import scanpy as sc

        sc.pp.filter_cells(adata, min_genes=self.config.min_genes)
        return adata

    def get_summary(self) -> dict:
        """Return a summary."""
        return {"min_genes": self.config.min_genes}


# Step 3: register the plugin.
from scLucid import AnalysisStepFactory

AnalysisStepFactory.register("strict_qc", StrictQCFilter)

# Step 4: use the plugin.
qc = AnalysisStepFactory.create("strict_qc")
qc.run(adata)
```

### Scenario 2: Compose Multiple Plugins

```python
steps = [
    "strict_qc",   # your custom QC
    "clustering",  # standard clustering
    "my_annotator",  # your custom annotator
]

configs = {
    "strict_qc": {"config": StrictQCConfig(min_genes=500)},
    "clustering": {"resolution": 0.8},
    "my_annotator": {...},
}

from scLucid.analysis import run_custom_analysis

adata = run_custom_analysis(adata, steps=steps, step_configs=configs)
```

---

## Factory Pattern

### Why Use a Factory?

Direct creation:
```python
qc1 = StrictQCFilter(config1)
qc2 = AnotherQC(config2)
# ❌ caller must know every concrete class name
```

Factory creation:
```python
qc1 = AnalysisStepFactory.create("strict_qc", config=config1)
qc2 = AnalysisStepFactory.create("another_qc", config=config2)
# ✅ only the plugin name is needed; more flexible
```

### Factory Benefits

1. **Decoupling** — callers do not need concrete class names.
2. **Dynamic** — selection can happen at runtime.
3. **Extensible** — adding a plugin does not require changing consumer code.
4. **Testable** — easy to mock and test.

---

## Benefits of Plugin Development

### For Developers

| Benefit | Description |
|---------|-------------|
| **Independent development** | Build in your own files without touching core code. |
| **Fast iteration** | No need to wait for core code review and merge. |
| **Version control** | Plugins can have their own versioning and release cycle. |
| **Opt-in usage** | Other users can decide whether to use your plugin. |

### For Users

| Benefit | Description |
|---------|-------------|
| **More choices** | Community-developed algorithms become available. |
| **Customization** | Lab-specific analysis pipelines are possible. |
| **Experimentation** | Experimental features can be tried without destabilizing core workflows. |

### For Maintainers

| Benefit | Description |
|---------|-------------|
| **Core stability** | Fewer changes to core code, lower regression risk. |
| **Lower burden** | Not every algorithm variant has to be maintained centrally. |
| **Community contributions** | Plugins can be shared and contributed externally. |

---

## Concrete Examples

### Example 1: Database Annotator

```python
class DatabaseAnnotator(AnalysisStep):
    """Annotate cell types from an external database."""

    def validate_input(self, adata):
        return "X_pca" in adata.obsm

    def run(self, adata, database_path, **kwargs):
        # connect to database
        # query the most similar cell type
        # store annotation in adata
        return adata


AnalysisStepFactory.register("db_annotator", DatabaseAnnotator)
```

### Example 2: Machine-Learning Scorer

```python
class MLScorer(AnalysisStep):
    """Score cells using a trained ML model."""

    def validate_input(self, adata):
        return adata.shape[1] > 0

    def run(self, adata, model_path, **kwargs):
        # load trained model
        # score each cell
        # store scores in adata.obs
        return adata


AnalysisStepFactory.register("ml_scorer", MLScorer)
```

### Example 3: Interactive Plotting Backend

```python
class InteractivePlotter(AnalysisStep):
    """Create interactive visualizations."""

    def validate_input(self, adata):
        return adata.n_obs > 0

    def run(self, adata, **kwargs):
        # create plotly figure
        # add interactive controls
        # return HTML
        return adata

    def get_summary(self):
        return {"type": "interactive"}


AnalysisStepFactory.register("interactive_plot", InteractivePlotter)
```

---

## Learning Resources

### Inside the Project
- `src/scLucid/base_interfaces.py` — abstract base class definitions
- `examples/plugin_development_example.py` — complete example
- `docs/dev/NAMING_CONVENTIONS.md` — naming conventions

### External Resources
- Python ABC: https://docs.python.org/3/library/abc.html
- Factory pattern: https://refactoring.guru/design-patterns/factory-method
- Plugin architecture: https://en.wikipedia.org/wiki/Plugin_(computing)

---

## FAQ

### Q1: Do I have to use plugins?

**A**: No. Plugins are optional. You can keep using existing functions and
classes. Plugins simply provide an extension mechanism.

### Q2: Do plugins affect performance?

**A**: Plugins are ordinary Python classes, so performance is comparable to
core functions. Factory creation overhead is negligible (<1 ms).

### Q3: How do I share my plugin?

**A**: You can:
- Publish it to PyPI as a standalone package.
- Share code snippets with users.
- Submit a PR to scLucid if the plugin is broadly useful.

### Q4: Can plugins break core code?

**A**: No. Plugins run independently and do not modify core code.

### Q5: Can I mix core functions and plugins?

**A**: Yes. You can combine core functions and custom plugins in the same
analysis pipeline.

---

## Quick Start

1. **Learn the ABCs**: read `src/scLucid/base_interfaces.py`.
2. **Run the example**: `python examples/plugin_development_example.py`.
3. **Create your first plugin**: inherit from `AnalysisStep` and implement the
   required methods.
4. **Register and use**: use `AnalysisStepFactory`.
5. **Share**: consider sharing useful plugins with the community.
