# Project Context

## Purpose
scLucid is a Python toolkit for end-to-end single-cell RNA-seq analysis. The project focuses on robust, reproducible workflows for quality control, preprocessing, clustering/annotation, differential expression, enrichment, proportion analysis, and specialized tumor and tool integrations.

## Tech Stack
- Python 3.10+ package using `src/` layout
- Core data model: `anndata.AnnData`
- Scientific stack: NumPy, Pandas, SciPy, scikit-learn, statsmodels
- Single-cell ecosystem: Scanpy (plus optional ecosystem backends)
- Configuration and validation: Pydantic v2
- Visualization: Matplotlib (primary), Plotly (optional)
- Testing: Pytest

## Project Conventions

### Code Style
- Follow PEP 8 with explicit, descriptive names.
- Add type hints for public APIs and key internal interfaces.
- Keep module-level side effects minimal; avoid expensive imports in `__init__.py`.
- Prefer explicit imports over wildcard imports.
- Use lowercase module names and keep public API names stable once released.

### Architecture Patterns
- Package modules are organized by domain: `qc`, `preprocess`, `analysis`, `plotting`, `tools`, `utils`, `tumor`, `web`.
- Config classes should inherit from shared base config types when available and use Pydantic validation.
- Public outputs should be stored in structured AnnData namespaces under `adata.uns['sclucid']`.
- Optional features MUST degrade gracefully when optional dependencies are unavailable.
- High-level `__init__.py` exports should be stable and defensive against missing optional backends.

### Testing Strategy
- Add unit tests for config validation and pure utility logic.
- Add smoke tests for package import surfaces (`scLucid`, `scLucid.qc`, `scLucid.analysis`, `scLucid.tools`, `scLucid.config`).
- Add lightweight integration tests with small synthetic AnnData objects for core workflows.
- New/changed public APIs require tests covering success and at least one failure/validation path.

### Git Workflow
- Use short-lived branches and focused commits.
- Keep changes scoped to one capability where possible.
- Prefer additive, backward-compatible changes unless a breaking change is intentional and documented.
- For behavior or API changes, update OpenSpec and corresponding docs/examples in the same change.

## Domain Context
- Input data are typically sparse count matrices from scRNA-seq with metadata in `obs` and feature metadata in `var`.
- Typical workflow order: QC -> preprocessing -> clustering/annotation -> downstream analysis.
- Reproducibility is important; random seeds and deterministic settings should be explicitly configurable.

## Important Constraints
- Datasets can be large; memory use and chunked processing options matter.
- Some dependencies are optional/heavy and may be absent in user environments.
- Imports should remain resilient so users can access available functionality without full optional stack installation.
- Backward compatibility for public API symbols is preferred across minor releases.

## External Dependencies
- Scanpy/AnnData ecosystem for core single-cell operations.
- Optional tool backends (for deconvolution, communication, trajectory, web) may vary by environment.
- Visualization backends may require system-level libraries (fonts, rendering, GUI/display stack).
