# Contributing to scLucid

Thank you for your interest in contributing to scLucid! This document provides guidelines and instructions for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Code Style Guidelines](#code-style-guidelines)
- [Testing Guidelines](#testing-guidelines)
- [Documentation Guidelines](#documentation-guidelines)
- [Submitting Changes](#submitting-changes)
- [Project Structure](#project-structure)

## Code of Conduct

- Be respectful and inclusive
- Provide constructive feedback
- Focus on what is best for the community
- Show empathy towards other community members

## Getting Started

### Prerequisites

- Python 3.9 or higher
- Git
- pip or conda/micromamba

### Development Setup

```bash
# 1. Fork and clone the repository
git clone https://github.com/yourusername/scLucid.git
cd scLucid

# 2. Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# 3. Install development dependencies
pip install -e ".[dev]"

# 4. Install pre-commit hooks (optional but recommended)
pre-commit install

# 5. Verify installation
pytest tests/ -v
```

### Development Dependencies

The `[dev]` extra includes:
- **Testing**: pytest, pytest-cov
- **Code quality**: black, flake8, mypy
- **Documentation**: sphinx, sphinx-rtd-theme
- **Pre-commit**: pre-commit hooks

## Development Workflow

### 1. Branch Naming

Use descriptive branch names:

```bash
# Feature branches
git checkout -b feature/add-new-clustering-method

# Bugfix branches
git checkout -b fix/qc-threshold-validation

# Documentation branches
git checkout -b docs/update-api-reference
```

### 2. Making Changes

```bash
# Make your changes
# ...

# Stage files
git add path/to/changed/files

# Commit with descriptive message
git commit -m "Add spectral clustering method

- Implement spectral_clustering() function
- Add SpectralClusteringConfig to config.py
- Update tests and documentation
- Closes #123"
```

### 3. Syncing with Upstream

```bash
# Add upstream (if not already added)
git remote add upstream https://github.com/yelu61/scLucid.git

# Fetch latest changes
git fetch upstream

# Rebase your branch on master
git rebase upstream/master
```

## Code Style Guidelines

### Python Style

scLucid follows **PEP 8** with these modifications:

- **Line length**: 100 characters (black default)
- **Indentation**: 4 spaces
- **Imports**: Grouped and sorted (isort style)
- **Docstrings**: Google style

### Formatting Code

```bash
# Format all code with black
black src/ tests/

# Check formatting without making changes
black --check src/ tests/

# Sort imports with isort
isort src/ tests/
```

### Linting

```bash
# Run flake8 linter
flake8 src/ tests/

# Run mypy type checker
mypy src/
```

### Naming Conventions

```python
# Modules: lowercase with underscores
# my_module.py

# Classes: CapWords (PascalCase)
class MyConfiguration:
    pass

# Functions and variables: lowercase_with_underscores
def my_function():
    my_variable = 1

# Constants: UPPER_CASE_WITH_UNDERSCORES
MY_CONSTANT = 42

# Private members: leading underscore
def _internal_helper():
    pass
```

### Configuration Classes

All configuration classes must use **Pydantic**:

```python
from pydantic import Field, field_validator
from ..base_config import SclucidBaseConfig

class MyConfig(SclucidBaseConfig):
    """Brief description of the configuration.

    Longer description with details about usage.
    """

    # Use Field() for all parameters
    threshold: float = Field(
        default=0.5,
        ge=0.0,
        le=1.0,
        description="Threshold parameter between 0 and 1"
    )

    @field_validator("threshold")
    @classmethod
    def validate_threshold(cls, v: float) -> float:
        """Custom validation logic."""
        if v < 0.1:
            logger.warning("Threshold is very small")
        return v
```

### Docstring Style

Use **Google style** docstrings:

```python
def calculate_metric(adata: anndata.AnnData, layer: str = "X") -> pd.DataFrame:
    """Calculate a QC metric for AnnData object.

    This function computes the specified metric for each cell in the
    AnnData object, using data from the specified layer.

    Args:
        adata: AnnData object containing the data.
        layer: Layer to use for calculation. Defaults to "X".

    Returns:
        DataFrame with metric values for each cell.

    Raises:
        ValueError: If specified layer does not exist.

    Examples:
        >>> metrics = calculate_metric(adata, layer="counts")
        >>> print(metrics.head())
    """
    pass
```

## Testing Guidelines

### Writing Tests

All new features must include tests:

```python
import pytest
import anndata
from scLucid.qc import calculate_qc_metric

def test_calculate_qc_metric_basic():
    """Test basic QC metric calculation."""
    # Create test data
    adata = anndata.AnnData(X=np.random.rand(100, 50))

    # Run function
    result = calculate_qc_metric(adata, species="human")

    # Assert results
    assert "n_genes_by_counts" in result.obs.columns
    assert "pct_counts_mt" in result.obs.columns
    assert result.obs["n_genes_by_counts"].min() >= 0

def test_calculate_qc_metric_invalid_species():
    """Test that invalid species raises error."""
    adata = anndata.AnnData(X=np.random.rand(100, 50))

    with pytest.raises(ValueError, match="Invalid species"):
        calculate_qc_metric(adata, species="invalid")
```

### Test Structure

```
tests/
├── fixtures/           # Test fixtures and data loaders
│   └── data_loader.py
├── qc/                 # QC module tests
│   ├── test_metrics.py
│   ├── test_doublet.py
│   └── test_workflow.py
├── preprocess/         # Preprocess module tests
│   └── test_normalize.py
├── analysis/           # Analysis module tests
│   └── test_clustering.py
└── conftest.py        # Shared pytest configuration
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/qc/test_metrics.py

# Run specific test function
pytest tests/qc/test_metrics.py::test_calculate_qc_metric_basic

# Run with coverage
pytest --cov=src/scLucid --cov-report=html

# Run only fast tests
pytest -m "not slow"

# Run only unit tests
pytest -m unit
```

### Test Markers

Use pytest markers to categorize tests:

```python
@pytest.mark.unit
def test_function_internal():
    pass

@pytest.mark.integration
def test_full_workflow():
    pass

@pytest.mark.slow
def test_large_dataset():
    pass
```

## Documentation Guidelines

### Docstring Requirements

- All public functions/classes must have docstrings
- Use Google style (see above)
- Include `Args`, `Returns`, `Raises`, `Examples` sections
- Keep descriptions concise but informative

### API Documentation

API docs are autogenerated from docstrings using Sphinx:

```bash
# Build documentation
cd docs
make html

# View documentation
open _build/html/index.html
```

### Updating Documentation

When adding new features:

1. **Update docstrings** in source code
2. **Add examples** to relevant notebooks in `docs/notebooks/`
3. **Update best practices** in `docs/source/best_practices.rst`
4. **Create examples** in `examples/` if appropriate
5. **Update OpenSpec specs** in `openspec/specs/` if breaking changes

### Documentation Files

- `README.md`: User-facing overview
- `CLAUDE.md`: AI assistant and developer guide
- `docs/source/`: Sphinx documentation (RST format)
- `docs/notebooks/`: Jupyter notebooks (tutorials)
- `examples/`: Runnable example scripts
- `openspec/`: Specifications and change proposals

## Submitting Changes

### Pull Request Process

1. **Update your branch** with latest master
2. **Run tests** and ensure they pass
3. **Run linters** and fix any issues
4. **Update documentation** if needed
5. **Commit changes** with clear messages
6. **Push to GitHub** and open PR

### Pull Request Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
- [ ] Tests added/updated
- [ ] All tests pass

## Documentation
- [ ] Docstrings updated
- [ ] README/docs updated
- [ ] Examples added

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Comments added to complex code
- [ ] Changes generate no new warnings
```

### Review Process

- **Automated checks**: CI runs tests, linting, type checking
- **Code review**: Maintainers review within 1-3 days
- **Feedback**: Address review comments promptly
- **Approval**: At least one maintainer approval required

## Project Structure

### Key Directories

```
scLucid/
├── src/scLucid/          # Main source code
│   ├── qc/              # Quality control module
│   ├── preprocess/      # Preprocessing module
│   ├── analysis/        # Analysis module
│   ├── tools/           # Specialized tools
│   ├── utils/           # Utilities
│   └── config.py        # Global configuration
├── tests/               # Test suite
├── docs/                # Documentation
├── examples/            # Example scripts
├── openspec/            # Specifications
└── pyproject.toml       # Project configuration
```

### Adding New Modules

When adding a new module (e.g., `new_module/`):

1. Create directory: `src/scLucid/new_module/`
2. Add `__init__.py` with exports
3. Add `config.py` with Pydantic configs
4. Implement core functionality
5. Add tests in `tests/new_module/`
6. Update main `__init__.py` to export module
7. Add API documentation in `docs/source/api/new_module.rst`
8. Update OpenSpec spec in `openspec/specs/new_module/`

### Configuration Files

- `pyproject.toml`: Project metadata, dependencies, build config
- `pytest.ini`: Pytest configuration
- `.pre-commit-config.yaml`: Pre-commit hooks
- `setup.cfg`: Package setup options

## Getting Help

### Resources

- **Documentation**: https://sclucid.readthedocs.io/
- **Issues**: https://github.com/yelu61/scLucid/issues
- **Discussions**: https://github.com/yelu61/scLucid/discussions

### Asking Questions

1. Check existing issues and discussions
2. Search documentation
3. Create new issue with:
   - Clear description
   - Minimal reproducible example
   - Error messages and traceback
   - Environment details (OS, Python version)

## Recognition

Contributors will be acknowledged in:
- `AUTHORS.md` file
- Release notes
- Documentation contributor list

Thank you for contributing to scLucid! 🎉

---

**Last Updated**: 2026-02-08
