"""
Pytest configuration and shared fixtures for scLucid tests.

This file is automatically discovered by pytest and provides:
- Synthetic data fixtures
- Test markers
- Shared utilities
"""

import os
import tempfile
import pytest
from anndata import AnnData
import numpy as np

# Import all fixtures from synthetic_data module
from tests.fixtures.synthetic_data import (
    synthetic_generator,
    minimal_adata,
    qc_test_adata,
    integration_test_adata,
    doublet_test_adata,
)

# Ensure writable runtime caches for matplotlib/numba in sandboxed CI environments.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "mplcfg"))
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "numba-cache"))


def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line("markers", "unit: Unit tests (fast, isolated)")
    config.addinivalue_line("markers", "integration: Integration tests (may use real data)")
    config.addinivalue_line("markers", "slow: Slow tests (skip in quick runs)")
    config.addinivalue_line("markers", "config: Configuration system tests")
    config.addinivalue_line("markers", "optional: Optional dependency/backend tests")
    config.addinivalue_line("markers", "smoke: Import/public-surface smoke tests")


@pytest.fixture(scope="session")
def test_data_dir():
    """Provide path to test data directory."""
    from pathlib import Path
    return Path(__file__).parent / "fixtures" / "data"


@pytest.fixture(scope="function")
def temp_output_dir(tmp_path):
    """Provide a temporary output directory for each test."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return str(output_dir)


@pytest.fixture(scope="function")
def tiny_adata():
    """A very small AnnData fixture for smoke/integration-style tests."""
    x = np.array(
        [
            [1, 0, 3, 4, 0, 2],
            [0, 2, 1, 3, 4, 0],
            [3, 1, 0, 2, 1, 5],
            [5, 2, 1, 0, 3, 1],
            [2, 0, 4, 1, 0, 2],
            [1, 3, 2, 4, 1, 0],
            [0, 1, 5, 2, 2, 1],
            [4, 0, 1, 3, 0, 3],
        ],
        dtype=float,
    )
    adata = AnnData(x.copy())
    adata.var_names = [f"gene_{i}" for i in range(adata.n_vars)]
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.layers["counts"] = x.copy()
    return adata
