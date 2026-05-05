"""
Tests for unified storage utilities.

Tests the standardized storage API for adata.uns['sclucid'].
"""

import numpy as np
import pytest
from anndata import AnnData

from scLucid.utils.storage import (
    STORAGE_ROOT,
    VALID_MODULES,
    clear_storage,
    get_storage,
    has_result,
    list_results,
    load_config,
    load_result,
    load_workflow_result,
    migrate_legacy_storage,
    save_result,
    save_workflow_result,
)


@pytest.fixture
def empty_adata():
    """Create an empty AnnData object for testing."""
    return AnnData(np.random.randn(10, 5))


@pytest.fixture
def adata_with_legacy_storage(empty_adata):
    """Create AnnData with legacy storage format."""
    adata = empty_adata.copy()
    # Simulate legacy top-level storage
    adata.uns["qc"] = {"workflow": "standard", "config_used": {"min_genes": 200}}
    return adata


class TestStorageConstants:
    """Test storage constants."""

    def test_storage_root_constant(self):
        """Test STORAGE_ROOT is correct."""
        assert STORAGE_ROOT == "sclucid"

    def test_valid_modules_contains_expected(self):
        """Test VALID_MODULES contains expected modules."""
        expected = {"qc", "preprocess", "analysis", "clustering", "annotation"}
        assert expected.issubset(VALID_MODULES)


class TestGetStorage:
    """Test get_storage function."""

    def test_get_storage_creates_structure(self, empty_adata):
        """Test that get_storage creates the storage structure."""
        storage = get_storage(empty_adata, "qc")

        assert STORAGE_ROOT in empty_adata.uns
        assert "qc" in empty_adata.uns[STORAGE_ROOT]
        assert isinstance(storage, dict)

    def test_get_storage_existing(self, empty_adata):
        """Test getting existing storage."""
        # Create storage first
        get_storage(empty_adata, "qc")
        empty_adata.uns[STORAGE_ROOT]["qc"]["test"] = "value"

        # Get again
        storage = get_storage(empty_adata, "qc")
        assert storage["test"] == "value"

    def test_get_storage_no_create(self, empty_adata):
        """Test get_storage with create=False."""
        storage = get_storage(empty_adata, "qc", create=False)

        assert STORAGE_ROOT not in empty_adata.uns
        assert storage == {}

    def test_get_storage_unknown_module(self, empty_adata):
        """Test get_storage with unknown module warns but works."""
        storage = get_storage(empty_adata, "unknown_module")

        assert isinstance(storage, dict)


class TestSaveAndLoadResult:
    """Test save_result and load_result functions."""

    def test_save_and_load_simple(self, empty_adata):
        """Test saving and loading a simple result."""
        data = {"n_cells": 100, "metrics": [1, 2, 3]}

        save_result(empty_adata, "qc", "metrics", data)
        loaded = load_result(empty_adata, "qc", "metrics")

        assert loaded == data

    def test_save_with_config(self, empty_adata):
        """Test saving result with config."""
        data = {"result": "value"}
        config = {"param1": 1, "param2": "test"}

        save_result(empty_adata, "qc", "analysis", data, config=config)

        # Check result
        assert load_result(empty_adata, "qc", "analysis") == data

        # Check config was saved
        saved_config = load_config(empty_adata, "qc", "analysis")
        assert saved_config["param1"] == 1
        assert saved_config["param2"] == "test"

    def test_load_result_default(self, empty_adata):
        """Test load_result with default value."""
        result = load_result(empty_adata, "qc", "nonexistent", default={"default": True})

        assert result == {"default": True}

    def test_save_no_overwrite(self, empty_adata):
        """Test save with overwrite=False."""
        save_result(empty_adata, "qc", "key", {"v1": 1})

        with pytest.raises(KeyError):
            save_result(empty_adata, "qc", "key", {"v2": 2}, overwrite=False)

    def test_save_overwrite_true(self, empty_adata):
        """Test save with overwrite=True."""
        save_result(empty_adata, "qc", "key", {"v1": 1})
        save_result(empty_adata, "qc", "key", {"v2": 2}, overwrite=True)

        result = load_result(empty_adata, "qc", "key")
        assert result == {"v2": 2}


class TestHasResult:
    """Test has_result function."""

    def test_has_result_true(self, empty_adata):
        """Test has_result returns True for existing key."""
        save_result(empty_adata, "qc", "exists", "value")

        assert has_result(empty_adata, "qc", "exists") is True

    def test_has_result_false(self, empty_adata):
        """Test has_result returns False for non-existing key."""
        assert has_result(empty_adata, "qc", "does_not_exist") is False


class TestListResults:
    """Test list_results function."""

    def test_list_all_results(self, empty_adata):
        """Test listing all results."""
        save_result(empty_adata, "qc", "metrics1", "data1")
        save_result(empty_adata, "qc", "metrics2", "data2")
        save_result(empty_adata, "preprocess", "normalized", "data3")

        results = list_results(empty_adata)

        assert "qc" in results
        assert "preprocess" in results
        assert "metrics1" in results["qc"]
        assert "metrics2" in results["qc"]
        assert "normalized" in results["preprocess"]

    def test_list_module_results(self, empty_adata):
        """Test listing results for specific module."""
        save_result(empty_adata, "qc", "metrics", "data")
        save_result(empty_adata, "preprocess", "normalized", "data")

        results = list_results(empty_adata, module="qc")

        assert "qc" in results
        assert "preprocess" not in results


class TestClearStorage:
    """Test clear_storage function."""

    def test_clear_specific_keys(self, empty_adata):
        """Test clearing specific keys."""
        save_result(empty_adata, "qc", "key1", "data1")
        save_result(empty_adata, "qc", "key2", "data2")

        result = clear_storage(empty_adata, module="qc", keys=["key1"])

        assert "key1" not in empty_adata.uns[STORAGE_ROOT]["qc"]
        assert "key2" in empty_adata.uns[STORAGE_ROOT]["qc"]
        # Key is prefixed with module name
        assert any("key1" in k for k in result["cleared"])

    def test_clear_module(self, empty_adata):
        """Test clearing entire module."""
        save_result(empty_adata, "qc", "key", "data")

        result = clear_storage(empty_adata, module="qc")

        assert "qc" not in empty_adata.uns[STORAGE_ROOT]
        assert "qc" in result["modules_cleared"]

    def test_clear_dry_run(self, empty_adata):
        """Test clear with dry_run."""
        save_result(empty_adata, "qc", "key", "data")

        result = clear_storage(empty_adata, module="qc", dry_run=True)

        # Should report but not delete
        assert "key" in empty_adata.uns[STORAGE_ROOT]["qc"]
        assert "qc" in result["modules_cleared"]


class TestMigrateLegacyStorage:
    """Test migrate_legacy_storage function."""

    def test_migrate_top_level_qc(self, adata_with_legacy_storage):
        """Test migrating top-level 'qc' key."""
        result = migrate_legacy_storage(adata_with_legacy_storage)

        # Check migration happened
        assert "qc (top-level -> sclucid.qc)" in result["migrated"]

        # Check data was moved
        assert STORAGE_ROOT in adata_with_legacy_storage.uns
        assert "qc" in adata_with_legacy_storage.uns[STORAGE_ROOT]

        # Check old key was removed
        assert "qc" not in adata_with_legacy_storage.uns

    def test_migrate_dry_run(self, adata_with_legacy_storage):
        """Test migrate with dry_run."""
        result = migrate_legacy_storage(adata_with_legacy_storage, dry_run=True)

        # Should report but not change
        assert "qc" in adata_with_legacy_storage.uns
        assert STORAGE_ROOT not in adata_with_legacy_storage.uns


class TestWorkflowResultHelpers:
    """Test workflow result helper functions."""

    def test_save_workflow_result(self, empty_adata):
        """Test save_workflow_result helper."""
        save_workflow_result(
            empty_adata,
            module="qc",
            workflow_name="standard",
            steps=["step1", "step2"],
            config={"param": "value"},
        )

        result = load_workflow_result(empty_adata, "qc", "standard")

        assert result["name"] == "standard"
        assert result["steps_executed"] == ["step1", "step2"]
        assert "completed_at" in result
        assert empty_adata.uns[STORAGE_ROOT]["qc"]["workflow_config"] == {"param": "value"}
        assert empty_adata.uns[STORAGE_ROOT]["qc"]["steps_executed"] == ["step1", "step2"]

    def test_save_canonical_review_summary_direct(self, empty_adata):
        """Canonical stable keys should be stored directly, not metadata-wrapped."""
        summary = {
            "schema_version": "1.0",
            "module": "qc",
            "workflow_name": "standard",
            "steps_executed": [],
            "data_shape": {"n_cells": 10, "n_genes": 5},
        }

        save_result(empty_adata, "qc", "review_summary", summary)

        assert empty_adata.uns[STORAGE_ROOT]["qc"]["review_summary"] == summary
        assert load_result(empty_adata, "qc", "review_summary") == summary

    def test_load_workflow_result_none(self, empty_adata):
        """Test load_workflow_result with no result."""
        result = load_workflow_result(empty_adata, "qc", "nonexistent")

        assert result is None


class TestStorageIntegration:
    """Integration tests for storage API."""

    def test_full_workflow_simulation(self, empty_adata):
        """Test simulating a full workflow with storage."""
        # QC step
        save_workflow_result(
            empty_adata, "qc", "standard", steps=["metrics", "filtering"], config={"min_genes": 200}
        )

        # Preprocess step
        save_workflow_result(
            empty_adata,
            "preprocess",
            "workflow",
            steps=["normalize", "hvg", "pca"],
            config={"n_top_genes": 2000},
        )

        # Analysis step
        save_result(empty_adata, "analysis", "clustering", {"n_clusters": 10, "method": "leiden"})

        # Verify all stored
        results = list_results(empty_adata)
        assert "qc" in results
        assert "preprocess" in results
        assert "analysis" in results

        # Verify can retrieve each
        assert load_workflow_result(empty_adata, "qc", "standard") is not None
        assert load_workflow_result(empty_adata, "preprocess", "workflow") is not None
        assert load_result(empty_adata, "analysis", "clustering") is not None
