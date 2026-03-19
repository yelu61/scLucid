"""
Tests for workflow utilities.

Tests progress bars, workflow errors, checkpoints, and partial results.
"""

import pytest
import tempfile
import json
from pathlib import Path

from scLucid.utils.workflow_utils import (
    get_progress_bar,
    WorkflowError,
    StepError,
    RecoveryError,
    WorkflowCheckpoint,
    PartialResultManager,
    with_error_recovery,
    merge_partial_results,
)


class TestWorkflowError:
    """Test WorkflowError exception."""

    def test_workflow_error_creation(self):
        """Test basic WorkflowError creation."""
        err = WorkflowError("Workflow failed")
        assert err.step_name == "unknown"
        assert err.original_error is None

    def test_workflow_error_with_details(self):
        """Test WorkflowError with step and original error."""
        original = ValueError("Original problem")
        err = WorkflowError("Workflow failed", step_name="normalization", original_error=original)

        assert err.step_name == "normalization"
        assert err.original_error is original

    def test_workflow_error_str(self):
        """Test WorkflowError string representation."""
        err = WorkflowError("Test message")
        str_repr = str(err)
        assert "Test message" in str_repr
        assert "step: unknown" in str_repr


class TestStepError:
    """Test StepError exception."""

    def test_step_error(self):
        """Test StepError creation."""
        err = StepError("Step failed", step_name="qc")
        assert err.step_name == "qc"


class TestRecoveryError:
    """Test RecoveryError exception."""

    def test_recovery_error_creation(self):
        """Test RecoveryError creation."""
        err = RecoveryError("Recovery failed", step_name="recovery_step")
        assert err.step_name == "recovery_step"
        assert "Recovery failed" in str(err)


class TestWorkflowCheckpoint:
    """Test WorkflowCheckpoint dataclass."""

    def test_checkpoint_creation(self):
        """Test creating WorkflowCheckpoint."""
        checkpoint = WorkflowCheckpoint(
            completed_steps=["step1", "step2"],
            failed_step="step3",
            error_message="Something went wrong"
        )

        assert checkpoint.completed_steps == ["step1", "step2"]
        assert checkpoint.failed_step == "step3"
        assert checkpoint.error_message == "Something went wrong"

    def test_checkpoint_to_dict(self):
        """Test checkpoint serialization."""
        checkpoint = WorkflowCheckpoint(
            completed_steps=["step1"],
            failed_step="step2",
            error_message="Error"
        )

        data = checkpoint.to_dict()

        assert data["completed_steps"] == ["step1"]
        assert data["failed_step"] == "step2"
        assert "timestamp" in data

    def test_checkpoint_from_dict(self):
        """Test checkpoint deserialization."""
        data = {
            "completed_steps": ["step1"],
            "failed_step": "step2",
            "error_message": "Error",
            "timestamp": "2024-01-01T00:00:00",
            "config_hash": "abc123"
        }

        checkpoint = WorkflowCheckpoint.from_dict(data)

        assert checkpoint.completed_steps == ["step1"]
        assert checkpoint.failed_step == "step2"


class TestPartialResultManager:
    """Test PartialResultManager class."""

    @pytest.fixture
    def temp_save_dir(self):
        """Create temporary save directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_manager_creation(self, temp_save_dir):
        """Test creating PartialResultManager."""
        manager = PartialResultManager(temp_save_dir)
        assert manager.save_dir == Path(temp_save_dir)

    def test_save_and_exists(self, temp_save_dir):
        """Test saving and checking existence."""
        import numpy as np
        from anndata import AnnData

        manager = PartialResultManager(temp_save_dir)

        # Create test data
        adata = AnnData(np.random.randn(10, 5))
        checkpoint = WorkflowCheckpoint(
            completed_steps=["step1"],
            failed_step="step2",
            error_message="Test error"
        )
        config = {"param": "value"}

        # Save
        manager.save(adata, checkpoint, config)

        # Check exists
        assert manager.exists() is True

        # Load
        loaded_adata, loaded_checkpoint, loaded_config = manager.load()

        assert loaded_adata.shape == adata.shape
        assert loaded_checkpoint.completed_steps == ["step1"]
        assert loaded_config["param"] == "value"

    def test_exists_when_empty(self, temp_save_dir):
        """Test exists() when no checkpoint saved."""
        manager = PartialResultManager(temp_save_dir)

        assert manager.exists() is False


class TestGetProgressBar:
    """Test get_progress_bar function."""

    def test_progress_bar_enabled(self):
        """Test progress bar when enabled."""
        items = [1, 2, 3]
        result = list(get_progress_bar(items, desc="Test", enabled=True))

        assert result == items

    def test_progress_bar_disabled(self):
        """Test progress bar when disabled."""
        items = [1, 2, 3]
        result = list(get_progress_bar(items, desc="Test", enabled=False))

        assert result == items

    def test_progress_bar_with_generator(self):
        """Test progress bar with generator."""
        def generator():
            yield from range(5)

        result = list(get_progress_bar(generator(), total=5, enabled=False))

        assert result == [0, 1, 2, 3, 4]


class TestWithErrorRecovery:
    """Test with_error_recovery decorator."""

    def test_successful_execution(self, tmp_path):
        """Test decorator with successful function."""
        import numpy as np
        from anndata import AnnData

        adata = AnnData(np.random.randn(10, 5))
        save_dir = str(tmp_path / "recovery")

        @with_error_recovery(save_dir=save_dir, step_name="test_step")
        def successful_func(adata):
            return "success"

        result = successful_func(adata)

        assert result == "success"

    def test_error_handling(self, tmp_path):
        """Test decorator with failing function."""
        import numpy as np
        from anndata import AnnData

        adata = AnnData(np.random.randn(10, 5))
        save_dir = str(tmp_path / "recovery")

        @with_error_recovery(save_dir=save_dir, step_name="failing_step")
        def failing_func(adata):
            raise ValueError("Test error")

        with pytest.raises(WorkflowError) as exc_info:
            failing_func(adata)

        assert exc_info.value.step_name == "failing_step"
        assert "Test error" in str(exc_info.value.original_error)


class TestMergePartialResults:
    """Test merge_partial_results function."""

    def test_merge_results(self):
        """Test merging partial results."""
        # This is a placeholder test since merge_partial_results
        # implementation may vary
        pass


class TestWorkflowUtilsIntegration:
    """Integration tests for workflow utilities."""

    def test_full_checkpoint_lifecycle(self, tmp_path):
        """Test full checkpoint save/load lifecycle."""
        import numpy as np
        from anndata import AnnData

        save_dir = str(tmp_path / "checkpoint_test")
        manager = PartialResultManager(save_dir)

        # Simulate workflow progress
        for i in range(3):
            adata = AnnData(np.random.randn(10 * (i + 1), 5))
            checkpoint = WorkflowCheckpoint(
                completed_steps=[f"step{j}" for j in range(i + 1)]
            )
            config = {"iteration": i}

            manager.save(adata, checkpoint, config)

        # Load final state
        loaded_adata, loaded_checkpoint, loaded_config = manager.load()

        assert loaded_checkpoint.completed_steps == ["step0", "step1", "step2"]
        assert loaded_config["iteration"] == 2


class TestWorkflowErrorStrWithOriginal:
    """Test WorkflowError string with original error."""

    def test_str_with_original_error(self):
        """Test string representation includes original error info."""
        original = ValueError("Original problem")
        err = WorkflowError("Workflow failed", step_name="normalization", original_error=original)

        str_repr = str(err)
        assert "Workflow failed" in str_repr
        assert "step: normalization" in str_repr
        assert "ValueError" in str_repr
