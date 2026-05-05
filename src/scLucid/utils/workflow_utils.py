"""
Generic workflow utilities for scLucid modules.

This module provides shared functionality for:
- Progress tracking (tqdm integration)
- Error recovery and partial results management
- Workflow step iteration
- Common workflow exceptions

Used by: qc, preprocess, and analysis modules.
"""

import logging
import pickle
from abc import ABC, abstractmethod
from collections.abc import Iterable
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import (
    Any,
    Callable,
    Dict,
    Generic,
    List,
    Optional,
    Tuple,
    TypeVar,
    Union,
)

from anndata import AnnData, read_h5ad

log = logging.getLogger(__name__)

T = TypeVar("T")
C = TypeVar("C")  # Config type


# ============================================================================
# Progress Tracking
# ============================================================================


def get_progress_bar(
    iterable: Optional[Iterable[T]] = None,
    desc: str = "Processing",
    enabled: bool = True,
    total: Optional[int] = None,
    unit: str = "step",
    **kwargs,
) -> Iterable[T]:
    """
    Wrap an iterable with a tqdm progress bar if enabled.

    Args:
        iterable: The iterable to wrap. If None, returns a tqdm instance.
        desc: Description for the progress bar.
        enabled: Whether to show the progress bar.
        total: Total number of items (optional).
        unit: Unit name displayed after the progress bar.
        **kwargs: Additional arguments passed to tqdm.

    Returns:
        The original iterable or a tqdm-wrapped iterable.

    Examples:
        >>> for item in get_progress_bar(items, desc="QC", enabled=True):
        ...     process(item)

        >>> # With total when iterable doesn't have __len__
        >>> for i in get_progress_bar(range(n), total=n, desc="Processing"):
        ...     process(i)
    """
    if not enabled:
        return iterable if iterable is not None else range(total or 0)

    try:
        from tqdm import tqdm

        if iterable is not None:
            return tqdm(iterable, desc=desc, total=total, unit=unit, **kwargs)
        else:
            return tqdm(total=total, desc=desc, unit=unit, **kwargs)
    except ImportError:
        log.debug("tqdm not available, running without progress bar")
        return iterable if iterable is not None else range(total or 0)


def progress_decorator(desc: str = "Processing", unit: str = "step"):
    """
    Decorator to add progress bar to a function.

    Args:
        desc: Description for the progress bar.
        unit: Unit name for the progress bar.

    Examples:
        >>> @progress_decorator(desc="Finding markers", unit="cluster")
        ... def find_markers(clusters, config):
        ...     for cluster in clusters:
        ...         yield process_cluster(cluster)
    """

    def decorator(func: Callable) -> Callable:
        def wrapper(*args, show_progress: bool = True, **kwargs):
            result = func(*args, **kwargs)
            if hasattr(result, "__iter__") and not isinstance(result, (str, bytes)):
                return get_progress_bar(result, desc=desc, enabled=show_progress, unit=unit)
            return result

        return wrapper

    return decorator


# ============================================================================
# Workflow Exceptions
# ============================================================================


class WorkflowError(Exception):
    """
    Generic workflow error with step context.

    Attributes:
        message: Error message
        step_name: Name of the step where error occurred
        original_error: The original exception that caused this error
    """

    def __init__(
        self,
        message: str,
        step_name: str = "unknown",
        original_error: Optional[Exception] = None,
        module: str = "unknown",
        hint: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None,
    ):
        super().__init__(message)
        self.step_name = step_name
        self.original_error = original_error
        self.module = module
        self.hint = hint
        self.context = context or {}

    def __str__(self) -> str:
        message = str(self.args[0]) if self.args else "Unknown error"
        prefix = f"[{self.module}:{self.step_name}] " if self.module != "unknown" else ""
        if self.original_error:
            message = (
                f"{prefix}{message} "
                f"(step: {self.step_name}, caused by: {type(self.original_error).__name__})"
            )
        else:
            message = f"{prefix}{message} (step: {self.step_name})"
        if self.hint:
            message = f"{message}. Hint: {self.hint}"
        return message

    def to_dict(self) -> Dict[str, Any]:
        """Return a structured, JSON-serializable error record."""
        return {
            "module": self.module,
            "step_name": self.step_name,
            "message": str(self.args[0]) if self.args else "Unknown error",
            "error_type": (
                type(self.original_error).__name__ if self.original_error else type(self).__name__
            ),
            "hint": self.hint,
            "context": self.context,
        }


class StepError(WorkflowError):
    """Error that occurred during a specific workflow step."""

    pass


class RecoveryError(WorkflowError):
    """Error during error recovery process."""

    pass


# ============================================================================
# Partial Results Management
# ============================================================================


@dataclass
class WorkflowCheckpoint:
    """
    Checkpoint data for workflow recovery.

    This class stores the state of a workflow at a specific point,
    allowing for recovery and resumption.
    """

    completed_steps: List[str] = field(default_factory=list)
    failed_step: Optional[str] = None
    error_message: Optional[str] = None
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert checkpoint to dictionary."""
        return {
            "completed_steps": self.completed_steps,
            "failed_step": self.failed_step,
            "error_message": self.error_message,
            "timestamp": self.timestamp,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "WorkflowCheckpoint":
        """Create checkpoint from dictionary."""
        return cls(
            completed_steps=data.get("completed_steps", []),
            failed_step=data.get("failed_step"),
            error_message=data.get("error_message"),
            timestamp=data.get("timestamp", datetime.now().isoformat()),
            metadata=data.get("metadata", {}),
        )


class PartialResultManager:
    """
    Manager for saving and loading partial workflow results.

    This class handles the persistence of workflow state, allowing
    workflows to be resumed after failure.
    """

    def __init__(self, save_dir: Union[str, Path]):
        self.save_dir = Path(save_dir)
        self.save_dir.mkdir(parents=True, exist_ok=True)

    def save(
        self,
        adata: AnnData,
        checkpoint: WorkflowCheckpoint,
        config: Optional[Any] = None,
        name: str = "partial",
    ) -> Path:
        """
        Save partial results to disk.

        Args:
            adata: AnnData object with processed data
            checkpoint: Workflow checkpoint information
            config: Workflow configuration (optional)
            name: Base name for saved files

        Returns:
            Path to the save directory
        """
        # Save AnnData
        adata_path = self.save_dir / f"{name}_result.h5ad"
        adata.write_h5ad(adata_path)

        # Save checkpoint
        checkpoint_path = self.save_dir / f"{name}_checkpoint.pkl"
        with open(checkpoint_path, "wb") as f:
            pickle.dump(checkpoint.to_dict(), f)

        # Save config if provided
        if config is not None:
            config_path = self.save_dir / f"{name}_config.pkl"
            with open(config_path, "wb") as f:
                pickle.dump(config, f)

        log.info(f"Partial results saved to: {self.save_dir}")
        return self.save_dir

    def load(self, name: str = "partial") -> Tuple[AnnData, WorkflowCheckpoint, Optional[Any]]:
        """
        Load partial results from disk.

        Args:
            name: Base name for saved files

        Returns:
            Tuple of (adata, checkpoint, config)
        """
        # Load AnnData
        adata_path = self.save_dir / f"{name}_result.h5ad"
        if not adata_path.exists():
            raise FileNotFoundError(f"No saved results found at: {adata_path}")
        adata = read_h5ad(adata_path)

        # Load checkpoint
        checkpoint_path = self.save_dir / f"{name}_checkpoint.pkl"
        with open(checkpoint_path, "rb") as f:
            checkpoint_data = pickle.load(f)
        checkpoint = WorkflowCheckpoint.from_dict(checkpoint_data)

        # Load config if exists
        config_path = self.save_dir / f"{name}_config.pkl"
        config = None
        if config_path.exists():
            with open(config_path, "rb") as f:
                config = pickle.load(f)

        return adata, checkpoint, config

    def exists(self, name: str = "partial") -> bool:
        """Check if saved results exist."""
        adata_path = self.save_dir / f"{name}_result.h5ad"
        checkpoint_path = self.save_dir / f"{name}_checkpoint.pkl"
        return adata_path.exists() and checkpoint_path.exists()

    def cleanup(self, name: str = "partial") -> None:
        """Remove saved partial results."""
        for file in self.save_dir.glob(f"{name}_*"):
            file.unlink()
        log.info(f"Cleaned up partial results in: {self.save_dir}")


# ============================================================================
# Generic Workflow Step Iterator
# ============================================================================


class WorkflowStepIterator:
    """
    Iterator for workflow steps with progress tracking and error handling.

    This class provides a standardized way to iterate through workflow steps
    with support for:
    - Progress bars
    - Error recovery
    - Step filtering (skip completed steps)
    - Error handling strategies
    """

    def __init__(
        self,
        steps: List[str],
        completed_steps: Optional[List[str]] = None,
        show_progress: bool = True,
        desc: str = "Workflow",
        on_error: str = "raise",  # "raise", "skip", "save"
        recovery_manager: Optional[PartialResultManager] = None,
    ):
        self.all_steps = steps
        self.completed_steps = set(completed_steps or [])
        self.show_progress = show_progress
        self.desc = desc
        self.on_error = on_error
        self.recovery_manager = recovery_manager

        # Filter steps to run
        self.steps_to_run = [s for s in steps if s not in self.completed_steps]
        self.current_step: Optional[str] = None
        self.successful_steps: List[str] = []

    def __iter__(self):
        """Iterate through steps with progress bar."""
        iterator = get_progress_bar(
            self.steps_to_run,
            desc=self.desc,
            enabled=self.show_progress,
            total=len(self.steps_to_run),
            unit="step",
        )
        for step in iterator:
            self.current_step = step
            yield step
            self.successful_steps.append(step)

    def get_checkpoint(self, error: Optional[Exception] = None) -> WorkflowCheckpoint:
        """Create a checkpoint from current state."""
        return WorkflowCheckpoint(
            completed_steps=self.successful_steps.copy(),
            failed_step=self.current_step if error else None,
            error_message=str(error) if error else None,
            metadata={
                "total_steps": len(self.all_steps),
                "remaining_steps": len(self.steps_to_run) - len(self.successful_steps),
            },
        )


# ============================================================================
# Base Workflow Class
# ============================================================================


class BaseWorkflow(ABC, Generic[C]):
    """
    Abstract base class for scLucid workflows.

    This class provides a common foundation for QC, preprocessing, and analysis
    workflows with standardized error handling, progress tracking, and recovery.

    Type Parameters:
        C: The configuration class type for this workflow
    """

    def __init__(
        self,
        config: Optional[C] = None,
        show_progress: bool = True,
        error_recovery: bool = False,
        recovery_save_dir: Optional[str] = None,
    ):
        self.config = config
        self.show_progress = show_progress
        self.error_recovery = error_recovery
        self.recovery_manager = (
            PartialResultManager(recovery_save_dir) if recovery_save_dir else None
        )
        self._step_iterator: Optional[WorkflowStepIterator] = None

    @abstractmethod
    def get_steps(self) -> List[str]:
        """Return list of workflow step names."""
        pass

    @abstractmethod
    def execute_step(self, adata: AnnData, step_name: str) -> AnnData:
        """Execute a single workflow step."""
        pass

    def run(
        self,
        adata: AnnData,
        resume_from: Optional[str] = None,
        on_error: str = "raise",
    ) -> AnnData:
        """
        Run the complete workflow.

        Args:
            adata: Input AnnData object
            resume_from: Directory to resume from (optional)
            on_error: Error handling strategy ("raise", "skip", "save")

        Returns:
            Processed AnnData object
        """
        # Handle resume
        completed_steps: List[str] = []
        if resume_from and self.recovery_manager is None:
            self.recovery_manager = PartialResultManager(resume_from)

        if resume_from and self.recovery_manager.exists():
            adata, checkpoint, _ = self.recovery_manager.load()
            completed_steps = checkpoint.completed_steps
            log.info(f"Resumed from checkpoint. Completed: {completed_steps}")

        # Create step iterator
        steps = self.get_steps()
        self._step_iterator = WorkflowStepIterator(
            steps=steps,
            completed_steps=completed_steps,
            show_progress=self.show_progress,
            desc=self.__class__.__name__,
            on_error=on_error,
            recovery_manager=self.recovery_manager,
        )

        # Execute steps
        try:
            for step_name in self._step_iterator:
                adata = self.execute_step(adata, step_name)

        except Exception as e:
            if self.error_recovery and self.recovery_manager:
                checkpoint = self._step_iterator.get_checkpoint(e)
                self.recovery_manager.save(adata, checkpoint, self.config, name="partial")
                log.warning(f"Workflow saved to: {self.recovery_manager.save_dir}")

                if on_error == "save":
                    return adata

            raise WorkflowError(
                f"Workflow failed at step '{self._step_iterator.current_step}'",
                step_name=self._step_iterator.current_step or "unknown",
                original_error=e,
            )

        return adata


# ============================================================================
# Convenience Functions
# ============================================================================


def with_error_recovery(
    save_dir: str,
    step_name: str = "step",
    adata_arg_name: str = "adata",
):
    """
    Decorator to execute a function with error recovery.

    Args:
        save_dir: Directory to save partial results on failure
        step_name: Name of the step for logging
        adata_arg_name: Name of the argument containing the AnnData object

    Returns:
        Decorated function

    Example:
        >>> @with_error_recovery(save_dir="/tmp/recovery", step_name="qc")
        ... def my_qc_step(adata):
        ...     return process_adata(adata)

    Raises:
        WorkflowError: If func fails and recovery is not possible
    """
    from functools import wraps

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            manager = PartialResultManager(save_dir)
            checkpoint = WorkflowCheckpoint(completed_steps=[], failed_step=None)

            # Get adata from args or kwargs
            adata = None
            if args:
                # Try to get first positional arg that looks like AnnData
                from anndata import AnnData as AnnDataClass

                for arg in args:
                    if isinstance(arg, AnnDataClass):
                        adata = arg
                        break
            if adata is None and adata_arg_name in kwargs:
                adata = kwargs[adata_arg_name]

            try:
                result = func(*args, **kwargs)
                return result
            except Exception as e:
                log.error(f"Step '{step_name}' failed: {e}")
                checkpoint.failed_step = step_name
                checkpoint.error_message = str(e)
                if adata is not None:
                    manager.save(adata, checkpoint, name="failed")
                raise WorkflowError(
                    f"Step '{step_name}' failed", step_name=step_name, original_error=e
                )

        return wrapper

    return decorator


def merge_partial_results(
    results: List[Tuple[str, Optional[AnnData]]],
    original_obs_names: List[str],
) -> Tuple[AnnData, List[str]]:
    """
    Merge partial results from parallel processing.

    Args:
        results: List of (name, adata) tuples. Failed items have None adata.
        original_obs_names: Original observation names to preserve order

    Returns:
        Tuple of (merged_adata, failed_names)
    """
    # Separate successful and failed
    successful = [(name, data) for name, data in results if data is not None]
    failed = [name for name, data in results if data is None]

    if not successful:
        raise RuntimeError("All samples failed processing")

    # Use first successful result as base
    merged = successful[0][1].copy()

    # Merge additional data if any
    obs_to_idx = {name: i for i, name in enumerate(original_obs_names)}

    for name, data in successful[1:]:
        for col in data.obs.columns:
            if col not in merged.obs:
                merged.obs[col] = np.nan
            for obs_name in data.obs_names:
                if obs_name in obs_to_idx:
                    merged.obs.loc[obs_name, col] = data.obs.loc[obs_name, col]

    return merged, failed


# ============================================================================
# Module Exports
# ============================================================================

__all__ = [
    # Progress tracking
    "get_progress_bar",
    "progress_decorator",
    # Exceptions
    "WorkflowError",
    "StepError",
    "RecoveryError",
    # Checkpoint management
    "WorkflowCheckpoint",
    "PartialResultManager",
    # Workflow iteration
    "WorkflowStepIterator",
    "BaseWorkflow",
    # Convenience functions
    "with_error_recovery",
    "merge_partial_results",
]
