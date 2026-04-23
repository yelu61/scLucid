"""
Workflow checkpoint and resumption system for scLucid.

Provides checkpoint/resume functionality for long-running workflows
to avoid recomputation after failures.

Example:
    # Run with automatic checkpointing
    checkpoint_manager = WorkflowCheckpoint(results_dir="./results")

    adata = checkpoint_manager.run_or_resume(
        adata,
        config,
        workflow_steps=["qc", "preprocess", "cluster", "annotate"]
    )

    # Or resume from a specific step
    adata = checkpoint_manager.run_or_resume(
        adata,
        config,
        resume_from="preprocess",  # Skip QC if already done
        workflow_steps=["qc", "preprocess", "cluster", "annotate"]
    )
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from anndata import AnnData

log = logging.getLogger(__name__)


class WorkflowStep(str, Enum):
    """Standard workflow steps."""

    QC = "qc"
    PREPROCESS = "preprocess"
    CLUSTER = "cluster"
    ANNOTATE = "annotate"
    DE = "differential_expression"
    ENRICHMENT = "enrichment"


class StepStatus(str, Enum):
    """Status of a workflow step."""

    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    SKIPPED = "skipped"


class CheckpointInfo:
    """Information about a saved checkpoint."""

    def __init__(
        self,
        step: str,
        timestamp: datetime,
        adata_path: Path,
        config_hash: str,
        metadata: Optional[Dict[str, Any]] = None,
    ):
        self.step = step
        self.timestamp = timestamp
        self.adata_path = adata_path
        self.config_hash = config_hash
        self.metadata = metadata or {}

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "step": self.step,
            "timestamp": self.timestamp.isoformat(),
            "adata_path": str(self.adata_path),
            "config_hash": self.config_hash,
            "metadata": self.metadata,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> CheckpointInfo:
        """Deserialize from dictionary."""
        return cls(
            step=data["step"],
            timestamp=datetime.fromisoformat(data["timestamp"]),
            adata_path=Path(data["adata_path"]),
            config_hash=data["config_hash"],
            metadata=data.get("metadata", {}),
        )


class WorkflowState:
    """Complete state of a workflow execution."""

    def __init__(
        self,
        workflow_id: str,
        steps: Dict[str, StepStatus] = None,
        checkpoints: Dict[str, CheckpointInfo] = None,
        current_step: Optional[str] = None,
        start_time: Optional[datetime] = None,
        end_time: Optional[datetime] = None,
    ):
        self.workflow_id = workflow_id
        self.steps = steps or {}
        self.checkpoints = checkpoints or {}
        self.current_step = current_step
        self.start_time = start_time or datetime.now()
        self.end_time = end_time

    def to_dict(self) -> Dict[str, Any]:
        """Serialize to dictionary."""
        return {
            "workflow_id": self.workflow_id,
            "steps": {k: v.value for k, v in self.steps.items()},
            "checkpoints": {k: v.to_dict() for k, v in self.checkpoints.items()},
            "current_step": self.current_step,
            "start_time": self.start_time.isoformat() if self.start_time else None,
            "end_time": self.end_time.isoformat() if self.end_time else None,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> WorkflowState:
        """Deserialize from dictionary."""
        return cls(
            workflow_id=data["workflow_id"],
            steps={k: StepStatus(v) for k, v in data.get("steps", {}).items()},
            checkpoints={
                k: CheckpointInfo.from_dict(v) for k, v in data.get("checkpoints", {}).items()
            },
            current_step=data.get("current_step"),
            start_time=(
                datetime.fromisoformat(data["start_time"]) if data.get("start_time") else None
            ),
            end_time=datetime.fromisoformat(data["end_time"]) if data.get("end_time") else None,
        )

    def get_last_completed_step(self) -> Optional[str]:
        """Get the last successfully completed step."""
        completed = [step for step, status in self.steps.items() if status == StepStatus.COMPLETED]
        return completed[-1] if completed else None

    def is_step_completed(self, step: str) -> bool:
        """Check if a step is completed."""
        return self.steps.get(step) == StepStatus.COMPLETED


class WorkflowCheckpoint:
    """
    Checkpoint manager for scLucid workflows.

    Manages saving and loading of workflow state to enable resumption
    after interruptions.

    Args:
        results_dir: Directory to save checkpoints and state
        workflow_id: Unique identifier for this workflow (default: timestamp)
        auto_save: Whether to auto-save after each step
        checkpoint_format: Format for checkpoint files ('h5ad' or 'zarr')
    """

    def __init__(
        self,
        results_dir: str,
        workflow_id: Optional[str] = None,
        auto_save: bool = True,
        checkpoint_format: str = "h5ad",
    ):
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.workflow_id = workflow_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        self.auto_save = auto_save
        self.checkpoint_format = checkpoint_format

        self.state = WorkflowState(workflow_id=self.workflow_id)
        self.checkpoint_dir = self.results_dir / "checkpoints"
        self.checkpoint_dir.mkdir(exist_ok=True)

        # Try to load existing state
        self._load_state()

    def _get_state_path(self) -> Path:
        """Get path to state file."""
        return self.results_dir / f"workflow_state_{self.workflow_id}.json"

    def _get_checkpoint_path(self, step: str) -> Path:
        """Get path for a step's checkpoint file."""
        ext = ".h5ad" if self.checkpoint_format == "h5ad" else ".zarr"
        return self.checkpoint_dir / f"{step}_checkpoint{ext}"

    def _load_state(self) -> None:
        """Load existing workflow state if available."""
        state_path = self._get_state_path()
        if state_path.exists():
            try:
                with open(state_path) as f:
                    data = json.load(f)
                self.state = WorkflowState.from_dict(data)
                log.info(f"Loaded existing workflow state: {self.workflow_id}")
            except Exception as e:
                log.warning(f"Could not load state: {e}")

    def _save_state(self) -> None:
        """Save current workflow state."""
        state_path = self._get_state_path()
        with open(state_path, "w") as f:
            json.dump(self.state.to_dict(), f, indent=2)

    def _compute_config_hash(self, config: Any) -> str:
        """Compute hash of configuration for change detection."""
        import hashlib

        # Convert config to string representation
        if hasattr(config, "to_dict"):
            config_str = json.dumps(config.to_dict(), sort_keys=True)
        elif hasattr(config, "__dict__"):
            config_str = json.dumps(config.__dict__, sort_keys=True, default=str)
        else:
            config_str = str(config)

        return hashlib.md5(config_str.encode()).hexdigest()[:16]

    def save_checkpoint(
        self,
        step: str,
        adata: AnnData,
        config: Any,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Path:
        """
        Save a checkpoint for the current step.

        Args:
            step: Name of the current step
            adata: AnnData object to save
            config: Configuration used for this step
            metadata: Optional metadata about the checkpoint

        Returns:
            Path to the saved checkpoint file
        """
        checkpoint_path = self._get_checkpoint_path(step)
        config_hash = self._compute_config_hash(config)

        # Save AnnData
        if self.checkpoint_format == "h5ad":
            adata.write_h5ad(checkpoint_path)
        else:
            adata.write_zarr(checkpoint_path)

        # Update state
        self.state.checkpoints[step] = CheckpointInfo(
            step=step,
            timestamp=datetime.now(),
            adata_path=checkpoint_path,
            config_hash=config_hash,
            metadata=metadata or {},
        )
        self.state.steps[step] = StepStatus.COMPLETED
        self._save_state()

        log.info(f"Checkpoint saved: {step} -> {checkpoint_path}")
        return checkpoint_path

    def load_checkpoint(
        self,
        step: str,
        config: Any,
        check_hash: bool = True,
    ) -> Optional[AnnData]:
        """
        Load a checkpoint if available and valid.

        Args:
            step: Name of the step to load
            config: Current configuration (for hash comparison)
            check_hash: Whether to verify config hasn't changed

        Returns:
            Loaded AnnData or None if checkpoint not available
        """
        if step not in self.state.checkpoints:
            return None

        checkpoint_info = self.state.checkpoints[step]

        # Verify config hash
        if check_hash:
            current_hash = self._compute_config_hash(config)
            if current_hash != checkpoint_info.config_hash:
                log.warning(f"Config changed since checkpoint '{step}'. " "Recomputing...")
                return None

        # Load AnnData
        try:
            if self.checkpoint_format == "h5ad":
                adata = AnnData.read_h5ad(checkpoint_info.adata_path)
            else:
                adata = AnnData.read_zarr(checkpoint_info.adata_path)

            log.info(f"Checkpoint loaded: {step}")
            return adata
        except Exception as e:
            log.error(f"Failed to load checkpoint '{step}': {e}")
            return None

    def run_step(
        self,
        step: str,
        adata: AnnData,
        config: Any,
        step_func: Callable[[AnnData, Any], AnnData],
        force_recompute: bool = False,
    ) -> AnnData:
        """
        Run a workflow step with checkpoint support.

        Args:
            step: Name of this step
            adata: Input AnnData
            config: Configuration for this step
            step_func: Function to execute (adata, config) -> adata
            force_recompute: Whether to recompute even if checkpoint exists

        Returns:
            Result AnnData (from checkpoint or computed)
        """
        # Check if we can use a checkpoint
        if not force_recompute and self.state.is_step_completed(step):
            cached = self.load_checkpoint(step, config)
            if cached is not None:
                log.info(f"Step '{step}' completed (from checkpoint)")
                return cached

        # Execute the step
        log.info(f"Running step: {step}")
        self.state.current_step = step
        self.state.steps[step] = StepStatus.RUNNING
        self._save_state()

        try:
            result = step_func(adata, config)

            # Save checkpoint
            if self.auto_save:
                self.save_checkpoint(step, result, config)

            return result

        except Exception as e:
            self.state.steps[step] = StepStatus.FAILED
            self._save_state()
            log.error(f"Step '{step}' failed: {e}")
            raise

    def get_resume_point(
        self,
        workflow_steps: List[str],
        resume_from: Optional[str] = None,
    ) -> tuple[int, Optional[AnnData]]:
        """
        Determine where to resume workflow.

        Args:
            workflow_steps: List of steps in order
            resume_from: Step name to resume from, or None to auto-detect

        Returns:
            Tuple of (step_index, checkpoint_adata or None)
        """
        if resume_from is None:
            # Auto-detect: find last completed step
            last_completed = self.state.get_last_completed_step()
            if last_completed is None:
                return 0, None

            if last_completed in workflow_steps:
                idx = workflow_steps.index(last_completed)
                checkpoint = self.load_checkpoint(last_completed, config=None, check_hash=False)
                return idx + 1, checkpoint

            return 0, None

        # Resume from specific step
        if resume_from in workflow_steps:
            idx = workflow_steps.index(resume_from)
            # Try to load checkpoint from previous step
            if idx > 0:
                prev_step = workflow_steps[idx - 1]
                checkpoint = self.load_checkpoint(prev_step, config=None, check_hash=False)
                return idx, checkpoint
            return idx, None

        log.warning(f"Unknown resume step: {resume_from}. Starting from beginning.")
        return 0, None

    def finalize(self, success: bool = True) -> None:
        """Mark workflow as complete."""
        self.state.end_time = datetime.now()
        self.state.current_step = None
        self._save_state()

        status = "completed" if success else "failed"
        log.info(f"Workflow {self.workflow_id} {status}")

    def list_checkpoints(self) -> Dict[str, CheckpointInfo]:
        """List all available checkpoints."""
        return self.state.checkpoints.copy()

    def clear_checkpoints(self, keep_last: int = 1) -> None:
        """
        Clear old checkpoints, keeping only the most recent.

        Args:
            keep_last: Number of recent checkpoints to keep
        """
        checkpoints = sorted(
            self.state.checkpoints.items(),
            key=lambda x: x[1].timestamp,
        )

        to_remove = checkpoints[:-keep_last] if keep_last > 0 else checkpoints

        for step, info in to_remove:
            if info.adata_path.exists():
                info.adata_path.unlink()
                log.debug(f"Removed checkpoint: {step}")
            del self.state.checkpoints[step]

        self._save_state()


def run_workflow_with_checkpoints(
    adata: AnnData,
    steps: List[tuple[str, Callable, Any]],
    results_dir: str,
    resume_from: Optional[str] = None,
    workflow_id: Optional[str] = None,
) -> AnnData:
    """
    Helper function to run a multi-step workflow with checkpoint support.

    Args:
        adata: Input AnnData
        steps: List of (step_name, step_function, config) tuples
        results_dir: Directory for checkpoints and results
        resume_from: Step name to resume from, or None
        workflow_id: Optional workflow identifier

    Returns:
        Final AnnData result

    Example:
        steps = [
            ("qc", run_qc, qc_config),
            ("preprocess", run_preprocessing, prep_config),
            ("cluster", run_clustering, cluster_config),
        ]

        adata = run_workflow_with_checkpoints(
            adata, steps, "./results", resume_from="preprocess"
        )
    """
    step_names = [s[0] for s in steps]
    checkpoint_mgr = WorkflowCheckpoint(
        results_dir=results_dir,
        workflow_id=workflow_id,
    )

    # Determine resume point
    start_idx, cached_adata = checkpoint_mgr.get_resume_point(step_names, resume_from)

    if cached_adata is not None:
        adata = cached_adata
        log.info(f"Resuming workflow from step {start_idx}: {step_names[start_idx]}")

    # Run remaining steps
    try:
        for step_name, step_func, config in steps[start_idx:]:
            adata = checkpoint_mgr.run_step(step_name, adata, config, step_func)

        checkpoint_mgr.finalize(success=True)

    except Exception:
        checkpoint_mgr.finalize(success=False)
        raise

    return adata
