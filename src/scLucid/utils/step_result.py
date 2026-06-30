"""Unified step result semantics for scLucid workflows.

Provides a lightweight, auditable record for every workflow step across all
modules. Step results are designed to be HDF5-serializable and human-readable,
supporting both machine audit (tumor readiness scoring) and manual review.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Literal, Optional, Union

from pydantic import ConfigDict, Field, field_validator

from ..base_config import SclucidBaseConfig

StepStatus = Literal["completed", "skipped", "failed", "degraded"]
EvidenceLevel = Literal["validated_core", "heuristic", "exploratory", "unavailable"]

STEP_STATUS_ORDER = ("failed", "degraded", "skipped", "completed")


class StepResult(SclucidBaseConfig):
    """A single auditable workflow step result.

    Attributes
    ----------
    name
        Stable step identifier (e.g., ``"malignancy_interpretation"``).
    status
        One of ``completed``, ``skipped``, ``failed``, ``degraded``.
    evidence_level
        Claim boundary for outputs produced by this step.
    outputs
        Compact summary of produced keys/artifacts (not the full data).
    warnings
        Human-readable warnings that do not block the step.
    error
        Error message when ``status`` is ``failed`` or ``degraded``.
    duration_seconds
        Step wall-clock duration, when available.
    """

    model_config = ConfigDict(extra="ignore")

    name: str = Field(description="Stable step identifier.")
    status: StepStatus = Field(default="completed")
    evidence_level: EvidenceLevel = Field(default="heuristic")
    outputs: Dict[str, Any] = Field(default_factory=dict)
    warnings: List[str] = Field(default_factory=list)
    error: Optional[str] = Field(default=None)
    duration_seconds: Optional[float] = Field(default=None, ge=0)
    timestamp: str = Field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    @field_validator("error")
    @classmethod
    def _error_requires_non_completed(cls, v: Optional[str], info) -> Optional[str]:
        values = info.data
        if v and values.get("status") == "completed":
            return f"{v} (note: status is completed)"
        return v

    def model_dump(self, **kwargs) -> Dict[str, Any]:
        # Ensure error/None stays JSON-safe
        data = super().model_dump(**kwargs)
        if "error" in data and data["error"] is None:
            data["error"] = ""
        return data

    def to_storage_dict(self) -> Dict[str, Any]:
        """Return a plain dict safe for adata.uns HDF5 storage."""
        return dict(self.model_dump(mode="json"))

    @classmethod
    def from_exception(
        cls,
        name: str,
        exc: Exception,
        evidence_level: EvidenceLevel = "unavailable",
        outputs: Optional[Dict[str, Any]] = None,
        warnings: Optional[List[str]] = None,
        degraded: bool = False,
    ) -> "StepResult":
        """Create a failed/degraded StepResult from an exception."""
        return cls(
            name=name,
            status="degraded" if degraded else "failed",
            evidence_level=evidence_level,
            outputs=outputs or {},
            warnings=warnings or [],
            error=f"{type(exc).__name__}: {exc}",
        )

    @classmethod
    def skipped(
        cls,
        name: str,
        reason: str,
        evidence_level: EvidenceLevel = "unavailable",
        outputs: Optional[Dict[str, Any]] = None,
    ) -> "StepResult":
        """Create a skipped StepResult with a structured reason."""
        return cls(
            name=name,
            status="skipped",
            evidence_level=evidence_level,
            outputs=outputs or {},
            warnings=[reason],
        )

    @classmethod
    def degraded(
        cls,
        name: str,
        reason: str,
        evidence_level: EvidenceLevel = "heuristic",
        outputs: Optional[Dict[str, Any]] = None,
    ) -> "StepResult":
        """Create a degraded StepResult with a structured reason."""
        return cls(
            name=name,
            status="degraded",
            evidence_level=evidence_level,
            outputs=outputs or {},
            warnings=[reason],
        )


def rollup_step_status(results: List[StepResult]) -> StepStatus:
    """Return the most severe status across a list of step results."""
    if not results:
        return "completed"
    statuses = {r.status for r in results}
    for status in STEP_STATUS_ORDER:
        if status in statuses:
            return status
    return "completed"


def summarize_step_results(results: List[StepResult]) -> Dict[str, Any]:
    """Produce a compact summary of step results for review summaries."""
    by_status: Dict[str, List[str]] = {
        "completed": [],
        "degraded": [],
        "skipped": [],
        "failed": [],
    }
    evidence_counts: Dict[str, int] = {
        "validated_core": 0,
        "heuristic": 0,
        "exploratory": 0,
        "unavailable": 0,
    }
    all_warnings: List[str] = []
    errors: List[str] = []
    for r in results:
        by_status[r.status].append(r.name)
        evidence_counts[r.evidence_level] = evidence_counts.get(r.evidence_level, 0) + 1
        all_warnings.extend(r.warnings)
        if r.error:
            errors.append(f"{r.name}: {r.error}")
    return {
        "n_steps": len(results),
        "overall_status": rollup_step_status(results),
        "by_status": {k: v for k, v in by_status.items() if v},
        "evidence_summary": evidence_counts,
        "warnings": list(dict.fromkeys(all_warnings)),
        "errors": errors,
    }


def step_results_to_storage(results: List[StepResult]) -> Dict[str, Dict[str, Any]]:
    """Serialize a list of StepResult objects for adata.uns storage."""
    return {str(i): r.to_storage_dict() for i, r in enumerate(results)}


def step_results_from_storage(
    data: Union[List[Dict[str, Any]], Dict[str, Dict[str, Any]]],
) -> List[StepResult]:
    """Deserialize StepResult objects from adata.uns storage."""
    if isinstance(data, dict):
        items = [
            data[key]
            for key in sorted(data, key=lambda value: int(value) if str(value).isdigit() else str(value))
        ]
        return [StepResult.model_validate(item) for item in items]
    return [StepResult.model_validate(item) for item in data]


__all__ = [
    "EvidenceLevel",
    "StepResult",
    "StepStatus",
    "rollup_step_status",
    "step_results_from_storage",
    "step_results_to_storage",
    "summarize_step_results",
]
