"""Shared evidence schema for scLucid review summaries.

The schema is intentionally small: modules keep their own domain-specific
review summaries, while this layer provides a common bundle for cross-stage
audit, reporting, and future evidence graph integration.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import ConfigDict, Field

from ..base_config import SclucidBaseConfig

EVIDENCE_SCHEMA_VERSION = "1.0"

EvidenceSource = Literal[
    "metric",
    "recommendation",
    "benchmark",
    "user_override",
    "context",
    "warning",
    "contract",
    "output_health",
    "downstream",
    "reproducibility",
]

ActionPriority = Literal["blocking", "required", "review", "optional"]


class EvidenceItem(SclucidBaseConfig):
    """A single piece of evidence supporting a module decision."""

    model_config = ConfigDict(extra="ignore")

    source: EvidenceSource = Field(description="Evidence source category.")
    name: str = Field(description="Stable evidence name.")
    value: Any = Field(default=None, description="Evidence value or compact summary.")
    confidence: float | None = Field(
        default=None,
        ge=0,
        le=1,
        description="Optional confidence score in [0, 1].",
    )
    rationale: str | None = Field(default=None, description="Why this evidence matters.")
    limitations: list[str] = Field(
        default_factory=list,
        description="Known caveats or situations requiring review.",
    )
    related_keys: list[str] = Field(
        default_factory=list,
        description="Pointers to source fields in the originating review summary.",
    )


class DecisionRecord(SclucidBaseConfig):
    """A common decision record linking recommended, applied, and evidence values."""

    model_config = ConfigDict(extra="ignore")

    parameter: str = Field(description="Parameter, label, or decision name.")
    recommended: Any = Field(default=None, description="Recommended value, if available.")
    applied: Any = Field(default=None, description="Actually applied value.")
    source: str = Field(default="unknown", description="Final decision source.")
    confidence: float | None = Field(default=None, ge=0, le=1)
    evidence: list[EvidenceItem] = Field(default_factory=list)
    user_override: bool = Field(default=False)
    downstream_impact: str | None = Field(default=None)


class ReviewAction(SclucidBaseConfig):
    """A prioritized human review action derived from evidence."""

    model_config = ConfigDict(extra="ignore")

    priority: ActionPriority = Field(description="Action priority.")
    action: str = Field(description="What the user should do.")
    rationale: str = Field(default="", description="Why this action is needed.")
    evidence_keys: list[str] = Field(default_factory=list)


class EvidenceBundle(SclucidBaseConfig):
    """Cross-module evidence bundle stored inside module review summaries."""

    model_config = ConfigDict(extra="ignore")

    schema_version: str = Field(default=EVIDENCE_SCHEMA_VERSION)
    module: str = Field(description="scLucid module name, e.g. 'qc'.")
    stage: str = Field(description="Workflow stage or entrypoint.")
    status: str = Field(default="unknown", description="Module-level readiness status.")
    confidence: float | None = Field(default=None, ge=0, le=1)
    context: dict[str, Any] = Field(default_factory=dict)
    decisions: list[DecisionRecord] = Field(default_factory=list)
    evidence_chain: list[EvidenceItem] = Field(default_factory=list)
    action_items: list[ReviewAction] = Field(default_factory=list)
    reproducibility: dict[str, Any] = Field(default_factory=dict)
    related_review_keys: list[str] = Field(default_factory=list)


def model_to_dict(model: SclucidBaseConfig) -> dict[str, Any]:
    """Serialize an evidence model to a plain JSON-compatible dictionary."""
    return model.model_dump(
        mode="json",
        fallback=lambda value: value.item() if hasattr(value, "item") else str(value),
    )


__all__ = [
    "ActionPriority",
    "DecisionRecord",
    "EVIDENCE_SCHEMA_VERSION",
    "EvidenceBundle",
    "EvidenceItem",
    "EvidenceSource",
    "ReviewAction",
    "model_to_dict",
]
