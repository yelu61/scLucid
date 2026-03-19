"""
Base configuration system for scLucid using Pydantic.

This module provides the foundation for all configuration classes in scLucid,
ensuring consistent validation, serialization, and documentation.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Literal, Optional, Self

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

logger = logging.getLogger(__name__)


class SclucidBaseConfig(BaseModel):
    """
    Base configuration class for all scLucid components.

    Provides:
    - Automatic validation via Pydantic
    - Serialization to/from dict and JSON
    - Path auto-creation for save_dir
    - Consistent configuration patterns

    Example:
        >>> config = MyConfig(save_dir="./results", n_jobs=4)
        >>> config.to_dict()  # Serialize
        >>> config.to_json_file("config.json")  # Save to file
        >>> config2 = MyConfig.from_json_file("config.json")  # Load from file
    """

    model_config = ConfigDict(
        # Allow extra fields for forward compatibility
        extra="ignore",
        # Validate during assignment
        validate_assignment=True,
        # Use field names in serialization, not aliases
        populate_by_name=True,
        # Allow arbitrary types (e.g., Callable)
        arbitrary_types_allowed=True,
    )

    # Common fields shared across most configs
    save_dir: Optional[str] = Field(
        default=None,
        description="Directory to save outputs. If None, nothing is saved.",
    )
    verbose: bool = Field(default=True, description="Whether to print progress messages.")
    plot: bool = Field(default=True, description="Whether to generate plots.")
    report: bool = Field(default=True, description="Whether to generate reports.")

    @field_validator("save_dir")
    @classmethod
    def _create_save_dir(cls, v: Optional[str]) -> Optional[str]:
        """Create save directory if specified."""
        if v is not None:
            Path(v).mkdir(parents=True, exist_ok=True)
        return v

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        return self.model_dump()

    def to_json(self, indent: int = 2) -> str:
        """Serialize configuration to JSON string."""
        return self.model_dump_json(indent=indent)

    def to_json_file(self, path: str) -> None:
        """Save configuration to JSON file."""
        Path(path).write_text(self.to_json())

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> Self:
        """Create configuration from dictionary."""
        return cls.model_validate(data)

    @classmethod
    def from_json(cls, json_str: str) -> Self:
        """Create configuration from JSON string."""
        return cls.model_validate_json(json_str)

    @classmethod
    def from_json_file(cls, path: str) -> Self:
        """Load configuration from JSON file."""
        return cls.from_json(Path(path).read_text())

    def validate(self) -> Self:
        """
        Explicit validation hook for subclasses.

        Called automatically during initialization, but can be called
        manually to re-validate after field changes.

        Returns:
            Self for method chaining
        """
        # Pydantic validates on creation, but subclasses can add custom logic
        return self


class WorkflowConfigBase(SclucidBaseConfig):
    """
    Base class for workflow-level configurations.

    Workflows are multi-step operations that process AnnData objects
    and may support checkpoint/resumption.
    """

    # Workflow-specific common fields
    n_jobs: int = Field(
        default=-1,
        ge=-1,
        description="Number of parallel jobs. -1 uses all CPUs.",
    )
    random_state: int = Field(default=42, description="Random seed for reproducibility.")

    @field_validator("n_jobs")
    @classmethod
    def _validate_n_jobs(cls, v: int) -> int:
        """Ensure n_jobs is valid."""
        if v == 0:
            raise ValueError("n_jobs cannot be 0. Use -1 for all CPUs or a positive integer.")
        return v


class ComputationConfig(SclucidBaseConfig):
    """Configuration for computational parameters."""

    backend: Literal["loky", "threading", "multiprocessing"] = Field(
        default="loky",
        description="Parallel backend for joblib.",
    )
    chunk_size: int = Field(
        default=1000,
        gt=0,
        description="Chunk size for processing large datasets.",
    )
    low_memory_mode: bool = Field(
        default=False,
        description="Enable memory-efficient processing (may be slower).",
    )


# Type alias for backward compatibility during migration
BaseConfig = SclucidBaseConfig

__all__ = [
    "SclucidBaseConfig",
    "BaseConfig",  # Backward compatibility
    "WorkflowConfigBase",
    "ComputationConfig",
]
