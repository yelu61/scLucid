"""Global configuration for scLucid."""

from __future__ import annotations

import logging
import threading
import warnings
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Literal, Optional

from pydantic import Field, model_validator

from .runtime import effective_n_jobs

try:
    from .base_config import SclucidBaseConfig
except ImportError:
    from pydantic import BaseModel

    class SclucidBaseConfig(BaseModel):
        """Fallback base config when base_config module is unavailable."""

        model_config = {"extra": "ignore"}


_config_lock = threading.Lock()


class GlobalConfig(SclucidBaseConfig):
    """Global configuration settings for scLucid."""

    model_config = {"extra": "ignore"}

    # Computational settings
    n_jobs: int = Field(default=-1, description="Number of parallel jobs (-1 for all cores)")
    random_state: int = Field(default=42, description="Random seed for reproducibility")
    backend: Literal["loky", "threading", "multiprocessing"] = Field(
        default="loky", description="Parallel backend for joblib"
    )

    # Logging settings
    verbosity: int = Field(default=1, description="Logging level (0: WARNING, 1: INFO, 2: DEBUG)")
    log_file: Optional[Path] = Field(default=None, description="Path to log file")

    # Cache settings
    cache_dir: Optional[Path] = Field(default=None, description="Directory for caching")
    use_cache: bool = Field(default=True, description="Whether to use caching")

    # Plotting settings
    plot_backend: Literal["matplotlib", "plotly"] = Field(
        default="matplotlib", description="Plotting backend"
    )
    figure_dpi: int = Field(default=100, description="Figure DPI")
    figure_format: str = Field(default="png", description="Figure format (png, pdf, svg)")
    color_palette: str = Field(default="tab20", description="Color palette name")
    plot_theme: str = Field(default="default", description="Plot theme")
    font_style: Optional[Literal["nature", "cell", "traditional"]] = Field(
        default=None, description="Academic font style"
    )

    # Memory settings
    chunk_size: int = Field(default=1000, description="Chunk size for large datasets")
    low_memory_mode: bool = Field(default=False, description="Enable memory-efficient mode")

    # Species-specific settings
    default_species: str = Field(default="human", description="Default species")
    default_dataset_type: Literal[
        "unknown",
        "pbmc_or_blood",
        "normal_tissue",
        "tumor_tissue",
        "cell_line",
        "organoid",
        "spatial",
    ] = Field(default="unknown", description="Default dataset context for workflows")

    # Resource paths
    marker_db_path: Optional[Path] = Field(default=None, description="Path to marker database")
    gene_set_path: Optional[Path] = Field(default=None, description="Path to gene sets")

    @model_validator(mode="before")
    @classmethod
    def setup_and_validate(cls, data: Any) -> Any:
        """Set up logging and validate settings before initialization."""
        if not isinstance(data, dict):
            return data

        # Convert string paths to Path objects
        path_fields = ["log_file", "marker_db_path", "gene_set_path", "cache_dir"]
        for field in path_fields:
            if field in data and isinstance(data[field], str):
                data[field] = Path(data[field])

        return data

    @model_validator(mode="after")
    def post_init_setup(self) -> GlobalConfig:
        """Set up logging after initialization."""
        self._setup_logging()
        self._validate_settings()
        return self

    def _setup_logging(self):
        """Configure logging based on verbosity."""
        logger = logging.getLogger("sclucid")
        log_levels = {0: logging.WARNING, 1: logging.INFO, 2: logging.DEBUG}

        level = log_levels.get(self.verbosity, logging.INFO)

        # Clear existing handlers
        logger.handlers.clear()

        # Configure root logger
        logger.setLevel(level)

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        formatter = logging.Formatter(
            "[%(asctime)s] %(levelname)s - %(name)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

        # File handler (if specified)
        if self.log_file:
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)

    def _validate_settings(self):
        """Validate configuration settings."""
        if self.n_jobs < -1 or self.n_jobs == 0:
            warnings.warn(
                f"n_jobs={self.n_jobs} is invalid. Setting to -1 (use all cores).", UserWarning
            )
            # Can't modify self during validation, so this warning is informational
            # The user should set n_jobs correctly when creating the config

        if self.verbosity not in [0, 1, 2]:
            warnings.warn(
                f"verbosity={self.verbosity} is invalid. Setting to 1 (INFO).", UserWarning
            )

        if effective_n_jobs(self.n_jobs) != self.n_jobs:
            logging.getLogger(__name__).debug(
                "Runtime will execute n_jobs=%s as n_jobs=%s in the current environment.",
                self.n_jobs,
                effective_n_jobs(self.n_jobs),
            )

    def set(self, **kwargs):
        """Update configuration settings."""
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise ValueError(f"Unknown configuration key: {key}")
        self._validate_settings()

        if "verbosity" in kwargs or "log_file" in kwargs:
            self._setup_logging()

        if any(k in kwargs for k in ["figure_dpi", "plot_theme", "font_style"]):
            from .settings import set_figure_params

            set_figure_params(
                dpi=self.figure_dpi, color_theme=self.plot_theme, font_style=self.font_style
            )

    def reset(self):
        """Reset to default configuration."""
        # Create a new default instance and copy its values
        default_config = GlobalConfig()
        for field_name in default_config.model_fields:
            setattr(self, field_name, getattr(default_config, field_name))


# Global instance
_config = GlobalConfig()


def get_config() -> GlobalConfig:
    """Get global configuration instance."""
    return _config


def set_config(**kwargs):
    """Set global configuration parameters."""
    with _config_lock:
        _config.set(**kwargs)


@contextmanager
def config_context(**kwargs):
    """Temporary configuration context manager."""
    with _config_lock:
        old_config = {k: getattr(_config, k) for k in kwargs}
        _config.set(**kwargs)
    try:
        yield _config
    finally:
        with _config_lock:
            _config.set(**old_config)


def reset_config():
    """Reset global configuration to defaults."""
    _config.reset()


__all__ = [
    "GlobalConfig",
    "get_config",
    "set_config",
    "config_context",
    "reset_config",
]
