"""Global configuration for scLucid."""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Dict, Any, Literal
import warnings


@dataclass
class GlobalConfig:
    """Global configuration settings for scLucid."""
    
    # Computational settings
    n_jobs: int = -1
    random_state: int = 42
    backend: Literal['loky', 'threading', 'multiprocessing'] = 'loky'
    
    # Logging settings
    verbosity: int = 1  # 0: WARNING, 1: INFO, 2: DEBUG
    log_file: Optional[Path] = None
    
    # Cache settings
    cache_dir: Optional[Path] = None
    use_cache: bool = True
    
    # Plotting settings
    plot_backend: Literal['matplotlib', 'plotly'] = 'matplotlib'
    figure_dpi: int = 100
    figure_format: str = 'png'
    color_palette: str = 'tab20'
    plot_theme: str = 'default'
    
    # Memory settings
    chunk_size: int = 1000
    low_memory_mode: bool = False
    
    # Species-specific settings
    default_species: str = 'human'
    
    # Resource paths
    marker_db_path: Optional[Path] = None
    gene_set_path: Optional[Path] = None
    
    def __post_init__(self):
        """Set up logging and validate settings."""
        self._setup_logging()
        self._validate_settings()
    
    def _setup_logging(self):
        """Configure logging based on verbosity."""
        log_levels = {
            0: logging.WARNING,
            1: logging.INFO,
            2: logging.DEBUG
        }
        
        level = log_levels.get(self.verbosity, logging.INFO)
        
        # Configure root logger
        logger = logging.getLogger('sclucid')
        logger.setLevel(level)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s - %(name)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        
        # File handler (if specified)
        if self.log_file:
            self.log_file = Path(self.log_file)
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(level)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
    
    def _validate_settings(self):
        """Validate configuration settings."""
        if self.n_jobs < -1 or self.n_jobs == 0:
            warnings.warn(
                f"n_jobs={self.n_jobs} is invalid. Setting to -1 (use all cores).",
                UserWarning
            )
            self.n_jobs = -1
        
        if self.verbosity not in [0, 1, 2]:
            warnings.warn(
                f"verbosity={self.verbosity} is invalid. Setting to 1 (INFO).",
                UserWarning
            )
            self.verbosity = 1
    
    def set(self, **kwargs):
        """Update configuration settings."""
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
                if key == 'verbosity':
                    self._setup_logging()
            else:
                raise ValueError(f"Unknown configuration key: {key}")
    
    def reset(self):
        """Reset to default configuration."""
        self.__init__()


# Global instance
_config = GlobalConfig()


def get_config() -> GlobalConfig:
    """Get global configuration instance."""
    return _config


def set_config(**kwargs):
    """Set global configuration parameters."""
    _config.set(**kwargs)


def reset_config():
    """Reset global configuration to defaults."""
    _config.reset()