"""Configuration module for single-cell RNA-seq data analysis"""

import logging
from typing import Dict, Any, Optional

log = logging.getLogger(__name__)

class Config:
    """
    Configuration class for the scRNA pipeline.
    Stores default parameters and allows for updates from user-provided dicts.
    """
    def __init__(self):
        # General parameters
        self.batch_key: str = "sampleID"     # Batch/sample identifier
        self.random_seed: int = 42          # Random seed
        
        # Layers
        self.layer_raw: str = "counts"       # Raw count layer
        self.layer_norm: str = "log1p_norm"  # Normalized layer
        self.layer_scale: str = "scaled"     # Scaled layer
        
        # QC
        self.min_genes: int = 300
        self.max_genes: int = 8000
        self.min_cells: int = 3
        self.max_mt_percent: float = 20.0
        
        # Preprocessing
        self.target_sum: float = 1e4           # Normalization target sum
        self.n_top_genes: int = 2000         # Number of HVGs
        
        # Analysis
        self.n_pcs: int = 50                 # Number of principal components
        self.n_neighbors: int = 30           # Number of neighbors for KNN
        self.resolution: float = 0.8
        
    def update(self, user_params: Dict[str, Any]):
        """Update configuration with user-provided parameters."""
        for key, value in user_params.items():
            if hasattr(self, key):
                setattr(self, key, value)
                log.info(f"Configuration updated: {key} = {value}")
            else:
                log.warning(f"Unknown configuration parameter: {key}")

settings = Config()