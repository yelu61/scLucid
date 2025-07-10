"""Configuration module for single-cell RNA-seq data analysis"""

import os
from typing import Dict, Any, Optional

# Default parameters
DEFAULT_PARAMS = {
    # General parameters
    "batch_key": "sampleID",     # Batch/sample identifier
    "random_seed": 42,           # Random seed
    
    # Layer names
    "layer_raw": "counts",       # Raw count layer
    "layer_norm": "log1p_norm",  # Normalized layer
    "layer_scale": "scaled",     # Scaled layer
    
    # Normalization parameters
    "target_sum": 1e4,           # Normalization target sum
    
    # Highly variable genes parameters
    "n_top_genes": 2000,         # Number of HVGs
    "min_mean": 0.0125,          # Minimum mean expression
    "max_mean": 3,               # Maximum mean expression
    "min_disp": 0.5,             # Minimum dispersion
    
    # Dimensionality reduction parameters
    "n_pcs": 50,                 # Number of principal components
    "n_neighbors": 30,           # Number of neighbors for KNN
    
    # Clustering parameters
    "resolution": 0.8,           # Clustering resolution
    
    # Quality control parameters
    "min_genes": 300,            # Minimum number of genes
    "max_genes": 6000,           # Maximum number of genes
    "min_cells": 3,              # Minimum number of cells
    "max_mt_percent": 20,        # Maximum mitochondrial percentage
}

def get_param(name: str, user_params: Optional[Dict[str, Any]] = None) -> Any:
    """Get parameter value, prioritizing user-specified values
    
    Args:
        name: Parameter name to retrieve
        user_params: Dictionary of user-specified parameters
        
    Returns:
        The parameter value
        
    Raises:
        ValueError: If the parameter name is not found
    """
    if user_params and name in user_params:
        return user_params[name]
    if name in DEFAULT_PARAMS:
        return DEFAULT_PARAMS[name]
    raise ValueError(f"Unknown parameter: {name}")

def load_config(config_file: str) -> Dict[str, Any]:
    """Load parameters from a configuration file
    
    Args:
        config_file: Path to the JSON configuration file
        
    Returns:
        Dictionary of user parameters
        
    Raises:
        FileNotFoundError: If the configuration file does not exist
        JSONDecodeError: If the configuration file contains invalid JSON
    """
    import json
    
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    
    try:
        with open(config_file, 'r') as f:
            user_params = json.load(f)
    except json.JSONDecodeError as e:
        raise json.JSONDecodeError(f"Invalid JSON in configuration file: {e.msg}", e.doc, e.pos)
    
    return user_params