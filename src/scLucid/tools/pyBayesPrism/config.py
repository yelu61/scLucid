"""
Configuration classes for BayesPrism (R-free)

BayesPrism configuration for deconvolution parameters and Gibbs sampling control.
"""

from __future__ import annotations

import logging
from typing import Optional, Dict, Any

from pydantic import Field, field_validator, model_validator

from ...base_config import SclucidBaseConfig

logger = logging.getLogger(__name__)


class PrismConfig(SclucidBaseConfig):
    """
    BayesPrism configuration parameters

    Parameters
    ----------
    n_iter : int
        Number of Gibbs sampling iterations (default: 100)
    n_chains : int
        Number of Markov chains for MCMC (default: 4)
    burnin : int
        Number of burn-in iterations (default: 50)
    update_bulk : bool
        Whether to update bulk expression estimates (default: True)
    pseudo_min : float
        Minimum pseudo-count to avoid zeros (default: 1e-8)
    key : Optional[str]
        Keyword for tumor samples (default: None)
    outlier_cutoff : float
        Outlier detection cutoff (default: 0.01)
    max_outlier : int
        Maximum number of outlier genes (default: 10)
    gibbs_control : Dict[str, Any]
        Additional Gibbs sampler control parameters

    Examples
    --------
    >>> config = PrismConfig(n_iter=200, n_chains=4, burnin=100)
    >>> config.gibbs_control['thinning'] = 5
    """

    model_config = {"extra": "ignore"}

    n_iter: int = Field(default=100, gt=0, description="Number of Gibbs sampling iterations")
    n_chains: int = Field(default=4, gt=0, description="Number of Markov chains for MCMC")
    burnin: int = Field(default=50, ge=0, description="Number of burn-in iterations")
    update_bulk: bool = Field(default=True, description="Whether to update bulk expression estimates")
    pseudo_min: float = Field(default=1e-8, gt=0, description="Minimum pseudo-count to avoid zeros")
    key: Optional[str] = Field(default=None, description="Keyword for tumor samples")
    outlier_cutoff: float = Field(default=0.01, gt=0, description="Outlier detection cutoff")
    max_outlier: int = Field(default=10, gt=0, description="Maximum number of outlier genes")
    gibbs_control: Dict[str, Any] = Field(default_factory=dict, description="Additional Gibbs sampler control parameters")

    @model_validator(mode="before")
    @classmethod
    def validate_gibbs_control(cls, data: Any) -> Any:
        """Initialize gibbs_control with defaults if not provided"""
        if not isinstance(data, dict):
            return data

        # Only set defaults if gibbs_control is empty or not provided
        if not data.get("gibbs_control"):
            n_iter = data.get("n_iter", 1000)
            burnin = data.get("burnin", 500)
            data["gibbs_control"] = {
                'chain_length': n_iter,
                'burn_in': burnin,
                'thinning': 1,
                'verbose': False,
            }
        return data

    @model_validator(mode="after")
    def validate_parameters(self) -> "PrismConfig":
        """Validate configuration parameters"""
        if self.burnin >= self.n_iter:
            raise ValueError(f"burnin ({self.burnin}) must be < n_iter ({self.n_iter})")
        return self


class ReferenceConfig(SclucidBaseConfig):
    """
    Configuration for reference data processing

    Parameters
    ----------
    input_type : str
        Type of input data ("count.matrix" or "GEP")
    pseudo_min : float
        Minimum pseudo-count
    min_cells_per_type : int
        Minimum cells required per cell type
    min_genes_per_cell : int
        Minimum genes expressed per cell
    """

    model_config = {"extra": "ignore"}

    input_type: str = Field(default="count.matrix", description="Type of input data")
    pseudo_min: float = Field(default=1e-8, gt=0, description="Minimum pseudo-count")
    min_cells_per_type: int = Field(default=10, gt=0, description="Minimum cells required per cell type")
    min_genes_per_cell: int = Field(default=200, gt=0, description="Minimum genes expressed per cell")

    @field_validator("input_type")
    @classmethod
    def validate_input_type(cls, v: str) -> str:
        """Validate input_type"""
        valid_types = ["count.matrix", "GEP"]
        if v not in valid_types:
            raise ValueError(f"input_type must be one of {valid_types}, got '{v}'")
        return v


class DeconvolutionConfig(SclucidBaseConfig):
    """
    Configuration for deconvolution execution

    Parameters
    ----------
    n_cores : int
        Number of CPU cores for parallel processing
    verbose : bool
        Whether to display progress
    return_samples : bool
        Whether to return full posterior samples
    confidence_level : float
        Confidence level for credible intervals
    """

    model_config = {"extra": "ignore"}

    n_cores: int = Field(default=1, ge=1, description="Number of CPU cores for parallel processing")
    verbose: bool = Field(default=True, description="Whether to display progress")
    return_samples: bool = Field(default=False, description="Whether to return full posterior samples")
    confidence_level: float = Field(default=0.95, gt=0, lt=1, description="Confidence level for credible intervals")


__all__ = [
    "PrismConfig",
    "ReferenceConfig",
    "DeconvolutionConfig",
]
