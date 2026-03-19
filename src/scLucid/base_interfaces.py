"""
Abstract base classes for scLucid analysis components.

This module provides abstract interfaces and base classes that define
the contract for various analysis components, ensuring consistency
and extensibility across the toolkit.

Key Abstract Classes:
- AnalysisStep: Base class for analysis steps
- QCFilter: Base class for QC filtering operations
- CellAnnotator: Base class for cell type annotation methods
- ScoringMethod: Base class for functional scoring methods
- PlottingBackend: Base class for visualization backends

Usage:
------
    from scLucid.base_interfaces import AnalysisStep, QCFilter

    class MyCustomAnalysis(AnalysisStep):
        def validate_input(self, adata):
            # Validate input data
            pass

        def run(self, adata, **kwargs):
            # Implement analysis
            pass
"""

from __future__ import annotations

import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Union

import matplotlib.pyplot as plt
from anndata import AnnData
import pandas as pd

from .base_config import SclucidBaseConfig

log = logging.getLogger(__name__)


class AnalysisStep(ABC):
    """
    Abstract base class for analysis steps.

    All analysis operations (QC, preprocessing, clustering, annotation, etc.)
    should inherit from this class to ensure consistent interface.

    Subclasses must implement:
    - validate_input(): Check if input data is valid
    - run(): Execute the analysis step
    - get_summary(): Return a summary of results

    Optional methods to override:
    - cleanup(): Clean up temporary resources
    """

    def __init__(self, config: Optional[SclucidBaseConfig] = None):
        """
        Initialize the analysis step.

        Parameters
        ----------
        config : SclucidBaseConfig, optional
            Configuration object for this step
        """
        self.config = config
        self._results: Optional[Dict[str, Any]] = None

    @abstractmethod
    def validate_input(self, adata: AnnData) -> bool:
        """
        Validate input data.

        Parameters
        ----------
        adata : AnnData
            Input data to validate

        Returns
        -------
        bool
            True if input is valid, False otherwise

        Raises
        ------
        ValueError
            If input is invalid and cannot be processed
        """
        pass

    @abstractmethod
    def run(self, adata: AnnData, **kwargs) -> AnnData:
        """
        Run the analysis step.

        Parameters
        ----------
        adata : AnnData
            Input data
        **kwargs
            Additional parameters

        Returns
        -------
        AnnData
            Processed data with results stored
        """
        pass

    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of the analysis results.

        Returns
        -------
        dict
            Dictionary containing summary statistics and metadata
        """
        if self._results is None:
            return {"status": "not_run"}
        return self._results

    def cleanup(self) -> None:
        """Clean up temporary resources."""
        pass

    def __repr__(self) -> str:
        config_info = f"config={self.config.__class__.__name__}" if self.config else "no config"
        return f"{self.__class__.__name__}({config_info})"


class QCFilter(ABC):
    """
    Abstract base class for QC filtering operations.

    QC filters should inherit from this class to ensure consistent
    interface for filtering cells or genes based on quality metrics.

    Subclasses must implement:
    - calculate_metric(): Compute the filtering metric
    - get_threshold(): Determine filtering threshold
    - apply_filter(): Apply the filter
    """

    @abstractmethod
    def calculate_metric(self, adata: AnnData) -> pd.Series:
        """
        Calculate the filtering metric for each cell/gene.

        Parameters
        ----------
        adata : AnnData
            Input data

        Returns
        -------
        pd.Series
            Series of metric values
        """
        pass

    @abstractmethod
    def get_threshold(self, metric: pd.Series) -> float:
        """
        Determine the filtering threshold.

        Parameters
        ----------
        metric : pd.Series
            Metric values

        Returns
        -------
        float
            Threshold value
        """
        pass

    @abstractmethod
    def apply_filter(self, adata: AnnData, threshold: float) -> AnnData:
        """
        Apply the filter to the data.

        Parameters
        ----------
        adata : AnnData
            Input data
        threshold : float
            Filtering threshold

        Returns
        -------
        AnnData
            Filtered data
        """
        pass


class CellAnnotator(ABC):
    """
    Abstract base class for cell type annotation methods.

    Cell annotation methods should inherit from this class to ensure
    consistent interface for assigning cell types.

    Subclasses must implement:
    - predict(): Predict cell types
    - get_confidence(): Get prediction confidence scores
    """

    @abstractmethod
    def predict(
        self,
        adata: AnnData,
        reference: Optional[AnnData] = None,
        **kwargs
    ) -> AnnData:
        """
        Predict cell types for the data.

        Parameters
        ----------
        adata : AnnData
            Data to annotate
        reference : AnnData, optional
            Reference data for supervised methods
        **kwargs
            Additional parameters

        Returns
        -------
        AnnData
            Data with cell type predictions stored in adata.obs
        """
        pass

    @abstractmethod
    def get_confidence(self, adata: AnnData) -> pd.Series:
        """
        Get confidence scores for predictions.

        Parameters
        ----------
        adata : AnnData
            Data with predictions

        Returns
        -------
        pd.Series
            Confidence scores
        """
        pass


class ScoringMethod(ABC):
    """
    Abstract base class for functional scoring methods.

    Scoring methods should inherit from this class to ensure
    consistent interface for scoring cells based on gene sets.

    Subclasses must implement:
    - score(): Calculate scores
    - normalize(): Normalize scores
    """

    @abstractmethod
    def score(
        self,
        adata: AnnData,
        gene_sets: Dict[str, List[str]],
        **kwargs
    ) -> AnnData:
        """
        Calculate functional scores for cells.

        Parameters
        ----------
        adata : AnnData
            Data to score
        gene_sets : dict
            Dictionary mapping gene set names to gene lists
        **kwargs
            Additional parameters

        Returns
        -------
        AnnData
            Data with scores stored in adata.obs
        """
        pass

    @abstractmethod
    def normalize(self, adata: AnnData, score_key: str) -> AnnData:
        """
        Normalize scores to a common scale.

        Parameters
        ----------
        adata : AnnData
            Data with scores
        score_key : str
            Key of score column in adata.obs

        Returns
        -------
        AnnData
            Data with normalized scores
        """
        pass


class PlottingBackend(ABC):
    """
    Abstract base class for plotting backends.

    Plotting backends should inherit from this class to ensure
    consistent interface for visualization.

    Subclasses must implement:
    - plot(): Create the plot
    - save(): Save the plot
    """

    @abstractmethod
    def plot(self, data: Any, **kwargs) -> plt.Figure:
        """
        Create a plot.

        Parameters
        ----------
        data : Any
            Data to plot
        **kwargs
            Plotting parameters

        Returns
        -------
        matplotlib.figure.Figure
            Created figure
        """
        pass

    @abstractmethod
    def save(self, fig: plt.Figure, path: str, **kwargs) -> None:
        """
        Save the plot to file.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure to save
        path : str
            Output path
        **kwargs
            Saving parameters (dpi, format, etc.)
        """
        pass


class ProportionMethod(ABC):
    """
    Abstract base class for cell proportion analysis methods.

    Proportion analysis methods should inherit from this class
    to ensure consistent interface.

    Subclasses must implement:
    - analyze(): Perform proportion analysis
    - test_significance(): Test statistical significance
    - plot_results(): Visualize results
    """

    @abstractmethod
    def analyze(
        self,
        adata: AnnData,
        groupby: str,
        **kwargs
    ) -> pd.DataFrame:
        """
        Analyze cell type proportions.

        Parameters
        ----------
        adata : AnnData
            Input data
        groupby : str
            Column to group by
        **kwargs
            Additional parameters

        Returns
        -------
        pd.DataFrame
            Proportion data
        """
        pass

    @abstractmethod
    def test_significance(
        self,
        proportions: pd.DataFrame,
        **kwargs
    ) -> pd.DataFrame:
        """
        Test statistical significance of proportion differences.

        Parameters
        ----------
        proportions : pd.DataFrame
            Proportion data
        **kwargs
            Test parameters

        Returns
        -------
        pd.DataFrame
            Statistical test results
        """
        pass

    @abstractmethod
    def plot_results(
        self,
        proportions: pd.DataFrame,
        **kwargs
    ) -> plt.Figure:
        """
        Visualize proportion analysis results.

        Parameters
        ----------
        proportions : pd.DataFrame
            Proportion data
        **kwargs
            Plotting parameters

        Returns
        -------
        matplotlib.figure.Figure
            Figure with plots
        """
        pass


# Factory function pattern
class AnalysisStepFactory:
    """
    Factory class for creating analysis steps.

    This factory allows dynamic registration and instantiation
    of analysis steps, enabling plugin-style extensibility.
    """

    _registered_steps: Dict[str, type[AnalysisStep]] = {}

    @classmethod
    def register(cls, name: str, step_class: type[AnalysisStep]) -> None:
        """
        Register a new analysis step.

        Parameters
        ----------
        name : str
            Name for the step
        step_class : type[AnalysisStep]
            Step class (must be subclass of AnalysisStep)
        """
        if not issubclass(step_class, AnalysisStep):
            raise TypeError(f"{step_class} must be a subclass of AnalysisStep")
        cls._registered_steps[name] = step_class
        log.info(f"Registered analysis step: {name}")

    @classmethod
    def create(cls, name: str, **kwargs) -> AnalysisStep:
        """
        Create an instance of a registered analysis step.

        Parameters
        ----------
        name : str
            Name of the step to create
        **kwargs
            Parameters to pass to the step constructor

        Returns
        -------
        AnalysisStep
            Instance of the requested step

        Raises
        ------
        KeyError
            If step name is not registered
        """
        if name not in cls._registered_steps:
            raise KeyError(f"Unknown analysis step: {name}")
        return cls._registered_steps[name](**kwargs)

    @classmethod
    def list_steps(cls) -> List[str]:
        """List all registered analysis steps."""
        return list(cls._registered_steps.keys())


__all__ = [
    # Abstract base classes
    "AnalysisStep",
    "QCFilter",
    "CellAnnotator",
    "ScoringMethod",
    "PlottingBackend",
    "ProportionMethod",
    # Factory
    "AnalysisStepFactory",
]
