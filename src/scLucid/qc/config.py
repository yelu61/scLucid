"""
Pydantic-based configuration classes for the Quality Control (QC) module.

This is v2 of the QC configuration system, migrating from dataclasses to Pydantic
for consistent validation, serialization, and documentation.
"""

from __future__ import annotations

import logging
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

from pydantic import ConfigDict, Field, field_validator, model_validator

from ..base_config import SclucidBaseConfig, WorkflowConfigBase

logger = logging.getLogger(__name__)


class MetricsReportingConfig(SclucidBaseConfig):
    """Configuration for reporting and visualization in the metrics calculation step."""

    model_config = ConfigDict(extra="ignore")

    include_standard_qc: bool = Field(default=True, description="Include standard QC metrics.")

    # Plotting Controls
    plot_violin: bool = Field(default=True, description="Plot violin plots for QC metrics.")
    plot_scatter: bool = Field(
        default=True, description="Plot scatter for total_counts vs n_genes."
    )
    plot_top_genes: bool = Field(
        default=True, description="Plot distribution of pct_counts_in_top_X_genes."
    )
    show_plots: bool = Field(default=True, description="Display plots interactively.")

    # File Export Controls
    export_stats: bool = Field(default=True, description="Export summary tables as CSV.")
    export_xlsx: bool = Field(default=False, description="Also export as Excel file.")

    # Logging Control
    print_stats: bool = Field(default=True, description="Print summary stats to log.")


class QCThresholds(SclucidBaseConfig):
    """Configuration for QC filtering thresholds."""

    model_config = ConfigDict(extra="ignore")

    # Gene-based filtering
    min_genes: Optional[int] = Field(default=200, ge=0, description="Minimum genes per cell.")
    max_genes: Optional[int] = Field(default=None, ge=0, description="Maximum genes per cell.")

    # Cell-based filtering
    min_counts: Optional[int] = Field(
        default=None, ge=0, description="Minimum UMI counts per cell."
    )
    max_counts: Optional[int] = Field(
        default=None, ge=0, description="Maximum UMI counts per cell."
    )

    # MAD-based outlier detection
    nmads: float = Field(default=5.0, gt=0, description="Number of MADs for outlier detection.")

    # Percentage-based filtering
    pc_mt: Optional[float] = Field(
        default=20.0, ge=0, le=100, description="Mitochondrial percentage threshold."
    )
    pc_hb: Optional[float] = Field(
        default=20.0, ge=0, le=100, description="Hemoglobin percentage threshold."
    )
    pc_top_genes: Dict[str, float] = Field(
        default_factory=dict, description="Top gene percentages, e.g., {'pc_top_20_genes': 50.0}"
    )
    use_fixed_top_gene_threshold: bool = Field(default=False)

    @model_validator(mode="after")
    def validate_thresholds(self) -> QCThresholds:
        """Validate threshold logical consistency."""
        if self.min_genes is not None and self.max_genes is not None:
            if self.min_genes > self.max_genes:
                raise ValueError("min_genes cannot be greater than max_genes")
        if self.min_counts is not None and self.max_counts is not None:
            if self.min_counts > self.max_counts:
                raise ValueError("min_counts cannot be greater than max_counts")

        # Validate pc_top_genes values
        for k, v in self.pc_top_genes.items():
            if not (0 <= v <= 100):
                raise ValueError(f"{k} in pc_top_genes should be between 0 and 100")

        return self

    def to_dict(self) -> Dict[str, Any]:
        """Convert thresholds to dictionary format."""
        return self.model_dump()


class MarkerConfig(SclucidBaseConfig):
    """Configuration for a single marker set in heuristic doublet detection."""

    model_config = ConfigDict(extra="ignore")

    genes: Union[List[str], str] = Field(description="List of genes or regex pattern")
    expression_threshold: float = Field(
        default=2.0, description="Minimum expression to be 'expressed'"
    )
    min_genes_required: int = Field(
        default=2, ge=1, description="Minimum genes required for lineage"
    )
    use_raw: bool = Field(default=True, description="Whether to use raw counts")
    is_regex: bool = Field(default=False, description="Determined from genes type")

    @model_validator(mode="before")
    @classmethod
    def set_is_regex(cls, data: Any) -> Any:
        """Set is_regex based on genes type before instantiation."""
        if isinstance(data, dict) and "genes" in data and "is_regex" not in data:
            data["is_regex"] = isinstance(data["genes"], str)
        return data


class DoubletConfig(SclucidBaseConfig):
    """High-level configuration for the main doublet detection workflow."""

    model_config = ConfigDict(extra="ignore")

    # --- Core Algorithm Parameters ---
    run_algorithm: bool = Field(
        default=True,
        description="If False, skips algorithmic detection. Useful for heuristic-only runs.",
    )
    method: Literal["scrublet", "solo", "doubletdetection"] = Field(
        default="scrublet", description="Algorithm for doublet score calculation"
    )

    # --- Scrublet Algorithm Specific Parameters ---
    scr_n_pcs: int = Field(default=30, gt=1, description="Number of PCs for scrublet")
    scr_plot_umap: bool = Field(default=False, description="Whether to plot UMAP")

    # --- SOLO Algorithm Specific Parameters ---
    solo_n_epochs: int = Field(default=400, ge=100, description="Training epochs for SOLO")
    solo_learning_rate: float = Field(default=1e-3, gt=0)
    solo_use_gpu: bool = Field(default=True)
    solo_use_raw: bool = Field(default=True)
    solo_clear_cache: bool = Field(default=True)

    # --- DoubletDetection Algorithm Specific Parameters ---
    dd_n_components: int = Field(default=30, gt=1, description="Number of PCs for DoubletDetection")
    dd_n_top_var_genes: int = Field(default=10000, ge=100)
    dd_p_thresh: float = Field(default=1e-16, gt=0, lt=1)
    dd_voter_thresh: float = Field(default=0.8, ge=0, le=1)
    dd_use_raw: bool = Field(default=True)

    # --- Heuristic Analysis ---
    use_heuristics: bool = Field(
        default=True, description="Master switch for marker-based heuristic analysis"
    )
    marker_species: str = Field(default="human", description="Species for marker loading")
    marker_tissue: Optional[str] = Field(default=None, description="Tissue context for markers")
    marker_configs: Optional[Dict[str, MarkerConfig]] = Field(
        default=None, description="Manual marker override (disables auto-loading)"
    )
    default_expression_threshold: float = Field(default=2.0)
    default_min_genes_required: int = Field(default=2, ge=1)
    default_use_raw: bool = Field(default=True)
    min_lineage_prevalence: float = Field(default=0.005, ge=0, le=1)
    min_lineages_for_doublet: int = Field(
        default=2, ge=2, description="Min lineages for doublet call"
    )
    ignore_coexpression_pairs: List[Tuple[str, str]] = Field(
        default_factory=lambda: [("Epithelial", "Mesenchymal")]
    )

    # --- Result Merging and Reporting ---
    merge_strategy: Literal["weighted_average", "max_score", "heuristic_boost"] = Field(
        default="weighted_average"
    )
    algorithm_weight: float = Field(default=0.7, ge=0, le=1)
    expected_doublet_rate: Optional[Union[float, Dict[str, float]]] = Field(default=None)
    score_threshold: Optional[float] = Field(
        default=None,
        description="Optional direct doublet-score threshold. When set, overrides expected_doublet_rate quantile-based threshold.",
    )
    random_state: int = Field(default=61)
    plot_summary: bool = Field(default=True)
    plot_bar: bool = Field(default=True)
    plot_scatter: bool = Field(default=True)
    plot_upset: bool = Field(default=True)
    export_stats: bool = Field(default=True)
    show_plots: bool = Field(default=True)

    @model_validator(mode="after")
    def validate_method_params(self) -> DoubletConfig:
        """Validate method-specific dependencies and parameters."""
        if self.method == "scrublet":
            try:
                import scrublet
            except ImportError:
                logger.warning("scrublet not found. Install with: pip install scrublet")

        elif self.method == "solo":
            try:
                import scvi
            except ImportError:
                logger.warning("scvi-tools not found. Install with: pip install scvi-tools")
            if self.solo_n_epochs < 100:
                logger.warning("solo_n_epochs < 100 may result in poor convergence")

        elif self.method == "doubletdetection":
            try:
                import doubletdetection
            except ImportError:
                logger.warning(
                    "doublet-detection not found. Install with: pip install doublet-detection"
                )

        return self

    @model_validator(mode="after")
    def validate_heuristics(self) -> DoubletConfig:
        """Validate heuristic-specific parameters."""
        if not self.use_heuristics:
            if self.marker_configs is not None or self.marker_tissue is not None:
                logger.warning(
                    "Heuristic-specific parameters are set, but use_heuristics is False. These will be ignored."
                )
        return self

    @field_validator("expected_doublet_rate")
    @classmethod
    def validate_doublet_rate(
        cls, v: Optional[Union[float, Dict[str, float]]]
    ) -> Optional[Union[float, Dict[str, float]]]:
        """Validate expected doublet rate."""
        if v is None:
            return v
        if isinstance(v, float):
            if not (0.0 < v < 1.0):
                raise ValueError("expected_doublet_rate must be between 0 and 1")
        elif isinstance(v, dict):
            for k, rate in v.items():
                if not (0.0 < rate < 1.0):
                    raise ValueError(
                        f"expected_doublet_rate for sample '{k}' must be between 0 and 1"
                    )
        return v


class MarkingConfig(SclucidBaseConfig):
    """Configuration for the cell marking workflow."""

    model_config = ConfigDict(extra="ignore")

    thresholds: QCThresholds = Field(default_factory=QCThresholds)
    qc_metrics_mad: List[Tuple[str, str]] = Field(
        default_factory=lambda: [
            ("log1p_total_counts", "both"),
            ("log1p_n_genes_by_counts", "both"),
            ("pct_counts_in_top_20_genes", "upper"),
        ]
    )
    custom_outlier_functions: Optional[Dict[str, Callable]] = Field(default=None)
    plot_outliers: bool = Field(default=True)
    cols_to_plot: Optional[List[str]] = Field(default=None)
    show_plots: bool = Field(default=True, description="Display plots interactively.")

    @field_validator("qc_metrics_mad")
    @classmethod
    def validate_qc_metrics(cls, v: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
        """Validate qc_metrics_mad format."""
        for item in v:
            if not (isinstance(item, tuple) and len(item) == 2):
                raise ValueError("qc_metrics_mad should be a list of (metric, direction) tuples")
            if item[1] not in ("both", "upper", "lower"):
                raise ValueError(f"Invalid direction {item[1]} in qc_metrics_mad")
        return v

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        d = self.model_dump()
        # Handle non-serializable functions
        if self.custom_outlier_functions:
            d["custom_outlier_functions"] = {
                name: func.__name__ for name, func in self.custom_outlier_functions.items()
            }
        return d


class FilterConfig(SclucidBaseConfig):
    """Configuration for cell filtering logic."""

    model_config = ConfigDict(extra="ignore")

    criteria_to_filter: List[str] = Field(
        default_factory=lambda: ["outlier_min_genes", "outlier_mt", "predicted_doublet"]
    )
    combination_logic: Literal["any", "all", "custom", "threshold"] = Field(default="any")
    custom_logic_expr: Optional[str] = Field(default=None)
    min_criteria_for_removal: int = Field(
        default=2, ge=1, description="Used when logic is 'threshold'"
    )
    metadata_filters: Optional[Dict[str, Any]] = Field(default=None)

    @model_validator(mode="after")
    def validate_filter_logic(self) -> FilterConfig:
        """Validate filter logic configuration."""
        if self.combination_logic == "custom" and not self.custom_logic_expr:
            raise ValueError(
                "combination_logic is 'custom', but 'custom_logic_expr' is not provided"
            )

        if self.combination_logic == "threshold" and self.min_criteria_for_removal < 1:
            raise ValueError("'min_criteria_for_removal' must be >= 1 for 'threshold' logic")

        if not self.criteria_to_filter:
            raise ValueError("criteria_to_filter must be a non-empty list")

        return self


class QCWorkflowConfig(WorkflowConfigBase):
    """A single configuration object for the entire QC workflow.

    Default-path semantics:
        - ``use_recommendations=True``: intelligent QC recommendations are applied
          only to fields the caller did not explicitly set.
        - ``threshold_mode="hierarchical"``: per-sample adaptive thresholds are
          computed when multiple samples are present; single-sample data falls
          back to pooled behavior automatically.
        - Tumor-aware adjustment is triggered when ``tissue_type`` contains
          "tumor" or "cancer" (e.g. ``outlier_mt`` is excluded from filtering).
        - When ``save_dir`` is set, the workflow writes a reviewer-facing summary
          (``qc_review_summary.json`` / ``qc_review_summary.md``) alongside the
          standard report.
    """

    model_config = ConfigDict(extra="ignore")

    # Sub-configurations for each step
    metrics_reporting_config: MetricsReportingConfig = Field(default_factory=MetricsReportingConfig)
    marking_config: MarkingConfig = Field(default_factory=MarkingConfig)
    doublet_config: DoubletConfig = Field(default_factory=DoubletConfig)
    filter_config: FilterConfig = Field(default_factory=FilterConfig)

    # Global parameters for the workflow
    sample_key: str = Field(default="sampleID", description="Column for sample identification")
    species: str = Field(default="human", description="Species for analysis")
    tissue: Optional[str] = Field(default=None, description="Tissue context for markers")
    tissue_type: Optional[str] = Field(
        default=None, description="Tissue type hint passed to recommender (e.g., 'lung_tumor')."
    )
    threshold_mode: Literal["hierarchical", "pooled", "independent"] = Field(
        default="hierarchical", description="Multi-sample threshold policy."
    )
    use_recommendations: bool = Field(
        default=True, description="Run intelligent QC recommendation and apply to thresholds."
    )
    # Note: save_dir is inherited from SclucidBaseConfig, used for results output

    # Performance parameters
    use_parallel: bool = Field(default=True, description="Use parallel processing for samples")

    # Backward compatibility alias
    @property
    def results_dir(self) -> str:
        """Backward compatibility alias for save_dir."""
        return self.save_dir if self.save_dir else "./qc_results"

    @results_dir.setter
    def results_dir(self, value: str):
        """Backward compatibility setter for save_dir."""
        self.save_dir = value

    def to_dict(self) -> Dict[str, Any]:
        """Convert workflow configuration to dictionary."""
        data = self.model_dump()
        # Handle nested configs with custom to_dict
        if hasattr(self.marking_config, "to_dict"):
            data["marking_config"] = self.marking_config.to_dict()
        return data

    @classmethod
    def from_simple_dict(cls, simple_config: Dict[str, Any]) -> QCWorkflowConfig:
        """
        Create QCWorkflowConfig from a simplified flat dictionary.

        This factory method allows users to create complex nested configurations
        from a simple flat dictionary, reducing the need to understand internal
        config structure.

        Args:
            simple_config: Flat dictionary with keys like:
                - thresholds_min_genes, thresholds_max_genes, thresholds_pc_mt
                - doublet_method, doublet_expected_rate
                - sample_key, species, tissue, results_dir
                - n_jobs, use_parallel

        Returns:
            QCWorkflowConfig: Fully configured workflow config

        Example:
            >>> config = QCWorkflowConfig.from_simple_dict({
            ...     "thresholds_min_genes": 200,
            ...     "thresholds_pc_mt": 20.0,
            ...     "doublet_method": "scrublet",
            ...     "species": "human",
            ...     "results_dir": "./qc_results"
            ... })
        """
        # Extract threshold parameters
        threshold_keys = [
            "min_genes",
            "max_genes",
            "min_counts",
            "max_counts",
            "pc_mt",
            "pc_hb",
            "nmads",
        ]
        thresholds = {}
        for key in threshold_keys:
            config_key = f"thresholds_{key}"
            if config_key in simple_config:
                thresholds[key] = simple_config.pop(config_key)

        # Extract doublet parameters
        doublet_keys = ["method", "expected_doublet_rate", "scr_n_pcs", "solo_n_epochs"]
        doublet_config = {}
        for key in doublet_keys:
            config_key = f"doublet_{key}"
            if config_key in simple_config:
                doublet_config[key] = simple_config.pop(config_key)

        # Build nested configs
        kwargs = simple_config.copy()

        # Handle backward compatibility: results_dir -> save_dir
        if "results_dir" in kwargs:
            kwargs["save_dir"] = kwargs.pop("results_dir")

        if thresholds:
            kwargs["marking_config"] = MarkingConfig(thresholds=QCThresholds(**thresholds))

        if doublet_config:
            kwargs["doublet_config"] = DoubletConfig(**doublet_config)

        return cls(**kwargs)

    @classmethod
    def quick(cls, min_genes: int = 200, pc_mt: float = 20.0, **kwargs) -> QCWorkflowConfig:
        """
        Quick configuration factory for common use cases.

        Args:
            min_genes: Minimum genes per cell threshold
            pc_mt: Mitochondrial percentage threshold
            **kwargs: Additional parameters (species, sample_key, etc.)

        Returns:
            QCWorkflowConfig: Pre-configured for standard analysis

        Example:
            >>> config = QCWorkflowConfig.quick(min_genes=300, pc_mt=15.0, species="mouse")
        """
        return cls(
            marking_config=MarkingConfig(thresholds=QCThresholds(min_genes=min_genes, pc_mt=pc_mt)),
            **kwargs,
        )


# Export all config classes
__all__ = [
    "MetricsReportingConfig",
    "QCThresholds",
    "MarkerConfig",
    "DoubletConfig",
    "MarkingConfig",
    "FilterConfig",
    "QCWorkflowConfig",
]
