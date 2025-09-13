"""
Configuration classes for the Quality Control (QC) module.

This file defines dataclass structures used to configure all aspects
of the QC pipeline, from threshold setting to doublet detection and filtering.
"""

import logging
from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Dict, List, Literal, Optional, Tuple, Union

import pandas as pd
from anndata import AnnData

__all__ = [
    "MetricsReportingConfig",
    "QCThresholds",
    "MarkerConfig",
    "DoubletConfig",
    "MarkingConfig",
    "FilterConfig",
    "QCWorkflowConfig",
]


@dataclass
class MetricsReportingConfig:
    """Configuration for reporting and visualization in the metrics calculation step."""

    include_standard_qc: bool = True

    # Plotting Controls
    plot_violin: bool = True  #: Whether to plot violin plots for QC metrics.
    plot_scatter: bool = (
        True  #: Whether to plot scatter plot for total_counts vs n_genes_by_counts.
    )
    plot_top_genes: bool = (
        True  #: Whether to plot distribution of pct_counts_in_top_X_genes.
    )
    show_plots: bool = True  #: Whether to display the plots interactively.

    # File Export Controls
    save_dir: Optional[str] = (
        None  #: Directory to save plots and statistics. If None, nothing is saved.
    )
    export_stats: bool = (
        True  #: Whether to export per-sample/global summary tables as CSV.
    )
    export_xlsx: bool = False  #: If export_stats is True, also export as an Excel file.

    # Logging Control
    print_stats: bool = True  #: Whether to print summary stats to the log.

    def validate(self):
        if not isinstance(self.include_standard_qc, bool):
            raise ValueError("include_standard_qc must be bool")
        for attr in [
            "plot_violin",
            "plot_scatter",
            "plot_top_genes",
            "show_plots",
            "export_stats",
            "export_xlsx",
            "print_stats",
        ]:
            if not isinstance(getattr(self, attr), bool):
                raise ValueError(f"{attr} must be bool")
        if self.save_dir is not None and not isinstance(self.save_dir, str):
            raise ValueError("save_dir must be None or str")


@dataclass
class QCThresholds:
    """Configuration for QC filtering thresholds."""

    # Gene-based filtering
    min_genes: Optional[int] = 200
    max_genes: Optional[int] = None
    # Cell-based filtering
    min_counts: Optional[int] = None
    max_counts: Optional[int] = None
    # MAD-based outlier detection
    nmads: float = 5.0  # Number of MADs for outlier detection
    # Percentage-based filtering
    pc_mt: Optional[float] = 20.0  # Mitochondrial percentage threshold
    pc_hb: Optional[float] = 20.0  # Hemoglobin percentage threshold
    pc_top_genes: Dict[str, float] = field(
        default_factory=dict
    )  # e.g. {'pc_top_20_genes': 50.0}
    use_fixed_top_gene_threshold: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert thresholds to dictionary format."""
        return self.__dict__

    def validate(self):
        """Validate threshold parameters for logical consistency."""
        if self.min_genes is not None and self.max_genes is not None:
            if self.min_genes > self.max_genes:
                raise ValueError("min_genes cannot be greater than max_genes")
        if self.min_counts is not None and self.max_counts is not None:
            if self.min_counts > self.max_counts:
                raise ValueError("min_counts cannot be greater than max_counts")
        if self.nmads is not None and self.nmads <= 0:
            raise ValueError("nmads must be positive")
        if self.pc_mt is not None and not (0 <= self.pc_mt <= 100):
            raise ValueError("pc_mt should be between 0 and 100")
        if self.pc_hb is not None and not (0 <= self.pc_hb <= 100):
            raise ValueError("pc_hb should be between 0 and 100")
        if self.pc_top_genes:
            for k, v in self.pc_top_genes.items():
                if v is not None and not (0 <= v <= 100):
                    raise ValueError(f"{k} in pc_top_genes should be between 0 and 100")
        # Add more as needed


@dataclass
class MarkerConfig:
    """Configuration for a single marker set in heuristic doublet detection."""

    genes: Union[List[str], str]  # List of genes or regex pattern
    expression_threshold: float = 2.0  # Minimum expression to be considered "expressed"
    min_genes_required: int = 2
    use_raw: bool = True  # Whether to use raw counts
    is_regex: bool = field(init=False)  # Determined post-init

    def __post_init__(self):
        self.is_regex = isinstance(self.genes, str)


@dataclass
class DoubletConfig:
    """
    High-level configuration for the main doublet detection workflow.
    """
    # --- Core Algorithm Parameters ---
    method: Literal["scrublet",] = "scrublet" #: The algorithm to use for doublet score calculation.
    n_pcs: int = 30 #: Number of principal components to use for the algorithm.
    plot_umap: bool = True #: Whether to plot the UMAP of the data.
    expected_doublet_rate: Optional[Union[float, Dict[str, float]]] = None #: Expected doublet rate. Can be a single float or a dict per sample.
    run_algorithm: bool = True
    """If False, skips the algorithmic detection step entirely. Useful for heuristic-only exploratory runs."""
    
    # --- Heuristic Analysis ---
    use_heuristics: bool = True #: Master switch to enable/disable marker-based heuristic analysis.
    """If True, enables the heuristic (marker-based) doublet detection workflow."""
    
    # Parameters for automatic marker loading (if use_heuristics is True)
    marker_species: str = "human" #: Species to load markers for ('human', 'mouse', etc.).
    marker_tissue: Optional[str] = None #: Tissue context for loading specific markers (e.g., 'Lung').
    # Manual override for markers. If provided, this dictionary will be used instead of automatic loading.
    marker_configs: Optional[Dict[str, MarkerConfig]] = None
    # Default evaluation parameters for heuristic analysis.
    default_expression_threshold: float = 2.0 #: Default expression value to consider a marker "on".
    default_min_genes_required: int = 2 #: Default number of markers required to define a lineage.
    default_use_raw: bool = True #: Default setting to use adata.raw for heuristics.
    min_lineage_prevalence: float = 0.005
    min_lineages_for_doublet: int = 2 #: How many distinct lineages must be co-expressed to be a doublet.
    ignore_coexpression_pairs: Optional[List[Tuple[str, str]]] = None #: Specific lineage pairs to ignore (e.g., [('T_cells', 'NK_cells')]).

    # --- Result Merging and Reporting ---
    merge_strategy: Literal['weighted_average', 'max_score', 'heuristic_boost'] = "weighted_average"
    algorithm_weight: float = 0.7
    random_state: int = 61 #: Random seed for reproducibility.
    plot_summary: bool = True #: Master switch to generate a summary plot at the end of the run.
    plot_bar: bool = True #: (If plot_summary=True) Include the bar plot.
    plot_scatter: bool = True #: (If plot_summary=True) Include the scatter plots.
    plot_upset: bool = True #: (If plot_summary=True) Include the UpSet plot.
    export_stats: bool = True #: Whether to export summary statistics to CSV files.
    save_dir: Optional[str] = None #: Directory to save plots and statistics.
    show_plots: bool = True #: Whether to display plots interactively.

    def validate(self):
        """Validate configuration parameters."""
        allowed_methods = ["scrublet"]
        if self.method not in allowed_methods:
            raise ValueError(f"method must be one of {allowed_methods}")
        if self.n_pcs <= 1:
            raise ValueError("n_pcs must be greater than 1")
        if (
            self.use_heuristics
            and (self.min_lineages_for_doublet is not None)
            and self.min_lineages_for_doublet < 2
        ):
            raise ValueError(
                "min_lineages_for_doublet must be at least 2 when using heuristics"
            )
        if isinstance(self.expected_doublet_rate, float):
            if not (0.0 < self.expected_doublet_rate < 1.0):
                raise ValueError("expected_doublet_rate must be between 0 and 1")
        if isinstance(self.expected_doublet_rate, dict):
            for k, v in self.expected_doublet_rate.items():
                if not (0.0 < v < 1.0):
                    raise ValueError(
                        f"expected_doublet_rate for sample '{k}' must be between 0 and 1"
                    )
        allowed_methods = ["scrublet", ]
        if self.method not in allowed_methods:
            raise ValueError(f"method must be one of {allowed_methods}")
        allowed_merge_strategies = [
            'weighted_average', 'max_score', 'heuristic_boost'
        ]
        if self.merge_strategy not in allowed_merge_strategies:
            raise ValueError(f"merge_strategy must be one of {allowed_merge_strategies}")
        
        if self.merge_strategy == "weighted_average" and not (0 <= self.algorithm_weight <= 1):
            raise ValueError("algorithm_weight must be between 0 and 1 for 'weighted_average' strategy.")

        if self.use_heuristics and self.min_lineages_for_doublet < 2:
            raise ValueError("min_lineages_for_doublet must be at least 2 when using heuristics.")
        
        if not self.use_heuristics and (
            self.marker_configs is not None or self.marker_tissue is not None
        ):
            logging.warning(
                "Heuristic-specific parameters (marker_configs, marker_tissue) are set, but use_heuristics is False. These will be ignored."
            )
        

@dataclass
class MarkingConfig:
    """Configuration for the cell marking workflow."""

    thresholds: QCThresholds = field(default_factory=QCThresholds)
    qc_metrics_mad: List[Tuple[str, str]] = field(
        default_factory=lambda: [
            ("log1p_total_counts", "both"),
            ("log1p_n_genes_by_counts", "both"),
            ("pct_counts_in_top_20_genes", "upper"),
        ]
    )
    custom_outlier_functions: Optional[Dict[str, Callable[[AnnData], pd.Series]]] = None
    plot_outliers: bool = True
    save_dir: Optional[str] = None
    show_plots: bool = True
    cols_to_plot: Optional[List[str]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert the configuration to a dictionary."""
        # Custom handling for non-serializable functions
        d = asdict(self)
        if self.custom_outlier_functions:
            d["custom_outlier_functions"] = {
                name: func.__name__
                for name, func in self.custom_outlier_functions.items()
            }
        return d

    def validate(self):
        if hasattr(self.thresholds, "validate"):
            self.thresholds.validate()
        if self.qc_metrics_mad:
            for item in self.qc_metrics_mad:
                if not (isinstance(item, tuple) and len(item) == 2):
                    raise ValueError(
                        "qc_metrics_mad should be a list of (metric, direction) tuples"
                    )
                if item[1] not in ("both", "upper", "lower"):
                    raise ValueError(f"Invalid direction {item[1]} in qc_metrics_mad")
        if self.custom_outlier_functions:
            for name, func in self.custom_outlier_functions.items():
                if not callable(func):
                    raise ValueError(
                        f"Custom outlier function '{name}' is not callable"
                    )


@dataclass
class FilterConfig:
    """Configuration for cell filtering logic."""

    criteria_to_filter: List[str] = field(
        default_factory=lambda: ["outlier_min_genes", "outlier_mt", "predicted_doublet"]
    )
    combination_logic: Literal["any", "all", "custom", "threshold"] = "any"
    custom_logic_expr: Optional[str] = None
    min_criteria_for_removal: int = 2  # Used when logic is 'threshold'
    metadata_filters: Optional[Dict[str, Any]] = None  # Additional sample/group filters

    def to_dict(self) -> Dict[str, Any]:
        """Convert the configuration to a dictionary."""
        return asdict(self)

    def validate(self):
        if self.combination_logic == "custom" and not self.custom_logic_expr:
            raise ValueError(
                "combination_logic is 'custom', but 'custom_logic_expr' is not provided."
            )
        if self.combination_logic == "threshold" and (
            self.min_criteria_for_removal is None or self.min_criteria_for_removal < 1
        ):
            raise ValueError(
                "'min_criteria_for_removal' must be >= 1 for 'threshold' logic."
            )
        if not self.criteria_to_filter or not isinstance(self.criteria_to_filter, list):
            raise ValueError("criteria_to_filter must be a non-empty list.")


@dataclass
class QCWorkflowConfig:
    """A single configuration object for the entire QC workflow."""

    # Sub-configurations for each step
    metrics_reporting_config: MetricsReportingConfig = field(
        default_factory=MetricsReportingConfig
    )
    marking_config: MarkingConfig = field(default_factory=MarkingConfig)
    doublet_config: DoubletConfig = field(default_factory=DoubletConfig)
    filter_config: FilterConfig = field(default_factory=FilterConfig)

    # Global parameters for the workflow
    sample_key: str = "sampleID"
    species: str = "human"
    tissue: Optional[str] = None
    results_dir: str = "./qc_results"

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the entire workflow configuration to a dictionary, handling custom sub-configs.
        """
        data = asdict(self)

        if hasattr(self.marking_config, "to_dict") and callable(
            self.marking_config.to_dict
        ):
            data["marking_config"] = self.marking_config.to_dict()

        return data

    def validate(self):
        self.metrics_reporting_config.validate()
        self.marking_config.validate()
        self.doublet_config.validate()
        self.filter_config.validate()
        if not isinstance(self.sample_key, str):
            raise ValueError("sample_key must be a string")
        if not isinstance(self.results_dir, str):
            raise ValueError("results_dir must be a string")
        if not isinstance(self.species, str):
            raise ValueError("species must be a string")
