"""Pydantic-based configuration classes for the tumor analysis workflow."""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Literal, Optional

from pydantic import ConfigDict, Field, field_validator

from ..base_config import SclucidBaseConfig, WorkflowConfigBase

logger = logging.getLogger(__name__)


class TumorAnalysisConfig(SclucidBaseConfig):
    """Configuration for tumor-specific analysis steps."""

    model_config = ConfigDict(extra="ignore")

    # Malignancy analysis
    run_malignancy: bool = Field(
        default=True, description="Run malignancy scoring and classification."
    )
    malignancy_method: Literal["cnv", "threshold", "ml"] = Field(
        default="cnv", description="Method for malignant cell classification."
    )
    malignancy_reference_key: Optional[str] = Field(
        default=None,
        description="Observation key identifying reference normal cells for threshold/ml methods.",
    )

    # TME analysis
    run_tme: bool = Field(
        default=True, description="Run tumor microenvironment deconvolution."
    )
    tme_cell_type_key: str = Field(
        default="cell_type_auto",
        description="Column in adata.obs containing cell type annotations for TME.",
    )

    # CNV analysis
    run_cnv: bool = Field(
        default=False,
        description="Run CNV inference. Default False because infercnvpy is optional.",
    )
    cnv_reference_key: Optional[str] = Field(
        default=None,
        description="Observation key identifying reference cells for CNV inference.",
    )

    # Therapy analysis
    run_therapy: bool = Field(
        default=False, description="Run therapy response prediction."
    )
    therapy_drugs: Optional[List[str]] = Field(
        default=None, description="List of drug names for resistance scoring."
    )

    @field_validator("malignancy_method")
    @classmethod
    def validate_malignancy_method(cls, v: str) -> str:
        """Warn about methods that require references."""
        if v in ("threshold", "ml"):
            logger.warning(
                f"malignancy_method='{v}' may require reference normal cells. "
                "Ensure malignancy_reference_key is set if needed."
            )
        return v


class TumorWorkflowConfig(WorkflowConfigBase):
    """Master configuration for the unified tumor analysis workflow."""

    model_config = ConfigDict(extra="ignore")

    # Sub-stage configs (expert override layer)
    qc_config: Optional[Any] = Field(
        default=None, description="Optional QCWorkflowConfig override."
    )
    preprocess_config: Optional[Any] = Field(
        default=None, description="Optional PreprocessingWorkflowConfig override."
    )
    analysis_config: Optional[Any] = Field(
        default=None, description="Optional AnalysisWorkflowConfig override."
    )
    tumor_config: Optional[TumorAnalysisConfig] = Field(
        default_factory=TumorAnalysisConfig,
        description="Tumor-specific analysis configuration.",
    )
    recommendation_config: Optional[Any] = Field(
        default=None, description="Optional RecommendationConfig override."
    )

    # Workflow behavior
    use_recommendations: bool = Field(
        default=True,
        description="Whether to run the recommendation engine and apply defaults.",
    )
    tissue_type: str = Field(
        default="tumor",
        description="Tissue type hint passed to QC and preprocessing recommenders.",
    )
    batch_key: Optional[str] = Field(
        default=None, description="Batch key for recommendation and integration."
    )
    cancer_type: Optional[str] = Field(
        default=None,
        description="Cancer type for loading cancer-specific markers (e.g., 'Lung Cancer').",
    )

    @classmethod
    def from_simple_dict(cls, simple_config: Dict[str, Any]) -> "TumorWorkflowConfig":
        """Create TumorWorkflowConfig from a simplified flat dictionary."""
        mapping = {
            "save_dir": "save_dir",
            "n_jobs": "n_jobs",
            "random_state": "random_state",
            "tissue_type": "tissue_type",
            "batch_key": "batch_key",
            "cancer_type": "cancer_type",
            "use_recommendations": "use_recommendations",
            "run_malignancy": ("tumor_config", "run_malignancy"),
            "malignancy_method": ("tumor_config", "malignancy_method"),
            "run_tme": ("tumor_config", "run_tme"),
            "run_cnv": ("tumor_config", "run_cnv"),
            "run_therapy": ("tumor_config", "run_therapy"),
        }

        nested: Dict[str, Any] = {}
        for key, value in simple_config.items():
            if key not in mapping:
                logger.warning(f"Unknown key in simple_config: {key}")
                continue
            target = mapping[key]
            if isinstance(target, tuple):
                section, field_name = target
                nested.setdefault(section, {})[field_name] = value
            else:
                nested[target] = value

        return cls.model_validate(nested)

    @classmethod
    def quick(
        cls,
        save_dir: Optional[str] = None,
        tissue_type: str = "tumor",
        batch_key: Optional[str] = None,
    ) -> "TumorWorkflowConfig":
        """Factory for a runnable default tumor workflow config."""
        from ..qc.config import (
            QCWorkflowConfig,
            MetricsReportingConfig,
            MarkingConfig,
            DoubletConfig,
            FilterConfig,
        )
        from ..preprocess.config import WorkflowConfig as PreprocessWorkflowConfig, IntegrationConfig
        from ..analysis.config import AnalysisWorkflowConfig

        from ..qc.config import QCThresholds

        qc_config = QCWorkflowConfig(
            metrics_reporting_config=MetricsReportingConfig(show_plots=False),
            marking_config=MarkingConfig(
                show_plots=False,
                thresholds=QCThresholds(min_genes=10, pc_mt=50.0),
            ),
            doublet_config=DoubletConfig(run_algorithm=False, use_heuristics=False),
            filter_config=FilterConfig(
                criteria_to_filter=["outlier_min_genes", "outlier_mt"]
            ),
        )

        return cls(
            save_dir=save_dir,
            tissue_type=tissue_type,
            batch_key=batch_key,
            use_recommendations=True,
            qc_config=qc_config,
            preprocess_config=PreprocessWorkflowConfig(
                integration=IntegrationConfig(method=None)
            ),
            analysis_config=AnalysisWorkflowConfig(),
            tumor_config=TumorAnalysisConfig(),
        )


__all__ = [
    "TumorAnalysisConfig",
    "TumorWorkflowConfig",
]
