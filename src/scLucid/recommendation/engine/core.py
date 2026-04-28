"""Core RecommendationEngine class — unified cross-stage parameter recommendation.

This is the main entry point, assembled from focused mix-in modules.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional, Union

import numpy as np
from anndata import AnnData

from ...analysis.clustering import find_resolution
from ...analysis.config import ResolutionSearchConfig
from ...preprocess.intelligent import (
    IntelligentPreprocessConfig,
    IntelligentPreprocessRecommender,
    PreprocessingStrategy,
)
from ...qc.intelligent_qc import IntelligentQCConfig, QCRecommendation, recommend_intelligent_qc
from ...utils.context import AnalysisContext, infer_analysis_context
from ..config import RecommendationConfig
from ..schema import WorkflowRecommendations
from ..tumor_adapter import adapt_tumor_recommendation
from .adapters import AdapterMixin
from .annotation import AnnotationMixin

log = logging.getLogger(__name__)


class RecommendationEngine(AnnotationMixin, AdapterMixin):
    """Single entry point for cross-stage parameter recommendations."""

    def __init__(self, config: Optional[RecommendationConfig] = None):
        self.config = config or RecommendationConfig()

    def recommend(
        self,
        adata: AnnData,
        *,
        context: Optional[Union[AnalysisContext, Dict[str, Any]]] = None,
        dataset_type: Optional[str] = None,
        batch_key: Optional[str] = None,
        tissue_type: str = "unknown",
        tissue: Optional[str] = None,
        species: str = "human",
        cancer_type: Optional[str] = None,
        plot: bool = False,
        save_dir: Optional[Path] = None,
    ) -> WorkflowRecommendations:
        """Run the selected recommenders and return a unified bundle."""
        sections: Dict[str, Any] = {}

        qc_raw: Optional[QCRecommendation] = None
        preprocess_raw: Optional[PreprocessingStrategy] = None
        analysis_context = infer_analysis_context(
            adata,
            context=context,
            dataset_type=dataset_type,
            species=species,
            tissue=tissue,
            tissue_type=tissue_type,
            cancer_type=cancer_type,
            batch_key=batch_key,
        )
        effective_tissue_type = analysis_context.qc_tissue_type
        effective_batch_key = batch_key or analysis_context.batch_key
        effective_cancer_type = cancer_type or analysis_context.cancer_type

        if "qc" in self.config.modules:
            qc_raw = self._recommend_qc(
                adata,
                tissue_type=effective_tissue_type,
                plot=plot,
                save_dir=save_dir,
            )
            sections["qc"] = self._adapt_qc(qc_raw)

        if "preprocess" in self.config.modules:
            preprocess_raw = self._recommend_preprocess(
                adata,
                batch_key=effective_batch_key,
                tissue_type=effective_tissue_type,
                plot=plot,
                save_dir=save_dir,
            )
            sections["preprocess"] = self._adapt_preprocess(preprocess_raw)

        if "clustering" in self.config.modules:
            if preprocess_raw is not None:
                sections["clustering"] = self._adapt_clustering_from_preprocess(preprocess_raw)
            else:
                sections["clustering"] = self._recommend_clustering_only(
                    adata,
                    plot=plot,
                    save_dir=save_dir,
                )

        if "annotation" in self.config.modules:
            sections["annotation"] = self._recommend_annotation(
                adata,
                context=analysis_context,
                clustering_section=sections.get("clustering"),
            )

        if "tumor" in self.config.modules:
            sections["tumor"] = self._recommend_tumor(
                adata,
                context=analysis_context,
                cancer_type=effective_cancer_type,
            )

        concerns = []
        for section in sections.values():
            concerns.extend(section.concerns)
        concerns = list(dict.fromkeys(concerns))

        context_dict = self._build_context(
            adata,
            analysis_context=analysis_context,
            batch_key=effective_batch_key,
        )
        overall_confidence = (
            float(np.mean([section.confidence for section in sections.values()]))
            if sections
            else 0.0
        )

        return WorkflowRecommendations(
            sections=sections,
            overall_confidence=overall_confidence,
            context=context_dict,
            concerns=concerns,
        )

    def _recommend_qc(
        self,
        adata: AnnData,
        *,
        tissue_type: str,
        plot: bool,
        save_dir: Optional[Path],
    ) -> QCRecommendation:
        qc_config = self.config.qc or IntelligentQCConfig()
        qc_overrides = {
            key: value
            for key, value in qc_config.model_dump().items()
            if key not in {"plot", "save_dir", "tissue_type", "sample_metadata", "verbose", "report"}
        }
        return recommend_intelligent_qc(
            adata,
            tissue_type=tissue_type,
            plot=plot,
            save_dir=save_dir,
            **qc_overrides,
        )

    def _recommend_preprocess(
        self,
        adata: AnnData,
        *,
        batch_key: Optional[str],
        tissue_type: str,
        plot: bool,
        save_dir: Optional[Path],
    ) -> PreprocessingStrategy:
        preprocess_config = self.config.preprocess or IntelligentPreprocessConfig()
        recommender = IntelligentPreprocessRecommender(config=preprocess_config)
        return recommender.recommend(
            adata,
            batch_key=batch_key,
            tissue_type=tissue_type,
            plot=plot,
            save_dir=save_dir,
        )

    def _recommend_clustering_only(
        self,
        adata: AnnData,
        *,
        plot: bool,
        save_dir: Optional[Path],
    ) -> Any:
        working = adata.copy()
        use_rep = self.config.clustering_use_rep

        if use_rep not in working.obsm:
            if not self.config.prepare_clustering_rep_if_missing:
                raise ValueError(
                    f"Representation '{use_rep}' not found and automatic preparation is disabled."
                )
            working = self._prepare_representation(working, use_rep=use_rep)

        resolution_config = self.config.resolution_search or ResolutionSearchConfig(
            method=self.config.clustering_method,
            use_rep=use_rep,
            plot=plot,
            save_dir=str(save_dir) if save_dir is not None else None,
        )
        resolution_config = resolution_config.model_copy(
            update={
                "method": self.config.clustering_method,
                "use_rep": use_rep,
                "plot": plot,
                "save_dir": str(save_dir) if save_dir is not None else None,
            }
        )

        eval_df, recommended = find_resolution(
            working,
            config=resolution_config,
            auto_select=True,
            selection_strategy=self.config.clustering_selection_strategy,
        )
        return self._adapt_clustering_search(eval_df, recommended)

    def _recommend_tumor(
        self,
        adata: AnnData,
        *,
        context: AnalysisContext,
        cancer_type: Optional[str] = None,
    ) -> Any:
        from ...tumor.config import TumorAnalysisConfig

        tumor_config = self.config.tumor or TumorAnalysisConfig()
        return adapt_tumor_recommendation(
            adata,
            config=tumor_config,
            context=context,
            cancer_type=cancer_type,
        )


def recommend_analysis_parameters(
    adata: AnnData,
    *,
    context: Optional[Union[AnalysisContext, Dict[str, Any]]] = None,
    dataset_type: Optional[str] = None,
    batch_key: Optional[str] = None,
    tissue_type: str = "unknown",
    tissue: Optional[str] = None,
    species: str = "human",
    cancer_type: Optional[str] = None,
    plot: bool = False,
    save_dir: Optional[Path] = None,
    config: Optional[RecommendationConfig] = None,
) -> WorkflowRecommendations:
    """Convenience entry point for the unified recommendation engine."""
    engine = RecommendationEngine(config=config)
    return engine.recommend(
        adata,
        context=context,
        dataset_type=dataset_type,
        batch_key=batch_key,
        tissue_type=tissue_type,
        tissue=tissue,
        species=species,
        cancer_type=cancer_type,
        plot=plot,
        save_dir=save_dir,
    )
