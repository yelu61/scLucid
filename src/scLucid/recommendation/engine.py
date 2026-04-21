"""Unified recommendation engine for QC, preprocessing, clustering, and annotation."""

from __future__ import annotations

import importlib.util
import logging
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ..analysis.clustering import find_resolution
from ..analysis.config import AnnotationConfig, ResolutionSearchConfig
from ..preprocess.intelligent import (
    IntelligentPreprocessConfig,
    IntelligentPreprocessRecommender,
    PreprocessingStrategy,
)
from ..qc.intelligent_qc import IntelligentQCConfig, QCRecommendation, recommend_intelligent_qc
from ..tumor.config import TumorAnalysisConfig
from ..utils import get_marker_manager
from .config import RecommendationConfig
from .schema import ParameterRecommendation, RecommendationSection, WorkflowRecommendations
from .tumor_adapter import adapt_tumor_recommendation

log = logging.getLogger(__name__)


class RecommendationEngine:
    """Single entry point for cross-stage parameter recommendations."""

    def __init__(self, config: Optional[RecommendationConfig] = None):
        self.config = config or RecommendationConfig()

    def recommend(
        self,
        adata: AnnData,
        *,
        batch_key: Optional[str] = None,
        tissue_type: str = "unknown",
        cancer_type: Optional[str] = None,
        plot: bool = False,
        save_dir: Optional[Path] = None,
    ) -> WorkflowRecommendations:
        """Run the selected recommenders and return a unified bundle."""
        sections: Dict[str, RecommendationSection] = {}

        qc_raw: Optional[QCRecommendation] = None
        preprocess_raw: Optional[PreprocessingStrategy] = None

        if "qc" in self.config.modules:
            qc_raw = self._recommend_qc(
                adata,
                tissue_type=tissue_type,
                plot=plot,
                save_dir=save_dir,
            )
            sections["qc"] = self._adapt_qc(qc_raw)

        if "preprocess" in self.config.modules:
            preprocess_raw = self._recommend_preprocess(
                adata,
                batch_key=batch_key,
                tissue_type=tissue_type,
                plot=plot,
                save_dir=save_dir,
            )
            sections["preprocess"] = self._adapt_preprocess(preprocess_raw)

        if "clustering" in self.config.modules:
            if preprocess_raw is not None:
                sections["clustering"] = self._adapt_clustering_from_preprocess(
                    preprocess_raw
                )
            else:
                sections["clustering"] = self._recommend_clustering_only(
                    adata,
                    plot=plot,
                    save_dir=save_dir,
                )

        if "annotation" in self.config.modules:
            sections["annotation"] = self._recommend_annotation(
                adata,
                clustering_section=sections.get("clustering"),
            )

        if "tumor" in self.config.modules:
            sections["tumor"] = self._recommend_tumor(adata, cancer_type=cancer_type)

        concerns = []
        for section in sections.values():
            concerns.extend(section.concerns)
        concerns = list(dict.fromkeys(concerns))

        context = self._build_context(
            adata,
            tissue_type=tissue_type,
            batch_key=batch_key,
        )
        overall_confidence = float(
            np.mean([section.confidence for section in sections.values()])
        ) if sections else 0.0

        return WorkflowRecommendations(
            sections=sections,
            overall_confidence=overall_confidence,
            context=context,
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
            if key not in {"plot", "save_dir", "tissue_type", "sample_metadata"}
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
    ) -> RecommendationSection:
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

    def _recommend_annotation(
        self,
        adata: AnnData,
        *,
        clustering_section: Optional[RecommendationSection] = None,
    ) -> RecommendationSection:
        annotation_config = self.config.annotation or AnnotationConfig()
        cluster_key = annotation_config.cluster_key

        expected_n_clusters = None
        if cluster_key in adata.obs.columns:
            expected_n_clusters = int(pd.Series(adata.obs[cluster_key]).nunique())
        elif clustering_section is not None:
            n_clusters_param = clustering_section.get_parameter("n_clusters")
            if n_clusters_param is not None:
                expected_n_clusters = int(n_clusters_param.value)

        annotation_evidence = self._assess_annotation_evidence(
            adata,
            cluster_key=cluster_key,
            marker_species=annotation_config.marker_species,
            marker_tissue=annotation_config.marker_tissue,
        )
        celltypist_available = self._celltypist_available()
        existing_celltypist = any(
            col.startswith("celltypist_") for col in adata.obs.columns
        )
        celltypist_evidence = self._evaluate_existing_celltypist_evidence(
            adata,
            cluster_key=cluster_key,
        )
        agreement_score = self._compute_annotation_agreement(
            annotation_evidence.get("cluster_best_labels", {}),
            celltypist_evidence.get("cluster_labels", {}),
        )
        strategy = self._select_annotation_strategy(
            marker_evidence=annotation_evidence,
            celltypist_evidence=celltypist_evidence,
            celltypist_available=celltypist_available,
            existing_celltypist=existing_celltypist,
            expected_n_clusters=expected_n_clusters,
        )
        marker_method = strategy["marker_method"]
        final_method = strategy["final_method"]
        run_celltypist = strategy["run_celltypist"]
        run_scoring = strategy["run_scoring"]
        confidence = float(strategy["confidence"])

        concerns = []
        if annotation_evidence["eligible_types"] == 0:
            concerns.append("Built-in marker sets have limited overlap with the current gene space.")
        if not (celltypist_available or existing_celltypist):
            concerns.append("CellTypist is unavailable; recommendation falls back to marker-based annotation.")
        if cluster_key not in adata.obs.columns:
            concerns.append(
                f"Cluster key '{cluster_key}' is not present yet; annotation recommendation uses expected cluster counts only."
            )
        if (
            existing_celltypist
            and agreement_score is not None
            and agreement_score < 0.4
        ):
            concerns.append(
                "Marker evidence and existing CellTypist labels disagree substantially across clusters."
            )
        if (
            existing_celltypist
            and celltypist_evidence["mean_confidence"] < annotation_config.celltypist_confidence_threshold
        ):
            concerns.append(
                "Existing CellTypist confidence is below the configured acceptance threshold."
            )

        notes = [
            f"Marker-supported cell types: {annotation_evidence['eligible_types']}/{annotation_evidence['total_types']}",
        ]
        if expected_n_clusters is not None:
            notes.append(f"Expected clusters for annotation: {expected_n_clusters}")
        if annotation_evidence["cluster_marker_signal"] is not None:
            notes.append(
                f"Cluster-level marker signal: {annotation_evidence['cluster_marker_signal']:.3f}"
            )
        if existing_celltypist:
            notes.append(
                f"Existing CellTypist mean confidence: {celltypist_evidence['mean_confidence']:.3f}"
            )
            notes.append(
                f"Existing CellTypist cluster purity: {celltypist_evidence['cluster_purity']:.3f}"
            )
        if agreement_score is not None:
            notes.append(f"Marker vs CellTypist agreement: {agreement_score:.3f}")

        recommended_config = annotation_config.model_copy(
            update={
                "run_celltypist": run_celltypist,
                "run_scoring": run_scoring,
                "final_method": final_method,
                "marker_method": marker_method,
            }
        )

        parameters = [
            ParameterRecommendation(
                name="final_method",
                value=final_method,
                method="annotation_evidence_integration",
                confidence=confidence,
                rationale="Primary annotation path selected from available evidence sources.",
                alternatives=[marker_method, "celltypist", "hybrid"],
                evidence={
                    "marker_supported_ratio": annotation_evidence["eligible_ratio"],
                    "cluster_marker_signal": annotation_evidence["cluster_marker_signal"],
                    "celltypist_mean_confidence": celltypist_evidence["mean_confidence"],
                    "celltypist_cluster_purity": celltypist_evidence["cluster_purity"],
                    "marker_celltypist_agreement": agreement_score,
                },
            ),
            ParameterRecommendation(
                name="marker_method",
                value=marker_method,
                method="marker_support_assessment",
                confidence=float(
                    annotation_evidence["cluster_marker_signal"]
                    if annotation_evidence["cluster_marker_signal"] is not None
                    else annotation_evidence["eligible_ratio"]
                ),
                rationale="Marker-based method recommended for the marker evidence branch.",
                evidence={
                    "cluster_best_labels": annotation_evidence.get("cluster_best_labels", {}),
                },
            ),
            ParameterRecommendation(
                name="run_celltypist",
                value=run_celltypist,
                method="dependency_and_evidence_check",
                confidence=0.9 if (celltypist_available or existing_celltypist) else 0.6,
                rationale="Whether CellTypist should be executed as part of annotation.",
            ),
            ParameterRecommendation(
                name="run_scoring",
                value=run_scoring,
                method="marker_support_assessment",
                confidence=float(annotation_evidence["eligible_ratio"]),
                rationale="Whether marker gene scoring should be part of the annotation workflow.",
            ),
            ParameterRecommendation(
                name="celltypist_model",
                value=annotation_config.celltypist_model,
                method="config_default",
                confidence=0.7,
                rationale="CellTypist model suggested for automated label transfer.",
            ),
            ParameterRecommendation(
                name="celltypist_confidence_threshold",
                value=float(annotation_config.celltypist_confidence_threshold),
                method="config_default",
                confidence=0.7,
                rationale="Minimum CellTypist confidence to accept direct labels.",
            ),
        ]

        return RecommendationSection(
            name="annotation",
            summary=(
                f"Recommend annotation strategy '{final_method}' "
                f"with marker branch '{marker_method}'."
            ),
            confidence=confidence,
            parameters=parameters,
            concerns=concerns,
            notes=notes,
            raw_result=recommended_config,
            metadata={
                "cluster_key": cluster_key,
                "marker_species": annotation_config.marker_species,
                "marker_tissue": annotation_config.marker_tissue,
                "key_added": annotation_config.key_added,
                "expected_n_clusters": expected_n_clusters,
                "marker_eligible_ratio": annotation_evidence["eligible_ratio"],
                "cluster_marker_signal": annotation_evidence["cluster_marker_signal"],
                "cluster_best_labels": annotation_evidence.get("cluster_best_labels", {}),
                "celltypist_mean_confidence": celltypist_evidence["mean_confidence"],
                "celltypist_cluster_purity": celltypist_evidence["cluster_purity"],
                "celltypist_cluster_labels": celltypist_evidence.get("cluster_labels", {}),
                "marker_celltypist_agreement": agreement_score,
            },
        )

    @staticmethod
    def _prepare_representation(adata: AnnData, *, use_rep: str) -> AnnData:
        """Build a basic PCA representation for standalone clustering recommendation."""
        if use_rep != "X_pca":
            raise ValueError(
                "Automatic representation preparation currently supports only 'X_pca'."
            )

        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()

        if "total_counts" not in adata.obs:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
        else:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

        n_top_genes = min(2000, max(200, adata.n_vars))
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=min(n_top_genes, adata.n_vars),
            flavor="seurat_v3",
            subset=False,
        )
        if "highly_variable" in adata.var and bool(adata.var["highly_variable"].sum()):
            adata = adata[:, adata.var["highly_variable"]].copy()

        sc.pp.scale(adata, max_value=10)
        n_comps = min(50, max(2, adata.n_vars - 1), max(2, adata.n_obs - 1))
        sc.tl.pca(adata, n_comps=n_comps)
        return adata

    @staticmethod
    def _adapt_qc(recommendation: QCRecommendation) -> RecommendationSection:
        parameters = [
            ParameterRecommendation(
                name="min_genes",
                value=int(recommendation.min_genes.threshold),
                ci_lower=int(recommendation.min_genes.ci_lower),
                ci_upper=int(recommendation.min_genes.ci_upper),
                method=recommendation.min_genes.method,
                confidence=recommendation.min_genes.confidence,
                rationale="Lower bound for retained genes per cell.",
                evidence=recommendation.min_genes.evidence,
            ),
            ParameterRecommendation(
                name="max_mt_percent",
                value=float(recommendation.max_mt_percent.threshold),
                ci_lower=float(recommendation.max_mt_percent.ci_lower),
                ci_upper=float(recommendation.max_mt_percent.ci_upper),
                method=recommendation.max_mt_percent.method,
                confidence=recommendation.max_mt_percent.confidence,
                rationale="Upper bound for mitochondrial fraction.",
                evidence=recommendation.max_mt_percent.evidence,
            ),
            ParameterRecommendation(
                name="doublet_threshold",
                value=float(recommendation.doublet_threshold.threshold),
                ci_lower=float(recommendation.doublet_threshold.ci_lower),
                ci_upper=float(recommendation.doublet_threshold.ci_upper),
                method=recommendation.doublet_threshold.method,
                confidence=recommendation.doublet_threshold.confidence,
                rationale="Score cutoff for flagging likely doublets.",
                evidence=recommendation.doublet_threshold.evidence,
            ),
            ParameterRecommendation(
                name="n_counts",
                value=float(recommendation.n_counts.threshold),
                ci_lower=float(recommendation.n_counts.ci_lower),
                ci_upper=float(recommendation.n_counts.ci_upper),
                method=recommendation.n_counts.method,
                confidence=recommendation.n_counts.confidence,
                rationale="Minimum count depth required per cell.",
                evidence=recommendation.n_counts.evidence,
            ),
        ]
        notes = [f"QC strategy: {recommendation.overall_strategy.value}"]
        notes.extend(recommendation.tumor_specific_considerations)
        return RecommendationSection(
            name="qc",
            summary=(
                f"QC strategy '{recommendation.overall_strategy.value}' with "
                f"data quality score {recommendation.data_quality_score:.1f}/100."
            ),
            confidence=recommendation.overall_confidence,
            parameters=parameters,
            concerns=recommendation.concerns,
            notes=notes,
            raw_result=recommendation,
            metadata={
                "data_quality_score": recommendation.data_quality_score,
                "strategy": recommendation.overall_strategy.value,
            },
        )

    @staticmethod
    def _adapt_preprocess(strategy: PreprocessingStrategy) -> RecommendationSection:
        parameters = [
            ParameterRecommendation(
                name="n_top_genes",
                value=int(strategy.hvg.n_top_genes),
                ci_lower=int(strategy.hvg.ci_lower),
                ci_upper=int(strategy.hvg.ci_upper),
                method=strategy.hvg.method,
                confidence=strategy.hvg.confidence,
                rationale="Number of highly variable genes to keep.",
                evidence=strategy.hvg.evidence,
            ),
            ParameterRecommendation(
                name="n_pcs",
                value=int(strategy.pca.n_pcs),
                ci_lower=int(strategy.pca.ci_lower),
                ci_upper=int(strategy.pca.ci_upper),
                method=strategy.pca.method,
                confidence=strategy.pca.confidence,
                rationale="Number of PCs to retain before graph construction.",
                evidence=strategy.pca.evidence,
            ),
            ParameterRecommendation(
                name="n_neighbors",
                value=int(strategy.neighbors.n_neighbors),
                ci_lower=int(strategy.neighbors.ci_lower_neighbors),
                ci_upper=int(strategy.neighbors.ci_upper_neighbors),
                method=strategy.neighbors.method,
                confidence=strategy.neighbors.confidence,
                rationale="Graph neighborhood size for manifold construction.",
                evidence=strategy.neighbors.evidence,
            ),
            ParameterRecommendation(
                name="graph_n_pcs",
                value=int(strategy.neighbors.n_pcs),
                ci_lower=int(strategy.neighbors.ci_lower_pcs),
                ci_upper=int(strategy.neighbors.ci_upper_pcs),
                method=strategy.neighbors.method,
                confidence=strategy.neighbors.confidence,
                rationale="PC count used for nearest-neighbor graph construction.",
                evidence=strategy.neighbors.evidence,
            ),
        ]
        notes = list(strategy.recommendations)
        metadata: Dict[str, Any] = {
            "strategy_type": strategy.data_profile.strategy_type,
            "data_quality_score": strategy.data_profile.data_quality_score,
        }
        if strategy.batch_correction is not None:
            parameters.append(
                ParameterRecommendation(
                    name="batch_correction_method",
                    value=strategy.batch_correction.recommended_method,
                    method="batch_effect_assessment",
                    confidence=strategy.batch_correction.confidence,
                    rationale="Recommended integration method based on batch severity.",
                    evidence=strategy.batch_correction.evidence,
                    alternatives=strategy.batch_correction.alternative_methods,
                )
            )
            metadata["batch_effect_severity"] = strategy.batch_correction.severity_score

        return RecommendationSection(
            name="preprocess",
            summary=(
                f"Preprocessing strategy '{strategy.data_profile.strategy_type}' "
                f"for {strategy.data_profile.n_cells} cells."
            ),
            confidence=strategy.overall_confidence,
            parameters=parameters,
            concerns=strategy.concerns,
            notes=notes,
            raw_result=strategy,
            metadata=metadata,
        )

    @staticmethod
    def _adapt_clustering_from_preprocess(
        strategy: PreprocessingStrategy,
    ) -> RecommendationSection:
        parameters = [
            ParameterRecommendation(
                name="resolution",
                value=float(strategy.resolution.resolution),
                ci_lower=float(strategy.resolution.ci_lower),
                ci_upper=float(strategy.resolution.ci_upper),
                method=strategy.resolution.method,
                confidence=strategy.resolution.confidence,
                rationale="Graph clustering resolution chosen from stability analysis.",
                evidence=strategy.resolution.evidence,
            ),
            ParameterRecommendation(
                name="n_clusters",
                value=int(strategy.resolution.n_clusters),
                method="derived_from_resolution",
                confidence=strategy.resolution.confidence,
                rationale="Expected cluster count at the recommended resolution.",
            ),
            ParameterRecommendation(
                name="n_neighbors",
                value=int(strategy.neighbors.n_neighbors),
                ci_lower=int(strategy.neighbors.ci_lower_neighbors),
                ci_upper=int(strategy.neighbors.ci_upper_neighbors),
                method=strategy.neighbors.method,
                confidence=strategy.neighbors.confidence,
                rationale="Neighborhood size that supports the recommended clustering.",
            ),
            ParameterRecommendation(
                name="n_pcs",
                value=int(strategy.neighbors.n_pcs),
                ci_lower=int(strategy.neighbors.ci_lower_pcs),
                ci_upper=int(strategy.neighbors.ci_upper_pcs),
                method=strategy.neighbors.method,
                confidence=strategy.neighbors.confidence,
                rationale="Embedding dimensionality used for clustering.",
            ),
        ]
        return RecommendationSection(
            name="clustering",
            summary=(
                f"Recommend resolution {strategy.resolution.resolution:.2f} "
                f"with ~{strategy.resolution.n_clusters} clusters."
            ),
            confidence=strategy.resolution.confidence,
            parameters=parameters,
            concerns=[],
            notes=[
                f"Stability score: {strategy.resolution.stability_score:.3f}",
            ],
            raw_result=strategy.resolution,
            metadata={
                "n_clusters": strategy.resolution.n_clusters,
                "stability_score": strategy.resolution.stability_score,
            },
        )

    @staticmethod
    def _adapt_clustering_search(
        eval_df: pd.DataFrame,
        recommended: Optional[float],
    ) -> RecommendationSection:
        if eval_df.empty or recommended is None:
            return RecommendationSection(
                name="clustering",
                summary="Clustering recommendation unavailable.",
                confidence=0.0,
                parameters=[],
                concerns=["Resolution search did not yield a valid recommendation."],
                raw_result=eval_df,
            )

        distances = (eval_df["resolution"] - float(recommended)).abs()
        idx = int(distances.idxmin())
        row = eval_df.loc[idx]
        ci_lower = float(eval_df.iloc[max(idx - 1, 0)]["resolution"])
        ci_upper = float(eval_df.iloc[min(idx + 1, len(eval_df) - 1)]["resolution"])

        confidence_parts = []
        if pd.notna(row.get("silhouette")):
            confidence_parts.append(float(np.clip((row["silhouette"] + 1.0) / 2.0, 0.0, 1.0)))
        if pd.notna(row.get("stability")):
            confidence_parts.append(float(np.clip(row["stability"], 0.0, 1.0)))
        if "marker_abundance" in eval_df.columns and eval_df["marker_abundance"].notna().any():
            denom = eval_df["marker_abundance"].max() - eval_df["marker_abundance"].min()
            if denom > 0 and pd.notna(row.get("marker_abundance")):
                normalized = (row["marker_abundance"] - eval_df["marker_abundance"].min()) / denom
                confidence_parts.append(float(np.clip(normalized, 0.0, 1.0)))
        confidence = float(np.mean(confidence_parts)) if confidence_parts else 0.5

        alternatives = (
            eval_df.assign(_distance=distances)
            .sort_values(by=["_distance", "resolution"])
            .head(3)["resolution"]
            .tolist()
        )
        evidence = row.to_dict()
        parameters = [
            ParameterRecommendation(
                name="resolution",
                value=float(recommended),
                ci_lower=ci_lower,
                ci_upper=ci_upper,
                method="resolution_grid_search",
                confidence=confidence,
                rationale="Best-scoring clustering resolution from grid search.",
                evidence=evidence,
                alternatives=alternatives,
            ),
            ParameterRecommendation(
                name="n_clusters",
                value=int(row["n_clusters"]),
                method="derived_from_resolution",
                confidence=confidence,
                rationale="Expected number of clusters at the selected resolution.",
            ),
        ]
        notes = []
        if pd.notna(row.get("silhouette")):
            notes.append(f"Silhouette: {float(row['silhouette']):.3f}")
        if pd.notna(row.get("stability")):
            notes.append(f"Stability: {float(row['stability']):.3f}")
        return RecommendationSection(
            name="clustering",
            summary=(
                f"Recommend resolution {float(recommended):.2f} "
                f"with {int(row['n_clusters'])} clusters."
            ),
            confidence=confidence,
            parameters=parameters,
            notes=notes,
            raw_result={"search": eval_df, "recommended_resolution": recommended},
            metadata=evidence,
        )

    @staticmethod
    def _build_context(
        adata: AnnData,
        *,
        tissue_type: str,
        batch_key: Optional[str],
    ) -> Dict[str, Any]:
        """Summarize input data characteristics."""
        context = {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "tissue_type": tissue_type,
        }
        if batch_key is not None and batch_key in adata.obs.columns:
            context["batch_key"] = batch_key
            context["n_batches"] = int(adata.obs[batch_key].nunique())
        return context

    @staticmethod
    def _celltypist_available() -> bool:
        """Check whether CellTypist is importable."""
        return importlib.util.find_spec("celltypist") is not None

    def _assess_annotation_evidence(
        self,
        adata: AnnData,
        *,
        cluster_key: str,
        marker_species: str,
        marker_tissue: Optional[str],
    ) -> Dict[str, Any]:
        """Assess marker availability and cluster-level marker signal."""
        marker_stats = self._assess_marker_support(
            adata,
            marker_species=marker_species,
            marker_tissue=marker_tissue,
        )
        marker_stats.setdefault("cluster_marker_signal", None)
        marker_stats.setdefault("cluster_best_labels", {})
        marker_stats.setdefault("cluster_signal_by_label", {})

        if cluster_key not in adata.obs.columns:
            return marker_stats

        try:
            mgr = get_marker_manager(species=marker_species, tissue=marker_tissue)
            mgr.intersect_with(adata.raw if adata.raw is not None else adata)

            base = adata.raw.to_adata() if adata.raw is not None else adata
            feature_names = pd.Index(base.var_names.astype(str))
            available_markers = sorted(
                {
                    marker
                    for cell in mgr.CELLS.values()
                    for marker in cell.markers
                    if marker in feature_names
                }
            )
            if not available_markers:
                return marker_stats

            marker_matrix = base[:, available_markers].X
            if hasattr(marker_matrix, "toarray"):
                marker_matrix = marker_matrix.toarray()
            marker_matrix = np.asarray(marker_matrix, dtype=float)
            if marker_matrix.size == 0:
                return marker_stats

            gene_mean = marker_matrix.mean(axis=0, keepdims=True)
            gene_std = marker_matrix.std(axis=0, keepdims=True) + 1e-8
            marker_matrix = (marker_matrix - gene_mean) / gene_std

            cluster_labels = adata.obs[cluster_key].astype(str).values
            expr_df = pd.DataFrame(
                marker_matrix,
                index=adata.obs_names,
                columns=available_markers,
            )
            cluster_means = expr_df.groupby(cluster_labels).mean()
            marker_lookup = {gene: idx for idx, gene in enumerate(available_markers)}

            cluster_scores: Dict[str, float] = {}
            cluster_best_labels: Dict[str, str] = {}
            cluster_signal_by_label: Dict[str, Dict[str, float]] = {}
            for cluster in cluster_means.index:
                label_scores: Dict[str, float] = {}
                for cell_type, cell in mgr.CELLS.items():
                    indices = [
                        marker_lookup[marker]
                        for marker in cell.markers
                        if marker in marker_lookup
                    ]
                    if len(indices) < 3:
                        continue
                    values = cluster_means.loc[cluster].iloc[indices]
                    label_scores[cell_type] = float(np.mean(values))

                if not label_scores:
                    continue

                sorted_scores = sorted(label_scores.items(), key=lambda item: item[1], reverse=True)
                best_label, best_score = sorted_scores[0]
                second_score = sorted_scores[1][1] if len(sorted_scores) > 1 else 0.0
                signal = 0.7 * self._sigmoid(best_score) + 0.3 * self._sigmoid(best_score - second_score)
                cluster_scores[str(cluster)] = float(signal)
                cluster_best_labels[str(cluster)] = best_label
                cluster_signal_by_label[str(cluster)] = label_scores

            if cluster_scores:
                marker_stats["cluster_marker_signal"] = float(
                    np.mean(list(cluster_scores.values()))
                )
                marker_stats["cluster_best_labels"] = cluster_best_labels
                marker_stats["cluster_signal_by_label"] = cluster_signal_by_label
                marker_stats["cluster_signal_by_cluster"] = cluster_scores
            return marker_stats
        except Exception as exc:
            log.warning(f"Annotation evidence assessment failed: {exc}")
            return marker_stats

    @staticmethod
    def _evaluate_existing_celltypist_evidence(
        adata: AnnData,
        *,
        cluster_key: str,
    ) -> Dict[str, Any]:
        """Summarize confidence and purity from existing CellTypist outputs."""
        majority_key = "celltypist_majority_voting"
        predicted_key = "celltypist_predicted_labels"
        conf_key = "celltypist_conf_score"

        if majority_key in adata.obs.columns:
            label_key = majority_key
        elif predicted_key in adata.obs.columns:
            label_key = predicted_key
        else:
            return {
                "mean_confidence": 0.0,
                "cluster_purity": 0.0,
                "cluster_labels": {},
                "label_source": None,
            }

        labels = adata.obs[label_key].astype(str)
        if conf_key in adata.obs.columns:
            confidences = pd.to_numeric(adata.obs[conf_key], errors="coerce").fillna(0.0)
        else:
            confidences = pd.Series(1.0, index=adata.obs.index, dtype=float)

        result = {
            "mean_confidence": float(confidences.mean()),
            "cluster_purity": 0.0,
            "cluster_labels": {},
            "label_source": label_key,
        }
        if cluster_key not in adata.obs.columns:
            return result

        cluster_series = adata.obs[cluster_key].astype(str)
        purities = []
        cluster_labels = {}
        for cluster in sorted(cluster_series.unique(), key=str):
            mask = cluster_series == str(cluster)
            cluster_pred = labels.loc[mask]
            if cluster_pred.empty:
                continue
            counts = cluster_pred.value_counts()
            top_label = str(counts.index[0])
            purity = float(counts.iloc[0] / counts.sum())
            cluster_labels[str(cluster)] = top_label
            purities.append(purity)

        if purities:
            result["cluster_purity"] = float(np.mean(purities))
            result["cluster_labels"] = cluster_labels
        return result

    @staticmethod
    def _compute_annotation_agreement(
        marker_labels: Dict[str, str],
        celltypist_labels: Dict[str, str],
    ) -> Optional[float]:
        """Compute agreement between marker-derived and CellTypist cluster labels."""
        shared = sorted(set(marker_labels) & set(celltypist_labels), key=str)
        if not shared:
            return None
        agreement = [
            1.0 if marker_labels[cluster] == celltypist_labels[cluster] else 0.0
            for cluster in shared
            if marker_labels[cluster] != "Unknown" and celltypist_labels[cluster] != "Unknown"
        ]
        if not agreement:
            return None
        return float(np.mean(agreement))

    def _select_annotation_strategy(
        self,
        *,
        marker_evidence: Dict[str, Any],
        celltypist_evidence: Dict[str, Any],
        celltypist_available: bool,
        existing_celltypist: bool,
        expected_n_clusters: Optional[int],
    ) -> Dict[str, Any]:
        """Select an annotation strategy from real evidence scores."""
        eligible_ratio = float(marker_evidence.get("eligible_ratio", 0.0))
        cluster_marker_signal = marker_evidence.get("cluster_marker_signal")
        marker_signal = (
            float(cluster_marker_signal)
            if cluster_marker_signal is not None
            else eligible_ratio
        )
        mean_confidence = float(celltypist_evidence.get("mean_confidence", 0.0))
        cluster_purity = float(celltypist_evidence.get("cluster_purity", 0.0))
        agreement = self._compute_annotation_agreement(
            marker_evidence.get("cluster_best_labels", {}),
            celltypist_evidence.get("cluster_labels", {}),
        )

        marker_method = "combined"
        if marker_signal < 0.25:
            marker_method = "enrichment"
        elif expected_n_clusters is not None and expected_n_clusters > 25:
            marker_method = "max_score"

        run_celltypist = bool(celltypist_available or existing_celltypist)
        if not run_celltypist:
            final_method = marker_method
            run_scoring = marker_method in {"max_score", "combined"}
        elif existing_celltypist:
            if mean_confidence >= 0.7 and cluster_purity >= 0.7 and (
                agreement is None or agreement >= 0.55
            ):
                final_method = "hybrid" if marker_signal >= 0.25 else "celltypist"
            elif marker_signal >= 0.45 and (agreement is not None and agreement < 0.4):
                final_method = marker_method
            elif mean_confidence >= 0.7 and marker_signal < 0.25:
                final_method = "celltypist"
            else:
                final_method = "hybrid" if marker_signal >= 0.35 else marker_method
            run_scoring = final_method in {"hybrid", "max_score", "combined"}
        else:
            if marker_signal >= 0.35:
                final_method = "hybrid"
                run_scoring = True
            elif eligible_ratio < 0.2:
                final_method = "celltypist"
                run_scoring = False
            else:
                final_method = "hybrid"
                run_scoring = marker_method in {"max_score", "combined"}

        confidence_parts = [marker_signal]
        if run_celltypist:
            confidence_parts.append(mean_confidence if existing_celltypist else 0.75)
            confidence_parts.append(cluster_purity if existing_celltypist else 0.7)
        else:
            confidence_parts.append(0.45)
        if agreement is not None:
            confidence_parts.append(agreement)
        if expected_n_clusters is not None:
            confidence_parts.append(0.75 if expected_n_clusters <= 20 else 0.6)

        return {
            "final_method": final_method,
            "marker_method": marker_method,
            "run_celltypist": run_celltypist,
            "run_scoring": run_scoring,
            "confidence": float(np.mean(confidence_parts)),
        }

    @staticmethod
    def _assess_marker_support(
        adata: AnnData,
        *,
        marker_species: str,
        marker_tissue: Optional[str],
    ) -> Dict[str, Any]:
        """Estimate how well built-in markers overlap the current dataset."""
        try:
            mgr = get_marker_manager(species=marker_species, tissue=marker_tissue)
            mgr.intersect_with(adata.raw if adata.raw is not None else adata)
            marker_lengths = [len(cell.markers) for cell in mgr.CELLS.values()]
            total_types = len(marker_lengths)
            eligible_types = sum(length >= 3 for length in marker_lengths)
            eligible_ratio = eligible_types / total_types if total_types else 0.0
            return {
                "total_types": total_types,
                "eligible_types": eligible_types,
                "eligible_ratio": float(eligible_ratio),
            }
        except Exception as exc:
            log.warning(f"Marker support assessment failed: {exc}")
            return {
                "total_types": 0,
                "eligible_types": 0,
                "eligible_ratio": 0.0,
            }

    @staticmethod
    def _sigmoid(value: float) -> float:
        """Numerically stable sigmoid helper for confidence scaling."""
        clipped = float(np.clip(value, -20.0, 20.0))
        return float(1.0 / (1.0 + np.exp(-clipped)))

    def _recommend_tumor(
        self,
        adata: AnnData,
        cancer_type: Optional[str] = None,
    ) -> RecommendationSection:
        """Generate tumor-specific analysis recommendations."""
        tumor_config = self.config.tumor or TumorAnalysisConfig()
        return adapt_tumor_recommendation(
            adata,
            config=tumor_config,
            cancer_type=cancer_type,
        )


def recommend_analysis_parameters(
    adata: AnnData,
    *,
    batch_key: Optional[str] = None,
    tissue_type: str = "unknown",
    cancer_type: Optional[str] = None,
    plot: bool = False,
    save_dir: Optional[Path] = None,
    config: Optional[RecommendationConfig] = None,
) -> WorkflowRecommendations:
    """Convenience entry point for the unified recommendation engine."""
    engine = RecommendationEngine(config=config)
    return engine.recommend(
        adata,
        batch_key=batch_key,
        tissue_type=tissue_type,
        cancer_type=cancer_type,
        plot=plot,
        save_dir=save_dir,
    )
