"""Adapter methods that convert recommendation outputs to structured sections.

Extracted from engine.py for maintainability.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ...qc.intelligent_qc import QCRecommendation
from ..schema import ParameterRecommendation, RecommendationSection


class AdapterMixin:
    """Mix-in providing adapter/converter methods for RecommendationEngine."""

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
    def _adapt_preprocess(strategy) -> RecommendationSection:
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
    def _adapt_clustering_from_preprocess(strategy) -> RecommendationSection:
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
        if pd.notna(row.get("interpretability_score")):
            confidence_parts.append(float(np.clip(row["interpretability_score"], 0.0, 1.0)))
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
                method="clustering_review",
                confidence=confidence,
                rationale="Best-supported practical clustering resolution from review evidence.",
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
    def _prepare_representation(adata: AnnData, *, use_rep: str) -> AnnData:
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
    def _build_context(
        adata: AnnData,
        *,
        analysis_context,
        batch_key: Optional[str],
    ) -> Dict[str, Any]:
        context = analysis_context.to_dict()
        context.update(
            {
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
            }
        )
        if batch_key is not None and batch_key in adata.obs.columns:
            context["batch_key"] = batch_key
            context["n_batches"] = int(adata.obs[batch_key].nunique())
        return context
