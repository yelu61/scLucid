"""Annotation evidence assessment and strategy selection methods.

Extracted from engine.py for maintainability.
"""

from __future__ import annotations

import importlib.util
import logging
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from ...analysis.config import AnnotationConfig
from ...utils import get_marker_manager
from ..schema import ParameterRecommendation, RecommendationSection

log = logging.getLogger(__name__)


class AnnotationMixin:
    """Mix-in providing annotation recommendation methods for RecommendationEngine."""

    def _recommend_annotation(
        self,
        adata: AnnData,
        *,
        context,
        clustering_section=None,
    ) -> RecommendationSection:
        from ...utils.context import AnalysisContext

        annotation_config = self.config.annotation or AnnotationConfig()
        if annotation_config.marker_tissue is None and context.tissue is not None:
            annotation_config = annotation_config.model_copy(
                update={"marker_tissue": context.tissue}
            )
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
        existing_celltypist = any(col.startswith("celltypist_") for col in adata.obs.columns)
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
            concerns.append(
                "Built-in marker sets have limited overlap with the current gene space."
            )
        if not (celltypist_available or existing_celltypist):
            concerns.append(
                "CellTypist is unavailable; recommendation falls back to marker-based annotation."
            )
        if cluster_key not in adata.obs.columns:
            concerns.append(
                f"Cluster key '{cluster_key}' is not present yet; annotation recommendation uses expected cluster counts only."
            )
        if existing_celltypist and agreement_score is not None and agreement_score < 0.4:
            concerns.append(
                "Marker evidence and existing CellTypist labels disagree substantially across clusters."
            )
        if (
            existing_celltypist
            and celltypist_evidence["mean_confidence"]
            < annotation_config.celltypist_confidence_threshold
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
                alternatives=["enrichment", "max_score", "combined"],
            ),
            ParameterRecommendation(
                name="run_celltypist",
                value=run_celltypist,
                method="dependency_and_evidence_check",
                confidence=0.9 if (celltypist_available or existing_celltypist) else 0.6,
                rationale="Whether CellTypist should be executed as part of annotation.",
                alternatives=[True, False],
            ),
            ParameterRecommendation(
                name="run_scoring",
                value=run_scoring,
                method="marker_support_assessment",
                confidence=float(annotation_evidence["eligible_ratio"]),
                rationale="Whether marker gene scoring should be part of the annotation workflow.",
                alternatives=[True, False],
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
                "dataset_type": context.dataset_type,
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

    def _assess_annotation_evidence(
        self,
        adata: AnnData,
        *,
        cluster_key: str,
        marker_species: str,
        marker_tissue: Optional[str],
    ) -> Dict[str, Any]:
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
                        marker_lookup[marker] for marker in cell.markers if marker in marker_lookup
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
                signal = 0.7 * self._sigmoid(best_score) + 0.3 * self._sigmoid(
                    best_score - second_score
                )
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
        eligible_ratio = float(marker_evidence.get("eligible_ratio", 0.0))
        cluster_marker_signal = marker_evidence.get("cluster_marker_signal")
        marker_signal = (
            float(cluster_marker_signal) if cluster_marker_signal is not None else eligible_ratio
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
            if (
                mean_confidence >= 0.7
                and cluster_purity >= 0.7
                and (agreement is None or agreement >= 0.55)
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
    def _celltypist_available() -> bool:
        return importlib.util.find_spec("celltypist") is not None

    @staticmethod
    def _sigmoid(value: float) -> float:
        clipped = float(np.clip(value, -20.0, 20.0))
        return float(1.0 / (1.0 + np.exp(-clipped)))
