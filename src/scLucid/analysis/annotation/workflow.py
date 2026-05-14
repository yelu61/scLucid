"""High-level annotation workflows and quality evaluation."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import logging
from importlib.metadata import version

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ...utils import Manager, get_marker_manager, sanitize_for_hdf5
from ..config import AnnotationConfig
from ..scoring import score_by_gene_sets
from .cluster import annotate_clusters
from .reference import (
    _build_celltypist_cluster_mapping,
    _combine_marker_and_celltypist_mappings,
    _get_celltypist_label_series,
    run_celltypist,
)
from .scoring import score_cell_types
from .utils import (
    _build_modular_annotation_label,
    _collect_state_signatures,
    _label_matches_target_lineage,
    _map_compartments,
    _rename_score_columns_for_manager,
    _resolve_annotation_manager,
    _state_applies_to_cell,
)

log = logging.getLogger(__name__)


def _sigmoid_score(value: float) -> float:
    """Map an unbounded evidence score onto [0, 1]."""
    clipped = float(np.clip(value, -20.0, 20.0))
    return float(1.0 / (1.0 + np.exp(-clipped)))


def _collect_marker_score_evidence(
    adata: AnnData,
    cluster_key: str,
    annotation_key: str,
) -> Dict[str, Dict[str, Any]]:
    """Summarize marker-score evidence per cluster when score columns exist."""
    score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
    if (
        not score_cols
        or cluster_key not in adata.obs.columns
        or annotation_key not in adata.obs.columns
    ):
        return {}

    means = adata.obs.groupby(cluster_key)[score_cols].mean()
    evidence: Dict[str, Dict[str, Any]] = {}
    for cluster in means.index:
        cluster_str = str(cluster)
        assigned_labels = adata.obs.loc[
            adata.obs[cluster_key].astype(str) == cluster_str, annotation_key
        ].astype(str)
        if assigned_labels.empty:
            continue
        assigned_label = assigned_labels.value_counts().index[0]
        assigned_score_col = f"{assigned_label}_score"
        row = means.loc[cluster]

        top_score_col = row.idxmax()
        top_label = top_score_col[:-6] if top_score_col.endswith("_score") else top_score_col
        top_score = float(row[top_score_col])
        assigned_score = (
            float(row[assigned_score_col]) if assigned_score_col in row.index else np.nan
        )

        if assigned_score_col in row.index:
            best_other = row.drop(labels=[assigned_score_col], errors="ignore").max()
            best_other = float(best_other) if pd.notna(best_other) else 0.0
            margin = assigned_score - best_other
            confidence = 0.7 * _sigmoid_score(assigned_score) + 0.3 * _sigmoid_score(margin)
        else:
            margin = np.nan
            confidence = np.nan

        evidence[cluster_str] = {
            "assigned_label": assigned_label,
            "top_marker_label": top_label,
            "top_marker_score": top_score,
            "assigned_marker_score": assigned_score,
            "marker_margin": float(margin) if pd.notna(margin) else np.nan,
            "marker_confidence": float(confidence) if pd.notna(confidence) else np.nan,
        }
    return evidence


def _build_annotation_evidence_summary(
    adata: AnnData,
    cluster_key: str,
    annotation_key: str,
    *,
    params: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """Build a unified cluster-level evidence summary for annotation outputs."""
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")
    if annotation_key not in adata.obs.columns:
        raise KeyError(f"'{annotation_key}' not found in adata.obs.")

    cluster_series = adata.obs[cluster_key].astype(str)
    annotation_series = adata.obs[annotation_key].astype(str)
    marker_evidence = _collect_marker_score_evidence(adata, cluster_key, annotation_key)

    celltypist_mapping: Dict[str, str] = {}
    celltypist_stats: Dict[str, Dict[str, float]] = {}
    if (
        "celltypist_majority_voting" in adata.obs.columns
        or "celltypist_predicted_labels" in adata.obs.columns
    ):
        try:
            celltypist_mapping, celltypist_stats = _build_celltypist_cluster_mapping(
                adata,
                cluster_key,
                min_confidence=0.0,
            )
        except Exception:
            celltypist_mapping, celltypist_stats = {}, {}

    hybrid_audit = {}
    if params is not None:
        hybrid_audit = {str(k): v for k, v in params.get("hybrid_audit", {}).items()}

    rows = []
    for cluster in sorted(cluster_series.unique(), key=str):
        mask = cluster_series == str(cluster)
        cluster_labels = annotation_series.loc[mask]
        label_counts = cluster_labels.value_counts()
        assigned_label = str(label_counts.index[0]) if not label_counts.empty else "Unknown"
        label_purity = (
            float(label_counts.iloc[0] / max(1, label_counts.sum()))
            if not label_counts.empty
            else 0.0
        )

        cluster_marker = marker_evidence.get(str(cluster), {})
        celltypist_label = celltypist_mapping.get(str(cluster), "Unknown")
        ct_stats = celltypist_stats.get(str(cluster), {})
        ct_mean_conf = float(ct_stats.get("mean_conf_score", np.nan))
        ct_majority_fraction = float(ct_stats.get("majority_fraction", np.nan))
        ct_agreement = (
            float(celltypist_label == assigned_label) if celltypist_label != "Unknown" else np.nan
        )

        confidence_parts = [label_purity]
        marker_conf = cluster_marker.get("marker_confidence", np.nan)
        if pd.notna(marker_conf):
            confidence_parts.append(float(marker_conf))

        if pd.notna(ct_mean_conf) or pd.notna(ct_majority_fraction):
            ct_component_parts = [
                value for value in [ct_mean_conf, ct_majority_fraction] if pd.notna(value)
            ]
            ct_component = float(np.mean(ct_component_parts)) if ct_component_parts else np.nan
            if pd.notna(ct_component):
                if pd.notna(ct_agreement):
                    ct_component = 0.7 * ct_component + 0.3 * float(ct_agreement)
                confidence_parts.append(float(ct_component))

        annotation_confidence = float(np.mean(confidence_parts)) if confidence_parts else 0.0
        evidence_sources = ["cluster_purity"]
        if pd.notna(marker_conf):
            evidence_sources.append("marker_scores")
        if pd.notna(ct_mean_conf) or pd.notna(ct_majority_fraction):
            evidence_sources.append("celltypist")

        hybrid_info = hybrid_audit.get(str(cluster), {})
        rows.append(
            {
                "cluster": str(cluster),
                "assigned_label": assigned_label,
                "n_cells": int(mask.sum()),
                "label_purity": label_purity,
                "marker_label": cluster_marker.get("top_marker_label"),
                "assigned_marker_score": cluster_marker.get("assigned_marker_score"),
                "marker_margin": cluster_marker.get("marker_margin"),
                "marker_confidence": marker_conf,
                "celltypist_label": None if celltypist_label == "Unknown" else celltypist_label,
                "celltypist_mean_confidence": ct_mean_conf,
                "celltypist_majority_fraction": ct_majority_fraction,
                "celltypist_agreement": ct_agreement,
                "hybrid_decision": hybrid_info.get("decision"),
                "annotation_confidence": annotation_confidence,
                "evidence_sources": ",".join(evidence_sources),
            }
        )

    return pd.DataFrame(rows)


def _store_annotation_evidence_summary(
    adata: AnnData,
    cluster_key: str,
    annotation_key: str,
    *,
    params: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    """Persist a standardized annotation evidence summary and cluster confidence."""
    summary_df = _build_annotation_evidence_summary(
        adata,
        cluster_key,
        annotation_key,
        params=params,
    )
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[f"{annotation_key}_evidence"] = summary_df
    annotation_ns[f"{annotation_key}_evidence_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "annotation_key": annotation_key,
            "n_clusters": int(summary_df.shape[0]),
            "columns": summary_df.columns.tolist(),
        }
    )

    cluster_conf_map = summary_df.set_index("cluster")["annotation_confidence"].to_dict()
    adata.obs[f"{annotation_key}_cluster_confidence"] = (
        adata.obs[cluster_key].astype(str).map(cluster_conf_map).astype(float)
    )
    if f"{annotation_key}_confidence" not in adata.obs.columns:
        adata.obs[f"{annotation_key}_confidence"] = adata.obs[
            f"{annotation_key}_cluster_confidence"
        ]
    return summary_df


def run_lineage_state_annotation(
    adata: AnnData,
    config: AnnotationConfig,
) -> AnnData:
    """
    Hierarchical annotation workflow:
    lineage -> optional subtype -> optional state/program -> modular display label.
    """
    use_raw = adata.raw is not None
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )

    lineage_mgr = _resolve_annotation_manager(
        species=config.marker_species,
        tissue=config.marker_tissue,
        states=config.marker_states,
        marker_config=config.lineage_marker_config,
    )
    lineage_mgr.intersect_with(adata.raw if use_raw else adata)

    if config.run_scoring:
        adata = score_cell_types(
            adata,
            marker_config=lineage_mgr,
            use_raw=use_raw,
            layer=None,
            score_name_suffix="_score",
        )

    adata = annotate_clusters(
        adata,
        cluster_key=config.cluster_key,
        marker_config=lineage_mgr,
        method=config.marker_method,
        key_added=config.lineage_key,
        min_score=config.min_confidence,
        use_raw=use_raw,
        plot=False,
    )
    _rename_score_columns_for_manager(adata, lineage_mgr, "_lineage_module")

    subtype_values = pd.Series("Not_applicable", index=adata.obs_names, dtype="object")
    target_mask = _label_matches_target_lineage(
        adata.obs[config.lineage_key], config.target_lineage
    )

    if config.subtype_marker_config:
        subtype_mgr = _resolve_annotation_manager(
            species=config.marker_species,
            tissue=config.marker_tissue,
            marker_config=config.subtype_marker_config,
        )
        subtype_mgr.intersect_with(adata.raw if use_raw else adata)
        if config.run_scoring:
            adata = score_cell_types(
                adata,
                marker_config=subtype_mgr,
                use_raw=use_raw,
                layer=None,
                score_name_suffix="_score",
            )
        adata = annotate_clusters(
            adata,
            cluster_key=config.cluster_key,
            marker_config=subtype_mgr,
            method=config.marker_method,
            key_added=config.subtype_key,
            min_score=config.min_confidence,
            use_raw=use_raw,
            plot=False,
        )
        subtype_values = adata.obs[config.subtype_key].astype(str)
        subtype_values = subtype_values.where(target_mask, "Not_applicable")
        adata.obs[config.subtype_key] = pd.Categorical(subtype_values)
        _rename_score_columns_for_manager(adata, subtype_mgr, "_subtype_module")
    else:
        adata.obs[config.subtype_key] = pd.Categorical(subtype_values)

    state_values = pd.Series("Not_applicable", index=adata.obs_names, dtype="object")
    state_confidence = pd.Series(np.nan, index=adata.obs_names, dtype=float)
    signatures, state_metadata = _collect_state_signatures(config)

    if signatures:
        adata = score_by_gene_sets(
            adata,
            signatures,
            use_raw=use_raw,
            layer=None,
            score_name_suffix=config.state_score_suffix,
            preserve_missing=True,
            min_genes_required=1,
        )
        state_score_cols = [
            f"{name}{config.state_score_suffix}"
            for name in signatures
            if f"{name}{config.state_score_suffix}" in adata.obs.columns
        ]
        if state_score_cols:
            state_score_df = adata.obs[state_score_cols].copy()
            lineage_labels = adata.obs[config.lineage_key].astype(str)
            subtype_labels = adata.obs[config.subtype_key].astype(str)

            scoped_score_df = state_score_df.copy()
            for state_name in signatures:
                score_col = f"{state_name}{config.state_score_suffix}"
                if score_col not in scoped_score_df.columns:
                    continue
                allowed_mask = [
                    _state_applies_to_cell(
                        state_metadata.get(state_name),
                        lineage_label=lin,
                        subtype_label=sub,
                    )
                    for lin, sub in zip(lineage_labels, subtype_labels)
                ]
                scoped_score_df.loc[
                    ~pd.Series(allowed_mask, index=scoped_score_df.index), score_col
                ] = -np.inf

            winner_cols = scoped_score_df.idxmax(axis=1)
            state_values = winner_cols.str.replace(config.state_score_suffix, "", regex=False)

            top1 = scoped_score_df.max(axis=1)
            top2 = scoped_score_df.apply(
                lambda row: (
                    row.replace(-np.inf, np.nan).nlargest(2).iloc[-1]
                    if row.replace(-np.inf, np.nan).dropna().shape[0] > 1
                    else row.replace(-np.inf, np.nan).fillna(0).iloc[0]
                ),
                axis=1,
            )
            state_confidence = (top1 - top2).astype(float)
            state_values = state_values.where(
                top1.replace(-np.inf, np.nan).notna(), "Not_applicable"
            )
            state_confidence = state_confidence.where(top1.replace(-np.inf, np.nan).notna(), np.nan)
            state_values = state_values.where(target_mask, "Not_applicable")
            state_confidence = state_confidence.where(target_mask, np.nan)

    adata.obs[config.state_key] = pd.Categorical(state_values)
    adata.obs[f"{config.state_key}_confidence"] = state_confidence

    if config.nomenclature_style == "modular":
        final_labels = _build_modular_annotation_label(
            adata.obs[config.lineage_key],
            adata.obs[config.subtype_key],
            adata.obs[config.state_key],
        )
    else:
        final_labels = adata.obs[config.subtype_key].astype(str)
        final_labels = final_labels.where(
            ~final_labels.isin(["Unknown", "Not_applicable"]),
            adata.obs[config.lineage_key].astype(str),
        )

    adata.obs[config.key_added] = pd.Categorical(final_labels)
    adata.obs[config.annotation_basis_key] = pd.Categorical(
        np.where(
            adata.obs[config.state_key].astype(str).isin(["Unknown", "Not_applicable"]),
            "lineage+subtype",
            "lineage+subtype+state",
        )
    )

    annotation_ns[f"{config.key_added}_hierarchical_params"] = sanitize_for_hdf5(
        {
            "cluster_key": config.cluster_key,
            "lineage_key": config.lineage_key,
            "subtype_key": config.subtype_key,
            "state_key": config.state_key,
            "target_lineage": config.target_lineage,
            "lineage_marker_config": config.lineage_marker_config,
            "subtype_marker_config": config.subtype_marker_config,
            "marker_states": config.marker_states,
            "state_marker_config": config.state_marker_config,
            "state_signature_names": config.state_signature_names,
            "state_signature_categories": config.state_signature_categories,
            "nomenclature_style": config.nomenclature_style,
        }
    )
    annotation_ns[f"{config.key_added}_state_signatures"] = sanitize_for_hdf5(
        {name: len(genes) for name, genes in signatures.items()}
    )
    annotation_ns[f"{config.key_added}_state_metadata"] = sanitize_for_hdf5(state_metadata)

    return adata


def evaluate_annotation(
    adata: AnnData,
    cluster_key: str,
    annotation_key: str,
    marker_config: Union[str, Manager],
    plot: bool = True,
    save_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Evaluate annotation quality: marker coverage, specificity, confidence.

    Enhancements:
    - Guard against empty DE results or missing columns.
    - Stable plotting and parameter trace saved to .uns.
    """
    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError("marker_config must be a file path or Manager instance.")
    if not mgr.CELLS:
        log.warning("Marker manager has no cell types. Evaluation may be uninformative.")

    # Ensure cluster_key exists/categorical
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")
    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    mgr.intersect_with(adata)
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        key_added=f"rank_genes_{cluster_key}",
    )
    results = []

    categories = list(adata.obs[cluster_key].cat.categories)
    for cluster in categories:
        de_genes = sc.get.rank_genes_groups_df(
            adata, key=f"rank_genes_{cluster_key}", group=cluster
        )
        if de_genes.empty or "pvals_adj" not in de_genes.columns or "names" not in de_genes.columns:
            log.warning(f"No DE genes or required columns missing for cluster '{cluster}'.")
            continue
        sig_genes = set(de_genes.loc[de_genes["pvals_adj"] < 0.05, "names"].astype(str))

        # Determine assigned type
        cluster_mask = adata.obs[cluster_key] == cluster
        if not cluster_mask.any() or annotation_key not in adata.obs.columns:
            continue
        assigned_series = adata.obs.loc[cluster_mask, annotation_key].astype(str)
        if assigned_series.empty:
            continue
        assigned_type = assigned_series.value_counts().index[0]
        if assigned_type == "Unknown" or assigned_type not in mgr.CELLS:
            continue

        expected_markers = set(mgr[assigned_type].markers)
        found_markers = sig_genes & expected_markers
        marker_coverage = len(found_markers) / len(expected_markers) if expected_markers else 0.0

        all_other_markers = set(
            m for t, c in mgr.CELLS.items() if t != assigned_type for m in c.markers
        )
        specificity = (
            1.0 - (len(found_markers & all_other_markers) / len(found_markers))
            if found_markers
            else 0.0
        )

        confidence = 0.6 * marker_coverage + 0.4 * specificity
        results.append(
            {
                "cluster": cluster,
                "cell_type": assigned_type,
                "marker_coverage": marker_coverage,
                "marker_specificity": specificity,
                "annotation_confidence": confidence,
                "found_markers": ", ".join(sorted(found_markers)) if found_markers else "",
                "expected_markers": len(expected_markers),
                "detected_markers": len(found_markers),
            }
        )

    results_df = pd.DataFrame(results)
    if plot:
        if not results_df.empty:
            results_df = results_df.sort_values("annotation_confidence")
            plt.figure(figsize=(12, 8))
            plt.barh(
                results_df["cluster"].astype(str) + " (" + results_df["cell_type"] + ")",
                results_df["annotation_confidence"],
                color="skyblue",
            )
            plt.xlabel("Annotation Confidence Score")
            plt.ylabel("Cluster (Cell Type)")
            plt.title("Annotation Confidence by Cluster")
            plt.xlim(0, 1.0)
            plt.tight_layout()
            if save_path:
                plt.savefig(f"{save_path}_confidence.png", dpi=300)
            plt.show()
        else:
            log.info("No evaluation results to plot (empty results).")

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    adata.uns["sclucid"]["analysis"]["annotation"][f"{annotation_key}_evaluation"] = results_df
    # Save params
    eval_params = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "annotation_key": annotation_key,
            "scanpy_version": version("scanpy"),
            "n_types_in_manager": len(getattr(mgr, "CELLS", {})),
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"][
        f"{annotation_key}_evaluation_params"
    ] = eval_params
    return results_df


def run_annotation(
    adata: AnnData,
    config: AnnotationConfig,
) -> AnnData:
    """
    Full annotation workflow: scoring -> auto annotation -> results in .obs/.uns.
    """
    if config.final_method == "hierarchical":
        adata = run_lineage_state_annotation(adata, config)
        _store_annotation_evidence_summary(
            adata,
            config.cluster_key,
            config.key_added,
            params=(
                adata.uns.get("sclucid", {})
                .get("analysis", {})
                .get("annotation", {})
                .get(f"{config.key_added}_hierarchical_params", {})
            ),
        )
        adata.obs["cell_compartment"] = _map_compartments(
            adata.obs[config.key_added],
            compartment_map=config.compartment_map,
        )
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
        adata.uns["sclucid"]["analysis"]["annotation"]["workflow_config"] = sanitize_for_hdf5(
            config.to_dict()
        )
        return adata

    use_raw = adata.raw is not None
    mgr = get_marker_manager(
        species=config.marker_species,
        tissue=config.marker_tissue,
        states=config.marker_states,
    )
    mgr.intersect_with(adata.raw if use_raw else adata)

    if config.run_celltypist or config.final_method in {"celltypist", "hybrid"}:
        adata = run_celltypist(adata, model=config.celltypist_model)

    if config.run_scoring and config.final_method != "celltypist":
        # Use raw for scoring by default
        adata = score_cell_types(
            adata, marker_config=mgr, use_raw=use_raw, layer=None if use_raw else None
        )

    params_for_summary: Dict[str, Any] = {}
    if config.final_method in {"max_score", "enrichment", "combined"}:
        adata = annotate_clusters(
            adata,
            cluster_key=config.cluster_key,
            marker_config=mgr,
            method=config.final_method,
            key_added=config.key_added,
            min_score=config.min_confidence,
            use_raw=use_raw,
            plot=config.plot,
        )
        params_for_summary = (
            adata.uns.get("sclucid", {})
            .get("analysis", {})
            .get("annotation", {})
            .get(f"{config.key_added}_params", {})
        )
    elif config.final_method == "celltypist":
        labels, confidences, label_key = _get_celltypist_label_series(adata)
        final_labels = labels.where(
            confidences >= config.celltypist_confidence_threshold, "Unknown"
        )
        adata.obs[config.key_added] = pd.Categorical(final_labels)
        adata.obs[f"{config.key_added}_confidence"] = confidences.astype(float)
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
        adata.uns["sclucid"]["analysis"]["annotation"][f"{config.key_added}_params"] = (
            sanitize_for_hdf5(
                {
                    "method": "celltypist",
                    "label_source": label_key,
                    "confidence_threshold": config.celltypist_confidence_threshold,
                }
            )
        )
        params_for_summary = (
            adata.uns.get("sclucid", {})
            .get("analysis", {})
            .get("annotation", {})
            .get(f"{config.key_added}_params", {})
        )
    elif config.final_method == "hybrid":
        marker_key = f"{config.key_added}_marker"
        adata = annotate_clusters(
            adata,
            cluster_key=config.cluster_key,
            marker_config=mgr,
            method=config.marker_method,
            key_added=marker_key,
            min_score=config.min_confidence,
            use_raw=use_raw,
            plot=False,
        )
        marker_params = (
            adata.uns.get("sclucid", {})
            .get("analysis", {})
            .get("annotation", {})
            .get(f"{marker_key}_params", {})
        )
        marker_mapping = {str(k): str(v) for k, v in marker_params.get("mapping", {}).items()}
        celltypist_mapping, celltypist_stats = _build_celltypist_cluster_mapping(
            adata,
            config.cluster_key,
            min_confidence=config.celltypist_confidence_threshold,
        )
        final_mapping, hybrid_audit = _combine_marker_and_celltypist_mappings(
            marker_mapping,
            celltypist_mapping,
            celltypist_stats,
            min_celltypist_confidence=config.celltypist_confidence_threshold,
        )
        cluster_codes = adata.obs[config.cluster_key].astype(str)
        adata.obs[config.key_added] = pd.Categorical(
            cluster_codes.map(final_mapping).fillna("Unknown")
        )
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
        adata.uns["sclucid"]["analysis"]["annotation"][f"{config.key_added}_params"] = (
            sanitize_for_hdf5(
                {
                    "method": "hybrid",
                    "marker_method": config.marker_method,
                    "confidence_threshold": config.celltypist_confidence_threshold,
                    "mapping": final_mapping,
                    "hybrid_audit": hybrid_audit,
                }
            )
        )
        params_for_summary = (
            adata.uns.get("sclucid", {})
            .get("analysis", {})
            .get("annotation", {})
            .get(f"{config.key_added}_params", {})
        )
        if config.plot:
            if "X_umap" not in adata.obsm:
                sc.tl.umap(adata)
            sc.pl.umap(adata, color=[config.cluster_key, marker_key, config.key_added], wspace=0.4)
    else:
        raise ValueError(f"Unknown annotation method: {config.final_method}")

    _store_annotation_evidence_summary(
        adata,
        config.cluster_key,
        config.key_added,
        params=params_for_summary,
    )

    # Map annotated labels to broad compartments for tumor-aware analysis
    adata.obs["cell_compartment"] = _map_compartments(
        adata.obs[config.key_added],
        compartment_map=config.compartment_map,
    )

    if config.report and config.save_dir:
        from ..plotting import export_annotation_report

        report_path = Path(config.save_dir) / f"{config.key_added}_annotation_report.png"
        try:
            export_annotation_report(
                adata,
                annotation_key=config.key_added,
                cluster_key=config.cluster_key,
                save=str(report_path),
                export_formats=("png", "pdf"),
                write_sidecars=True,
                show=False,
            )
        except Exception as exc:
            log.warning(f"Failed to export annotation report: {exc}")

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    adata.uns["sclucid"]["analysis"]["annotation"]["workflow_config"] = sanitize_for_hdf5(
        config.to_dict()
    )
    return adata
