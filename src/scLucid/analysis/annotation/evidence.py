"""Evidence-table helpers for annotation review and final label application."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ...utils import Manager, sanitize_for_hdf5
from .utils import _classify_annotation_marker, _map_compartments

__all__ = [
    "ANNOTATION_REVIEW_SCHEMA",
    "ANALYSIS_REVIEW_SUMMARY_SCHEMA",
    "standardize_cluster_marker_table",
    "run_marker_annotation_evidence",
    "build_llm_annotation_bundle",
    "merge_annotation_evidence",
    "apply_final_annotation",
]

ANNOTATION_REVIEW_SCHEMA = [
    "cluster",
    "n_cells",
    "pct_cells",
    "reference_label",
    "reference_confidence",
    "marker_label",
    "marker_confidence",
    "llm_label",
    "llm_confidence",
    "final_label",
    "annotation_confidence",
    "decision",
    "conflicts",
    "warnings",
    "needs_review",
    "top_markers",
    "top_terms",
]

ANALYSIS_REVIEW_SUMMARY_SCHEMA = [
    "module",
    "workflow_name",
    "steps_executed",
    "clustering",
    "annotation",
    "warnings",
    "artifacts",
]


def _format_top_distribution(values: pd.Series, n: int = 3) -> str:
    """Return a compact top-category distribution string."""
    if values.empty:
        return ""
    counts = values.astype(str).value_counts(normalize=True).head(n)
    return ", ".join(f"{name}:{frac:.2f}" for name, frac in counts.items())


def standardize_cluster_marker_table(
    markers_df: pd.DataFrame,
    *,
    cluster_col: str = "group",
    gene_col: str = "names",
    score_col: Optional[str] = "scores",
    logfc_col: Optional[str] = "logfoldchanges",
    padj_col: Optional[str] = "pvals_adj",
    pct_in_col: Optional[str] = "pct_nz_group",
    pct_out_col: Optional[str] = "pct_nz_reference",
    keep_top_n_per_cluster: Optional[int] = None,
    drop_noise: bool = False,
    key_added: Optional[str] = None,
    adata: Optional[AnnData] = None,
) -> pd.DataFrame:
    """
    Normalize marker-DE output into the annotation evidence schema.

    The standardized table preserves the original DE statistics when present
    and adds `noise_category`, `is_annotation_informative`, and `marker_rank`.
    """
    if cluster_col not in markers_df.columns:
        raise KeyError(f"'{cluster_col}' not found in markers_df.")
    if gene_col not in markers_df.columns:
        raise KeyError(f"'{gene_col}' not found in markers_df.")

    standardized = pd.DataFrame(
        {
            "cluster": markers_df[cluster_col].astype(str),
            "gene": markers_df[gene_col].astype(str),
        }
    )
    optional_cols = {
        "score": score_col,
        "logfoldchanges": logfc_col,
        "pvals_adj": padj_col,
        "pct_in": pct_in_col,
        "pct_out": pct_out_col,
    }
    for target, source in optional_cols.items():
        if source and source in markers_df.columns:
            standardized[target] = pd.to_numeric(markers_df[source], errors="coerce")
        else:
            standardized[target] = np.nan

    sort_cols = ["cluster"]
    ascending = [True]
    if "score" in standardized.columns and standardized["score"].notna().any():
        sort_cols.append("score")
        ascending.append(False)
    elif standardized["logfoldchanges"].notna().any():
        sort_cols.append("logfoldchanges")
        ascending.append(False)
    standardized = standardized.sort_values(sort_cols, ascending=ascending).copy()
    standardized["marker_rank"] = standardized.groupby("cluster").cumcount() + 1
    standardized["noise_category"] = standardized["gene"].map(_classify_annotation_marker)
    standardized["is_annotation_informative"] = standardized["noise_category"].isna()

    if drop_noise:
        standardized = standardized[standardized["is_annotation_informative"]].copy()
    if keep_top_n_per_cluster is not None:
        standardized = (
            standardized.groupby("cluster", group_keys=False)
            .head(int(keep_top_n_per_cluster))
            .copy()
        )

    standardized = standardized.reset_index(drop=True)
    if adata is not None:
        annotation_ns = (
            adata.uns.setdefault("sclucid", {})
            .setdefault("analysis", {})
            .setdefault("annotation", {})
        )
        target_key = key_added or "cluster_marker_table"
        annotation_ns[target_key] = standardized
        annotation_ns[f"{target_key}_params"] = sanitize_for_hdf5(
            {
                "cluster_col": cluster_col,
                "gene_col": gene_col,
                "score_col": score_col,
                "logfc_col": logfc_col,
                "padj_col": padj_col,
                "pct_in_col": pct_in_col,
                "pct_out_col": pct_out_col,
                "keep_top_n_per_cluster": keep_top_n_per_cluster,
                "drop_noise": drop_noise,
            }
        )
    return standardized


def _resolve_marker_cell_metadata(cell: Any) -> Dict[str, Any]:
    """Return marker manager metadata when available."""
    metadata = getattr(cell, "metadata", {}) or {}
    return dict(metadata) if isinstance(metadata, dict) else {}


def run_marker_annotation_evidence(
    adata: AnnData,
    cluster_key: str,
    marker_config: Union[str, Manager],
    *,
    markers_df: Optional[pd.DataFrame] = None,
    marker_gene_col: str = "names",
    marker_group_col: str = "group",
    top_n_markers: int = 50,
    min_overlap: int = 1,
    min_confidence: float = 0.05,
    use_raw: bool = True,
    key_added: str = "marker_annotation_evidence",
) -> pd.DataFrame:
    """
    Score cluster labels from marker-manager overlap evidence.

    This function is intentionally cluster-level and auditable: it reports the
    winning label, confidence, matched genes, and runner-up for each cluster.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")
    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError("marker_config must be a file path or Manager instance.")
    mgr.intersect_with(adata.raw if use_raw and adata.raw is not None else adata)

    if markers_df is None:
        rank_key = f"rank_genes_{cluster_key}"
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw and adata.raw is not None,
            pts=True,
            key_added=rank_key,
        )
        markers_df = sc.get.rank_genes_groups_df(adata, key=rank_key, group=None)

    marker_table = standardize_cluster_marker_table(
        markers_df,
        cluster_col=marker_group_col,
        gene_col=marker_gene_col,
        keep_top_n_per_cluster=top_n_markers,
        drop_noise=True,
    )

    rows: List[Dict[str, Any]] = []
    cluster_series = adata.obs[cluster_key].astype(str)
    for cluster in cluster_series.drop_duplicates().tolist():
        cluster_genes = marker_table.loc[marker_table["cluster"] == str(cluster), "gene"].tolist()
        cluster_gene_set = {gene.upper() for gene in cluster_genes}
        scored_labels: List[Dict[str, Any]] = []
        for label, cell in mgr.CELLS.items():
            markers = [str(g) for g in getattr(cell, "markers", [])]
            marker_set = {gene.upper() for gene in markers}
            if not marker_set:
                continue
            matched = [gene for gene in cluster_genes if gene.upper() in marker_set]
            overlap = len(matched)
            recall = overlap / max(1, len(marker_set))
            precision = overlap / max(1, len(cluster_gene_set))
            confidence = 0.65 * recall + 0.35 * precision
            if overlap >= min_overlap:
                scored_labels.append(
                    {
                        "label": str(label),
                        "overlap": int(overlap),
                        "recall": float(recall),
                        "precision": float(precision),
                        "confidence": float(confidence),
                        "matched_markers": matched,
                        "metadata": _resolve_marker_cell_metadata(cell),
                    }
                )

        scored_labels = sorted(
            scored_labels,
            key=lambda row: (row["confidence"], row["overlap"], row["label"]),
            reverse=True,
        )
        winner = scored_labels[0] if scored_labels else None
        runner_up = scored_labels[1] if len(scored_labels) > 1 else None
        if winner is None or float(winner["confidence"]) < min_confidence:
            marker_label = "Unknown"
            confidence = 0.0
            matched_markers = ""
        else:
            marker_label = str(winner["label"])
            confidence = float(winner["confidence"])
            matched_markers = ", ".join(winner["matched_markers"][:12])

        rows.append(
            {
                "cluster": str(cluster),
                "n_cells": int((cluster_series == str(cluster)).sum()),
                "marker_label": marker_label,
                "marker_confidence": confidence,
                "marker_overlap": int(winner["overlap"]) if winner else 0,
                "marker_recall": float(winner["recall"]) if winner else 0.0,
                "marker_precision": float(winner["precision"]) if winner else 0.0,
                "matched_markers": matched_markers,
                "runner_up_marker_label": runner_up["label"] if runner_up else None,
                "runner_up_marker_confidence": (
                    float(runner_up["confidence"]) if runner_up else np.nan
                ),
                "top_informative_markers": ", ".join(cluster_genes[:12]),
            }
        )

    evidence_df = pd.DataFrame(rows)
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[key_added] = evidence_df
    annotation_ns[f"{key_added}_marker_table"] = marker_table
    annotation_ns[f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "top_n_markers": int(top_n_markers),
            "min_overlap": int(min_overlap),
            "min_confidence": float(min_confidence),
            "n_marker_labels": len(getattr(mgr, "CELLS", {})),
        }
    )
    return evidence_df


def build_llm_annotation_bundle(
    adata: AnnData,
    cluster_key: str,
    *,
    markers_df: Optional[pd.DataFrame] = None,
    enrichment_dict: Optional[Dict[str, pd.DataFrame]] = None,
    marker_evidence: Optional[pd.DataFrame] = None,
    reference_key: Optional[str] = None,
    top_n_markers: int = 15,
    top_n_terms: int = 5,
    sample_col: Optional[str] = None,
    group_col: Optional[str] = None,
    key_added: str = "llm_annotation_bundle",
) -> Dict[str, Any]:
    """
    Build a compact, auditable input bundle for data-driven LLM annotation.

    The bundle is stored in `.uns` and deliberately does not call an LLM or
    mutate final labels.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")

    marker_table = None
    if markers_df is not None:
        gene_col = "gene" if "gene" in markers_df.columns else "names"
        cluster_col = "cluster" if "cluster" in markers_df.columns else "group"
        marker_table = standardize_cluster_marker_table(
            markers_df,
            cluster_col=cluster_col,
            gene_col=gene_col,
            keep_top_n_per_cluster=top_n_markers,
            drop_noise=False,
        )
    marker_evidence_by_cluster = (
        marker_evidence.set_index("cluster").to_dict(orient="index")
        if isinstance(marker_evidence, pd.DataFrame)
        and not marker_evidence.empty
        and "cluster" in marker_evidence.columns
        else {}
    )

    cluster_series = adata.obs[cluster_key].astype(str)
    clusters: Dict[str, Any] = {}
    for cluster in cluster_series.drop_duplicates().tolist():
        mask = cluster_series == str(cluster)
        obs_subset = adata.obs.loc[mask]

        marker_rows = (
            marker_table[marker_table["cluster"] == str(cluster)].copy()
            if marker_table is not None
            else pd.DataFrame()
        )
        informative = (
            marker_rows[marker_rows["is_annotation_informative"]]
            if not marker_rows.empty
            else pd.DataFrame()
        )
        noisy = (
            marker_rows[~marker_rows["is_annotation_informative"]]
            if not marker_rows.empty
            else pd.DataFrame()
        )

        terms: List[str] = []
        if enrichment_dict:
            term_df = enrichment_dict.get(cluster) or enrichment_dict.get(str(cluster))
            if isinstance(term_df, pd.DataFrame) and not term_df.empty:
                sort_col = "Adjusted P-value" if "Adjusted P-value" in term_df.columns else None
                term_df = term_df.sort_values(sort_col) if sort_col else term_df
                if "Term" in term_df.columns:
                    terms = term_df["Term"].astype(str).head(top_n_terms).tolist()

        reference_summary = None
        if reference_key and reference_key in obs_subset.columns:
            ref_counts = obs_subset[reference_key].astype(str).value_counts(normalize=True).head(5)
            reference_summary = [
                {"label": str(label), "fraction": float(frac)} for label, frac in ref_counts.items()
            ]

        clusters[str(cluster)] = {
            "n_cells": int(mask.sum()),
            "pct_cells": float(mask.mean()),
            "top_informative_markers": (
                informative["gene"].astype(str).head(top_n_markers).tolist()
                if not informative.empty
                else []
            ),
            "top_noisy_markers": (
                noisy["gene"].astype(str).head(8).tolist() if not noisy.empty else []
            ),
            "top_pathways": terms,
            "marker_manager_evidence": marker_evidence_by_cluster.get(str(cluster), {}),
            "reference_annotation": reference_summary,
            "top_samples": (
                _format_top_distribution(obs_subset[sample_col])
                if sample_col and sample_col in obs_subset.columns
                else ""
            ),
            "group_distribution": (
                _format_top_distribution(obs_subset[group_col])
                if group_col and group_col in obs_subset.columns
                else ""
            ),
        }

    instructions = (
        "Assign one conservative cell lineage label per cluster. Prefer broad lineage labels "
        "when evidence is mixed. Do not use ribosomal, mitochondrial, stress, or housekeeping "
        "genes as primary lineage evidence. Return a table with cluster, llm_label, "
        "llm_confidence from 0 to 1, evidence, conflicts, and needs_review."
    )
    bundle = {
        "schema_version": "analysis_annotation_bundle_v1",
        "cluster_key": cluster_key,
        "instructions": instructions,
        "clusters": sanitize_for_hdf5(clusters),
    }
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[key_added] = bundle
    annotation_ns[f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "top_n_markers": int(top_n_markers),
            "top_n_terms": int(top_n_terms),
            "reference_key": reference_key,
            "sample_col": sample_col,
            "group_col": group_col,
        }
    )
    return bundle


def _cluster_label_from_obs(
    adata: AnnData,
    cluster_key: str,
    label_key: Optional[str],
    confidence_key: Optional[str] = None,
) -> pd.DataFrame:
    """Aggregate a cell-level label column into cluster-level evidence."""
    if not label_key or label_key not in adata.obs.columns:
        return pd.DataFrame(columns=["cluster", "label", "confidence"])
    cluster_series = adata.obs[cluster_key].astype(str)
    rows = []
    for cluster in cluster_series.drop_duplicates().tolist():
        mask = cluster_series == str(cluster)
        labels = adata.obs.loc[mask, label_key].astype(str)
        labels = labels[labels.notna() & ~labels.isin(["nan", "None"])]
        if labels.empty:
            rows.append({"cluster": str(cluster), "label": "Unknown", "confidence": 0.0})
            continue
        counts = labels.value_counts()
        label = str(counts.index[0])
        majority = float(counts.iloc[0] / max(1, counts.sum()))
        if confidence_key and confidence_key in adata.obs.columns:
            conf = pd.to_numeric(adata.obs.loc[mask, confidence_key], errors="coerce").mean()
            confidence = float(np.nanmean([majority, conf])) if pd.notna(conf) else majority
        else:
            confidence = majority
        rows.append({"cluster": str(cluster), "label": label, "confidence": confidence})
    return pd.DataFrame(rows)


def _coerce_llm_evidence(
    llm_annotations: Optional[Union[pd.DataFrame, Dict[str, Any]]],
) -> pd.DataFrame:
    """Normalize optional LLM annotations into cluster-level evidence columns."""
    if llm_annotations is None:
        return pd.DataFrame(columns=["cluster", "llm_label", "llm_confidence"])
    if isinstance(llm_annotations, pd.DataFrame):
        df = llm_annotations.copy()
    elif isinstance(llm_annotations, dict):
        rows = []
        for cluster, value in llm_annotations.items():
            if isinstance(value, dict):
                rows.append({"cluster": str(cluster), **value})
            else:
                rows.append({"cluster": str(cluster), "llm_label": str(value)})
        df = pd.DataFrame(rows)
    else:
        raise TypeError("llm_annotations must be a DataFrame, dict, or None.")
    if "cluster" not in df.columns:
        raise KeyError("llm_annotations must contain a 'cluster' column.")
    if "llm_label" not in df.columns and "label" in df.columns:
        df = df.rename(columns={"label": "llm_label"})
    if "llm_confidence" not in df.columns and "confidence" in df.columns:
        df = df.rename(columns={"confidence": "llm_confidence"})
    if "llm_label" not in df.columns:
        df["llm_label"] = "Unknown"
    if "llm_confidence" not in df.columns:
        df["llm_confidence"] = np.nan
    df["cluster"] = df["cluster"].astype(str)
    return df


def merge_annotation_evidence(
    adata: AnnData,
    cluster_key: str,
    *,
    marker_evidence: Optional[pd.DataFrame] = None,
    reference_key: Optional[str] = None,
    reference_confidence_key: Optional[str] = None,
    llm_annotations: Optional[Union[pd.DataFrame, Dict[str, Any]]] = None,
    review_table: Optional[pd.DataFrame] = None,
    min_final_confidence: float = 0.2,
    prefer_llm_when_confident: bool = True,
    key_added: str = "annotation_review_table",
) -> pd.DataFrame:
    """
    Merge reference, marker-manager, and LLM evidence into final cluster labels.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")

    cluster_series = adata.obs[cluster_key].astype(str)
    base = pd.DataFrame(
        {
            "cluster": cluster_series.drop_duplicates().astype(str).tolist(),
        }
    )
    base["n_cells"] = base["cluster"].map(cluster_series.value_counts().to_dict()).astype(int)
    base["pct_cells"] = base["n_cells"] / max(1, adata.n_obs)

    reference_df = _cluster_label_from_obs(
        adata,
        cluster_key,
        reference_key,
        confidence_key=reference_confidence_key,
    ).rename(columns={"label": "reference_label", "confidence": "reference_confidence"})
    marker_df = (
        marker_evidence.copy()
        if isinstance(marker_evidence, pd.DataFrame) and not marker_evidence.empty
        else pd.DataFrame(columns=["cluster", "marker_label", "marker_confidence"])
    )
    llm_df = _coerce_llm_evidence(llm_annotations)

    merged = base.merge(reference_df, on="cluster", how="left")
    merged = merged.merge(marker_df, on="cluster", how="left")
    merged = merged.merge(llm_df, on="cluster", how="left")

    if review_table is not None and not review_table.empty and "cluster" in review_table.columns:
        optional_cols = [
            c for c in ["cluster", "top_markers", "top_terms"] if c in review_table.columns
        ]
        merged = merged.merge(review_table[optional_cols], on="cluster", how="left")
    else:
        merged["top_markers"] = ""
        merged["top_terms"] = ""

    for col in ["reference_label", "marker_label", "llm_label"]:
        if col not in merged.columns:
            merged[col] = "Unknown"
        merged[col] = merged[col].fillna("Unknown").astype(str)
    for col in ["reference_confidence", "marker_confidence", "llm_confidence"]:
        if col not in merged.columns:
            merged[col] = np.nan
        merged[col] = pd.to_numeric(merged[col], errors="coerce")

    final_rows = []
    for _, row in merged.iterrows():
        candidates = [
            ("marker", row["marker_label"], row["marker_confidence"]),
            ("reference", row["reference_label"], row["reference_confidence"]),
            ("llm", row["llm_label"], row["llm_confidence"]),
        ]
        usable = [
            (source, label, float(conf) if pd.notna(conf) else 0.0)
            for source, label, conf in candidates
            if label not in {"Unknown", "nan", "None", ""}
        ]
        votes: Dict[str, List[tuple[str, float]]] = {}
        for source, label, conf in usable:
            votes.setdefault(label, []).append((source, conf))

        conflicts = []
        final_label = "Unknown"
        decision = "insufficient_evidence"
        confidence = 0.0
        if votes:
            if len(votes) > 1:
                conflicts = sorted(votes)
            agreement = {
                label: {
                    "sources": [s for s, _ in entries],
                    "mean_conf": float(np.mean([c for _, c in entries])),
                    "n_sources": len(entries),
                }
                for label, entries in votes.items()
            }
            agreed = [(label, info) for label, info in agreement.items() if info["n_sources"] >= 2]
            if agreed:
                final_label, info = sorted(
                    agreed,
                    key=lambda item: (item[1]["n_sources"], item[1]["mean_conf"], item[0]),
                    reverse=True,
                )[0]
                confidence = min(1.0, 0.15 + 0.85 * float(info["mean_conf"]))
                decision = "multi_source_agreement"
            elif (
                prefer_llm_when_confident
                and row["llm_label"] not in {"Unknown", "nan", "None", ""}
                and pd.notna(row["llm_confidence"])
                and float(row["llm_confidence"]) >= 0.75
            ):
                final_label = row["llm_label"]
                confidence = float(row["llm_confidence"])
                decision = "llm_high_confidence"
            else:
                source, final_label, confidence = sorted(
                    usable,
                    key=lambda item: (item[2], item[0] == "marker", item[1]),
                    reverse=True,
                )[0]
                decision = f"{source}_best_available"

        if confidence < min_final_confidence:
            final_label = "Unknown"
            decision = "below_confidence_threshold"

        warnings = []
        if conflicts:
            warnings.append("evidence_conflict")
        if row.get("marker_label") == "Unknown" and row.get("reference_label") == "Unknown":
            warnings.append("weak_marker_reference_evidence")
        needs_review = bool(conflicts or final_label == "Unknown" or confidence < 0.5)
        final_rows.append(
            {
                "final_label": final_label,
                "annotation_confidence": float(confidence),
                "decision": decision,
                "conflicts": ",".join(conflicts),
                "warnings": ",".join(warnings),
                "needs_review": needs_review,
            }
        )

    final_df = pd.concat([merged.reset_index(drop=True), pd.DataFrame(final_rows)], axis=1)
    for col in ANNOTATION_REVIEW_SCHEMA:
        if col not in final_df.columns:
            final_df[col] = "" if col not in {"needs_review", "n_cells", "pct_cells"} else np.nan
    final_df = final_df[
        ANNOTATION_REVIEW_SCHEMA
        + [c for c in final_df.columns if c not in ANNOTATION_REVIEW_SCHEMA]
    ]

    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[key_added] = final_df
    annotation_ns[f"{key_added}_schema"] = list(ANNOTATION_REVIEW_SCHEMA)
    annotation_ns[f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "reference_key": reference_key,
            "reference_confidence_key": reference_confidence_key,
            "min_final_confidence": float(min_final_confidence),
            "prefer_llm_when_confident": bool(prefer_llm_when_confident),
        }
    )
    return final_df


def apply_final_annotation(
    adata: AnnData,
    cluster_key: str,
    annotation_review_table: Optional[pd.DataFrame] = None,
    *,
    label_col: str = "final_label",
    confidence_col: str = "annotation_confidence",
    status_col: str = "needs_review",
    key_added: str = "cell_type_final",
) -> AnnData:
    """
    Apply final cluster-level annotation labels back to cell-level `.obs`.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")
    if annotation_review_table is None:
        annotation_review_table = (
            adata.uns.get("sclucid", {})
            .get("analysis", {})
            .get("annotation", {})
            .get("annotation_review_table")
        )
    if not isinstance(annotation_review_table, pd.DataFrame):
        raise ValueError("annotation_review_table must be provided or stored in adata.uns.")
    if "cluster" not in annotation_review_table.columns:
        raise KeyError("annotation_review_table must contain a 'cluster' column.")
    if label_col not in annotation_review_table.columns:
        raise KeyError(f"'{label_col}' not found in annotation_review_table.")

    cluster_series = adata.obs[cluster_key].astype(str)
    label_map = annotation_review_table.set_index("cluster")[label_col].astype(str).to_dict()
    adata.obs[key_added] = pd.Categorical(cluster_series.map(label_map).fillna("Unknown"))

    if confidence_col in annotation_review_table.columns:
        conf_map = annotation_review_table.set_index("cluster")[confidence_col].to_dict()
        adata.obs[f"{key_added}_confidence"] = pd.to_numeric(
            cluster_series.map(conf_map), errors="coerce"
        )
    if status_col in annotation_review_table.columns:
        status_map = annotation_review_table.set_index("cluster")[status_col].to_dict()
        needs_review = cluster_series.map(status_map).fillna(True).astype(bool)
        adata.obs[f"{key_added}_status"] = pd.Categorical(
            np.where(needs_review, "needs_review", "accepted")
        )

    adata.obs["cell_compartment"] = _map_compartments(adata.obs[key_added])
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[f"{key_added}_apply_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "label_col": label_col,
            "confidence_col": confidence_col,
            "status_col": status_col,
            "n_labels": int(adata.obs[key_added].nunique()),
        }
    )
    return adata
