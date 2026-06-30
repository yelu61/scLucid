"""Annotation evidence report generation."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Literal, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from ...utils import sanitize_for_hdf5
from ..config import AnnotationConfig
from ..trace import build_annotation_cluster_evidence_table


def _format_distribution(values: pd.Series, n: int = 3) -> str:
    """Return a compact top-category distribution string."""
    if values.empty:
        return ""
    counts = values.astype(str).value_counts(normalize=True).head(n)
    return ", ".join(f"{name}:{frac:.2f}" for name, frac in counts.items())


def _resolve_cluster_and_annotation_keys(
    adata: AnnData,
    cluster_key: Optional[str],
    annotation_key: Optional[str],
) -> tuple[str, str]:
    """Resolve sensible defaults for cluster and annotation columns."""
    if cluster_key is None:
        for candidate in ("leiden_clusters", "leiden", "cluster"):
            if candidate in adata.obs.columns:
                cluster_key = candidate
                break
        if cluster_key is None:
            raise ValueError("No cluster column found; provide cluster_key.")

    if annotation_key is None:
        for candidate in ("cell_type", "cell_type_auto"):
            if candidate in adata.obs.columns:
                annotation_key = candidate
                break
        if annotation_key is None:
            raise ValueError("No annotation column found; provide annotation_key.")

    return str(cluster_key), str(annotation_key)


def _build_per_cluster_report(
    adata: AnnData,
    cluster_key: str,
    annotation_key: str,
    review_table: Optional[pd.DataFrame],
) -> List[Dict[str, Any]]:
    """Build per-cluster evidence rows using existing trace helpers."""
    # Use the existing cluster-evidence builder; pass a minimal config-like object.
    config = AnnotationConfig(key_added=annotation_key)
    cluster_evidence = build_annotation_cluster_evidence_table(adata, config=config)

    # Pull suspect-cluster reasons if they have been computed.
    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    suspect_key = f"{cluster_key}_suspect_flags"
    suspect_df = annotation_ns.get(suspect_key)
    if isinstance(suspect_df, pd.DataFrame) and not suspect_df.empty:
        suspect_map = suspect_df.set_index("cluster")["suspect_reasons"].to_dict()
        flag_map = suspect_df.set_index("cluster")["suspect_flag"].to_dict()
    else:
        suspect_map = {}
        flag_map = {}

    rows: List[Dict[str, Any]] = []
    for entry in cluster_evidence:
        cluster = str(entry.get("cluster", ""))
        row = {
            "cluster": cluster,
            "top_marker_label": entry.get("marker_label"),
            "celltypist_label": entry.get("reference_model_label"),
            "celltypist_confidence": entry.get("reference_model_confidence"),
            "annotation_confidence": entry.get("annotation_confidence"),
            "confidence_level": entry.get("confidence_level"),
            "evidence_status": entry.get("evidence_status"),
            "conflicting_labels": entry.get("contradictory_labels", []),
            "conflicts": entry.get("conflicts", []),
            "warnings": entry.get("warnings", []),
            "requires_manual_review": entry.get("requires_manual_review", True),
            "suspect_flag": flag_map.get(cluster, "clean"),
            "suspect_reasons": suspect_map.get(cluster, ""),
            "top_markers": entry.get("positive_marker_support", []),
            "n_cells": entry.get("n_cells"),
        }
        rows.append(row)
    return rows


def _build_per_celltype_report(
    adata: AnnData,
    annotation_key: str,
    cluster_key: str,
    sample_col: Optional[str],
) -> Dict[str, Dict[str, Any]]:
    """Build per-cell-type summary statistics."""
    confidence_key = f"{annotation_key}_confidence"
    has_confidence = confidence_key in adata.obs.columns

    report: Dict[str, Dict[str, Any]] = {}
    for cell_type in sorted(adata.obs[annotation_key].astype(str).unique()):
        mask = adata.obs[annotation_key].astype(str) == cell_type
        ct_obs = adata.obs.loc[mask]

        cluster_counts = ct_obs[cluster_key].astype(str).value_counts()
        purity = (
            float(cluster_counts.iloc[0] / max(1, cluster_counts.sum()))
            if not cluster_counts.empty
            else np.nan
        )

        median_confidence = np.nan
        if has_confidence:
            conf_values = pd.to_numeric(ct_obs[confidence_key], errors="coerce")
            if conf_values.notna().any():
                median_confidence = float(conf_values.median())

        report[cell_type] = {
            "n_cells": int(mask.sum()),
            "median_annotation_confidence": median_confidence,
            "cluster_purity": purity,
            "clusters_present": sorted(cluster_counts.index.tolist()),
            "sample_distribution": _format_distribution(ct_obs[sample_col])
            if sample_col and sample_col in ct_obs.columns
            else "",
        }
    return report


def _build_global_report(
    adata: AnnData,
    cluster_key: str,
    annotation_key: str,
    confidence_threshold: float,
    config: Optional[AnnotationConfig],
    per_cluster_rows: List[Dict[str, Any]],
    sample_col: Optional[str],
) -> Dict[str, Any]:
    """Build global report section and action items."""
    confidence_key = f"{annotation_key}_confidence"
    low_conf_cells: Optional[int] = None
    if confidence_key in adata.obs.columns:
        conf_values = pd.to_numeric(adata.obs[confidence_key], errors="coerce")
        low_conf_cells = int((conf_values < confidence_threshold).sum())

    n_review_clusters = sum(row.get("requires_manual_review", True) for row in per_cluster_rows)
    n_suspect_clusters = sum(
        row.get("suspect_flag", "clean") != "clean" for row in per_cluster_rows
    )
    n_conflicting_clusters = sum(
        bool(row.get("conflicting_labels") or row.get("conflicts")) for row in per_cluster_rows
    )

    action_items: List[str] = []
    if low_conf_cells:
        action_items.append(
            f"{low_conf_cells} cells have annotation confidence below {confidence_threshold}; "
            "review before treating labels as final."
        )
    if n_conflicting_clusters:
        action_items.append(
            f"{n_conflicting_clusters} cluster(s) show conflicting annotation evidence; "
            "inspect marker/reference support."
        )
    if n_suspect_clusters:
        action_items.append(
            f"{n_suspect_clusters} cluster(s) are flagged as suspect; review QC context."
        )
    if n_review_clusters:
        action_items.append(
            f"{n_review_clusters} cluster(s) require manual review based on evidence thresholds."
        )
    if not action_items:
        action_items.append("No critical annotation evidence issues detected.")

    reference_model = None
    if config is not None:
        reference_model = config.celltypist_model if config.run_celltypist else None

    return {
        "annotation_method": config.final_method if config is not None else "unknown",
        "reference_model": reference_model,
        "confidence_threshold": float(confidence_threshold),
        "n_low_confidence_cells": low_conf_cells,
        "n_cells": int(adata.n_obs),
        "n_clusters": int(adata.obs[cluster_key].nunique()),
        "n_cell_types": int(adata.obs[annotation_key].nunique()),
        "sample_col": sample_col,
        "n_review_clusters": int(n_review_clusters),
        "n_suspect_clusters": int(n_suspect_clusters),
        "n_conflicting_clusters": int(n_conflicting_clusters),
        "action_items": action_items,
    }


def build_annotation_evidence_report(
    adata: AnnData,
    *,
    cluster_key: Optional[str] = None,
    annotation_key: Optional[str] = "cell_type",
    confidence_threshold: float = 0.5,
    sample_col: Optional[str] = None,
    config: Optional[AnnotationConfig] = None,
    key_added: str = "annotation_evidence_report",
) -> Dict[str, Any]:
    """
    Build a structured annotation evidence report after annotation.

    The report combines existing annotation evidence tables into:
    - per-cluster: top marker label, CellTypist/reference label, confidence scores,
      conflict flags, and suspect-cluster reasons
    - per-cell-type: cell counts, median confidence, cluster purity, sample distribution
    - global: annotation method, reference model, confidence threshold,
      number of low-confidence cells, and action items

    Results are stored at
    ``adata.uns["sclucid"]["analysis"]["annotation_evidence_report"]``.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix after annotation.
    cluster_key : str, optional
        Cluster column. Defaults to ``leiden_clusters`` / ``leiden`` / ``cluster``.
    annotation_key : str, optional
        Final annotation column. Defaults to ``cell_type`` / ``cell_type_auto``.
    confidence_threshold : float, default=0.5
        Threshold for flagging low-confidence cells.
    sample_col : str, optional
        Column with sample identifiers for per-cell-type distribution.
    config : AnnotationConfig, optional
        Annotation configuration for method/model metadata.
    key_added : str, default="annotation_evidence_report"
        Storage key under ``adata.uns["sclucid"]["analysis"]``.

    Returns:
    -------
    dict
        The evidence report dictionary.
    """
    cluster_key, annotation_key = _resolve_cluster_and_annotation_keys(
        adata, cluster_key, annotation_key
    )

    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    review_table = annotation_ns.get("annotation_review_table")
    if isinstance(review_table, dict):
        # Possible HDF5 round-trip; attempt to recover
        review_table = pd.DataFrame.from_dict(review_table, orient="columns")

    per_cluster_rows = _build_per_cluster_report(
        adata, cluster_key, annotation_key, review_table
    )
    per_celltype = _build_per_celltype_report(
        adata, annotation_key, cluster_key, sample_col
    )
    global_report = _build_global_report(
        adata,
        cluster_key,
        annotation_key,
        confidence_threshold,
        config,
        per_cluster_rows,
        sample_col,
    )

    report: Dict[str, Any] = {
        "schema_version": "annotation_evidence_report_v1",
        "cluster_key": cluster_key,
        "annotation_key": annotation_key,
        "per_cluster": per_cluster_rows,
        "per_cell_type": per_celltype,
        "global": global_report,
    }

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {})[key_added] = sanitize_for_hdf5(
        report
    )
    return report


def _cell_str(value: Any) -> str:
    """Safely render a table cell that may be a list or missing."""
    if isinstance(value, (list, tuple)):
        return ", ".join(str(v) for v in value)
    if isinstance(value, np.ndarray):
        return ", ".join(str(v) for v in value.tolist())
    try:
        if pd.isna(value):
            return ""
    except (ValueError, TypeError):
        pass
    return str(value)


def _dataframe_to_markdown(df: pd.DataFrame, title: str) -> str:
    """Render a small DataFrame as a markdown table."""
    if df.empty:
        return f"### {title}\n\nNo data.\n"
    lines = [f"### {title}"]
    lines.append("| " + " | ".join(str(c) for c in df.columns) + " |")
    lines.append("| " + " | ".join(["---"] * len(df.columns)) + " |")
    for _, row in df.iterrows():
        lines.append(
            "| " + " | ".join(_cell_str(row[c]) for c in df.columns) + " |"
        )
    return "\n".join(lines) + "\n"


def export_annotation_evidence_report(
    adata: AnnData,
    path: str,
    *,
    format: Literal["markdown", "html"] = "markdown",
    cluster_key: Optional[str] = None,
    annotation_key: Optional[str] = "cell_type",
    confidence_threshold: float = 0.5,
    sample_col: Optional[str] = None,
    config: Optional[AnnotationConfig] = None,
) -> str:
    """
    Export the annotation evidence report to a markdown or HTML file.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    path : str
        Output file path.
    format : {"markdown", "html"}, default="markdown"
        Export format. HTML is a minimal wrapper around the markdown content.
    cluster_key, annotation_key, confidence_threshold, sample_col, config
        Passed to ``build_annotation_evidence_report`` if no report is present.

    Returns:
    -------
    str
        The path to the exported file.
    """
    analysis_ns = adata.uns.get("sclucid", {}).get("analysis", {})
    report = analysis_ns.get("annotation_evidence_report")
    if report is None:
        report = build_annotation_evidence_report(
            adata,
            cluster_key=cluster_key,
            annotation_key=annotation_key,
            confidence_threshold=confidence_threshold,
            sample_col=sample_col,
            config=config,
        )

    per_cluster_df = pd.DataFrame(report.get("per_cluster", []))
    per_celltype_df = pd.DataFrame.from_dict(
        report.get("per_cell_type", {}), orient="index"
    ).reset_index().rename(columns={"index": "cell_type"})
    global_section = report.get("global", {})

    lines: List[str] = []
    lines.append("# Annotation Evidence Report")
    lines.append(f"**Schema version:** {report.get('schema_version', 'unknown')}")
    lines.append(f"**Cluster key:** {report.get('cluster_key')}")
    lines.append(f"**Annotation key:** {report.get('annotation_key')}")
    lines.append("")

    lines.append("## Global Summary")
    for key, value in global_section.items():
        if key == "action_items":
            continue
        lines.append(f"- **{key}:** {value}")
    lines.append("")
    lines.append("### Action Items")
    for item in global_section.get("action_items", []):
        lines.append(f"- {item}")
    lines.append("")

    lines.append(_dataframe_to_markdown(per_cluster_df, "Per-Cluster Evidence"))
    lines.append("")
    lines.append(_dataframe_to_markdown(per_celltype_df, "Per-Cell-Type Summary"))

    content = "\n".join(lines)

    if format == "html":
        content = (
            "<!DOCTYPE html>\n<html>\n<head>"
            "<meta charset=\"utf-8\"><title>Annotation Evidence Report</title>"
            "</head>\n<body>\n"
            + content.replace("\n", "<br>\n")
            + "\n</body>\n</html>"
        )
        suffix = ".html"
    else:
        suffix = ".md"

    out_path = Path(path)
    if not out_path.suffix:
        out_path = out_path.with_suffix(suffix)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(content, encoding="utf-8")
    return str(out_path)


__all__ = [
    "build_annotation_evidence_report",
    "export_annotation_evidence_report",
]
