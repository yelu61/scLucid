"""Plotting helpers for annotation evidence, confidence review, and reports."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

log = logging.getLogger(__name__)

_PRIORITY_RANK = {"Review Now": 0, "Optional Check": 1, "Safe To Proceed": 2}
_PRIORITY_COLORS = {
    "Review Now": "#c93c37",
    "Optional Check": "#d98c10",
    "Safe To Proceed": "#6c757d",
}


def _load_annotation_evidence(
    adata: AnnData,
    annotation_key: str,
    evidence_key: Optional[str] = None,
) -> pd.DataFrame:
    """Load and validate annotation evidence table from ``adata.uns``."""
    annotation_ns = adata.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
    evidence_key = evidence_key or f"{annotation_key}_evidence"
    if evidence_key not in annotation_ns:
        raise KeyError(
            f"Annotation evidence '{evidence_key}' not found. Run `run_annotation()` first."
        )

    evidence = annotation_ns[evidence_key]
    if not isinstance(evidence, pd.DataFrame):
        evidence = pd.DataFrame(evidence)
    if evidence.empty:
        raise ValueError(f"Annotation evidence '{evidence_key}' is empty.")
    return evidence.copy()


def _resolve_annotation_palette(
    df: pd.DataFrame, palette: Optional[Dict[str, str]] = None
) -> Dict[str, str]:
    """Build a label palette for annotation plots."""
    if palette is not None:
        return palette
    labels = df["assigned_label"].dropna().astype(str).unique().tolist()
    colors = sns.color_palette("tab20", n_colors=max(1, len(labels)))
    return {label: colors[i] for i, label in enumerate(labels)}


def _compute_top_markers_summary(
    adata: AnnData,
    cluster_key: str,
    *,
    n_markers: int = 3,
    max_clusters: int = 8,
) -> Dict[str, List[str]]:
    """Compute a compact top-marker summary for clusters."""
    if cluster_key not in adata.obs.columns:
        return {}

    use_raw = adata.raw is not None
    tmp_key = f"__annotation_report_markers_{cluster_key}"
    markers: Dict[str, List[str]] = {}
    try:
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw,
            key_added=tmp_key,
            n_genes=max(n_markers, 10),
        )
        cluster_values = pd.Series(adata.obs[cluster_key]).astype(str)
        cluster_order = cluster_values.value_counts().index.astype(str).tolist()[:max_clusters]
        for cluster in cluster_order:
            marker_df = sc.get.rank_genes_groups_df(adata, key=tmp_key, group=cluster)
            if marker_df.empty or "names" not in marker_df.columns:
                markers[str(cluster)] = []
                continue
            top_genes = marker_df["names"].astype(str).dropna().head(n_markers).tolist()
            markers[str(cluster)] = top_genes
    except Exception as exc:
        log.warning(f"Could not compute top marker summary for report: {exc}")
    finally:
        if tmp_key in adata.uns:
            del adata.uns[tmp_key]
    return markers


def _collect_annotation_warnings(
    df: pd.DataFrame,
    *,
    low_confidence_threshold: float = 0.55,
    low_purity_threshold: float = 0.7,
    low_agreement_threshold: float = 0.5,
    max_warnings: int = 8,
) -> List[str]:
    """Collect cluster warnings that deserve manual review."""
    warnings = _build_warning_records(
        df,
        low_confidence_threshold=low_confidence_threshold,
        low_purity_threshold=low_purity_threshold,
        low_agreement_threshold=low_agreement_threshold,
        max_warnings=max_warnings,
    )
    return _format_warning_records(warnings)


def _build_warning_records(
    df: pd.DataFrame,
    *,
    low_confidence_threshold: float = 0.55,
    low_purity_threshold: float = 0.7,
    low_agreement_threshold: float = 0.5,
    max_warnings: int = 8,
) -> List[Dict[str, Any]]:
    """Build structured warning records used by report text and sidecars."""
    warning_records = []
    ranked = df.sort_values(by="annotation_confidence", ascending=True)
    for _, row in ranked.iterrows():
        cluster = str(row["cluster"])
        label = str(row["assigned_label"])
        reasons = []
        confidence = float(row.get("annotation_confidence", np.nan))
        purity = float(row.get("label_purity", np.nan))
        agreement = row.get("celltypist_agreement", np.nan)
        hybrid_decision = row.get("hybrid_decision")
        action_hints = []
        priority = "Optional Check"

        if pd.notna(confidence) and confidence < low_confidence_threshold:
            reasons.append(f"low confidence {confidence:.2f}")
            action_hints.append("review top markers and cluster-level label separation")
            priority = "Review Now"
        if pd.notna(purity) and purity < low_purity_threshold:
            reasons.append(f"label purity {purity:.2f}")
            action_hints.append("consider re-clustering or subsetting this cluster")
            priority = "Review Now"
        if pd.notna(agreement) and float(agreement) < low_agreement_threshold:
            reasons.append(f"marker/CellTypist disagreement {float(agreement):.2f}")
            action_hints.append("check marker overlap and reduce CellTypist weight if needed")
            priority = "Review Now"
        if pd.notna(hybrid_decision) and str(hybrid_decision) in {
            "marker_fallback",
            "celltypist_majority_fallback",
            "insufficient_evidence",
        }:
            reasons.append(f"decision {hybrid_decision}")
            if str(hybrid_decision) == "marker_fallback":
                action_hints.append("validate CellTypist model choice or confidence threshold")
                priority = "Review Now"
            elif str(hybrid_decision) == "celltypist_majority_fallback":
                action_hints.append("verify marker coverage before accepting majority-vote labels")
                priority = "Optional Check"
            elif str(hybrid_decision) == "insufficient_evidence":
                action_hints.append(
                    "flag for manual annotation and inspect raw expression patterns"
                )
                priority = "Review Now"

        if reasons:
            unique_hints = list(dict.fromkeys(action_hints))
            warning_records.append(
                {
                    "priority": priority,
                    "priority_rank": _PRIORITY_RANK.get(priority, 9),
                    "priority_color": _PRIORITY_COLORS.get(priority, "#6c757d"),
                    "cluster": cluster,
                    "label": label,
                    "annotation_confidence": (
                        round(float(confidence), 4) if pd.notna(confidence) else None
                    ),
                    "reasons": reasons,
                    "actions": unique_hints,
                }
            )

    warning_records = sorted(
        warning_records,
        key=lambda item: (
            item["priority_rank"],
            item["annotation_confidence"] if item["annotation_confidence"] is not None else 1.0,
        ),
    )[:max_warnings]

    if not warning_records:
        warning_records.append(
            {
                "priority": "Safe To Proceed",
                "priority_rank": _PRIORITY_RANK["Safe To Proceed"],
                "priority_color": _PRIORITY_COLORS["Safe To Proceed"],
                "cluster": "all",
                "label": "consistent",
                "annotation_confidence": None,
                "reasons": ["No major warnings. Annotation evidence is internally consistent."],
                "actions": ["proceed to manual biological review of top markers only"],
            }
        )
    return warning_records


def _format_warning_records(warning_records: List[Dict[str, Any]]) -> List[str]:
    """Format structured warning records into report-friendly text lines."""
    lines: List[str] = []
    for record in warning_records:
        action_text = ""
        if record.get("actions"):
            action_text = " | Action: " + " / ".join(record["actions"])
        lines.append(
            f"[{record['priority']}] Cluster {record['cluster']} ({record['label']}): "
            + "; ".join(record.get("reasons", []))
            + action_text
        )
    return lines


def _build_risk_summary(
    df: pd.DataFrame,
    warning_records: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Aggregate high-level risk counts for report overview and sidecars."""
    priority_counts = dict.fromkeys(_PRIORITY_RANK, 0)
    flagged_clusters = []
    for record in warning_records:
        priority = record.get("priority", "Safe To Proceed")
        priority_counts[priority] = priority_counts.get(priority, 0) + 1
        if record.get("cluster") != "all":
            flagged_clusters.append(str(record["cluster"]))

    confidence = pd.to_numeric(df.get("annotation_confidence"), errors="coerce")
    purity = pd.to_numeric(df.get("label_purity"), errors="coerce")
    agreement = pd.to_numeric(df.get("celltypist_agreement"), errors="coerce")
    unique_labels = df.get("assigned_label", pd.Series(dtype=str)).astype(str).nunique()
    return {
        "n_clusters": int(df["cluster"].astype(str).nunique()),
        "n_labels": int(unique_labels),
        "mean_annotation_confidence": (
            round(float(confidence.mean()), 4) if confidence.notna().any() else None
        ),
        "median_annotation_confidence": (
            round(float(confidence.median()), 4) if confidence.notna().any() else None
        ),
        "mean_label_purity": round(float(purity.mean()), 4) if purity.notna().any() else None,
        "mean_celltypist_agreement": (
            round(float(agreement.mean()), 4) if agreement.notna().any() else None
        ),
        "priority_counts": priority_counts,
        "flagged_clusters": flagged_clusters,
    }


def _render_risk_overview(ax: plt.Axes, summary: Dict[str, Any]) -> None:
    """Render a compact risk overview with color-coded counts."""
    ax.axis("off")
    ax.set_title("Risk Overview", pad=10)
    ax.add_patch(
        plt.Rectangle(
            (0.0, 0.0),
            1.0,
            1.0,
            transform=ax.transAxes,
            facecolor="#f8f9fb",
            edgecolor="#cfd6df",
            linewidth=1.0,
            zorder=0,
        )
    )
    count_y = 0.82
    for priority in ["Review Now", "Optional Check", "Safe To Proceed"]:
        ax.text(
            0.06,
            count_y,
            priority,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=9,
            color=_PRIORITY_COLORS[priority],
            fontweight="bold",
        )
        ax.text(
            0.92,
            count_y,
            str(summary["priority_counts"].get(priority, 0)),
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=12,
            color=_PRIORITY_COLORS[priority],
            fontweight="bold",
        )
        ax.add_line(
            plt.Line2D(
                [0.06, 0.92],
                [count_y - 0.08, count_y - 0.08],
                transform=ax.transAxes,
                color="#e5e9f0",
                linewidth=0.8,
            )
        )
        count_y -= 0.18

    stats_lines = [
        f"Clusters: {summary['n_clusters']}",
        f"Labels: {summary['n_labels']}",
        "Mean confidence: "
        + (
            f"{summary['mean_annotation_confidence']:.2f}"
            if summary["mean_annotation_confidence"] is not None
            else "n/a"
        ),
        "Mean purity: "
        + (
            f"{summary['mean_label_purity']:.2f}"
            if summary["mean_label_purity"] is not None
            else "n/a"
        ),
    ]
    ax.text(
        0.06,
        0.28,
        "\n".join(stats_lines),
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        family="monospace",
    )


def _resolve_export_base_path(save: str | Path) -> Path:
    """Resolve a stable base path used to write multiple formats and sidecars."""
    save_path = Path(save)
    if save_path.suffix.lower() in {".png", ".pdf", ".svg"}:
        return save_path.with_suffix("")
    return save_path


def _write_report_sidecars(
    base_path: Path,
    *,
    annotation_key: str,
    cluster_key: str,
    summary: Dict[str, Any],
    warning_records: List[Dict[str, Any]],
    top_markers: Dict[str, List[str]],
) -> None:
    """Write machine-readable and reviewer-friendly sidecar summaries."""
    payload = {
        "annotation_key": annotation_key,
        "cluster_key": cluster_key,
        "summary": summary,
        "warnings": warning_records,
        "top_markers": top_markers,
    }
    json_path = base_path.with_suffix(".json")
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n")

    markdown_lines = [
        f"# {annotation_key.replace('_', ' ').title()} Annotation Report",
        "",
        f"- Cluster key: `{cluster_key}`",
        f"- Total clusters: `{summary['n_clusters']}`",
        f"- Total labels: `{summary['n_labels']}`",
        "- Mean annotation confidence: `"
        + (
            f"{summary['mean_annotation_confidence']:.2f}`"
            if summary["mean_annotation_confidence"] is not None
            else "n/a`"
        ),
        "- Mean label purity: `"
        + (
            f"{summary['mean_label_purity']:.2f}`"
            if summary["mean_label_purity"] is not None
            else "n/a`"
        ),
        "",
        "## Risk Summary",
        "",
    ]
    for priority in ["Review Now", "Optional Check", "Safe To Proceed"]:
        markdown_lines.append(f"- {priority}: `{summary['priority_counts'].get(priority, 0)}`")
    markdown_lines.extend(["", "## Warnings", ""])
    for line in _format_warning_records(warning_records):
        markdown_lines.append(f"- {line}")
    markdown_lines.extend(["", "## Top Markers", ""])
    for cluster, markers in top_markers.items():
        marker_text = ", ".join(markers) if markers else "markers unavailable"
        markdown_lines.append(f"- Cluster {cluster}: {marker_text}")
    md_path = base_path.with_suffix(".md")
    md_path.write_text("\n".join(markdown_lines).rstrip() + "\n")


def _render_text_box(
    ax: plt.Axes,
    title: str,
    lines: List[str],
    *,
    facecolor: str = "#f7f7f7",
    edgecolor: str = "#d0d0d0",
) -> None:
    """Render a simple titled text box inside an axis."""
    ax.axis("off")
    ax.set_title(title, pad=10)
    ax.add_patch(
        plt.Rectangle(
            (0.0, 0.0),
            1.0,
            1.0,
            transform=ax.transAxes,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=1.0,
            zorder=0,
        )
    )
    content = "\n".join(lines)
    ax.text(
        0.03,
        0.97,
        content,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        family="monospace",
    )


def _draw_annotation_evidence_panel(
    fig: plt.Figure,
    df: pd.DataFrame,
    palette: Dict[str, str],
) -> None:
    """Draw the multi-panel annotation evidence layout onto a figure-like object."""
    gs = fig.add_gridspec(
        nrows=2,
        ncols=2,
        width_ratios=[1.3, 1.0],
        height_ratios=[1.0, 1.0],
        wspace=0.28,
        hspace=0.28,
    )

    # Panel 1: confidence bars
    ax_bar = fig.add_subplot(gs[:, 0])
    bar_colors = [palette.get(label, "#999999") for label in df["assigned_label"]]
    ax_bar.barh(df["cluster"], df["annotation_confidence"], color=bar_colors, alpha=0.9)
    for idx, row in df.iterrows():
        ax_bar.text(
            min(float(row["annotation_confidence"]) + 0.015, 0.98),
            idx,
            row["assigned_label"],
            va="center",
            fontsize=9,
        )
    ax_bar.set_xlim(0, 1.0)
    ax_bar.set_xlabel("Annotation Confidence")
    ax_bar.set_ylabel("Cluster")
    ax_bar.set_title("Cluster Confidence")
    ax_bar.grid(axis="x", linestyle="--", alpha=0.3)

    # Panel 2: marker vs CellTypist evidence
    ax_scatter = fig.add_subplot(gs[0, 1])
    scatter_df = df.copy()
    scatter_df["marker_confidence"] = pd.to_numeric(
        scatter_df.get("marker_confidence", np.nan), errors="coerce"
    )
    scatter_df["celltypist_mean_confidence"] = pd.to_numeric(
        scatter_df.get("celltypist_mean_confidence", np.nan), errors="coerce"
    )
    scatter_df["n_cells"] = pd.to_numeric(scatter_df.get("n_cells", 0), errors="coerce").fillna(0)
    scatter_df["celltypist_mean_confidence"] = scatter_df["celltypist_mean_confidence"].fillna(0.0)
    scatter_df["marker_confidence"] = scatter_df["marker_confidence"].fillna(0.0)
    size_scale = np.clip(scatter_df["n_cells"].to_numpy(dtype=float), 1, None)
    size_scale = 40 + 160 * (size_scale / size_scale.max()) if len(size_scale) else 60

    for idx, row in scatter_df.iterrows():
        ax_scatter.scatter(
            row["marker_confidence"],
            row["celltypist_mean_confidence"],
            s=size_scale[idx],
            color=palette.get(row["assigned_label"], "#999999"),
            alpha=0.85,
            edgecolor="white",
            linewidth=0.8,
        )
        ax_scatter.text(
            row["marker_confidence"] + 0.01,
            row["celltypist_mean_confidence"] + 0.01,
            row["cluster"],
            fontsize=8,
        )
    ax_scatter.plot([0, 1], [0, 1], linestyle="--", color="#bbbbbb", linewidth=1)
    ax_scatter.set_xlim(0, 1.02)
    ax_scatter.set_ylim(0, 1.02)
    ax_scatter.set_xlabel("Marker Confidence")
    ax_scatter.set_ylabel("CellTypist Confidence")
    ax_scatter.set_title("Evidence Concordance")
    ax_scatter.grid(linestyle="--", alpha=0.25)

    # Panel 3: compact evidence heatmap
    ax_heat = fig.add_subplot(gs[1, 1])
    heat_cols = [
        col
        for col in [
            "label_purity",
            "marker_confidence",
            "celltypist_mean_confidence",
            "celltypist_agreement",
            "annotation_confidence",
        ]
        if col in df.columns
    ]
    heat_df = df.set_index("cluster")[heat_cols].apply(pd.to_numeric, errors="coerce")
    heat_df = heat_df.fillna(0.0)
    sns.heatmap(
        heat_df,
        ax=ax_heat,
        cmap="RdYlBu_r",
        vmin=0,
        vmax=1,
        cbar=True,
        linewidths=0.5,
        linecolor="white",
    )
    ax_heat.set_title("Evidence Matrix")
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("Cluster")

    # Decision strip appended to the right side of heatmap panel.
    decision_series = (
        df.set_index("cluster").get("hybrid_decision") if "hybrid_decision" in df.columns else None
    )
    if decision_series is not None and decision_series.notna().any():
        unique_decisions = [str(x) for x in decision_series.dropna().unique().tolist()]
        decision_palette = {
            decision: color
            for decision, color in zip(
                unique_decisions,
                sns.color_palette("Set2", n_colors=max(1, len(unique_decisions))),
            )
        }
        for i, cluster in enumerate(heat_df.index.tolist()):
            decision = decision_series.get(cluster)
            if pd.isna(decision):
                continue
            ax_heat.add_patch(
                plt.Rectangle(
                    (len(heat_cols) + 0.05, i),
                    0.35,
                    1.0,
                    facecolor=decision_palette[str(decision)],
                    edgecolor="white",
                    linewidth=0.5,
                    clip_on=False,
                )
            )
            ax_heat.text(
                len(heat_cols) + 0.48,
                i + 0.5,
                str(decision),
                va="center",
                fontsize=8,
                clip_on=False,
            )
        ax_heat.set_xlim(0, len(heat_cols) + 2.2)
        ax_heat.text(
            len(heat_cols) + 0.05,
            -0.55,
            "Hybrid Decision",
            fontsize=9,
            fontweight="bold",
            ha="left",
            clip_on=False,
        )

    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=label,
            markerfacecolor=color,
            markersize=8,
        )
        for label, color in palette.items()
    ]
    if handles:
        ax_bar.legend(
            handles=handles,
            title="Assigned Label",
            loc="lower right",
            frameon=False,
            fontsize=8,
            title_fontsize=9,
        )


def plot_annotation_evidence_panel(
    adata: AnnData,
    annotation_key: str = "cell_type_auto",
    evidence_key: Optional[str] = None,
    sort_by: str = "annotation_confidence",
    figsize: Tuple[float, float] = (15, 9),
    palette: Optional[Dict[str, str]] = None,
    save: Optional[str] = None,
    dpi: int = 300,
    show: bool = True,
) -> plt.Figure:
    """
    Plot a multi-panel summary of annotation evidence.

    Parameters
    ----------
    adata
        AnnData object with annotation evidence stored in
        ``adata.uns['sclucid']['analysis']['annotation']``.
    annotation_key
        Annotation column prefix used by ``run_annotation``.
    evidence_key
        Optional override for the evidence table key in ``.uns``.
    sort_by
        Column used to order clusters in the panel.
    figsize
        Figure size.
    palette
        Optional mapping from assigned labels to colors.
    save
        Optional output path.
    dpi
        Resolution for saved figure.
    show
        Whether to display the figure.
    """
    df = _load_annotation_evidence(adata, annotation_key, evidence_key=evidence_key)
    if sort_by not in df.columns:
        raise KeyError(f"Column '{sort_by}' not found in annotation evidence.")
    df["cluster"] = df["cluster"].astype(str)
    df["assigned_label"] = df["assigned_label"].astype(str)
    df = df.sort_values(by=sort_by, ascending=True).reset_index(drop=True)
    palette = _resolve_annotation_palette(df, palette=palette)

    fig = plt.figure(figsize=figsize)
    _draw_annotation_evidence_panel(fig, df, palette)

    title = annotation_key.replace("_", " ").title()
    fig.suptitle(f"{title} Evidence Panel", fontsize=15, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
        log.info(f"Saved annotation evidence panel to {save}")

    if show:
        plt.show()

    return fig


def export_annotation_report(
    adata: AnnData,
    annotation_key: str = "cell_type_auto",
    cluster_key: str = "leiden_clusters",
    evidence_key: Optional[str] = None,
    figsize: Tuple[float, float] = (15, 14),
    palette: Optional[Dict[str, str]] = None,
    top_marker_count: int = 3,
    max_marker_clusters: int = 8,
    save: Optional[str] = None,
    export_formats: Tuple[str, ...] = ("png", "pdf"),
    write_sidecars: bool = True,
    dpi: int = 300,
    show: bool = True,
) -> plt.Figure:
    """
    Export a compact annotation review report.

    The report contains:
    - UMAP colored by final annotation (or a placeholder if unavailable)
    - Cluster summary table
    - Full annotation evidence panel
    """
    df = _load_annotation_evidence(adata, annotation_key, evidence_key=evidence_key)
    df["cluster"] = df["cluster"].astype(str)
    df["assigned_label"] = df["assigned_label"].astype(str)
    df = df.sort_values(by="annotation_confidence", ascending=False).reset_index(drop=True)
    palette = _resolve_annotation_palette(df, palette=palette)
    warning_records = _build_warning_records(df)
    warning_lines = _format_warning_records(warning_records)
    risk_summary = _build_risk_summary(df, warning_records)
    top_markers = _compute_top_markers_summary(
        adata,
        cluster_key,
        n_markers=top_marker_count,
        max_clusters=max_marker_clusters,
    )
    marker_lines = []
    for _, row in df.head(max_marker_clusters).iterrows():
        cluster = str(row["cluster"])
        label = str(row["assigned_label"])
        genes = top_markers.get(cluster, [])
        marker_lines.append(
            f"C{cluster} {label}: " + (", ".join(genes) if genes else "markers unavailable")
        )
    if not marker_lines:
        marker_lines.append("No marker summary available.")

    fig = plt.figure(figsize=figsize)
    subfigs = fig.subfigures(2, 1, height_ratios=[0.42, 0.58], hspace=0.06)
    top_subfig = subfigs[0]
    bottom_subfig = subfigs[1]

    top_gs = top_subfig.add_gridspec(
        nrows=2,
        ncols=4,
        width_ratios=[1.2, 1.0, 0.95, 1.05],
        height_ratios=[1.0, 1.0],
        wspace=0.25,
        hspace=0.22,
    )
    ax_umap = top_subfig.add_subplot(top_gs[:, 0])
    ax_table = top_subfig.add_subplot(top_gs[:, 1])
    ax_overview = top_subfig.add_subplot(top_gs[0, 2])
    ax_warnings = top_subfig.add_subplot(top_gs[1, 2])
    ax_markers = top_subfig.add_subplot(top_gs[:, 3])

    if "X_umap" in adata.obsm and annotation_key in adata.obs.columns:
        sc.pl.umap(
            adata,
            color=annotation_key,
            ax=ax_umap,
            show=False,
            palette=palette,
            frameon=False,
            legend_loc="on data" if adata.obs[annotation_key].nunique() <= 12 else "right margin",
        )
        ax_umap.set_title("Annotation UMAP")
    else:
        ax_umap.axis("off")
        message = "UMAP unavailable"
        if "X_umap" not in adata.obsm:
            message = "UMAP embedding unavailable"
        elif annotation_key not in adata.obs.columns:
            message = f"'{annotation_key}' not found"
        ax_umap.text(
            0.5,
            0.5,
            message,
            ha="center",
            va="center",
            fontsize=12,
        )

    ax_table.axis("off")
    summary_cols = [
        col
        for col in [
            "cluster",
            "assigned_label",
            "n_cells",
            "annotation_confidence",
            "hybrid_decision",
        ]
        if col in df.columns
    ]
    summary_df = df[summary_cols].copy()
    if "annotation_confidence" in summary_df.columns:
        summary_df["annotation_confidence"] = summary_df["annotation_confidence"].map(
            lambda x: f"{float(x):.2f}"
        )
    if "n_cells" in summary_df.columns:
        summary_df["n_cells"] = summary_df["n_cells"].astype(int).astype(str)
    if "hybrid_decision" in summary_df.columns:
        summary_df["hybrid_decision"] = summary_df["hybrid_decision"].fillna("-")
    summary_df = summary_df.head(12)
    table = ax_table.table(
        cellText=summary_df.values,
        colLabels=[col.replace("_", " ").title() for col in summary_df.columns],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1.0, 1.35)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold")
            cell.set_facecolor("#ececec")
        elif col == 1 and row > 0:
            label = summary_df.iloc[row - 1, col]
            color = palette.get(label)
            if color is not None:
                cell.set_facecolor(color)
                luminance = 0.299 * color[0] + 0.587 * color[1] + 0.114 * color[2]
                cell.set_text_props(color="black" if luminance > 0.6 else "white")
    ax_table.set_title("Cluster Summary", pad=12)
    _render_risk_overview(ax_overview, risk_summary)
    _render_text_box(
        ax_warnings,
        "Warnings",
        warning_lines,
        facecolor="#fff5f0",
        edgecolor="#f0c2ad",
    )
    _render_text_box(
        ax_markers,
        "Top Markers",
        marker_lines,
        facecolor="#f5f7fb",
        edgecolor="#cad5e6",
    )

    _draw_annotation_evidence_panel(
        bottom_subfig, df.sort_values(by="annotation_confidence", ascending=True), palette
    )
    bottom_subfig.suptitle("Evidence Review", fontsize=13, y=0.99)

    fig.suptitle(
        f"{annotation_key.replace('_', ' ').title()} Annotation Report",
        fontsize=16,
        y=0.995,
    )

    if save:
        base_path = _resolve_export_base_path(save)
        base_path.parent.mkdir(parents=True, exist_ok=True)
        normalized_formats = tuple(
            dict.fromkeys(fmt.lower().lstrip(".") for fmt in export_formats if fmt)
        )
        if not normalized_formats:
            normalized_formats = ("png",)
        for fmt in normalized_formats:
            output_path = base_path.with_suffix(f".{fmt}")
            fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
            log.info(f"Saved annotation report to {output_path}")
        if write_sidecars:
            _write_report_sidecars(
                base_path,
                annotation_key=annotation_key,
                cluster_key=cluster_key,
                summary=risk_summary,
                warning_records=warning_records,
                top_markers=top_markers,
            )

    if show:
        plt.show()

    return fig
