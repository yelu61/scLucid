"""Annotation review tables and marker filtering."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

from ...utils import sanitize_for_hdf5
from .utils import _classify_annotation_marker, _resolve_score_columns


def summarize_annotation_evidence(
    adata: AnnData,
    markers_df: pd.DataFrame,
    enrichment_dict: dict,
    groupby: str = "leiden",
    n_markers: int = 10,
    n_terms: int = 5,
    out_file: Optional[str] = None,
) -> Dict[str, str]:
    """
    Export per-cluster marker+enrichment summary (markdown) for AI/manual annotation.
    """
    summary = {}
    for g in markers_df["group"].unique():
        top_genes = (
            markers_df[markers_df["group"] == g]
            .sort_values("logfoldchanges", ascending=False)
            .head(n_markers)["names"]
            .tolist()
        )
        top_terms = (
            enrichment_dict.get(g, pd.DataFrame())
            .sort_values("Adjusted P-value")
            .head(n_terms)["Term"]
            .tolist()
            if g in enrichment_dict and not enrichment_dict[g].empty
            else []
        )
        s = f"### Cluster {g}\nTop markers: {', '.join(top_genes)}\nTop pathways: {', '.join(top_terms)}"
        summary[g] = s
    if out_file:
        with open(out_file, "w") as f:
            f.write("\n\n".join(summary.values()))
    return summary


def filter_marker_table_for_annotation(
    markers_df: pd.DataFrame,
    *,
    gene_col: str = "names",
    group_col: str = "group",
    score_col: Optional[str] = "logfoldchanges",
    keep_categories: Optional[Iterable[str]] = None,
    keep_top_n_per_group: Optional[int] = None,
    min_score: Optional[float] = None,
    custom_drop_genes: Optional[Iterable[str]] = None,
    drop_noise: bool = True,
) -> pd.DataFrame:
    """
    Filter marker tables for manual annotation by removing common noisy genes.

    The returned DataFrame always includes:
    - `annotation_noise_category`
    - `is_annotation_informative`
    """
    if gene_col not in markers_df.columns:
        raise KeyError(f"'{gene_col}' not found in markers_df.")
    if group_col not in markers_df.columns:
        raise KeyError(f"'{group_col}' not found in markers_df.")
    if score_col is not None and score_col not in markers_df.columns:
        raise KeyError(f"'{score_col}' not found in markers_df.")

    filtered = markers_df.copy()
    filtered["annotation_noise_category"] = filtered[gene_col].map(_classify_annotation_marker)

    if custom_drop_genes is not None:
        custom_drop = {str(g).upper() for g in custom_drop_genes}
        custom_mask = filtered[gene_col].astype(str).str.upper().isin(custom_drop)
        filtered.loc[custom_mask, "annotation_noise_category"] = "custom_drop"

    informative_mask = filtered["annotation_noise_category"].isna()
    if keep_categories is not None:
        allowed = set(keep_categories)
        informative_mask = informative_mask | filtered["annotation_noise_category"].isin(allowed)

    filtered["is_annotation_informative"] = informative_mask.astype(bool)

    if min_score is not None and score_col is not None:
        filtered = filtered[filtered[score_col] >= float(min_score)].copy()

    if keep_top_n_per_group is not None:
        sort_col = score_col if score_col is not None else gene_col
        filtered = (
            filtered.sort_values([group_col, sort_col], ascending=[True, False])
            .groupby(group_col, group_keys=False)
            .head(int(keep_top_n_per_group))
            .copy()
        )

    if drop_noise:
        filtered = filtered[filtered["is_annotation_informative"]].copy()

    return filtered.reset_index(drop=True)


def flag_suspect_clusters(
    adata: AnnData,
    cluster_key: str,
    *,
    markers_df: Optional[pd.DataFrame] = None,
    marker_gene_col: str = "names",
    marker_group_col: str = "group",
    marker_score_col: str = "logfoldchanges",
    top_n_markers: int = 15,
    doublet_flag_cols: Sequence[str] = (
        "predicted_doublet",
        "scrublet_predicted",
        "doubletdetection_predicted",
    ),
    mt_col: Optional[str] = "pct_counts_mt",
    n_genes_col: Optional[str] = "n_genes_by_counts",
    score_cols: Optional[Sequence[str]] = None,
    ribosomal_fraction_threshold: float = 0.5,
    stress_fraction_threshold: float = 0.4,
    doublet_fraction_threshold: float = 0.2,
    mt_mean_threshold: float = 15.0,
    min_informative_markers: int = 3,
    key_added: Optional[str] = None,
) -> pd.DataFrame:
    """
    Flag clusters dominated by low-information, stress, or doublet-like signals.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")

    available_scores = _resolve_score_columns(adata, score_cols)
    cluster_series = adata.obs[cluster_key].astype(str)
    cluster_order = cluster_series.drop_duplicates().tolist()
    rows: List[Dict[str, Any]] = []

    marker_subset = None
    if markers_df is not None:
        if marker_gene_col not in markers_df.columns:
            raise KeyError(f"'{marker_gene_col}' not found in markers_df.")
        if marker_group_col not in markers_df.columns:
            raise KeyError(f"'{marker_group_col}' not found in markers_df.")
        if marker_score_col not in markers_df.columns:
            raise KeyError(f"'{marker_score_col}' not found in markers_df.")

        marker_subset = (
            markers_df.copy()
            .sort_values([marker_group_col, marker_score_col], ascending=[True, False])
            .groupby(marker_group_col, group_keys=False)
            .head(int(top_n_markers))
            .copy()
        )
        marker_subset["annotation_noise_category"] = marker_subset[marker_gene_col].map(
            _classify_annotation_marker
        )

    for cluster in cluster_order:
        mask = cluster_series == cluster
        cluster_obs = adata.obs.loc[mask]
        row: Dict[str, Any] = {
            "cluster": cluster,
            "n_cells": int(mask.sum()),
            "cell_fraction": float(mask.mean()),
        }

        if mt_col and mt_col in cluster_obs.columns:
            row["mean_pct_counts_mt"] = float(cluster_obs[mt_col].mean())
        else:
            row["mean_pct_counts_mt"] = np.nan

        if n_genes_col and n_genes_col in cluster_obs.columns:
            row["mean_n_genes_by_counts"] = float(cluster_obs[n_genes_col].mean())
        else:
            row["mean_n_genes_by_counts"] = np.nan

        present_doublet_cols = [
            col
            for col in doublet_flag_cols
            if col in cluster_obs.columns and pd.api.types.is_bool_dtype(cluster_obs[col])
        ]
        if present_doublet_cols:
            doublet_mask = cluster_obs[present_doublet_cols].any(axis=1)
            row["doublet_fraction"] = float(doublet_mask.mean())
        else:
            row["doublet_fraction"] = np.nan

        for score_col in available_scores:
            row[f"mean_{score_col}"] = float(cluster_obs[score_col].mean())

        suspect_reasons: List[str] = []
        primary_flag = "clean"

        if marker_subset is not None:
            cluster_markers = marker_subset[
                marker_subset[marker_group_col].astype(str) == str(cluster)
            ].copy()
            noise_counts = cluster_markers["annotation_noise_category"].value_counts()
            n_markers = int(cluster_markers.shape[0])
            n_informative = int(cluster_markers["annotation_noise_category"].isna().sum())
            row["n_top_markers"] = n_markers
            row["n_informative_markers"] = n_informative
            row["ribosomal_marker_fraction"] = (
                float(noise_counts.get("ribosomal", 0) / n_markers) if n_markers else np.nan
            )
            row["stress_marker_fraction"] = (
                float(noise_counts.get("stress", 0) / n_markers) if n_markers else np.nan
            )
            row["mitochondrial_marker_fraction"] = (
                float(noise_counts.get("mitochondrial", 0) / n_markers) if n_markers else np.nan
            )
            row["housekeeping_marker_fraction"] = (
                float(noise_counts.get("housekeeping", 0) / n_markers) if n_markers else np.nan
            )
            row["top_marker_preview"] = ", ".join(
                cluster_markers[marker_gene_col].astype(str).head(6).tolist()
            )

            if n_markers and row["ribosomal_marker_fraction"] >= ribosomal_fraction_threshold:
                suspect_reasons.append("ribosomal_dominant")
            if n_markers and row["stress_marker_fraction"] >= stress_fraction_threshold:
                suspect_reasons.append("stress_high")
            if n_markers and n_informative < min_informative_markers:
                suspect_reasons.append("low_information")
        else:
            row["n_top_markers"] = np.nan
            row["n_informative_markers"] = np.nan
            row["ribosomal_marker_fraction"] = np.nan
            row["stress_marker_fraction"] = np.nan
            row["mitochondrial_marker_fraction"] = np.nan
            row["housekeeping_marker_fraction"] = np.nan
            row["top_marker_preview"] = ""

        if (
            pd.notna(row["doublet_fraction"])
            and row["doublet_fraction"] >= doublet_fraction_threshold
        ):
            suspect_reasons.append("doublet_suspect")
        if pd.notna(row["mean_pct_counts_mt"]) and row["mean_pct_counts_mt"] >= mt_mean_threshold:
            suspect_reasons.append("mt_high")

        for severity in (
            "doublet_suspect",
            "stress_high",
            "ribosomal_dominant",
            "low_information",
            "mt_high",
        ):
            if severity in suspect_reasons:
                primary_flag = severity
                break

        row["suspect_reasons"] = ",".join(suspect_reasons)
        row["suspect_flag"] = primary_flag
        rows.append(row)

    summary_df = pd.DataFrame(rows)
    target_key = key_added or f"{cluster_key}_suspect_flags"
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[target_key] = summary_df
    annotation_ns[f"{target_key}_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "top_n_markers": int(top_n_markers),
            "doublet_flag_cols": list(doublet_flag_cols),
            "mt_col": mt_col,
            "n_genes_col": n_genes_col,
            "score_cols": available_scores,
        }
    )

    flag_map = summary_df.set_index("cluster")["suspect_flag"].to_dict()
    reason_map = summary_df.set_index("cluster")["suspect_reasons"].to_dict()
    adata.obs[f"{target_key}_flag"] = cluster_series.map(flag_map).astype("category")
    adata.obs[f"{target_key}_reasons"] = cluster_series.map(reason_map).astype(str)
    return summary_df


def build_annotation_review_table(
    adata: AnnData,
    cluster_key: str,
    *,
    markers_df: Optional[pd.DataFrame] = None,
    marker_gene_col: str = "names",
    marker_group_col: str = "group",
    marker_score_col: Optional[str] = "logfoldchanges",
    enrichment_dict: Optional[Dict[str, pd.DataFrame]] = None,
    annotation_key: Optional[str] = None,
    sample_col: Optional[str] = None,
    group_col: Optional[str] = None,
    time_col: Optional[str] = None,
    score_cols: Optional[Sequence[str]] = None,
    top_n_markers: int = 8,
    top_n_terms: int = 3,
    key_added: Optional[str] = None,
) -> pd.DataFrame:
    """
    Build a per-cluster review table for manual annotation and notebook reporting.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")

    cluster_series = adata.obs[cluster_key].astype(str)
    total_cells = max(adata.n_obs, 1)
    score_cols = _resolve_score_columns(adata, score_cols)

    filtered_markers = None
    if markers_df is not None:
        filtered_markers = filter_marker_table_for_annotation(
            markers_df,
            gene_col=marker_gene_col,
            group_col=marker_group_col,
            score_col=marker_score_col if marker_score_col in markers_df.columns else None,
            keep_top_n_per_group=top_n_markers,
            drop_noise=True,
        )

    rows: List[Dict[str, Any]] = []
    for cluster in cluster_series.drop_duplicates().tolist():
        mask = cluster_series == cluster
        cluster_obs = adata.obs.loc[mask]
        row: Dict[str, Any] = {
            "cluster": cluster,
            "n_cells": int(mask.sum()),
            "pct_cells": float(mask.sum() / total_cells),
        }

        if annotation_key and annotation_key in cluster_obs.columns:
            row["annotation"] = str(cluster_obs[annotation_key].mode(dropna=True).iloc[0])
        else:
            row["annotation"] = None

        if filtered_markers is not None and not filtered_markers.empty:
            top_markers = (
                filtered_markers[filtered_markers[marker_group_col].astype(str) == str(cluster)][
                    marker_gene_col
                ]
                .astype(str)
                .head(top_n_markers)
                .tolist()
            )
            row["top_markers"] = ", ".join(top_markers)
        else:
            row["top_markers"] = ""

        if enrichment_dict:
            cluster_terms = enrichment_dict.get(cluster)
            if cluster_terms is None:
                cluster_terms = enrichment_dict.get(str(cluster))
            if isinstance(cluster_terms, pd.DataFrame) and not cluster_terms.empty:
                sort_col = (
                    "Adjusted P-value"
                    if "Adjusted P-value" in cluster_terms.columns
                    else cluster_terms.columns[0]
                )
                row["top_terms"] = (
                    ", ".join(
                        cluster_terms.sort_values(sort_col)
                        .head(top_n_terms)["Term"]
                        .astype(str)
                        .tolist()
                    )
                    if "Term" in cluster_terms.columns
                    else ""
                )
            else:
                row["top_terms"] = ""
        else:
            row["top_terms"] = ""

        if sample_col and sample_col in cluster_obs.columns:
            sample_counts = cluster_obs[sample_col].astype(str).value_counts(normalize=True).head(3)
            row["top_samples"] = ", ".join(
                f"{sample}:{frac:.2f}" for sample, frac in sample_counts.items()
            )
        else:
            row["top_samples"] = ""

        if group_col and group_col in cluster_obs.columns:
            group_counts = cluster_obs[group_col].astype(str).value_counts(normalize=True).head(3)
            row["group_distribution"] = ", ".join(
                f"{name}:{frac:.2f}" for name, frac in group_counts.items()
            )
        else:
            row["group_distribution"] = ""

        if time_col and time_col in cluster_obs.columns:
            time_counts = cluster_obs[time_col].astype(str).value_counts(normalize=True).head(3)
            row["time_distribution"] = ", ".join(
                f"{name}:{frac:.2f}" for name, frac in time_counts.items()
            )
        else:
            row["time_distribution"] = ""

        row["mean_scores"] = ", ".join(f"{col}:{cluster_obs[col].mean():.3f}" for col in score_cols)
        rows.append(row)

    review_df = pd.DataFrame(rows)
    target_key = key_added or f"{cluster_key}_review_table"
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[target_key] = review_df
    annotation_ns[f"{target_key}_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "marker_gene_col": marker_gene_col,
            "marker_group_col": marker_group_col,
            "marker_score_col": marker_score_col,
            "annotation_key": annotation_key,
            "sample_col": sample_col,
            "group_col": group_col,
            "time_col": time_col,
            "score_cols": score_cols,
            "top_n_markers": int(top_n_markers),
            "top_n_terms": int(top_n_terms),
        }
    )
    return review_df
