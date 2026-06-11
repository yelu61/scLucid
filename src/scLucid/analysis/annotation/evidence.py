"""Evidence-table helpers for annotation review and final label application."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ...utils import Manager, get_marker_manager, sanitize_for_hdf5
from .utils import _classify_annotation_marker, _map_compartments

__all__ = [
    "ANNOTATION_REVIEW_SCHEMA",
    "ANALYSIS_REVIEW_SUMMARY_SCHEMA",
    "standardize_cluster_marker_table",
    "build_hierarchical_annotation_plan",
    "run_subset_annotation_refinement",
    "build_subset_annotation_reconciliation",
    "apply_subset_annotation_reconciliation",
    "run_marker_annotation_evidence",
    "run_program_annotation_evidence",
    "run_annotation_evidence",
    "build_llm_annotation_bundle",
    "merge_annotation_evidence",
    "build_annotation_consensus",
    "apply_final_annotation",
    "evaluate_annotation_benchmark",
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
    "lineage_label",
    "lineage_confidence",
    "subtype_label",
    "subtype_confidence",
    "state_label",
    "state_confidence",
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


def build_hierarchical_annotation_plan(
    adata: AnnData,
    *,
    cluster_key: str,
    lineage_key: str,
    min_cells_per_lineage: int = 50,
    min_clusters_per_lineage: int = 2,
    min_lineage_purity: float = 0.7,
    target_lineages: Optional[List[str]] = None,
    key_added: str = "hierarchical_annotation_plan",
) -> pd.DataFrame:
    """
    Build a lineage-gated refinement plan for hierarchical annotation.

    The plan keeps the annotation workflow conservative: major lineages are
    established first, and subtype/state interpretation is recommended only for
    lineages with enough cells and cluster support to justify refinement.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")
    if lineage_key not in adata.obs.columns:
        raise KeyError(f"'{lineage_key}' not found in adata.obs.")
    if min_cells_per_lineage < 1:
        raise ValueError("min_cells_per_lineage must be >= 1.")
    if min_clusters_per_lineage < 1:
        raise ValueError("min_clusters_per_lineage must be >= 1.")

    cluster_series = adata.obs[cluster_key].astype(str)
    lineage_series = adata.obs[lineage_key].astype(str)
    target_set = {str(item) for item in target_lineages} if target_lineages else None

    rows: List[Dict[str, Any]] = []
    for lineage in sorted(lineage_series.dropna().unique(), key=str):
        if lineage in {"Unknown", "nan", "None", "Not_applicable", ""}:
            continue
        if target_set is not None and lineage not in target_set:
            continue

        lineage_mask = lineage_series == lineage
        n_cells = int(lineage_mask.sum())
        cluster_counts = cluster_series.loc[lineage_mask].value_counts()
        n_clusters = int(cluster_counts.shape[0])
        dominant_fraction = (
            float(cluster_counts.iloc[0] / max(1, n_cells)) if n_cells else 0.0
        )
        cluster_distribution = ", ".join(
            f"{cluster}:{count}" for cluster, count in cluster_counts.head(8).items()
        )

        if n_cells < min_cells_per_lineage:
            action = "keep_lineage"
            reason = "too_few_cells_for_reclustering"
        elif n_clusters < min_clusters_per_lineage:
            action = "keep_lineage"
            reason = "single_cluster_or_low_cluster_count"
        elif dominant_fraction < min_lineage_purity:
            action = "review_lineage_gate"
            reason = "lineage_spread_across_multiple_clusters"
        else:
            action = "subset_recluster_for_subtype"
            reason = "sufficient_cells_and_cluster_support"

        rows.append(
            {
                "lineage_label": lineage,
                "n_cells": n_cells,
                "n_clusters": n_clusters,
                "dominant_cluster_fraction": dominant_fraction,
                "cluster_distribution": cluster_distribution,
                "recommended_action": action,
                "reason": reason,
            }
        )

    plan = pd.DataFrame(rows)
    if not plan.empty:
        plan = plan.sort_values(
            ["recommended_action", "n_cells", "lineage_label"],
            ascending=[True, False, True],
        ).reset_index(drop=True)

    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[key_added] = sanitize_for_hdf5(plan)
    annotation_ns[f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "lineage_key": lineage_key,
            "min_cells_per_lineage": int(min_cells_per_lineage),
            "min_clusters_per_lineage": int(min_clusters_per_lineage),
            "min_lineage_purity": float(min_lineage_purity),
            "target_lineages": target_lineages or [],
        }
    )
    return plan


def _safe_obs_key(value: str) -> str:
    """Return a compact obs-key-safe token."""
    return (
        str(value)
        .strip()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("+", "pos")
        .replace("-", "_")
        .replace("|", "_")
    )


def _extract_subset_for_reprocessing(
    adata: AnnData,
    mask: pd.Series,
    *,
    counts_layer: str = "counts",
    prefer_raw: bool = False,
) -> tuple[AnnData, str]:
    """Extract a subset and reset ``X`` to raw-count-like input when available."""
    subset = adata[mask.to_numpy()].copy()
    input_mode = "current_X_fallback"

    if counts_layer and counts_layer in subset.layers:
        subset.X = subset.layers[counts_layer].copy()
        input_mode = f"layer:{counts_layer}"
    elif prefer_raw and adata.raw is not None:
        raw_subset = adata.raw.to_adata()[mask.to_numpy()].copy()
        raw_subset.obs = subset.obs.copy()
        input_mode = "adata.raw"
        subset = raw_subset

    subset.layers["subset_annotation_input"] = subset.X.copy()
    return subset, input_mode


def _is_missing_label(value: Any) -> bool:
    """Return True for labels that should not be treated as real annotations."""
    if pd.isna(value):
        return True
    return str(value) in {"", "nan", "None", "Unknown", "Not_applicable", "unassigned"}


def _as_bool_review_flag(value: Any) -> bool:
    """Interpret common manual-review exclusion values as booleans."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {
        "1",
        "true",
        "yes",
        "y",
        "exclude",
        "drop",
        "remove",
        "discard",
        "low_quality",
    }


def _infer_subset_cluster_key(subset: AnnData, key_added: str) -> Optional[str]:
    """Infer the cluster key created by subset refinement."""
    annotation_ns = (
        subset.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
        if isinstance(subset.uns, dict)
        else {}
    )
    metadata = annotation_ns.get(key_added, {})
    if isinstance(metadata, dict) and metadata.get("cluster_key") in subset.obs.columns:
        return str(metadata["cluster_key"])
    candidates = [col for col in subset.obs.columns if str(col).endswith("_clusters")]
    return str(candidates[0]) if candidates else None


def _infer_subset_label_series(
    subset: AnnData,
    *,
    subset_cluster_key: Optional[str],
    subset_label_key: Optional[str],
    key_added: str,
) -> pd.Series:
    """Infer per-cell subset labels from a direct obs column or marker-evidence table."""
    if subset_label_key and subset_label_key in subset.obs.columns:
        return subset.obs[subset_label_key].astype(str)

    labels = pd.Series("Unknown", index=subset.obs_names, dtype="object")
    if not subset_cluster_key or subset_cluster_key not in subset.obs.columns:
        return labels

    annotation_ns = (
        subset.uns.get("sclucid", {}).get("analysis", {}).get("annotation", {})
        if isinstance(subset.uns, dict)
        else {}
    )
    evidence_tables = [
        value
        for name, value in annotation_ns.items()
        if str(name).startswith(f"{key_added}_") and str(name).endswith("_marker_evidence")
    ]
    for table in evidence_tables:
        if isinstance(table, pd.DataFrame) and {"cluster", "marker_label"}.issubset(table.columns):
            label_map = table.set_index("cluster")["marker_label"].astype(str).to_dict()
            labels = subset.obs[subset_cluster_key].astype(str).map(label_map).fillna("Unknown")
            break
    return labels.astype(str)


def build_subset_annotation_reconciliation(
    adata: AnnData,
    subset_results: Dict[str, AnnData],
    *,
    lineage_key: str,
    global_subtype_key: Optional[str] = None,
    subset_label_key: Optional[str] = None,
    subset_cluster_key: Optional[str] = None,
    subset_exclude_key: str = "subset_review_exclude",
    key_added: str = "subset_annotation_refinement",
) -> pd.DataFrame:
    """
    Build a cell-level reconciliation table for subset annotation results.

    The table is the review boundary between lineage-gated subset annotation and
    global ``adata.obs``. It records subtype conflicts and manual exclusion flags
    without modifying global labels or dropping cells.
    """
    if lineage_key not in adata.obs.columns:
        raise KeyError(f"'{lineage_key}' not found in adata.obs.")

    rows: List[Dict[str, Any]] = []
    global_obs_names = set(adata.obs_names.astype(str))
    for lineage_label, subset in subset_results.items():
        cluster_key = subset_cluster_key or _infer_subset_cluster_key(subset, key_added)
        subset_labels = _infer_subset_label_series(
            subset,
            subset_cluster_key=cluster_key,
            subset_label_key=subset_label_key,
            key_added=key_added,
        )
        subset_clusters = (
            subset.obs[cluster_key].astype(str)
            if cluster_key and cluster_key in subset.obs.columns
            else pd.Series("Unknown", index=subset.obs_names, dtype="object")
        )
        exclude_flags = (
            subset.obs[subset_exclude_key].map(_as_bool_review_flag)
            if subset_exclude_key in subset.obs.columns
            else pd.Series(False, index=subset.obs_names, dtype=bool)
        )

        for obs_name in subset.obs_names.astype(str):
            if obs_name not in global_obs_names:
                continue
            global_lineage = adata.obs.at[obs_name, lineage_key]
            global_subtype = (
                adata.obs.at[obs_name, global_subtype_key]
                if global_subtype_key and global_subtype_key in adata.obs.columns
                else "Not_applicable"
            )
            subset_label = subset_labels.loc[obs_name]
            excluded = bool(exclude_flags.loc[obs_name])
            lineage_conflict = str(global_lineage) != str(lineage_label)
            subtype_conflict = (
                not _is_missing_label(global_subtype)
                and not _is_missing_label(subset_label)
                and str(global_subtype) != str(subset_label)
            )
            if excluded:
                action = "exclude_from_global_review"
            elif lineage_conflict:
                action = "review_lineage_conflict"
            elif subtype_conflict:
                action = "review_subtype_conflict"
            elif not _is_missing_label(subset_label):
                action = "accept_subset_label"
            else:
                action = "keep_global_label"

            rows.append(
                {
                    "obs_name": obs_name,
                    "lineage_label": str(lineage_label),
                    "global_lineage_label": str(global_lineage),
                    "global_subtype_label": str(global_subtype),
                    "subset_cluster": str(subset_clusters.loc[obs_name]),
                    "subset_label": str(subset_label),
                    "lineage_conflict": bool(lineage_conflict),
                    "subtype_conflict": bool(subtype_conflict),
                    "exclude_from_global": excluded,
                    "recommended_action": action,
                }
            )

    reconciliation = pd.DataFrame(rows)
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[f"{key_added}_reconciliation"] = reconciliation
    return reconciliation


def apply_subset_annotation_reconciliation(
    adata: AnnData,
    reconciliation: pd.DataFrame,
    *,
    target_key: str = "subtype_refined",
    global_subtype_key: Optional[str] = None,
    action_col: str = "recommended_action",
    subset_label_col: str = "subset_label",
    exclusion_key: str = "subset_refinement_exclude",
    action_key: str = "subset_refinement_action",
    apply_actions: Optional[List[str]] = None,
    overwrite: bool = False,
) -> AnnData:
    """
    Apply reviewed subset labels conservatively to global ``adata.obs``.

    This function never removes cells. Cells marked as exclusions are flagged in
    ``exclusion_key`` so downstream QC/filtering can make an explicit decision.
    """
    required = {"obs_name", action_col, subset_label_col, "exclude_from_global"}
    missing = required.difference(reconciliation.columns)
    if missing:
        raise KeyError(f"reconciliation is missing required columns: {sorted(missing)}")

    allowed_actions = set(apply_actions or ["accept_subset_label"])
    if global_subtype_key and global_subtype_key in adata.obs.columns:
        initial = adata.obs[global_subtype_key].astype(str)
    else:
        initial = pd.Series("Not_applicable", index=adata.obs_names, dtype="object")
    if target_key not in adata.obs.columns:
        adata.obs[target_key] = initial.astype("object")
    else:
        adata.obs[target_key] = adata.obs[target_key].astype("object")
    adata.obs[exclusion_key] = False
    adata.obs[action_key] = "not_in_subset_reconciliation"

    indexed = reconciliation.set_index("obs_name", drop=False)
    common = indexed.index.intersection(adata.obs_names.astype(str))
    for obs_name in common:
        row = indexed.loc[obs_name]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        action = str(row[action_col])
        label = row[subset_label_col]
        adata.obs.at[obs_name, action_key] = action
        adata.obs.at[obs_name, exclusion_key] = bool(row["exclude_from_global"])
        current = adata.obs.at[obs_name, target_key]
        can_write = action in allowed_actions and not _is_missing_label(label)
        if can_write and (overwrite or _is_missing_label(current)):
            adata.obs.at[obs_name, target_key] = str(label)

    adata.obs[target_key] = pd.Categorical(adata.obs[target_key].astype(str))
    adata.obs[action_key] = pd.Categorical(adata.obs[action_key].astype(str))
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[f"{target_key}_reconciliation_apply"] = sanitize_for_hdf5(
        {
            "target_key": target_key,
            "global_subtype_key": global_subtype_key,
            "action_col": action_col,
            "subset_label_col": subset_label_col,
            "exclusion_key": exclusion_key,
            "action_key": action_key,
            "apply_actions": sorted(allowed_actions),
            "overwrite": bool(overwrite),
            "n_rows": int(reconciliation.shape[0]),
            "n_applied": int(
                reconciliation[action_col].astype(str).isin(allowed_actions).sum()
            ),
            "n_excluded": int(reconciliation["exclude_from_global"].astype(bool).sum()),
        }
    )
    return adata


def run_subset_annotation_refinement(
    adata: AnnData,
    *,
    lineage_key: str,
    lineages: Optional[List[str]] = None,
    plan: Optional[pd.DataFrame] = None,
    cluster_resolution: float = 0.5,
    counts_layer: str = "counts",
    prefer_raw: bool = False,
    n_top_hvgs: int = 2000,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    min_cells: int = 50,
    marker_config: Optional[Union[str, Manager]] = None,
    global_subtype_key: Optional[str] = None,
    subset_exclude_key: str = "subset_review_exclude",
    write_back: bool = False,
    key_added: str = "subset_annotation_refinement",
) -> Dict[str, AnnData]:
    """
    Reprocess and recluster lineage subsets for subtype/state annotation.

    The default output is an artifact dictionary of subset AnnData objects plus a
    summary and reconciliation table in ``adata.uns``. Writing labels back to
    global ``adata.obs`` is optional because subset refinement often needs manual
    review first. Cells judged problematic in a subset should be marked in the
    subset with ``subset_exclude_key`` and reconciled before global application.
    """
    if lineage_key not in adata.obs.columns:
        raise KeyError(f"'{lineage_key}' not found in adata.obs.")
    if min_cells < 1:
        raise ValueError("min_cells must be >= 1.")
    if n_top_hvgs < 1:
        raise ValueError("n_top_hvgs must be >= 1.")

    lineage_series = adata.obs[lineage_key].astype(str)
    if lineages is None:
        if plan is not None and not plan.empty and "lineage_label" in plan.columns:
            if "recommended_action" in plan.columns:
                selected = plan.loc[
                    plan["recommended_action"].astype(str)
                    == "subset_recluster_for_subtype",
                    "lineage_label",
                ]
                lineages = selected.astype(str).tolist()
            else:
                lineages = plan["lineage_label"].astype(str).tolist()
        else:
            lineages = [
                item
                for item in sorted(lineage_series.dropna().unique(), key=str)
                if item not in {"Unknown", "nan", "None", "Not_applicable", ""}
            ]

    subset_results: Dict[str, AnnData] = {}
    summary_rows: List[Dict[str, Any]] = []
    marker_mgr = marker_config
    if isinstance(marker_config, str):
        marker_mgr = Manager(marker_config)

    if write_back:
        adata.obs[f"{key_added}_cluster"] = pd.Series(
            "Not_applicable", index=adata.obs_names, dtype="object"
        )
        if marker_config is not None:
            adata.obs[f"{key_added}_marker_label"] = pd.Series(
                "Not_applicable", index=adata.obs_names, dtype="object"
            )

    for lineage in lineages or []:
        mask = lineage_series == str(lineage)
        n_cells = int(mask.sum())
        result_key = _safe_obs_key(str(lineage))
        if n_cells < min_cells:
            summary_rows.append(
                {
                    "lineage_label": str(lineage),
                    "n_cells": n_cells,
                    "status": "skipped",
                    "reason": "too_few_cells",
                    "input_mode": "",
                    "cluster_key": "",
                    "n_clusters": 0,
                }
            )
            continue

        subset, input_mode = _extract_subset_for_reprocessing(
            adata,
            mask,
            counts_layer=counts_layer,
            prefer_raw=prefer_raw,
        )
        sc.pp.normalize_total(subset, target_sum=1e4)
        sc.pp.log1p(subset)
        n_hvgs = min(int(n_top_hvgs), max(1, subset.n_vars - 1))
        sc.pp.highly_variable_genes(subset, n_top_genes=n_hvgs, flavor="seurat")
        if "highly_variable" in subset.var.columns and subset.var["highly_variable"].any():
            subset = subset[:, subset.var["highly_variable"].to_numpy()].copy()
        sc.pp.scale(subset, max_value=10)
        n_comps = min(int(n_pcs), max(1, subset.n_obs - 1), max(1, subset.n_vars - 1))
        sc.tl.pca(subset, n_comps=n_comps)
        sc.pp.neighbors(
            subset,
            n_neighbors=min(int(n_neighbors), max(2, subset.n_obs - 1)),
            n_pcs=n_comps,
        )
        cluster_key = f"{key_added}_{result_key}_clusters"
        sc.tl.leiden(subset, resolution=float(cluster_resolution), key_added=cluster_key)

        marker_evidence = pd.DataFrame()
        if marker_mgr is not None:
            marker_evidence = run_marker_annotation_evidence(
                subset,
                cluster_key,
                marker_mgr,
                top_n_markers=50,
                use_raw=False,
                key_added=f"{key_added}_{result_key}_marker_evidence",
            )

        subset.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
            "annotation", {}
        )[key_added] = sanitize_for_hdf5(
            {
                "lineage_label": str(lineage),
                "input_mode": input_mode,
                "cluster_key": cluster_key,
                "write_back": bool(write_back),
            }
        )
        subset_results[str(lineage)] = subset

        n_clusters = int(subset.obs[cluster_key].nunique())
        summary_rows.append(
            {
                "lineage_label": str(lineage),
                "n_cells": n_cells,
                "status": "completed",
                "reason": "",
                "input_mode": input_mode,
                "cluster_key": cluster_key,
                "n_clusters": n_clusters,
                "n_hvgs": int(subset.n_vars),
                "marker_evidence_rows": int(marker_evidence.shape[0]),
            }
        )

        if write_back:
            cluster_values = subset.obs[cluster_key].astype(str)
            adata.obs.loc[mask, f"{key_added}_cluster"] = [
                f"{result_key}:{value}" for value in cluster_values
            ]
            if marker_config is not None and not marker_evidence.empty:
                label_map = marker_evidence.set_index("cluster")["marker_label"].astype(str).to_dict()
                labels = cluster_values.map(label_map).fillna("Unknown")
                adata.obs.loc[mask, f"{key_added}_marker_label"] = labels.to_numpy()

    if write_back:
        adata.obs[f"{key_added}_cluster"] = pd.Categorical(
            adata.obs[f"{key_added}_cluster"].astype(str)
        )
        if marker_config is not None:
            adata.obs[f"{key_added}_marker_label"] = pd.Categorical(
                adata.obs[f"{key_added}_marker_label"].astype(str)
            )

    if subset_results:
        build_subset_annotation_reconciliation(
            adata,
            subset_results,
            lineage_key=lineage_key,
            global_subtype_key=global_subtype_key,
            subset_exclude_key=subset_exclude_key,
            key_added=key_added,
        )

    summary = pd.DataFrame(summary_rows)
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[key_added] = sanitize_for_hdf5(summary)
    annotation_ns[f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "lineage_key": lineage_key,
            "lineages": lineages or [],
            "cluster_resolution": float(cluster_resolution),
            "counts_layer": counts_layer,
            "prefer_raw": bool(prefer_raw),
            "n_top_hvgs": int(n_top_hvgs),
            "n_pcs": int(n_pcs),
            "n_neighbors": int(n_neighbors),
            "min_cells": int(min_cells),
            "global_subtype_key": global_subtype_key,
            "subset_exclude_key": subset_exclude_key,
            "write_back": bool(write_back),
            "has_marker_config": marker_config is not None,
        }
    )
    return subset_results


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

            # Negative marker conflict detection
            negative_markers = [str(g) for g in getattr(cell, "negative_markers", [])]
            negative_hits: List[str] = []
            if negative_markers:
                negative_set = {gene.upper() for gene in negative_markers}
                negative_hits = [gene for gene in cluster_genes if gene.upper() in negative_set]

            confidence = 0.65 * recall + 0.35 * precision
            # Penalize confidence if negative markers are hit
            if negative_hits:
                penalty = min(0.5, len(negative_hits) * 0.15)
                confidence *= (1.0 - penalty)

            if overlap >= min_overlap:
                scored_labels.append(
                    {
                        "label": str(label),
                        "overlap": int(overlap),
                        "recall": float(recall),
                        "precision": float(precision),
                        "confidence": float(confidence),
                        "matched_markers": matched,
                        "negative_markers_hit": negative_hits,
                        "n_negative_hits": len(negative_hits),
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

        negative_hit_str = (
            ", ".join(winner["negative_markers_hit"][:8])
            if winner and winner.get("negative_markers_hit")
            else ""
        )
        n_negative_hits = int(winner["n_negative_hits"]) if winner else 0
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
                "negative_markers_hit": negative_hit_str,
                "n_negative_hits": n_negative_hits,
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
    annotation_ns[key_added] = sanitize_for_hdf5(evidence_df)
    annotation_ns[f"{key_added}_marker_table"] = sanitize_for_hdf5(marker_table)
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


def run_program_annotation_evidence(
    adata: AnnData,
    cluster_key: str,
    program_config: Optional[Union[str, Manager]] = None,
    *,
    species: str = "human",
    tissue: Optional[str] = None,
    use_raw: bool = False,
    min_genes: int = 3,
    top_n_programs: int = 5,
    score_suffix: str = "_program_score",
    key_added: str = "program_annotation_evidence",
) -> pd.DataFrame:
    """
    Score geneset/program evidence by cluster for annotation interpretation.

    Program evidence answers "what is this cluster doing?" rather than "what
    cell identity is this cluster?". The returned table is intended for review,
    visualization, and LLM annotation context, not direct final-label voting.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")
    if min_genes < 1:
        raise ValueError("min_genes must be >= 1.")
    if top_n_programs < 1:
        raise ValueError("top_n_programs must be >= 1.")

    if program_config is None:
        mgr = get_marker_manager(species=species, tissue=tissue, view="program_scoring")
    elif isinstance(program_config, str):
        mgr = Manager(program_config)
    elif isinstance(program_config, Manager):
        mgr = program_config
    else:
        raise TypeError("program_config must be None, a file path, or Manager instance.")

    source_var_names = (
        set(adata.raw.var_names.astype(str))
        if use_raw and adata.raw is not None
        else set(adata.var_names.astype(str))
    )
    score_columns: Dict[str, str] = {}
    program_metadata: Dict[str, Dict[str, Any]] = {}
    skipped_rows: List[Dict[str, Any]] = []
    for program_name, cell in mgr.CELLS.items():
        genes = [str(g) for g in getattr(cell, "markers", []) if isinstance(g, str)]
        matched = [gene for gene in genes if gene in source_var_names]
        metadata = _resolve_marker_cell_metadata(cell)
        program_metadata[str(program_name)] = metadata
        if len(matched) < min_genes:
            skipped_rows.append(
                {
                    "program": str(program_name),
                    "n_genes": len(genes),
                    "n_genes_used": len(matched),
                    "reason": "too_few_matched_genes",
                }
            )
            continue
        score_col = f"{_safe_obs_key(str(program_name))}{score_suffix}"
        # Ensure repeated calls remain deterministic without clobbering unrelated columns.
        if score_col in adata.obs.columns:
            adata.obs = adata.obs.drop(columns=[score_col])
        sc.tl.score_genes(
            adata,
            gene_list=matched,
            score_name=score_col,
            use_raw=use_raw and adata.raw is not None,
        )
        score_columns[str(program_name)] = score_col

    cluster_series = adata.obs[cluster_key].astype(str)
    rows: List[Dict[str, Any]] = []
    for cluster in cluster_series.drop_duplicates().tolist():
        mask = cluster_series == str(cluster)
        for program_name, score_col in score_columns.items():
            values = pd.to_numeric(adata.obs.loc[mask, score_col], errors="coerce")
            metadata = program_metadata.get(program_name, {})
            rows.append(
                {
                    "cluster": str(cluster),
                    "program": str(program_name),
                    "program_score_mean": float(values.mean()) if values.notna().any() else np.nan,
                    "program_score_median": (
                        float(values.median()) if values.notna().any() else np.nan
                    ),
                    "program_score_fraction_positive": (
                        float((values > 0).mean()) if values.notna().any() else np.nan
                    ),
                    "score_column": score_col,
                    "n_cells": int(mask.sum()),
                    "n_genes_used": int(
                        len([g for g in mgr.CELLS[program_name].markers if str(g) in source_var_names])
                    ),
                    "category": metadata.get("category", ""),
                    "source_collection": metadata.get("source_collection", ""),
                    "source_ids": ",".join(map(str, metadata.get("source_ids", [])))
                    if isinstance(metadata.get("source_ids", []), list)
                    else str(metadata.get("source_ids", "")),
                    "not_for": ",".join(map(str, metadata.get("not_for", [])))
                    if isinstance(metadata.get("not_for", []), list)
                    else str(metadata.get("not_for", "")),
                }
            )

    evidence_df = pd.DataFrame(rows)
    if not evidence_df.empty:
        evidence_df = evidence_df.sort_values(
            ["cluster", "program_score_mean", "program"],
            ascending=[True, False, True],
        ).copy()
        evidence_df["program_rank"] = evidence_df.groupby("cluster").cumcount() + 1
        evidence_df["is_top_program"] = evidence_df["program_rank"] <= int(top_n_programs)
    else:
        evidence_df = pd.DataFrame(
            columns=[
                "cluster",
                "program",
                "program_score_mean",
                "program_rank",
                "is_top_program",
            ]
        )

    top_summary_rows = []
    if not evidence_df.empty:
        for cluster, sub in evidence_df.loc[evidence_df["is_top_program"]].groupby(
            "cluster", observed=False
        ):
            top_summary_rows.append(
                {
                    "cluster": str(cluster),
                    "top_programs": "; ".join(
                        f"{row.program}:{row.program_score_mean:.3f}"
                        for row in sub.itertuples(index=False)
                    ),
                }
            )
    top_summary = pd.DataFrame(top_summary_rows, columns=["cluster", "top_programs"])
    if not top_summary.empty:
        evidence_df = evidence_df.merge(top_summary, on="cluster", how="left")
    else:
        evidence_df["top_programs"] = ""

    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[key_added] = sanitize_for_hdf5(evidence_df)
    annotation_ns[f"{key_added}_skipped"] = sanitize_for_hdf5(pd.DataFrame(skipped_rows))
    annotation_ns[f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "species": species,
            "tissue": tissue,
            "use_raw": bool(use_raw and adata.raw is not None),
            "min_genes": int(min_genes),
            "top_n_programs": int(top_n_programs),
            "score_suffix": score_suffix,
            "n_programs_scored": int(len(score_columns)),
            "n_programs_skipped": int(len(skipped_rows)),
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
    program_evidence: Optional[pd.DataFrame] = None,
    reference_key: Optional[str] = None,
    lineage_key: Optional[str] = None,
    subtype_key: Optional[str] = None,
    state_key: Optional[str] = None,
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
    program_evidence_by_cluster: Dict[str, Any] = {}
    if (
        isinstance(program_evidence, pd.DataFrame)
        and not program_evidence.empty
        and {"cluster", "program", "program_score_mean"}.issubset(program_evidence.columns)
    ):
        top_programs = program_evidence.copy()
        if "is_top_program" in top_programs.columns:
            top_programs = top_programs[top_programs["is_top_program"].astype(bool)]
        for cluster, sub in top_programs.groupby("cluster", observed=False):
            program_evidence_by_cluster[str(cluster)] = [
                {
                    "program": str(row["program"]),
                    "score_mean": float(row["program_score_mean"])
                    if pd.notna(row["program_score_mean"])
                    else None,
                    "category": str(row.get("category", "")),
                    "source_collection": str(row.get("source_collection", "")),
                }
                for _, row in sub.head(8).iterrows()
            ]

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

        lineage_summary = None
        if lineage_key and lineage_key in obs_subset.columns:
            lineage_counts = (
                obs_subset[lineage_key].astype(str).value_counts(normalize=True).head(5)
            )
            lineage_summary = [
                {"label": str(label), "fraction": float(frac)}
                for label, frac in lineage_counts.items()
            ]

        subtype_summary = None
        if subtype_key and subtype_key in obs_subset.columns:
            subtype_counts = (
                obs_subset[subtype_key].astype(str).value_counts(normalize=True).head(5)
            )
            subtype_summary = [
                {"label": str(label), "fraction": float(frac)}
                for label, frac in subtype_counts.items()
            ]

        state_summary = None
        if state_key and state_key in obs_subset.columns:
            state_counts = obs_subset[state_key].astype(str).value_counts(normalize=True).head(5)
            state_summary = [
                {"label": str(label), "fraction": float(frac)}
                for label, frac in state_counts.items()
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
            "program_evidence": program_evidence_by_cluster.get(str(cluster), []),
            "lineage_annotation": lineage_summary,
            "subtype_annotation": subtype_summary,
            "state_annotation": state_summary,
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
        "Assign conservative labels using hierarchical annotation logic. First respect the "
        "major lineage evidence; only assign subtype/state labels when the lineage gate and "
        "marker evidence agree. Use program/geneset evidence to interpret cell state, function, "
        "tumor context, and artifacts, but do not use program names as direct cell identity "
        "labels. Prefer broad lineage labels when evidence is mixed. Do not use "
        "ribosomal, mitochondrial, stress, or housekeeping genes as primary lineage evidence. "
        "Return a table with cluster, llm_label, llm_confidence from 0 to 1, evidence, "
        "conflicts, and needs_review."
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
            "program_evidence": isinstance(program_evidence, pd.DataFrame)
            and not program_evidence.empty,
            "lineage_key": lineage_key,
            "subtype_key": subtype_key,
            "state_key": state_key,
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


def evaluate_annotation_benchmark(
    adata: AnnData,
    *,
    label_key: str = "cell_type_final",
    truth_key: Optional[str] = None,
    cluster_key: Optional[str] = None,
    review_table: Optional[pd.DataFrame] = None,
    evidence_label_cols: Optional[List[str]] = None,
    key_added: str = "annotation_benchmark",
) -> Dict[str, Any]:
    """Build benchmark and disagreement summaries for annotation review.

    ``truth_key`` is optional. When absent, the function still reports
    disagreement among reference/marker/LLM/final labels and conservative label
    categories, which is useful before a gold-standard benchmark is available.
    """
    if label_key not in adata.obs.columns and review_table is None:
        raise KeyError(f"'{label_key}' not found in adata.obs and no review_table supplied.")

    if review_table is None:
        if cluster_key and cluster_key in adata.obs.columns:
            cluster_series = adata.obs[cluster_key].astype(str)
            rows = []
            for cluster in cluster_series.drop_duplicates().tolist():
                mask = cluster_series == str(cluster)
                labels = adata.obs.loc[mask, label_key].astype(str)
                rows.append(
                    {
                        "cluster": str(cluster),
                        "final_label": labels.value_counts().index[0],
                        "n_cells": int(mask.sum()),
                    }
                )
            review = pd.DataFrame(rows)
        else:
            review = pd.DataFrame(
                {
                    "cluster": adata.obs.index.astype(str),
                    "final_label": adata.obs[label_key].astype(str).to_numpy(),
                    "n_cells": 1,
                }
            )
    else:
        review = review_table.copy()

    if "final_label" not in review.columns:
        if label_key in review.columns:
            review["final_label"] = review[label_key]
        else:
            review["final_label"] = "Unknown"
    review["final_label"] = review["final_label"].fillna("Unknown").astype(str)

    default_evidence_cols = [
        "reference_label",
        "marker_label",
        "llm_label",
        "lineage_label",
        "subtype_label",
        "state_label",
    ]
    evidence_cols = evidence_label_cols or [c for c in default_evidence_cols if c in review.columns]

    disagreement_rows = []
    cols_for_matrix = evidence_cols + ["final_label"]
    for left in cols_for_matrix:
        for right in cols_for_matrix:
            if left >= right:
                continue
            valid = review[left].notna() & review[right].notna()
            valid &= ~review[left].astype(str).isin(["Unknown", "nan", "None"])
            valid &= ~review[right].astype(str).isin(["Unknown", "nan", "None"])
            n = int(valid.sum())
            disagreement = (
                float((review.loc[valid, left].astype(str) != review.loc[valid, right].astype(str)).mean())
                if n
                else np.nan
            )
            disagreement_rows.append(
                {
                    "source_a": left,
                    "source_b": right,
                    "n_compared": n,
                    "disagreement_rate": disagreement,
                }
            )
    disagreement_df = pd.DataFrame(disagreement_rows)

    def _category(row: pd.Series) -> str:
        final = str(row.get("final_label", "Unknown"))
        if final in {"Unknown", "nan", "None", ""}:
            return "Unknown"
        labels = [
            str(row[c])
            for c in evidence_cols
            if c in row and pd.notna(row[c]) and str(row[c]) not in {"Unknown", "nan", "None", ""}
        ]
        unique = set(labels)
        if len(unique) >= 2 and final not in unique:
            return "Mixed"
        if len(unique) >= 2:
            return "Ambiguous"
        if bool(row.get("needs_review", False)):
            return "Ambiguous"
        return "Resolved"

    review["conservative_label_category"] = review.apply(_category, axis=1)
    category_counts = review["conservative_label_category"].value_counts().to_dict()

    confusion = pd.DataFrame()
    accuracy = np.nan
    if truth_key and truth_key in adata.obs.columns and label_key in adata.obs.columns:
        truth = adata.obs[truth_key].astype(str)
        pred = adata.obs[label_key].astype(str)
        confusion = pd.crosstab(truth, pred, rownames=["truth"], colnames=["predicted"])
        accuracy = float((truth == pred).mean()) if len(truth) else np.nan

    provenance_cols = [
        c
        for c in review.columns
        if c.endswith("_source") or "provenance" in c or c in {"marker_database", "reference_model"}
    ]
    provenance_summary = {
        col: review[col].astype(str).value_counts().head(10).to_dict() for col in provenance_cols
    }

    result = {
        "schema_version": "annotation_benchmark_v1",
        "label_key": label_key,
        "truth_key": truth_key,
        "n_units": int(len(review)),
        "conservative_label_counts": sanitize_for_hdf5(category_counts),
        "disagreement_matrix": sanitize_for_hdf5(disagreement_df),
        "confusion_matrix": sanitize_for_hdf5(confusion),
        "accuracy": accuracy,
        "marker_provenance_summary": sanitize_for_hdf5(provenance_summary),
        "review_table": sanitize_for_hdf5(review),
        "recommendation": (
            "Use Unknown/Ambiguous/Mixed labels for clusters with disagreement before publication."
        ),
    }
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})[
        key_added
    ] = result
    return result


def merge_annotation_evidence(
    adata: AnnData,
    cluster_key: str,
    *,
    marker_evidence: Optional[pd.DataFrame] = None,
    reference_key: Optional[str] = None,
    reference_confidence_key: Optional[str] = None,
    llm_annotations: Optional[Union[pd.DataFrame, Dict[str, Any]]] = None,
    review_table: Optional[pd.DataFrame] = None,
    suspect_flags: Optional[pd.DataFrame] = None,
    program_evidence: Optional[pd.DataFrame] = None,
    lineage_evidence: Optional[pd.DataFrame] = None,
    subtype_evidence: Optional[pd.DataFrame] = None,
    state_evidence: Optional[pd.DataFrame] = None,
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

    if (
        program_evidence is not None
        and not program_evidence.empty
        and {"cluster", "top_programs"}.issubset(program_evidence.columns)
    ):
        program_summary = (
            program_evidence[["cluster", "top_programs"]]
            .dropna()
            .drop_duplicates(subset=["cluster"], keep="first")
            .copy()
        )
        merged = merged.merge(program_summary, on="cluster", how="left")
    else:
        merged["top_programs"] = ""

    # Merge suspect cluster flags if available
    if suspect_flags is not None and not suspect_flags.empty and "cluster" in suspect_flags.columns:
        flag_cols = [c for c in ["cluster", "suspect_flag", "suspect_reasons"] if c in suspect_flags.columns]
        # Drop existing flag columns from merged to avoid duplicates
        for col in flag_cols:
            if col in merged.columns and col != "cluster":
                merged = merged.drop(columns=[col])
        merged = merged.merge(suspect_flags[flag_cols], on="cluster", how="left")
    else:
        merged["suspect_flag"] = "clean"
        merged["suspect_reasons"] = ""

    # Merge hierarchical annotation evidence (lineage / subtype / state)
    for evidence_df, label_col, conf_col in [
        (lineage_evidence, "lineage_label", "lineage_confidence"),
        (subtype_evidence, "subtype_label", "subtype_confidence"),
        (state_evidence, "state_label", "state_confidence"),
    ]:
        if evidence_df is not None and not evidence_df.empty and "cluster" in evidence_df.columns:
            sub = evidence_df[["cluster", "marker_label", "marker_confidence"]].copy()
            sub = sub.rename(columns={"marker_label": label_col, "marker_confidence": conf_col})
            merged = merged.merge(sub, on="cluster", how="left")
        else:
            merged[label_col] = "Unknown"
            merged[conf_col] = np.nan

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

        # Apply suspect-flag penalties before threshold check
        suspect_flag = str(row.get("suspect_flag", "clean"))
        suspect_penalties = {
            "doublet_suspect": 0.3,
            "stress_high": 0.6,
            "ribosomal_dominant": 0.6,
            "low_information": 0.7,
            "mt_high": 0.8,
        }
        if suspect_flag in suspect_penalties:
            confidence *= suspect_penalties[suspect_flag]
            decision = f"{decision}_with_{suspect_flag}"

        if confidence < min_final_confidence:
            final_label = "Unknown"
            decision = "below_confidence_threshold"

        warnings = []
        if conflicts:
            warnings.append("evidence_conflict")
        if row.get("marker_label") == "Unknown" and row.get("reference_label") == "Unknown":
            warnings.append("weak_marker_reference_evidence")
        if suspect_flag != "clean":
            warnings.append(suspect_flag)
        needs_review = bool(
            conflicts or final_label == "Unknown" or confidence < 0.5 or suspect_flag != "clean"
        )
        final_rows.append(
            {
                "final_label": final_label,
                "lineage_label": str(row.get("lineage_label", "Unknown")),
                "lineage_confidence": float(row.get("lineage_confidence", 0.0)) if pd.notna(row.get("lineage_confidence")) else 0.0,
                "subtype_label": str(row.get("subtype_label", "Unknown")),
                "subtype_confidence": float(row.get("subtype_confidence", 0.0)) if pd.notna(row.get("subtype_confidence")) else 0.0,
                "state_label": str(row.get("state_label", "Unknown")),
                "state_confidence": float(row.get("state_confidence", 0.0)) if pd.notna(row.get("state_confidence")) else 0.0,
                "annotation_confidence": float(confidence),
                "decision": decision,
                "conflicts": ",".join(conflicts),
                "warnings": ",".join(warnings),
                "needs_review": needs_review,
                "suspect_flag": suspect_flag,
            }
        )

    final_df = pd.concat([merged.reset_index(drop=True), pd.DataFrame(final_rows)], axis=1)
    for col in ANNOTATION_REVIEW_SCHEMA:
        if col not in final_df.columns:
            final_df[col] = "" if col not in {"needs_review", "n_cells", "pct_cells"} else np.nan
    # Drop duplicate columns before storing
    final_df = final_df.loc[:, ~final_df.columns.duplicated()]
    final_df = final_df[
        ANNOTATION_REVIEW_SCHEMA
        + [c for c in final_df.columns if c not in ANNOTATION_REVIEW_SCHEMA]
    ]

    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[key_added] = sanitize_for_hdf5(final_df)
    annotation_ns[f"{key_added}_schema"] = list(ANNOTATION_REVIEW_SCHEMA)
    annotation_ns[f"{key_added}_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "reference_key": reference_key,
            "reference_confidence_key": reference_confidence_key,
            "min_final_confidence": float(min_final_confidence),
            "prefer_llm_when_confident": bool(prefer_llm_when_confident),
            "has_program_evidence": bool(
                program_evidence is not None and not program_evidence.empty
            ),
        }
    )

    # P3: Post-hoc doublet evidence backflow to QC namespace
    posthoc_rows = []
    for _, row in final_df.iterrows():
        warnings_str = str(row.get("warnings", ""))
        suspect = str(row.get("suspect_flag", ""))
        if "doublet_suspect" in warnings_str or "doublet_suspect" in suspect:
            posthoc_rows.append(
                {
                    "cluster": str(row["cluster"]),
                    "evidence_source": "annotation_suspect_flag",
                    "details": f"suspect_flag={suspect}, warnings={warnings_str}",
                    "recommendation": "inspect_doublet_fraction_in_qc",
                }
            )
        n_neg = row.get("n_negative_hits")
        if n_neg is not None and int(n_neg) > 0:
            posthoc_rows.append(
                {
                    "cluster": str(row["cluster"]),
                    "evidence_source": "negative_marker_conflict",
                    "details": f"n_negative_hits={int(n_neg)}, hits={row.get('negative_markers_hit', '')}",
                    "recommendation": "review_cluster_annotation_and_qc",
                }
            )
    if posthoc_rows:
        adata.uns.setdefault("sclucid", {}).setdefault("qc", {})[
            "posthoc_doublet_evidence"
        ] = sanitize_for_hdf5(pd.DataFrame(posthoc_rows).to_dict(orient="records"))

    return final_df


def run_annotation_evidence(
    adata: AnnData,
    cluster_key: str,
    *,
    markers_df: Optional[pd.DataFrame] = None,
    methods: tuple[str, ...] = ("reference", "marker_manager", "data_driven"),
    marker_config: Optional[Union[str, Manager]] = None,
    reference_key: Optional[str] = None,
    reference_confidence_key: Optional[str] = None,
    llm_annotations: Optional[Union[pd.DataFrame, Dict[str, Any]]] = None,
    llm_annotator: Optional[Any] = None,
    enrichment_dict: Optional[Dict[str, pd.DataFrame]] = None,
    sample_col: Optional[str] = None,
    group_col: Optional[str] = None,
    top_n_markers: int = 15,
    hierarchical: bool = False,
    species: str = "human",
    tissue: Optional[str] = None,
    key_added: str = "annotation_review_table",
) -> pd.DataFrame:
    """
    Build cluster-level annotation evidence without treating any path as truth.

    The workflow can combine reference labels already present in ``obs``,
    marker-manager overlap evidence, and data-driven LLM suggestions. LLM use is
    callback-based: pass either precomputed ``llm_annotations`` or a callable that
    accepts one cluster payload and returns a dict with ``llm_label`` and optional
    ``llm_confidence``.

    When ``hierarchical=True`` and ``marker_config`` is None, the function uses
    built-in marker manager views (``lineage_annotation``, ``subtype_annotation``,
    ``state_annotation``) to produce layered evidence rather than a single label.
    """
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")

    active_methods = set(methods or ())
    marker_evidence = None
    lineage_evidence = None
    subtype_evidence = None
    state_evidence = None
    program_evidence = None

    if "marker_manager" in active_methods:
        if hierarchical and marker_config is None:
            lineage_mgr = get_marker_manager(species=species, tissue=tissue, view="lineage_annotation")
            lineage_evidence = run_marker_annotation_evidence(
                adata,
                cluster_key,
                lineage_mgr,
                markers_df=markers_df,
                top_n_markers=max(top_n_markers, 20),
                key_added="lineage_annotation_evidence",
            )
            subtype_mgr = get_marker_manager(species=species, tissue=tissue, view="subtype_annotation")
            subtype_evidence = run_marker_annotation_evidence(
                adata,
                cluster_key,
                subtype_mgr,
                markers_df=markers_df,
                top_n_markers=max(top_n_markers, 20),
                key_added="subtype_annotation_evidence",
            )
            state_mgr = get_marker_manager(species=species, tissue=tissue, view="state_annotation")
            state_evidence = run_marker_annotation_evidence(
                adata,
                cluster_key,
                state_mgr,
                markers_df=markers_df,
                top_n_markers=max(top_n_markers, 20),
                key_added="state_annotation_evidence",
            )
        elif marker_config is not None:
            marker_evidence = run_marker_annotation_evidence(
                adata,
                cluster_key,
                marker_config,
                markers_df=markers_df,
                top_n_markers=max(top_n_markers, 20),
            )

    if "program" in active_methods or "program_scoring" in active_methods:
        program_evidence = run_program_annotation_evidence(
            adata,
            cluster_key,
            species=species,
            tissue=tissue,
            top_n_programs=5,
            key_added="program_annotation_evidence",
        )

    bundle = None
    if "data_driven" in active_methods:
        bundle = build_llm_annotation_bundle(
            adata,
            cluster_key,
            markers_df=markers_df,
            enrichment_dict=enrichment_dict,
            marker_evidence=marker_evidence,
            program_evidence=program_evidence,
            reference_key=reference_key,
            top_n_markers=top_n_markers,
            sample_col=sample_col,
            group_col=group_col,
        )
        if llm_annotations is None and llm_annotator is not None:
            llm_rows = []
            for cluster, payload in bundle.get("clusters", {}).items():
                result = llm_annotator({"cluster": cluster, **payload})
                if result is None:
                    continue
                if not isinstance(result, dict):
                    result = {"llm_label": str(result)}
                row = {"cluster": str(cluster), **result}
                if "label" in row and "llm_label" not in row:
                    row["llm_label"] = row.pop("label")
                if "confidence" in row and "llm_confidence" not in row:
                    row["llm_confidence"] = row.pop("confidence")
                llm_rows.append(row)
            llm_annotations = pd.DataFrame(llm_rows)

    review = merge_annotation_evidence(
        adata,
        cluster_key,
        marker_evidence=marker_evidence,
        reference_key=reference_key
        if {"reference", "celltypist"} & active_methods
        else None,
        reference_confidence_key=reference_confidence_key,
        llm_annotations=llm_annotations if "data_driven" in active_methods else None,
        program_evidence=program_evidence,
        lineage_evidence=lineage_evidence,
        subtype_evidence=subtype_evidence,
        state_evidence=state_evidence,
        min_final_confidence=0.2,
        key_added=key_added,
    )
    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns[f"{key_added}_evidence_methods"] = sanitize_for_hdf5(sorted(active_methods))
    return review


def build_annotation_consensus(
    adata: AnnData,
    cluster_key: str,
    annotation_review_table: Optional[pd.DataFrame] = None,
    *,
    key_added: str = "cell_type",
    lineage_key: Optional[str] = "lineage",
    label_col: str = "final_label",
) -> pd.DataFrame:
    """
    Apply consensus annotation labels and return the review table.

    This is a semantic wrapper around ``apply_final_annotation`` for the standard
    analysis workflow. It keeps final labels cell-level while preserving the
    cluster-level evidence table in ``uns``.
    """
    if annotation_review_table is None:
        annotation_review_table = (
            adata.uns.get("sclucid", {})
            .get("analysis", {})
            .get("annotation", {})
            .get("annotation_review_table")
        )
    if not isinstance(annotation_review_table, pd.DataFrame):
        raise ValueError("annotation_review_table must be provided or stored in adata.uns.")

    apply_final_annotation(
        adata,
        cluster_key,
        annotation_review_table,
        label_col=label_col,
        key_added=key_added,
    )
    # Backward-compatible lineage alias: expose the caller's requested key even
    # when apply_final_annotation produced the canonical ``{key_added}_lineage``.
    lineage_obs_key = f"{key_added}_lineage"
    if lineage_key and lineage_key not in adata.obs.columns:
        source_key = lineage_obs_key if lineage_obs_key in adata.obs.columns else key_added
        adata.obs[lineage_key] = adata.obs[source_key].astype(str)

    annotation_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annotation_ns["annotation_consensus_table"] = sanitize_for_hdf5(annotation_review_table)
    annotation_ns["annotation_consensus_params"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "key_added": key_added,
            "lineage_key": lineage_key,
            "label_col": label_col,
        }
    )
    return annotation_review_table


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

    # Write hierarchical labels if available in review table
    def _review_column(name: str) -> pd.Series:
        values = annotation_review_table.loc[:, name]
        if isinstance(values, pd.DataFrame):
            values = values.iloc[:, 0]
        return values

    for layer_col, layer_conf_col in [
        ("lineage_label", "lineage_confidence"),
        ("subtype_label", "subtype_confidence"),
        ("state_label", "state_confidence"),
    ]:
        if layer_col in annotation_review_table.columns:
            layer_map = pd.Series(
                _review_column(layer_col).astype(str).to_numpy(),
                index=annotation_review_table["cluster"],
            ).to_dict()
            layer_obs_key = f"{key_added}_{layer_col.replace('_label', '')}"
            adata.obs[layer_obs_key] = pd.Categorical(
                cluster_series.map(layer_map).fillna("Unknown")
            )
            if layer_conf_col in annotation_review_table.columns:
                conf_map = pd.Series(
                    _review_column(layer_conf_col).to_numpy(),
                    index=annotation_review_table["cluster"],
                ).to_dict()
                adata.obs[f"{layer_obs_key}_confidence"] = pd.to_numeric(
                    cluster_series.map(conf_map), errors="coerce"
                )

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
