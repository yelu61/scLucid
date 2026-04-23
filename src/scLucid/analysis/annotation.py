"""
Cell type annotation functions for single-cell RNA-seq data.
Supports: score-based, enrichment-based, combined, CellTypist, reference transfer,
manual/AI mapping import, and annotation evaluation.
"""

import logging
from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ..utils import Manager, get_marker_manager, sanitize_for_hdf5, use_layer_as_X
from .config import AnnotationConfig
from .scoring import FunctionalSignatureManager, score_by_gene_sets

log = logging.getLogger(__name__)


__all__ = [
    "score_cell_types",
    "annotate_clusters",
    "run_celltypist",
    "transfer_labels",
    "evaluate_annotation",
    "summarize_annotation_evidence",
    "run_lineage_state_annotation",
    "filter_marker_table_for_annotation",
    "flag_suspect_clusters",
    "build_annotation_review_table",
    "apply_annotation_mapping",
    "remap_labels",
    "run_annotation",
]


_ANNOTATION_NOISE_EXACT_GENES = {
    "MALAT1": "housekeeping",
    "NEAT1": "housekeeping",
    "XIST": "housekeeping",
}

_ANNOTATION_NOISE_PREFIXES = {
    "ribosomal": ("RPL", "RPS", "MRPL", "MRPS"),
    "mitochondrial": ("MT-",),
    "stress": (
        "HSP",
        "HSPA",
        "HSPB",
        "HSPC",
        "HSPD",
        "HSPE",
        "HSPH",
        "DNAJ",
        "FOS",
        "JUN",
        "ATF3",
        "IER",
        "DDIT3",
        "PPP1R15A",
    ),
}


# --------------------- Helper Functions -----------------------


def _get_default_compartment_map() -> Dict[str, str]:
    """Return a default compartment mapper for common cell type labels."""
    return {
        # Tumor / malignant
        "Malignant": "tumor",
        "Tumor": "tumor",
        "Cancer": "tumor",
        "Carcinoma": "tumor",
        "Epithelial": "tumor",
        "Pan-Cancer": "tumor",
        "LUAD": "tumor",
        "LUSC": "tumor",
        "SCLC": "tumor",
        # Immune
        "T cells": "immune",
        "B cells": "immune",
        "NK cells": "immune",
        "Macrophage": "immune",
        "Monocyte": "immune",
        "DC": "immune",
        "Neutrophil": "immune",
        "Mast": "immune",
        "Plasma": "immune",
        # Stromal
        "Fibroblast": "stromal",
        "Endothelial": "stromal",
        "Pericyte": "stromal",
        "Smooth muscle": "stromal",
        "Stromal": "stromal",
        # Uncertain
        "Unknown": "uncertain",
        "Unassigned": "uncertain",
    }


def _map_compartments(
    labels: pd.Series,
    compartment_map: Optional[Dict[str, str]] = None,
) -> pd.Series:
    """Map cell type labels to broad compartments (tumor/immune/stromal/uncertain/mixed)."""
    if compartment_map is None:
        compartment_map = _get_default_compartment_map()

    def _map_label(label: str) -> str:
        if pd.isna(label):
            return "uncertain"
        label_str = str(label)
        # Exact match first
        if label_str in compartment_map:
            return compartment_map[label_str]
        # Partial match fallback
        for key, compartment in compartment_map.items():
            if key.lower() in label_str.lower():
                return compartment
        return "uncertain"

    mapped = labels.apply(_map_label)
    return mapped


def _read_table_file(path: str) -> pd.DataFrame:
    """
    Read a tabular mapping file with robust handling:
    - Supports .xlsx/.xls (preferred) and .csv
    - For CSV, tries common encodings, then chardet probing as fallback
    Returns a DataFrame.
    """
    import os

    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        # Excel is robust against encoding issues
        return pd.read_excel(path)
    elif ext == ".csv":
        # Try multiple encodings
        last_err = None
        encodings = ["utf-8", "utf-8-sig", "gbk", "gb2312", "big5", "cp1252", "latin1"]
        for enc in encodings:
            try:
                return pd.read_csv(path, encoding=enc)
            except Exception as e:
                last_err = e
        # chardet fallback
        try:
            import chardet

            with open(path, "rb") as f:
                raw = f.read(200000)
            guess = chardet.detect(raw).get("encoding") or "utf-8"
            return pd.read_csv(path, encoding=guess, errors="replace")
        except Exception:
            pass
        # final fallback
        try:
            return pd.read_csv(path, encoding="latin1", errors="replace")
        except Exception:
            raise (
                last_err
                if last_err
                else RuntimeError(f"Failed to read CSV with multiple encodings: {path}")
            )
    else:
        raise ValueError("Unsupported table format. Use .xlsx, .xls, or .csv")


def _read_json_file(path: str) -> dict:
    """
    Read JSON with robust encoding handling.
    """
    import json

    last_err = None
    for enc in ["utf-8", "utf-8-sig", "gbk", "cp1252", "latin1"]:
        try:
            with open(path, encoding=enc) as f:
                return json.load(f)
        except Exception as e:
            last_err = e
    try:
        import chardet

        with open(path, "rb") as fb:
            raw = fb.read(200000)
        guess = chardet.detect(raw).get("encoding") or "utf-8"
        with open(path, encoding=guess, errors="replace") as f:
            return json.load(f)
    except Exception:
        pass
    raise (last_err if last_err else RuntimeError("Failed to read JSON with multiple encodings."))


def _classify_annotation_marker(gene: Any) -> Optional[str]:
    """Classify common non-informative marker genes seen during manual annotation."""
    if pd.isna(gene):
        return None

    gene_upper = str(gene).upper()
    if gene_upper in _ANNOTATION_NOISE_EXACT_GENES:
        return _ANNOTATION_NOISE_EXACT_GENES[gene_upper]

    for category, prefixes in _ANNOTATION_NOISE_PREFIXES.items():
        if any(gene_upper.startswith(prefix) for prefix in prefixes):
            return category
    return None


def _resolve_score_columns(
    adata: AnnData,
    score_cols: Optional[Sequence[str]] = None,
) -> List[str]:
    """Resolve module score columns from obs."""
    if score_cols is not None:
        return [col for col in score_cols if col in adata.obs.columns]
    return [
        col
        for col in adata.obs.columns
        if col.endswith("_score") and pd.api.types.is_numeric_dtype(adata.obs[col])
    ]


def _resolve_annotation_manager(
    *,
    species: str,
    tissue: Optional[str],
    states: Optional[List[str]] = None,
    marker_config: Optional[str] = None,
):
    """Resolve either a custom marker file or the built-in combined marker manager."""
    if marker_config:
        return Manager(marker_config, case_sensitive=True)
    return get_marker_manager(species=species, tissue=tissue, states=states)


def _label_matches_target_lineage(labels: pd.Series, target_lineage: Optional[str]) -> pd.Series:
    """Return a boolean mask for cells matching the requested target lineage."""
    if not target_lineage:
        return pd.Series(True, index=labels.index)

    target = str(target_lineage).strip().lower()
    values = labels.astype(str).str.lower()
    return values.eq(target) | values.str.contains(target, regex=False)


def _build_modular_annotation_label(
    lineage: pd.Series,
    subtype: pd.Series,
    state: pd.Series,
) -> pd.Series:
    """Construct a compact modular display label from lineage/subtype/state columns."""
    result = []
    for lin, sub, st in zip(lineage.astype(str), subtype.astype(str), state.astype(str)):
        label = lin
        if sub not in {"Unknown", "Not_applicable", "nan", ""} and sub != lin:
            label = sub
        if st not in {"Unknown", "Not_applicable", "nan", ""}:
            label = f"{label} | {st}"
        result.append(label)
    return pd.Series(result, index=lineage.index, dtype="object")


def _rename_score_columns_for_manager(
    adata: AnnData,
    manager: Manager,
    suffix: str,
) -> None:
    """Rename manager-derived *_score columns to a scoped suffix."""
    rename_map = {}
    for cell_type in manager.CELLS:
        source = f"{cell_type}_score"
        if source in adata.obs.columns:
            rename_map[source] = f"{cell_type}{suffix}"
    if rename_map:
        adata.obs.rename(columns=rename_map, inplace=True)


def _collect_state_signatures(
    config: AnnotationConfig,
) -> tuple[Dict[str, List[str]], Dict[str, Dict[str, object]]]:
    """Collect state/program signatures plus optional scope metadata."""
    signatures: Dict[str, List[str]] = {}
    metadata: Dict[str, Dict[str, object]] = {}

    if config.marker_states:
        state_mgr = _resolve_annotation_manager(
            species=config.marker_species,
            tissue=None,
            marker_config=config.state_marker_config or f"cell_state_{config.marker_species}",
        )
        selected = state_mgr.select_cells(config.marker_states, include_children=True)
        for name, cell in selected.CELLS.items():
            if cell.markers:
                signatures[name] = list(cell.markers)
                metadata[name] = dict(cell.metadata)

    if config.custom_state_signatures:
        signatures.update(config.custom_state_signatures)
        if config.custom_state_metadata:
            metadata.update(config.custom_state_metadata)

    if config.state_signature_names or config.state_signature_categories:
        manager = FunctionalSignatureManager(species=config.marker_species)
        for category in config.state_signature_categories:
            signatures.update(manager.get_category(category))
        for name in config.state_signature_names:
            signatures[name] = manager.get_signature(name)

    return signatures, metadata


def _state_applies_to_cell(
    state_meta: Optional[Dict[str, object]],
    lineage_label: str,
    subtype_label: str,
) -> bool:
    """Check whether a state program is valid for the current lineage/subtype context."""
    if not state_meta:
        return True

    scope = str(state_meta.get("scope", "all")).lower()
    applies_to = state_meta.get("applies_to", ["all"])
    if isinstance(applies_to, str):
        applies = [applies_to]
    else:
        applies = [str(x) for x in applies_to]

    applies_lower = [x.lower() for x in applies]
    if "all" in applies_lower or scope == "all":
        return True

    lineage_lower = str(lineage_label).lower()
    subtype_lower = str(subtype_label).lower()

    for allowed in applies_lower:
        if allowed and (allowed == lineage_lower or allowed in lineage_lower):
            return True
        if allowed and subtype_lower not in {"not_applicable", "unknown", "nan"}:
            if allowed == subtype_lower or allowed in subtype_lower:
                return True
    return False


# --------------------- Scoring Function -----------------------


def score_cell_types(
    adata: AnnData,
    marker_config: Union[str, Manager],
    layer: Optional[str] = "normalized",
    use_raw: bool = True,
    min_genes: int = 3,
    ctrl_size: int = 50,
    score_name_suffix: str = "_score",
    copy: bool = False,
) -> AnnData:
    """
    Score cells for cell type marker gene sets.
    Adds score columns to adata.obs.
    """
    if copy:
        adata = adata.copy()
    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError("marker_config must be a file path or Manager instance.")
    mgr.intersect_with(adata.raw if use_raw and adata.raw is not None else adata)
    n_scored, n_skipped = 0, 0
    if use_raw:
        for cell_type, cell in mgr.CELLS.items():
            if len(cell.markers) >= min_genes:
                sc.tl.score_genes(
                    adata,
                    cell.markers,
                    score_name=f"{cell_type}{score_name_suffix}",
                    use_raw=True,
                    ctrl_size=ctrl_size,
                )
                n_scored += 1
            else:
                n_skipped += 1
    else:
        with use_layer_as_X(adata, layer):
            for cell_type, cell in mgr.CELLS.items():
                if len(cell.markers) >= min_genes:
                    sc.tl.score_genes(
                        adata,
                        cell.markers,
                        score_name=f"{cell_type}{score_name_suffix}",
                        ctrl_size=ctrl_size,
                    )
                    n_scored += 1
                else:
                    n_skipped += 1
    log.info(f"Scored {n_scored} cell types ({n_skipped} skipped).")
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    scoring_params = sanitize_for_hdf5(
        {
            "use_raw": use_raw,
            "layer": layer,
            "min_genes": min_genes,
            "ctrl_size": ctrl_size,
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"]["scoring_params"] = scoring_params
    return adata


# --------------------- Core annotation function --------------------------


def annotate_clusters(
    adata: AnnData,
    cluster_key: str,
    marker_config: Union[str, Manager],
    method: Literal["max_score", "enrichment", "combined"] = "max_score",
    use_raw: bool = False,
    key_added: Optional[str] = None,
    min_confidence: float = 0.3,
    confidence_key: Optional[str] = None,
    min_score: float = 0.1,
    n_genes: int = 100,
    score_weight: float = 0.6,
    enrichment_weight: float = 0.4,
    plot: bool = False,
    copy: bool = False,
) -> AnnData:
    """
    Assign cell type labels to clusters using various evidence.

    Enhancements:
    - Robust to missing score columns and non-categorical cluster keys.
    - Stable Unknown handling and category order preservation.
    - Parameter trace including scanpy version and marker stats.
    """
    if copy:
        adata = adata.copy()
    if key_added is None:
        key_added = f"{cluster_key}_annotated"

    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError("marker_config must be a file path or Manager instance.")
    mgr.intersect_with(adata.raw if use_raw and adata.raw is not None else adata)

    # Ensure categorical
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")
    if not pd.api.types.is_categorical_dtype(adata.obs[cluster_key]):
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

    # 1. Score-based
    def annotate_by_max_score():
        score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
        if not score_cols:
            raise RuntimeError("No *_score columns found. Please run score_cell_types first.")
        means = adata.obs.groupby(cluster_key)[score_cols].mean()
        result = {}
        for cluster in means.index:
            best = means.loc[cluster].idxmax()
            best_score = float(means.loc[cluster, best])
            cell_type = best[:-6] if best.endswith("_score") else best
            result[str(cluster)] = cell_type if best_score >= min_score else "Unknown"
        return result

    # 2. Enrichment-based
    def annotate_by_enrichment():
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw,
            key_added=f"rank_genes_{cluster_key}",
        )
        markers_df = sc.get.rank_genes_groups_df(adata, key=f"rank_genes_{cluster_key}")
        result = {}
        categories = list(adata.obs[cluster_key].cat.categories)
        for cluster in categories:
            genes = (
                markers_df.loc[markers_df["group"] == cluster, "names"]
                .head(n_genes)
                .astype(str)
                .tolist()
            )
            best_score, best_type = -1.0, "Unknown"
            for cell_type, cell in mgr.CELLS.items():
                if not cell.markers:
                    continue
                denom = max(1, len(cell.markers))
                overlap = len(set(genes) & set(cell.markers)) / denom
                if overlap > best_score:
                    best_score, best_type = overlap, cell_type
            result[str(cluster)] = best_type if best_score >= min_score else "Unknown"
        return result

    # 3. Combined
    def annotate_by_combined():
        score_cols = [col for col in adata.obs.columns if col.endswith("_score")]
        if not score_cols:
            raise RuntimeError("No *_score columns found. Please run score_cell_types first.")
        means = adata.obs.groupby(cluster_key)[score_cols].mean()
        sc.tl.rank_genes_groups(
            adata,
            groupby=cluster_key,
            method="wilcoxon",
            use_raw=use_raw,
            key_added=f"rank_genes_{cluster_key}",
        )
        markers_df = sc.get.rank_genes_groups_df(adata, key=f"rank_genes_{cluster_key}")
        categories = list(adata.obs[cluster_key].cat.categories)
        result = {}
        for cluster in categories:
            genes = (
                markers_df.loc[markers_df["group"] == cluster, "names"]
                .head(n_genes)
                .astype(str)
                .tolist()
            )
            combined_scores = {}
            for cell_type, cell in mgr.CELLS.items():
                score_col = f"{cell_type}_score"
                score_val = (
                    float(means.loc[cluster, score_col]) if score_col in means.columns else 0.0
                )
                denom = max(1, len(cell.markers))
                overlap_val = (len(set(genes) & set(cell.markers)) / denom) if cell.markers else 0.0
                combined_scores[cell_type] = (
                    score_weight * score_val + enrichment_weight * overlap_val
                )
            best_type = max(combined_scores, key=combined_scores.get)
            best_score = combined_scores[best_type]
            result[str(cluster)] = best_type if best_score >= min_score else "Unknown"
        return result

    # Select method
    if method == "max_score":
        mapping = annotate_by_max_score()
    elif method == "enrichment":
        mapping = annotate_by_enrichment()
    elif method == "combined":
        mapping = annotate_by_combined()
    else:
        raise ValueError(f"Unknown annotation method: {method}")

    # Assign labels
    cluster_codes = adata.obs[cluster_key].astype(str)
    assigned = cluster_codes.map(mapping)
    assigned = assigned.fillna("Unknown")
    adata.obs[key_added] = pd.Categorical(assigned)

    # Save and optional plot
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    params_dict = sanitize_for_hdf5(
        {
            "method": method,
            "min_score": min_score,
            "score_weight": score_weight,
            "enrichment_weight": enrichment_weight,
            "mapping": mapping,
            "scanpy_version": getattr(sc, "__version__", "unknown"),
            "n_markers": {k: len(v.markers) for k, v in getattr(mgr, "CELLS", {}).items()},
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_params"] = params_dict
    if plot:
        if "X_umap" not in adata.obsm:
            sc.tl.umap(adata)
        sc.pl.umap(adata, color=[cluster_key, key_added], wspace=0.4)
    return adata


def _get_celltypist_label_series(
    adata: AnnData,
    key_prefix: str = "celltypist",
) -> tuple[pd.Series, pd.Series, str]:
    """Return the preferred CellTypist label series and confidence scores."""
    majority_key = f"{key_prefix}_majority_voting"
    predicted_key = f"{key_prefix}_predicted_labels"
    conf_key = f"{key_prefix}_conf_score"

    if majority_key in adata.obs.columns:
        label_key = majority_key
    elif predicted_key in adata.obs.columns:
        label_key = predicted_key
    else:
        raise KeyError(
            f"CellTypist results not found for prefix '{key_prefix}'. "
            f"Expected '{majority_key}' or '{predicted_key}'."
        )

    labels = adata.obs[label_key].astype(str)
    confidences = (
        pd.to_numeric(adata.obs[conf_key], errors="coerce")
        if conf_key in adata.obs.columns
        else pd.Series(1.0, index=adata.obs.index, dtype=float)
    )
    return labels, confidences, label_key


def _build_celltypist_cluster_mapping(
    adata: AnnData,
    cluster_key: str,
    *,
    key_prefix: str = "celltypist",
    min_confidence: float = 0.5,
) -> tuple[Dict[str, str], Dict[str, Dict[str, float]]]:
    """Aggregate CellTypist predictions to cluster-level labels."""
    labels, confidences, label_key = _get_celltypist_label_series(adata, key_prefix=key_prefix)
    mapping: Dict[str, str] = {}
    stats: Dict[str, Dict[str, float]] = {}

    cluster_labels = adata.obs[cluster_key].astype(str)
    for cluster in sorted(cluster_labels.unique(), key=str):
        mask = cluster_labels == str(cluster)
        cluster_pred = labels.loc[mask]
        cluster_conf = confidences.loc[mask]
        cluster_pred = cluster_pred[cluster_pred.notna() & (cluster_pred != "nan")]

        if cluster_pred.empty:
            mapping[str(cluster)] = "Unknown"
            stats[str(cluster)] = {
                "mean_conf_score": 0.0,
                "majority_fraction": 0.0,
                "label_source": label_key,
            }
            continue

        counts = Counter(cluster_pred.tolist())
        best_label, best_count = counts.most_common(1)[0]
        majority_fraction = best_count / max(1, len(cluster_pred))
        mean_conf = float(cluster_conf.mean()) if len(cluster_conf) else majority_fraction

        mapping[str(cluster)] = best_label if mean_conf >= min_confidence else "Unknown"
        stats[str(cluster)] = {
            "mean_conf_score": mean_conf,
            "majority_fraction": float(majority_fraction),
            "label_source": label_key,
        }

    return mapping, stats


def _combine_marker_and_celltypist_mappings(
    marker_mapping: Dict[str, str],
    celltypist_mapping: Dict[str, str],
    celltypist_stats: Dict[str, Dict[str, float]],
    *,
    min_celltypist_confidence: float = 0.7,
) -> tuple[Dict[str, str], Dict[str, Dict[str, Union[str, float]]]]:
    """Resolve final labels from marker-based and CellTypist evidence."""
    final_mapping: Dict[str, str] = {}
    audit: Dict[str, Dict[str, Union[str, float]]] = {}

    all_clusters = sorted(set(marker_mapping) | set(celltypist_mapping), key=str)
    for cluster in all_clusters:
        marker_label = marker_mapping.get(cluster, "Unknown")
        celltypist_label = celltypist_mapping.get(cluster, "Unknown")
        stats = celltypist_stats.get(cluster, {})
        mean_conf = float(stats.get("mean_conf_score", 0.0))
        majority_fraction = float(stats.get("majority_fraction", 0.0))

        if (
            marker_label != "Unknown"
            and celltypist_label != "Unknown"
            and marker_label == celltypist_label
        ):
            final_label = marker_label
            decision = "agreement"
        elif celltypist_label != "Unknown" and mean_conf >= min_celltypist_confidence:
            final_label = celltypist_label
            decision = "celltypist_high_confidence"
        elif marker_label != "Unknown":
            final_label = marker_label
            decision = "marker_fallback"
        elif celltypist_label != "Unknown" and majority_fraction >= 0.6:
            final_label = celltypist_label
            decision = "celltypist_majority_fallback"
        else:
            final_label = "Unknown"
            decision = "insufficient_evidence"

        final_mapping[cluster] = final_label
        audit[cluster] = {
            "marker_label": marker_label,
            "celltypist_label": celltypist_label,
            "final_label": final_label,
            "decision": decision,
            "mean_conf_score": mean_conf,
            "majority_fraction": majority_fraction,
        }

    return final_mapping, audit


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


# --------------------- CellTypist/Transfer --------------------


def run_celltypist(
    adata: AnnData,
    model: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
    use_raw: bool = True,
    key_added: str = "celltypist",
) -> AnnData:
    """
    Run CellTypist for automated annotation. Stores results in adata.obs.
    """
    import celltypist

    input_adata = adata.raw.to_adata() if use_raw and adata.raw is not None else adata
    pred = celltypist.annotate(input_adata, model=model, majority_voting=majority_voting)
    adata.obs[f"{key_added}_predicted_labels"] = pred.predicted_labels["predicted_labels"].reindex(
        adata.obs.index
    )
    adata.obs[f"{key_added}_conf_score"] = pred.predicted_labels["conf_score"].reindex(
        adata.obs.index
    )
    if majority_voting:
        adata.obs[f"{key_added}_majority_voting"] = pred.predicted_labels[
            "majority_voting"
        ].reindex(adata.obs.index)
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    celltypist_info = {
        "model": model,
        "majority_voting": majority_voting,
        "timestamp": str(pd.Timestamp.now()),
    }
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_results"] = sanitize_for_hdf5(
        celltypist_info
    )
    return adata


def transfer_labels(
    adata: AnnData,
    ref_adata: AnnData,
    ref_label_key: str,
    n_neighbors: int = 30,
    use_rep: str = "X_pca",
    key_added: Optional[str] = "predicted_labels",
    normalize_weights: bool = True,
    confidence_threshold: float = 0.7,
    copy: bool = False,
) -> AnnData:
    """
    Transfer cell type labels from a reference dataset (kNN label transfer).
    """
    from sklearn.neighbors import NearestNeighbors

    if copy:
        adata = adata.copy()
    query, ref = adata.obsm[use_rep], ref_adata.obsm[use_rep]
    ref_labels = ref_adata.obs[ref_label_key].values
    nn = NearestNeighbors(n_neighbors=n_neighbors).fit(ref)
    distances, indices = nn.kneighbors(query)
    result, confidences = [], []
    for i in range(len(adata)):
        neighbor_labels = ref_labels[indices[i]]
        if normalize_weights:
            weights = 1.0 / (distances[i] + 1e-10)
            label_votes = {}
            for l in np.unique(neighbor_labels):
                label_votes[l] = np.sum(weights[neighbor_labels == l])
            best_label, confidence = (
                max(label_votes, key=label_votes.get),
                max(label_votes.values()) / sum(label_votes.values()),
            )
        else:
            from collections import Counter

            c = Counter(neighbor_labels)
            best_label, confidence = c.most_common(1)[0][0], c[best_label] / n_neighbors
        result.append(best_label if confidence >= confidence_threshold else "Unknown")
        confidences.append(confidence)
    adata.obs[key_added] = pd.Categorical(result)
    adata.obs[f"{key_added}_confidence"] = confidences
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    transfer_params = sanitize_for_hdf5(
        {
            "method": "knn_transfer",
            "n_neighbors": n_neighbors,
            "use_rep": use_rep,
            "normalize_weights": normalize_weights,
            "confidence_threshold": confidence_threshold,
            "reference_label_key": ref_label_key,
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_params"] = transfer_params
    return adata


# --------------------- AI/Manual annotation -------------------------


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


def apply_annotation_mapping(
    adata: AnnData,
    cluster_key: str,
    mapping: Union[Dict[str, str], str],
    key_added: str = "cell_type",
) -> AnnData:
    """
    Apply AI/manual cluster-to-celltype mapping from dict, Excel (.xlsx/.xls), CSV, or JSON.
    Robust to type mismatches and common header conventions.

    Supported table schemas:
    - Columns named ['cluster','cell_type'] (preferred)
    - Or the first two columns are treated as [cluster, cell_type]
    """
    import datetime
    import os

    # 0) Validate cluster_key
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")

    # 1) Load mapping
    if isinstance(mapping, str):
        ext = os.path.splitext(mapping)[1].lower()
        if ext in [".xlsx", ".xls", ".csv"]:
            df = _read_table_file(mapping)
            if {"cluster", "cell_type"}.issubset(df.columns):
                mapping_dict = pd.Series(df["cell_type"].values, index=df["cluster"]).to_dict()
            elif len(df.columns) >= 2:
                mapping_dict = pd.Series(df.iloc[:, 1].values, index=df.iloc[:, 0]).to_dict()
                log.info(
                    "Mapping file has no standard headers; used first two columns as [cluster, cell_type]."
                )
            else:
                raise ValueError(
                    "Mapping table must have either ['cluster','cell_type'] columns or at least 2 columns."
                )
            mapping_source = {"type": "file", "path": mapping}
        elif ext == ".json":
            mapping_dict = _read_json_file(mapping)
            mapping_source = {"type": "file", "path": mapping}
        else:
            raise ValueError("Unsupported mapping file format. Use .xlsx, .xls, .csv or .json.")
    elif isinstance(mapping, dict):
        mapping_dict = mapping
        mapping_source = {"type": "dict"}
    else:
        raise TypeError("mapping must be a dictionary or a file path string.")

    # 2) Ensure robust type matching (convert keys to string)
    source_clusters = adata.obs[cluster_key].astype(str)
    mapping_dict_str_keys = {str(k): v for k, v in mapping_dict.items()}

    # 3) Apply the mapping
    new_series = source_clusters.map(mapping_dict_str_keys)
    adata.obs[key_added] = pd.Categorical(new_series)

    # 4) Unmapped handling
    unmapped_mask = adata.obs[key_added].isnull()
    if unmapped_mask.any():
        unmapped_ids = source_clusters[unmapped_mask].unique().tolist()
        log.warning(
            f"{unmapped_mask.sum()} cells could not be mapped. "
            f"Missing cluster IDs in mapping: {unmapped_ids}"
        )
        adata.obs[key_added] = (
            adata.obs[key_added].cat.add_categories("Unmapped").fillna("Unmapped")
        )

    # 5) Store metadata snapshot
    annot_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annot_ns[f"{key_added}_mapping"] = sanitize_for_hdf5(mapping_dict)
    annot_ns[f"{key_added}_mapping_meta"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "created_at": datetime.datetime.now().isoformat(timespec="seconds"),
            **mapping_source,
        }
    )

    log.info(
        f"Applied mapping to column '{key_added}'. Categories: {list(adata.obs[key_added].cat.categories)}"
    )
    return adata


def remap_labels(
    adata: AnnData,
    column: str,
    mapping: Optional[Dict[str, str]] = None,
    where: Optional[pd.Series] = None,
    to: Optional[str] = None,
    in_place: bool = True,
    key_added: Optional[str] = None,
    tidy_categories: bool = True,
) -> AnnData:
    """
    Partially remap or correct labels in an existing categorical column of adata.obs.

    Modes:
    1) Dictionary-based rename:
       - mapping={'OldA':'NewA', 'OldB':'NewB'}
    2) Condition-based assign:
       - where: boolean Series aligned to adata.obs (True = replace), to='NewLabel'

    Args:
        column: Existing obs column to modify (e.g., 'celltype').
        mapping: Old -> New label mapping.
        where: Boolean mask of cells to set to `to`.
        to: Target label for cells where mask is True.
        in_place: If False, write to `key_added` and keep original intact.
        key_added: Required if in_place=False.
        tidy_categories: Drop unused categories after remap (clean legend).
    """
    if column not in adata.obs.columns:
        raise KeyError(f"Column '{column}' not found in adata.obs.")
    if not in_place and not key_added:
        raise ValueError("When in_place=False, you must provide key_added for the new column name.")

    src = adata.obs[column].astype(str)
    s = src.copy()

    changed_by_mapping = 0
    changed_by_where = 0

    if mapping is not None:
        before = s.copy()
        s = s.replace(mapping)
        changed_by_mapping = int((s != before).sum())

    if where is not None:
        if to is None:
            raise ValueError("When providing 'where', you must also provide 'to' label.")
        if not isinstance(where, pd.Series) or not where.index.equals(adata.obs.index):
            raise ValueError("'where' must be a boolean Series aligned to adata.obs index.")
        before = s.copy()
        s.loc[where] = to
        changed_by_where = int((s != before).sum())

    # Cast back to categorical
    cat = pd.Categorical(s)
    target_col = column if in_place else key_added
    adata.obs[target_col] = cat

    # Optionally tidy categories (drop unused)
    if tidy_categories:
        adata.obs[target_col] = adata.obs[target_col].cat.remove_unused_categories()

    # Audit
    annot_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    audit = annot_ns.setdefault("remap_audit", [])
    audit_entry = sanitize_for_hdf5(
        {
            "source_column": column,
            "target_column": target_col,
            "mode": (
                "mapping+where"
                if (mapping is not None and where is not None)
                else ("mapping" if mapping is not None else "where")
            ),
            "mapping": mapping if mapping is not None else None,
            "to": to if where is not None else None,
            "changed_by_mapping": changed_by_mapping,
            "changed_by_where": changed_by_where,
            "final_categories": list(adata.obs[target_col].cat.categories),
        }
    )
    audit.append(audit_entry)

    log.info(
        f"Remapped labels in '{target_col}': +{changed_by_mapping} (mapping), +{changed_by_where} (where). Categories: {list(adata.obs[target_col].cat.categories)}"
    )
    return adata


# --------------------- Evaluation & Main Function --------------------------


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
            "scanpy_version": getattr(sc, "__version__", "unknown"),
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
