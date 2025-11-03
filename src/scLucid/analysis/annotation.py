"""
Cell type annotation functions for single-cell RNA-seq data.
Supports: score-based, enrichment-based, combined, CellTypist, reference transfer,
manual/AI mapping import, and annotation evaluation.
"""

import logging
from typing import Dict, Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from ..utils import sanitize_for_hdf5, use_layer_as_X
from ..markers import Manager, get_marker_manager
from .config import AnnotationConfig

log = logging.getLogger(__name__)


__all__ = [
    "score_cell_types",
    "annotate_clusters",
    "run_celltypist",
    "transfer_labels",
    "evaluate_annotation",
    "summarize_annotation_evidence",
    "apply_annotation_mapping",
    "remap_labels",
    "run_annotation",
]


# --------------------- Helper Functions -----------------------
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
            with open(path, "r", encoding=enc) as f:
                return json.load(f)
        except Exception as e:
            last_err = e
    try:
        import chardet

        with open(path, "rb") as fb:
            raw = fb.read(200000)
        guess = chardet.detect(raw).get("encoding") or "utf-8"
        with open(path, "r", encoding=guess, errors="replace") as f:
            return json.load(f)
    except Exception:
        pass
    raise (
        last_err
        if last_err
        else RuntimeError("Failed to read JSON with multiple encodings.")
    )


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
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
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
            raise RuntimeError(
                "No *_score columns found. Please run score_cell_types first."
            )
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
            raise RuntimeError(
                "No *_score columns found. Please run score_cell_types first."
            )
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
                    float(means.loc[cluster, score_col])
                    if score_col in means.columns
                    else 0.0
                )
                denom = max(1, len(cell.markers))
                overlap_val = (
                    (len(set(genes) & set(cell.markers)) / denom)
                    if cell.markers
                    else 0.0
                )
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
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
    params_dict = sanitize_for_hdf5(
        {
            "method": method,
            "min_score": min_score,
            "score_weight": score_weight,
            "enrichment_weight": enrichment_weight,
            "mapping": mapping,
            "scanpy_version": getattr(sc, "__version__", "unknown"),
            "n_markers": {
                k: len(v.markers) for k, v in getattr(mgr, "CELLS", {}).items()
            },
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_params"] = params_dict
    if plot:
        if "X_umap" not in adata.obsm:
            sc.tl.umap(adata)
        sc.pl.umap(adata, color=[cluster_key, key_added], wspace=0.4)
    return adata


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
    pred = celltypist.annotate(
        input_adata, model=model, majority_voting=majority_voting
    )
    adata.obs[f"{key_added}_predicted_labels"] = pred.predicted_labels[
        "predicted_labels"
    ].reindex(adata.obs.index)
    adata.obs[f"{key_added}_conf_score"] = pred.predicted_labels["conf_score"].reindex(
        adata.obs.index
    )
    if majority_voting:
        adata.obs[f"{key_added}_majority_voting"] = pred.predicted_labels[
            "majority_voting"
        ].reindex(adata.obs.index)
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
    celltypist_info = {
        "model": model,
        "majority_voting": majority_voting,
        "timestamp": str(pd.Timestamp.now()),
    }
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_results"] = (
        sanitize_for_hdf5(celltypist_info)
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
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
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
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_params"] = (
        transfer_params
    )
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
                mapping_dict = pd.Series(
                    df["cell_type"].values, index=df["cluster"]
                ).to_dict()
            elif len(df.columns) >= 2:
                mapping_dict = pd.Series(
                    df.iloc[:, 1].values, index=df.iloc[:, 0]
                ).to_dict()
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
            raise ValueError(
                "Unsupported mapping file format. Use .xlsx, .xls, .csv or .json."
            )
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
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("annotation", {})
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
        raise ValueError(
            "When in_place=False, you must provide key_added for the new column name."
        )

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
            raise ValueError(
                "When providing 'where', you must also provide 'to' label."
            )
        if not isinstance(where, pd.Series) or not where.index.equals(adata.obs.index):
            raise ValueError(
                "'where' must be a boolean Series aligned to adata.obs index."
            )
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
        adata.uns.setdefault("sclucid", {})
        .setdefault("analysis", {})
        .setdefault("annotation", {})
    )
    audit = annot_ns.setdefault("remap_audit", [])
    audit_entry = sanitize_for_hdf5(
        {
            "source_column": column,
            "target_column": target_col,
            "mode": "mapping+where"
            if (mapping is not None and where is not None)
            else ("mapping" if mapping is not None else "where"),
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
        log.warning(
            "Marker manager has no cell types. Evaluation may be uninformative."
        )

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
        if (
            de_genes.empty
            or "pvals_adj" not in de_genes.columns
            or "names" not in de_genes.columns
        ):
            log.warning(
                f"No DE genes or required columns missing for cluster '{cluster}'."
            )
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
        marker_coverage = (
            len(found_markers) / len(expected_markers) if expected_markers else 0.0
        )

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
                "found_markers": ", ".join(sorted(found_markers))
                if found_markers
                else "",
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
                results_df["cluster"].astype(str)
                + " ("
                + results_df["cell_type"]
                + ")",
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

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
    adata.uns["sclucid"]["analysis"]["annotation"][f"{annotation_key}_evaluation"] = (
        results_df
    )
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
    mgr = get_marker_manager(species=config.marker_species, tissue=config.marker_tissue)
    mgr.intersect_with(adata.raw if adata.raw is not None else adata)

    if config.run_celltypist:
        adata = run_celltypist(adata, model=config.celltypist_model)

    if config.run_scoring:
        # Use raw for scoring by default
        adata = score_cell_types(adata, marker_config=mgr, use_raw=True)

    adata = annotate_clusters(
        adata,
        cluster_key=config.cluster_key,
        marker_config=mgr,
        method=config.final_method,
        key_added=config.key_added,
        min_score=config.min_confidence,
        use_raw=True,
        plot=config.plot,
    )

    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
    adata.uns["sclucid"]["analysis"]["annotation"]["workflow_config"] = (
        sanitize_for_hdf5(config.to_dict())
    )
    return adata
