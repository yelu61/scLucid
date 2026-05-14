"""Reference-based annotation helpers."""

from __future__ import annotations

from collections import Counter
from typing import Dict, Union

import logging
import pandas as pd
from anndata import AnnData

from ...utils import sanitize_for_hdf5

log = logging.getLogger(__name__)


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
