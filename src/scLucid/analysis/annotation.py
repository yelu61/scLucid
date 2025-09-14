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

from ..utils import use_layer_as_X
from ..utils.marker_manager import Manager, get_marker_manager
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
    "run_annotation",
]

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
    adata.uns["sclucid"]["analysis"]["annotation"]["scoring_params"] = {
        "use_raw": use_raw,
        "layer": layer,
        "min_genes": min_genes,
        "ctrl_size": ctrl_size,
    }
    return adata


# --------------------- Core annotation function --------------------------


def annotate_clusters(
    adata: AnnData,
    cluster_key: str,
    marker_config: Union[str, Manager],
    method: Literal["max_score", "enrichment", "combined"] = "max_score",
    use_raw: bool = False,
    key_added: Optional[str] = None,
    min_score: float = 0.1,
    n_genes: int = 100,
    score_weight: float = 0.6,
    enrichment_weight: float = 0.4,
    plot: bool = False,
    copy: bool = False,
) -> AnnData:
    """
    Assign cell type labels to clusters using various evidence.
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
            best_score = means.loc[cluster, best]
            cell_type = best.replace("_score", "")
            result[cluster] = cell_type if best_score >= min_score else "Unknown"
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
        for cluster in adata.obs[cluster_key].cat.categories:
            genes = (
                markers_df[markers_df["group"] == cluster]["names"]
                .head(n_genes)
                .tolist()
            )
            best_score, best_type = -1, "Unknown"
            for cell_type, cell in mgr.CELLS.items():
                if not cell.markers:
                    continue
                overlap = len(set(genes) & set(cell.markers)) / len(cell.markers)
                if overlap > best_score:
                    best_score, best_type = overlap, cell_type
            result[cluster] = best_type if best_score >= min_score else "Unknown"
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
        result = {}
        for cluster in adata.obs[cluster_key].cat.categories:
            genes = (
                markers_df[markers_df["group"] == cluster]["names"]
                .head(n_genes)
                .tolist()
            )
            combined_scores = {}
            for cell_type, cell in mgr.CELLS.items():
                score_val = (
                    means.loc[cluster, f"{cell_type}_score"]
                    if f"{cell_type}_score" in means.columns
                    else 0
                )
                overlap_val = (
                    len(set(genes) & set(cell.markers)) / len(cell.markers)
                    if cell.markers
                    else 0
                )
                combined_scores[cell_type] = (
                    score_weight * score_val + enrichment_weight * overlap_val
                )
            best_type = max(combined_scores, key=combined_scores.get)
            best_score = combined_scores[best_type]
            result[cluster] = best_type if best_score >= min_score else "Unknown"
        return result

    # Select annotation method
    if method == "max_score":
        mapping = annotate_by_max_score()
    elif method == "enrichment":
        mapping = annotate_by_enrichment()
    elif method == "combined":
        mapping = annotate_by_combined()
    else:
        raise ValueError(f"Unknown annotation method: {method}")
    adata.obs[key_added] = adata.obs[cluster_key].map(mapping).astype("category")
    # Save and plot
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_params"] = {
        "method": method,
        "min_score": min_score,
        "score_weight": score_weight,
        "enrichment_weight": enrichment_weight,
        "mapping": mapping,
    }
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
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_results"] = pred
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
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_params"] = {
        "method": "knn_transfer",
        "n_neighbors": n_neighbors,
        "use_rep": use_rep,
        "normalize_weights": normalize_weights,
        "confidence_threshold": confidence_threshold,
        "reference_label_key": ref_label_key,
    }
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
    Apply AI/manual cluster-to-celltype mapping from dict or csv/json file.
    This robust version handles data type mismatches between clusters and mapping keys.
    """
    # 1. --- Load mapping from file path if necessary ---
    if isinstance(mapping, str):
        if mapping.endswith(".csv"):
            df = pd.read_csv(mapping)
            if len(df.columns) < 2:
                raise ValueError("Mapping CSV must have at least two columns (cluster, cell_type)")
            # Use the first column for keys and the second for values
            mapping_dict = pd.Series(df.iloc[:, 1].values, index=df.iloc[:, 0]).to_dict()
        elif mapping.endswith(".json"):
            import json
            with open(mapping) as f:
                mapping_dict = json.load(f)
        else:
            raise ValueError("Unsupported mapping file format. Use .csv or .json")
    elif isinstance(mapping, dict):
        mapping_dict = mapping
    else:
        raise TypeError("mapping must be a dictionary or a file path string.")

    # 2. --- [CRITICAL FIX] Ensure robust type matching by converting to string ---
    # Convert cluster IDs in adata to string
    source_clusters = adata.obs[cluster_key].astype(str)
    # Convert keys in the mapping dictionary to string
    mapping_dict_str_keys = {str(k): v for k, v in mapping_dict.items()}

    # 3. --- Apply the mapping ---
    adata.obs[key_added] = source_clusters.map(mapping_dict_str_keys).astype("category")

    # 4. --- Check for unmapped clusters and provide a helpful warning ---
    unmapped_mask = adata.obs[key_added].isnull()
    if unmapped_mask.any():
        unmapped_ids = source_clusters[unmapped_mask].unique().tolist()
        log.warning(
            f"{unmapped_mask.sum()} cells could not be mapped. "
            f"This is likely because the following cluster IDs were not found in your mapping file: {unmapped_ids}"
        )
        # Optionally, fill with a default value like 'Unmapped'
        adata.obs[key_added] = adata.obs[key_added].cat.add_categories("Unmapped").fillna("Unmapped")


    # 5. --- Store results and return ---
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    adata.uns["sclucid"]["analysis"]["annotation"][f"{key_added}_mapping"] = mapping_dict
    
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
    """
    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError("marker_config must be a file path or Manager instance.")
    mgr.intersect_with(adata)
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
        key_added=f"rank_genes_{cluster_key}",
    )
    results = []
    for cluster in adata.obs[cluster_key].cat.categories:
        de_genes = sc.get.rank_genes_groups_df(
            adata, key=f"rank_genes_{cluster_key}", group=cluster
        )
        sig_genes = set(de_genes[de_genes["pvals_adj"] < 0.05]["names"])
        cell_types = adata.obs.loc[adata.obs[cluster_key] == cluster, annotation_key]
        assigned_type = cell_types.value_counts().index[0]
        if assigned_type == "Unknown" or assigned_type not in mgr.CELLS:
            continue
        expected_markers = set(mgr[assigned_type].markers)
        found_markers = sig_genes & expected_markers
        marker_coverage = (
            len(found_markers) / len(expected_markers) if expected_markers else 0
        )
        all_other_markers = set(
            m for t, c in mgr.CELLS.items() if t != assigned_type for m in c.markers
        )
        specificity = (
            1.0 - (len(found_markers & all_other_markers) / len(found_markers))
            if found_markers
            else 0
        )
        confidence = 0.6 * marker_coverage + 0.4 * specificity
        results.append(
            {
                "cluster": cluster,
                "cell_type": assigned_type,
                "marker_coverage": marker_coverage,
                "marker_specificity": specificity,
                "annotation_confidence": confidence,
                "found_markers": ", ".join(found_markers),
                "expected_markers": len(expected_markers),
                "detected_markers": len(found_markers),
            }
        )
    results_df = pd.DataFrame(results)
    if plot and not results_df.empty:
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
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
    adata.uns["sclucid"]["analysis"]["annotation"][f"{annotation_key}_evaluation"] = (
        results_df
    )
    return results_df


def run_annotation(
    adata: AnnData,
    config: AnnotationConfig,
) -> AnnData:
    """
    Full annotation workflow: scoring → auto annotation → results in .obs/.uns.
    """
    
    mgr = get_marker_manager(species=config.marker_species, tissue=config.marker_tissue)
    mgr.intersect_with(adata.raw if adata.raw is not None else adata)
    if config.run_celltypist:
        adata = run_celltypist(adata, model=config.celltypist_model)
    if config.run_scoring:
        adata = score_cell_types(adata, marker_config=mgr, use_raw=True)
    adata = annotate_clusters(
        adata,
        cluster_key=config.cluster_key,
        marker_config=mgr,
        method=config.final_method,
        key_added=config.key_added,
        min_score=config.min_confidence,
        use_raw=True,
    )
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
        "annotation", {}
    )
    adata.uns["sclucid"]["analysis"]["annotation"]["workflow_config"] = config.to_dict()
    return adata
