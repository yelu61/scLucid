#!/usr/bin/env python3
"""Pilot PCA / graph handoff stability benchmark on real datasets.

This benchmark now uses Leiden clustering (via scanpy) for seed-stability
assessment, matching the production preprocessing workflow. The Leiden results
are reported as ``leiden_seed`` stability to distinguish them from the earlier
KMeans-based proxy.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.dataset_registry import DATASETS


def _matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _subset(adata: ad.AnnData, max_cells: int | None, seed: int) -> ad.AnnData:
    if max_cells is None or max_cells <= 0 or adata.n_obs <= max_cells:
        return adata
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(adata.n_obs, size=max_cells, replace=False))
    return adata[keep].copy()


def _log_normalized_matrix(adata: ad.AnnData, target_sum: float = 1e4):
    X = _matrix(adata).astype(float)
    totals = np.asarray(X.sum(axis=1)).ravel()
    scale = np.divide(target_sum, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
    if sp.issparse(X):
        Xn = X.multiply(scale[:, None]).tocsr(copy=False)
        Xn.data = np.log1p(Xn.data)
        return Xn
    return np.log1p(np.asarray(X) * scale[:, None])


def _top_variable_subset(X, var_names, n_genes: int = 2000):
    if X.shape[1] <= n_genes:
        return X, list(map(str, var_names))
    if sp.issparse(X):
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
        variance = np.maximum(mean_sq - mean * mean, 0.0)
    else:
        variance = np.asarray(X).var(axis=0)
    idx = np.argsort(-variance, kind="mergesort")[:n_genes]
    return X[:, idx], [str(var_names[i]) for i in idx]


def _pca_embedding(X, n_pcs: int, seed: int):
    n_components = max(2, min(n_pcs, X.shape[0] - 1, X.shape[1] - 1))
    svd = TruncatedSVD(n_components=n_components, random_state=seed)
    embedding = svd.fit_transform(X)
    return embedding, svd.explained_variance_ratio_


def _cluster_leiden(embedding: np.ndarray, seed: int, resolution: float = 1.0) -> np.ndarray:
    """Cluster an embedding with Leiden via a temporary AnnData object.

    This mirrors the production workflow (neighbors + Leiden) and does not
    require a pre-specified number of clusters.
    """
    adata = ad.AnnData(X=np.zeros((embedding.shape[0], 2)))
    adata.obsm["X_emb"] = embedding
    # Use a small number of neighbors for the stability assessment; this is a
    # proxy metric, not a final clustering.
    n_neighbors = min(15, embedding.shape[0] - 1)
    sc.pp.neighbors(adata, use_rep="X_emb", n_neighbors=max(2, n_neighbors))
    # Use igraph backend and two iterations to match scanpy future defaults and
    # avoid the leidenalg deprecation warning.
    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=seed,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    return adata.obs["leiden"].to_numpy()


def _neighbor_overlap(embedding_a: np.ndarray, embedding_b: np.ndarray, n_neighbors: int) -> float:
    n_neighbors = max(2, min(n_neighbors, embedding_a.shape[0] - 1, embedding_b.shape[0] - 1))
    neigh_a = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(embedding_a).kneighbors(return_distance=False)[:, 1:]
    neigh_b = NearestNeighbors(n_neighbors=n_neighbors + 1).fit(embedding_b).kneighbors(return_distance=False)[:, 1:]
    overlaps = []
    for left, right in zip(neigh_a, neigh_b):
        overlaps.append(len(set(left) & set(right)) / n_neighbors)
    return float(np.mean(overlaps)) if overlaps else float("nan")


def _rare_population_rows(dataset: str, labels: pd.Series, clusters: np.ndarray, n_pcs: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if labels.empty:
        rows.append(
            {
                "dataset": dataset,
                "n_pcs": n_pcs,
                "label_key": labels.name if labels.name else "unknown",
                "rare_population": "__missing__",
                "n_cells": 0,
                "fraction": float("nan"),
                "max_cluster_fraction": float("nan"),
                "preservation_proxy": "missing_annotation",
                "review_required": False,
            }
        )
        return rows
    counts = labels.astype(str).value_counts(dropna=False)
    rare = counts[(counts >= 10) & (counts / len(labels) <= 0.05)]
    cluster_series = pd.Series(clusters, index=labels.index)
    for label, count in rare.items():
        mask = labels.astype(str) == str(label)
        cluster_counts = cluster_series[mask].value_counts()
        max_cluster_fraction = float(cluster_counts.max() / max(int(mask.sum()), 1))
        rows.append(
            {
                "dataset": dataset,
                "n_pcs": n_pcs,
                "label_key": labels.name,
                "rare_population": str(label),
                "n_cells": int(count),
                "fraction": float(count / len(labels)),
                "max_cluster_fraction": max_cluster_fraction,
                "preservation_proxy": "concentrated" if max_cluster_fraction >= 0.5 else "diffuse",
                "review_required": max_cluster_fraction < 0.5,
            }
        )
    return rows


def run(
    output_dir: Path,
    datasets: set[str] | None,
    max_cells: int | None,
    seed: int,
    n_pcs_values: list[int],
    n_neighbors_values: list[int],
    leiden_resolution: float = 1.0,
) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    pca_rows: list[dict[str, Any]] = []
    stability_rows: list[dict[str, Any]] = []
    rare_rows: list[dict[str, Any]] = []
    figure_rows: list[dict[str, Any]] = []

    for spec in DATASETS:
        if datasets and spec.key not in datasets:
            continue
        if not spec.preprocess_roles or not spec.path.exists() or spec.key == "cellbender_tiny":
            continue
        adata = _subset(ad.read_h5ad(spec.path), max_cells=max_cells, seed=seed)
        X, _ = _top_variable_subset(_log_normalized_matrix(adata), adata.var_names, n_genes=2000)
        embeddings: dict[int, np.ndarray] = {}
        cluster_by_pcs_seed: dict[tuple[int, int], np.ndarray] = {}
        for n_pcs in n_pcs_values:
            embedding, variance_ratio = _pca_embedding(X, n_pcs=n_pcs, seed=seed)
            embeddings[n_pcs] = embedding
            pca_rows.append(
                {
                    "dataset": spec.key,
                    "n_pcs": int(embedding.shape[1]),
                    "variance_explained_total": float(np.sum(variance_ratio)),
                    "variance_explained_top3": json.dumps([float(v) for v in variance_ratio[:3]]),
                    "review_required": bool(np.sum(variance_ratio) < 0.2),
                }
            )
            for cluster_seed in (seed, seed + 1, seed + 2):
                cluster_by_pcs_seed[(n_pcs, cluster_seed)] = _cluster_leiden(
                    embedding, seed=cluster_seed, resolution=leiden_resolution
                )
            base = cluster_by_pcs_seed[(n_pcs, seed)]
            for cluster_seed in (seed + 1, seed + 2):
                ari = adjusted_rand_score(base, cluster_by_pcs_seed[(n_pcs, cluster_seed)])
                stability_rows.append(
                    {
                        "dataset": spec.key,
                        "stability_type": "leiden_seed",
                        "n_pcs": n_pcs,
                        "n_neighbors": "",
                        "seed_a": seed,
                        "seed_b": cluster_seed,
                        "score": float(ari),
                        "review_required": bool(ari < 0.5),
                    }
                )
            if "cell_type" in adata.obs:
                rare_rows.extend(_rare_population_rows(spec.key, adata.obs["cell_type"], base, n_pcs))
            else:
                rare_rows.extend(
                    _rare_population_rows(
                        spec.key,
                        pd.Series(name="cell_type"),
                        base,
                        n_pcs,
                    )
                )
        sorted_pcs = sorted(embeddings)
        if len(sorted_pcs) >= 2:
            low, high = sorted_pcs[0], sorted_pcs[-1]
            for n_neighbors in n_neighbors_values:
                overlap = _neighbor_overlap(embeddings[low], embeddings[high], n_neighbors=n_neighbors)
                stability_rows.append(
                    {
                        "dataset": spec.key,
                        "stability_type": "neighbor_npcs",
                        "n_pcs": f"{low}_vs_{high}",
                        "n_neighbors": n_neighbors,
                        "seed_a": seed,
                        "seed_b": seed,
                        "score": overlap,
                        "review_required": bool(overlap < 0.5),
                    }
                )
                figure_rows.append(
                    {
                        "figure_panel": "3D",
                        "dataset": spec.key,
                        "metric": "neighbor_overlap_low_vs_high_pcs",
                        "value": overlap,
                        "context": json.dumps({"low_n_pcs": low, "high_n_pcs": high, "n_neighbors": n_neighbors}),
                    }
                )

    paths = {
        "pca": output_dir / "pca_neighbors_stability.tsv",
        "clustering": output_dir / "leiden_seed_stability.tsv",
        "rare": output_dir / "rare_population_preservation.tsv",
        "figure3": output_dir / "figure3_graph_data.tsv",
    }
    pd.DataFrame(pca_rows).to_csv(paths["pca"], sep="\t", index=False)
    pd.DataFrame(stability_rows).to_csv(paths["clustering"], sep="\t", index=False)
    pd.DataFrame(rare_rows).to_csv(paths["rare"], sep="\t", index=False)
    pd.DataFrame(figure_rows).to_csv(paths["figure3"], sep="\t", index=False)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/work/preprocess_graph_stability"),
    )
    parser.add_argument("--datasets", nargs="*", help="Optional dataset keys to include.")
    parser.add_argument("--max-cells", type=int, default=2000, help="Pilot subset size. Use 0 for full data.")
    parser.add_argument("--seed", type=int, default=19)
    parser.add_argument("--n-pcs", nargs="*", type=int, default=[20, 40])
    parser.add_argument("--n-neighbors", nargs="*", type=int, default=[15, 30])
    parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=1.0,
        help="Resolution passed to scanpy.tl.leiden for stability assessment.",
    )
    args = parser.parse_args()
    paths = run(
        args.output_dir,
        datasets=set(args.datasets) if args.datasets else None,
        max_cells=args.max_cells,
        seed=args.seed,
        n_pcs_values=args.n_pcs,
        n_neighbors_values=args.n_neighbors,
        leiden_resolution=args.leiden_resolution,
    )
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
