#!/usr/bin/env python3
"""Run the controlled mixology held-out preprocessing and integration benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sparse
from scipy.spatial.distance import jensenshannon
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import (
    adjusted_rand_score,
    balanced_accuracy_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import scLucid as scl
from validation.qc_preprocess.locked_acceptance import evaluate_preprocess_policy_acceptance

CANDIDATES = {
    "standard_unintegrated": ("log1p", "batch_aware_hvg"),
    "pearson_residuals": ("pearson", "batch_aware_hvg"),
    "multinomial_deviance": ("log1p", "deviance"),
    "pearson_residuals_deviance": ("pearson", "deviance"),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _counts(adata: ad.AnnData) -> sparse.csr_matrix:
    matrix = adata.layers.get("counts", adata.X)
    matrix = sparse.csr_matrix(matrix)
    if matrix.data.size == 0 or np.any(matrix.data < 0):
        raise ValueError("A non-empty, non-negative count matrix is required.")
    if not np.allclose(matrix.data, np.rint(matrix.data), atol=1e-6):
        raise ValueError("The mixology benchmark requires integer-like counts.")
    return matrix


def _log_normalize(counts: sparse.csr_matrix, target_sum: float = 1e4) -> sparse.csr_matrix:
    totals = np.asarray(counts.sum(axis=1)).ravel().astype(float)
    scale = np.divide(target_sum, totals, out=np.zeros_like(totals), where=totals > 0)
    result = counts.astype(float).multiply(scale[:, None]).tocsr()
    result.data = np.log1p(result.data)
    return result


def _batch_aware_hvg(
    counts: sparse.csr_matrix,
    protocols: pd.Series,
    n_top_genes: int,
) -> np.ndarray:
    normalized = _log_normalize(counts)
    work = ad.AnnData(X=normalized, obs=pd.DataFrame({"protocol": protocols.to_numpy()}))
    n_top = min(int(n_top_genes), work.n_vars)
    sc.pp.highly_variable_genes(
        work,
        n_top_genes=n_top,
        flavor="seurat",
        batch_key="protocol",
        subset=False,
        inplace=True,
    )
    selected = np.flatnonzero(work.var["highly_variable"].to_numpy(bool))
    if len(selected) < 2:
        raise RuntimeError("Batch-aware HVG selection returned fewer than two genes.")
    return selected


def _multinomial_deviance(counts: sparse.csr_matrix, n_top_genes: int) -> np.ndarray:
    matrix = counts.tocsc().astype(float)
    cell_totals = np.asarray(matrix.sum(axis=1)).ravel()
    gene_totals = np.asarray(matrix.sum(axis=0)).ravel()
    grand_total = float(gene_totals.sum())
    if grand_total <= 0:
        raise ValueError("Cannot calculate deviance from an empty count matrix.")
    probabilities = gene_totals / grand_total
    deviance = np.zeros(matrix.shape[1], dtype=float)
    for gene in range(matrix.shape[1]):
        start, stop = matrix.indptr[gene], matrix.indptr[gene + 1]
        values = matrix.data[start:stop]
        rows = matrix.indices[start:stop]
        expected = cell_totals[rows] * probabilities[gene]
        valid = (values > 0) & (expected > 0)
        if valid.any():
            deviance[gene] = 2.0 * np.sum(values[valid] * np.log(values[valid] / expected[valid]))
    n_top = min(int(n_top_genes), matrix.shape[1])
    return np.argsort(-deviance, kind="mergesort")[:n_top]


def _transform_pair(
    train_counts: sparse.csr_matrix,
    test_counts: sparse.csr_matrix,
    features: np.ndarray,
    method: str,
) -> tuple[np.ndarray, np.ndarray]:
    if method == "log1p":
        train = _log_normalize(train_counts)[:, features].toarray()
        test = _log_normalize(test_counts)[:, features].toarray()
        return train, test
    if method != "pearson":
        raise ValueError(f"Unknown transformation: {method}")

    train_total = np.asarray(train_counts.sum(axis=1)).ravel().astype(float)
    test_total = np.asarray(test_counts.sum(axis=1)).ravel().astype(float)
    grand_total = max(float(train_total.sum()), 1.0)
    gene_probability = np.asarray(train_counts[:, features].sum(axis=0)).ravel() / grand_total
    theta = 100.0
    clip = np.sqrt(max(train_counts.shape[0], 1))

    def residuals(matrix, totals):
        observed = matrix[:, features].toarray().astype(float)
        expected = totals[:, None] * gene_probability[None, :]
        denominator = np.sqrt(expected + expected * expected / theta)
        values = np.divide(
            observed - expected,
            denominator,
            out=np.zeros_like(expected),
            where=denominator > 0,
        )
        return np.clip(values, -clip, clip)

    return residuals(train_counts, train_total), residuals(test_counts, test_total)


def _fit_pca(
    train: np.ndarray,
    test: np.ndarray,
    *,
    n_pcs: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    scaler = StandardScaler()
    train_scaled = scaler.fit_transform(train)
    test_scaled = scaler.transform(test)
    components = max(2, min(int(n_pcs), train.shape[0] - 1, train.shape[1] - 1))
    model = PCA(n_components=components, svd_solver="randomized", random_state=seed)
    return model.fit_transform(train_scaled), model.transform(test_scaled)


def _heldout_accuracy(
    train_embedding: np.ndarray,
    test_embedding: np.ndarray,
    train_labels: pd.Series,
    test_labels: pd.Series,
) -> float:
    neighbors = min(15, len(train_labels))
    classifier = KNeighborsClassifier(n_neighbors=neighbors, weights="distance")
    classifier.fit(train_embedding, train_labels.astype(str))
    predicted = classifier.predict(test_embedding)
    return float(balanced_accuracy_score(test_labels.astype(str), predicted))


def _neighbor_indices(embedding: np.ndarray, n_neighbors: int = 30) -> np.ndarray:
    neighbors = min(int(n_neighbors), len(embedding) - 1)
    model = NearestNeighbors(n_neighbors=neighbors + 1).fit(embedding)
    return model.kneighbors(embedding, return_distance=False)[:, 1:]


def _cross_protocol_identity_purity(
    embedding: np.ndarray,
    protocols: pd.Series,
    identities: pd.Series,
    n_neighbors: int = 15,
) -> float:
    protocol_values = protocols.astype(str).to_numpy()
    identity_values = identities.astype(str).to_numpy()
    scores: list[float] = []
    for protocol in sorted(pd.unique(protocol_values)):
        query = np.flatnonzero(protocol_values == protocol)
        reference = np.flatnonzero(protocol_values != protocol)
        neighbors = min(int(n_neighbors), len(reference))
        index = NearestNeighbors(n_neighbors=neighbors).fit(embedding[reference])
        found = index.kneighbors(embedding[query], return_distance=False)
        matched = identity_values[reference][found] == identity_values[query, None]
        scores.extend(matched.mean(axis=1).tolist())
    return float(np.mean(scores))


def _batch_mixing_js_similarity(
    embedding: np.ndarray,
    protocols: pd.Series,
    n_neighbors: int = 30,
) -> float:
    protocol_values = protocols.astype(str).to_numpy()
    levels = sorted(pd.unique(protocol_values))
    global_counts = pd.Series(protocol_values).value_counts().reindex(levels, fill_value=0)
    global_distribution = global_counts.to_numpy(dtype=float) / len(protocol_values)
    neighbor_indices = _neighbor_indices(embedding, n_neighbors=n_neighbors)
    similarities: list[float] = []
    for neighbors in neighbor_indices:
        local_counts = pd.Series(protocol_values[neighbors]).value_counts().reindex(levels, fill_value=0)
        local_distribution = local_counts.to_numpy(dtype=float) / len(neighbors)
        divergence = float(jensenshannon(local_distribution, global_distribution, base=2.0) ** 2)
        similarities.append(1.0 - divergence)
    return float(np.mean(similarities))


def _graph_overlap(left: np.ndarray, right: np.ndarray, n_neighbors: int = 30) -> float:
    left_neighbors = _neighbor_indices(left, n_neighbors=n_neighbors)
    right_neighbors = _neighbor_indices(right, n_neighbors=n_neighbors)
    return float(
        np.mean(
            [
                len(set(a).intersection(b)) / len(a)
                for a, b in zip(left_neighbors, right_neighbors)
            ]
        )
    )


def _embedding_metrics(
    embedding: np.ndarray,
    identities: pd.Series,
    protocols: pd.Series,
    *,
    seed: int,
) -> dict[str, float]:
    labels = identities.astype(str).to_numpy()
    clusters = KMeans(n_clusters=len(pd.unique(labels)), n_init=20, random_state=seed).fit_predict(
        embedding
    )
    neighbors = _neighbor_indices(embedding)
    purity = float(np.mean(labels[neighbors] == labels[:, None]))
    return {
        "identity_neighbor_purity": purity,
        "cross_protocol_identity_purity": _cross_protocol_identity_purity(
            embedding, protocols, identities
        ),
        "identity_silhouette": float(silhouette_score(embedding, labels)),
        "kmeans_identity_ari": float(adjusted_rand_score(labels, clusters)),
        "kmeans_identity_nmi": float(normalized_mutual_info_score(labels, clusters)),
        "batch_mixing_js_similarity": _batch_mixing_js_similarity(embedding, protocols),
    }


def _cramers_v(left: pd.Series, right: pd.Series) -> float:
    table = pd.crosstab(left.astype(str), right.astype(str)).to_numpy(dtype=float)
    total = table.sum()
    if total <= 0 or min(table.shape) < 2:
        return 0.0
    expected = table.sum(axis=1, keepdims=True) @ table.sum(axis=0, keepdims=True) / total
    chi2 = np.divide(
        (table - expected) ** 2,
        expected,
        out=np.zeros_like(table),
        where=expected > 0,
    ).sum()
    return float(np.sqrt((chi2 / total) / max(min(table.shape) - 1, 1)))


def _harmony_embedding(embedding: np.ndarray, obs: pd.DataFrame) -> np.ndarray:
    import scanpy.external as sce

    work = ad.AnnData(X=np.zeros((len(obs), 1), dtype=np.float32), obs=obs.copy())
    work.obsm["X_pca"] = embedding.copy()
    sce.pp.harmony_integrate(
        work,
        key="protocol",
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        max_iter_harmony=50,
    )
    return np.asarray(work.obsm["X_pca_harmony"])


def _run_product_policy(adata: ad.AnnData) -> tuple[str, dict[str, Any]]:
    context = scl.ProjectContext(
        dataset_type="cell_line",
        assay="scrna",
        sample_key="protocol",
        batch_key="protocol",
        condition_key="condition" if "condition" in adata.obs else None,
        cell_type_key="mixology_identity",
        is_multi_sample=True,
        input_provenance="filtered_counts",
    )
    card = scl.recommend_preprocess_policy(adata, context, consumer="exploration")
    policy = card.policy
    mapping = {
        ("standard", "scanpy", False): "standard_unintegrated",
        ("pearson_residuals", "scanpy", False): "pearson_residuals",
        ("standard", "deviance", False): "multinomial_deviance",
        ("pearson_residuals", "deviance", False): "pearson_residuals_deviance",
    }
    selected = mapping.get(
        (
            policy.normalization_method,
            policy.feature_selection_method,
            bool(policy.run_integration),
        ),
        "UNMAPPED_POLICY",
    )
    evidence: dict[str, Any] = {
        "review_status": card.status,
        "recommended": card.recommended,
        "reason": card.reason,
        "selected_candidate": selected,
        "integration_review": card.details.get("integration_review", {}),
        "execution_status": "NOT_RUN",
        "representation_contract": {},
    }
    if card.status != "BLOCKED":
        execution = scl.apply_preprocess_policy(adata, policy)
        result = execution.adata
        evidence["execution_status"] = execution.status
        evidence["representation_contract"] = {
            "counts_preserved": "counts" in result.layers,
            "normalized_full_present": "normalized_full" in result.layers,
            "discovery_feature_present": "discovery_feature" in result.var,
            "discovery_rep_present": "X_pca" in result.obsm,
            "integrated_rep_selected": bool(policy.run_integration),
        }
    return selected, evidence


def _integration_pareto(
    baseline: dict[str, float],
    candidate: dict[str, float],
    *,
    tolerances: dict[str, float] | None = None,
    improvements: dict[str, float] | None = None,
) -> dict[str, Any]:
    tolerances = tolerances or {
        "cross_protocol_identity_purity": 0.02,
        "identity_silhouette": 0.02,
        "batch_mixing_js_similarity": 0.01,
        "graph_seed_stability": 0.02,
    }
    improvements = improvements or {
        "cross_protocol_identity_purity": 0.01,
        "identity_silhouette": 0.01,
        "batch_mixing_js_similarity": 0.02,
        "graph_seed_stability": 0.01,
    }
    deltas = {key: float(candidate[key] - baseline[key]) for key in tolerances}
    no_material_loss = all(deltas[key] >= -tolerances[key] for key in tolerances)
    meaningful_gain = any(deltas[key] >= improvements[key] for key in improvements)
    return {
        "pareto_dominates_baseline": bool(no_material_loss and meaningful_gain),
        "no_material_loss": no_material_loss,
        "meaningful_gain": meaningful_gain,
        "deltas": deltas,
        "loss_tolerances": tolerances,
        "gain_thresholds": improvements,
    }


def run(
    input_path: Path,
    output_dir: Path,
    *,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    seed: int = 17,
    run_harmony: bool = True,
    contract_path: Path | None = None,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_path)
    required = {"protocol", "mixology_identity"}
    missing = sorted(required - set(adata.obs.columns))
    if missing:
        raise ValueError(f"Missing controlled mixology metadata: {missing}")
    counts = _counts(adata)
    protocols = adata.obs["protocol"].astype(str).reset_index(drop=True)
    identities = adata.obs["mixology_identity"].astype(str).reset_index(drop=True)
    if protocols.nunique() < 3 or identities.nunique() < 2:
        raise ValueError("At least three protocols and two identities are required.")
    selected_candidate, product_policy = _run_product_policy(adata)

    fold_rows: list[dict[str, Any]] = []
    for holdout in sorted(protocols.unique()):
        test_mask = protocols.to_numpy() == holdout
        train_mask = ~test_mask
        train_counts = counts[train_mask]
        test_counts = counts[test_mask]
        train_protocols = protocols[train_mask].reset_index(drop=True)
        features = {
            "batch_aware_hvg": _batch_aware_hvg(
                train_counts, train_protocols, n_top_genes=n_top_genes
            ),
            "deviance": _multinomial_deviance(train_counts, n_top_genes=n_top_genes),
        }
        for candidate, (normalization, feature_method) in CANDIDATES.items():
            train, test = _transform_pair(
                train_counts,
                test_counts,
                features[feature_method],
                normalization,
            )
            train_embedding, test_embedding = _fit_pca(
                train, test, n_pcs=n_pcs, seed=seed
            )
            fold_rows.append(
                {
                    "dataset": "public_mixology",
                    "candidate": candidate,
                    "holdout_protocol": holdout,
                    "n_train": int(train_mask.sum()),
                    "n_test": int(test_mask.sum()),
                    "n_features": int(len(features[feature_method])),
                    "identity_balanced_accuracy": _heldout_accuracy(
                        train_embedding,
                        test_embedding,
                        identities[train_mask],
                        identities[test_mask],
                    ),
                }
            )

    full_features = {
        "batch_aware_hvg": _batch_aware_hvg(counts, protocols, n_top_genes=n_top_genes),
        "deviance": _multinomial_deviance(counts, n_top_genes=n_top_genes),
    }
    full_embeddings: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    candidate_rows: list[dict[str, Any]] = []
    fold_frame = pd.DataFrame(fold_rows)
    for candidate, (normalization, feature_method) in CANDIDATES.items():
        transformed, _ = _transform_pair(
            counts,
            counts,
            full_features[feature_method],
            normalization,
        )
        embedding_a, _ = _fit_pca(transformed, transformed, n_pcs=n_pcs, seed=seed)
        embedding_b, _ = _fit_pca(transformed, transformed, n_pcs=n_pcs, seed=seed + 1)
        full_embeddings[candidate] = (embedding_a, embedding_b)
        metrics = _embedding_metrics(embedding_a, identities, protocols, seed=seed)
        metrics["graph_seed_stability"] = _graph_overlap(embedding_a, embedding_b)
        heldout = fold_frame[fold_frame["candidate"] == candidate]
        candidate_rows.append(
            {
                "dataset": "public_mixology",
                "candidate": candidate,
                "selected": candidate == selected_candidate,
                "preregistered_task_utility": float(
                    heldout["identity_balanced_accuracy"].mean()
                ),
                "heldout_protocol_sd": float(heldout["identity_balanced_accuracy"].std(ddof=1)),
                "n_features": int(len(full_features[feature_method])),
                **metrics,
            }
        )

    candidate_frame = pd.DataFrame(candidate_rows)
    baseline_purity = float(
        candidate_frame.loc[
            candidate_frame["candidate"] == "standard_unintegrated",
            "cross_protocol_identity_purity",
        ].iloc[0]
    )
    candidate_frame["biology_loss"] = np.maximum(
        0.0,
        (baseline_purity - candidate_frame["cross_protocol_identity_purity"])
        / max(abs(baseline_purity), 1e-12),
    )
    contract = json.loads(contract_path.read_text()) if contract_path else {}
    primary_contract = contract.get("preprocess_primary_endpoint", {})
    acceptance = evaluate_preprocess_policy_acceptance(
        candidate_frame,
        max_regret=float(primary_contract.get("max_held_out_regret", 0.05)),
        max_biology_loss=float(primary_contract.get("max_biology_loss", 0.02)),
    )

    baseline_a, baseline_b = full_embeddings["standard_unintegrated"]
    baseline_metrics = _embedding_metrics(baseline_a, identities, protocols, seed=seed)
    baseline_metrics["graph_seed_stability"] = _graph_overlap(baseline_a, baseline_b)
    integration_rows = [
        {
            "method": "unintegrated",
            "status": "EVALUATED",
            "selected": True,
            **baseline_metrics,
        }
    ]
    confounding = _cramers_v(protocols, identities)
    integration_review: dict[str, Any] = {
        "status": "NOT_RUN",
        "batch_identity_cramers_v": confounding,
        "reason": "Harmony comparison was disabled.",
    }
    if confounding >= 0.7:
        integration_review = {
            "status": "BLOCKED",
            "batch_identity_cramers_v": confounding,
            "reason": "Protocol is strongly confounded with cell-line identity.",
        }
    elif run_harmony:
        try:
            harmony_a = _harmony_embedding(baseline_a, adata.obs)
            harmony_b = _harmony_embedding(baseline_b, adata.obs)
            harmony_metrics = _embedding_metrics(harmony_a, identities, protocols, seed=seed)
            harmony_metrics["graph_seed_stability"] = _graph_overlap(harmony_a, harmony_b)
            design_contract = contract.get("preprocess_validation_design", {})
            pareto = _integration_pareto(
                baseline_metrics,
                harmony_metrics,
                tolerances=design_contract.get("integration_absolute_loss_tolerances"),
                improvements=design_contract.get("integration_minimum_gain_thresholds"),
            )
            integration_rows.append(
                {
                    "method": "harmony",
                    "status": "EVALUATED",
                    "selected": False,
                    **harmony_metrics,
                    **{f"delta_{key}": value for key, value in pareto["deltas"].items()},
                    "pareto_dominates_baseline": pareto["pareto_dominates_baseline"],
                }
            )
            integration_review = {
                "status": "REVIEW" if pareto["pareto_dominates_baseline"] else "PASS_BASELINE",
                "batch_identity_cramers_v": confounding,
                **pareto,
                "reason": (
                    "Harmony shows a controlled-dataset Pareto gain; selector update requires external confirmation."
                    if pareto["pareto_dominates_baseline"]
                    else "Harmony did not Pareto-dominate the unintegrated baseline."
                ),
            }
        except Exception as exc:
            integration_review = {
                "status": "NOT_EVALUABLE",
                "batch_identity_cramers_v": confounding,
                "reason": f"Harmony comparison failed: {type(exc).__name__}: {exc}",
            }

    fold_path = output_dir / "mixology_heldout_protocol_metrics.tsv"
    candidate_path = output_dir / "mixology_candidate_metrics.tsv"
    integration_path = output_dir / "mixology_integration_pareto.tsv"
    fold_frame.to_csv(fold_path, sep="\t", index=False)
    candidate_frame.to_csv(candidate_path, sep="\t", index=False)
    pd.DataFrame(integration_rows).to_csv(integration_path, sep="\t", index=False)

    report = {
        "schema_version": "sclucid_mixology_preprocess_benchmark_v1",
        "status": "REVIEW",
        "dataset": "public_mixology",
        "dataset_registry_id": "scmixology_gse118767",
        "preregistered_endpoint_ids": [
            "pp_selector_regret",
            "pp_integration_pareto",
            "pp_identity_preservation",
        ],
        "source_path": str(input_path.resolve()),
        "source_sha256": _sha256(input_path),
        "acceptance_contract": (
            {"path": str(contract_path.resolve()), "sha256": _sha256(contract_path)}
            if contract_path
            else None
        ),
        "estimand": "Mean balanced cell-line identity accuracy across held-out protocols.",
        "experimental_unit": "protocol",
        "n_protocols": int(protocols.nunique()),
        "n_cells": int(adata.n_obs),
        "candidate_acceptance": acceptance,
        "product_policy": product_policy,
        "integration_review": integration_review,
        "release_gate": {
            "status": "BLOCKED",
            "reason": "Controlled mixology evidence cannot replace external real-project validation.",
        },
        "claim_boundary": {
            "supported": [
                "Candidate representations were compared using leave-one-protocol-out identity truth."
            ],
            "exploratory": [
                "Performance and integration Pareto results are specific to this controlled mixture."
            ],
            "unsupported": [
                "Universal preprocessing superiority or tumor-project benefit."
            ],
        },
        "artifacts": {
            "fold_metrics": str(fold_path),
            "candidate_metrics": str(candidate_path),
            "integration_pareto": str(integration_path),
        },
        "next_action": (
            "Investigate held-out regret before changing the simple baseline."
            if acceptance["status"] != "PASS"
            else "Keep the simple baseline provisionally and test the same contract on external projects."
        ),
    }
    json_path = output_dir / "mixology_preprocess_benchmark.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n")
    lines = [
        "# Mixology preprocessing benchmark",
        "",
        f"Status: **{report['status']}**",
        "",
        f"Candidate acceptance: **{acceptance['status']}**",
        f"Integration review: **{integration_review['status']}**",
        "",
        "The experimental unit is protocol; cell-level observations are not treated as independent replicates.",
        "",
        "## Next action",
        "",
        report["next_action"],
        "",
    ]
    (output_dir / "mixology_preprocess_benchmark.md").write_text("\n".join(lines))
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/public_mixology.h5ad"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--no-harmony", action="store_true")
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("validation/qc_preprocess/acceptance_contract.json"),
    )
    args = parser.parse_args()
    report = run(
        args.input,
        args.output_dir,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        seed=args.seed,
        run_harmony=not args.no_harmony,
        contract_path=args.contract,
    )
    print(json.dumps({"status": report["status"], "output_dir": str(args.output_dir)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
