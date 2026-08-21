#!/usr/bin/env python3
"""Exercise the public preprocessing policy API against its execution contract.

This benchmark intentionally answers an engineering question: does the public
review/apply boundary preserve the promised expression spaces and fail closed?
It does not estimate scientific performance or replace held-out mixology and
real-project validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy import sparse

import scLucid as scl

DEFAULT_OUTPUT = Path(
    "validation_outputs/current/preprocess_contract/controlled_preprocess_contract_benchmark.json"
)


def _matrix_digest(matrix: Any) -> str:
    digest = hashlib.sha256()
    if sparse.issparse(matrix):
        csr = matrix.tocsr()
        digest.update(str(csr.shape).encode("utf-8"))
        digest.update(str(csr.dtype).encode("utf-8"))
        digest.update(np.ascontiguousarray(csr.data).tobytes())
        digest.update(np.ascontiguousarray(csr.indices).tobytes())
        digest.update(np.ascontiguousarray(csr.indptr).tobytes())
    else:
        values = np.ascontiguousarray(matrix)
        digest.update(str(values.shape).encode("utf-8"))
        digest.update(str(values.dtype).encode("utf-8"))
        digest.update(values.tobytes())
    return digest.hexdigest()


def _frame_digest(frame: pd.DataFrame) -> str:
    values = pd.util.hash_pandas_object(frame, index=True).to_numpy(np.uint64)
    return hashlib.sha256(values.tobytes()).hexdigest()


def _snapshot(adata: AnnData) -> dict[str, Any]:
    """Capture all mutable slots populated in the controlled fixture."""
    return {
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "X": _matrix_digest(adata.X),
        "obs": _frame_digest(adata.obs),
        "var": _frame_digest(adata.var),
        "layers": {str(key): _matrix_digest(value) for key, value in sorted(adata.layers.items())},
        "obsm_keys": sorted(map(str, adata.obsm.keys())),
        "varm_keys": sorted(map(str, adata.varm.keys())),
        "obsp_keys": sorted(map(str, adata.obsp.keys())),
        "uns_repr": repr(adata.uns),
        "raw_present": adata.raw is not None,
    }


def _matrix_equal(left: Any, right: Any, *, atol: float = 0.0) -> bool:
    if left.shape != right.shape:
        return False
    if sparse.issparse(left) or sparse.issparse(right):
        left_csr = sparse.csr_matrix(left)
        right_csr = sparse.csr_matrix(right)
        difference = left_csr - right_csr
        if difference.nnz == 0:
            return True
        return bool(np.max(np.abs(difference.data)) <= atol)
    return bool(np.allclose(np.asarray(left), np.asarray(right), atol=atol, rtol=0.0))


def _controlled_umi_adata(seed: int = 29) -> AnnData:
    """Create a small deterministic UMI object with crossed sample biology."""
    rng = np.random.default_rng(seed)
    n_cells = 120
    n_genes = 160
    sample = np.repeat(["S1", "S2"], n_cells // 2)
    lineage = np.tile(np.repeat(["T", "B", "Myeloid"], 20), 2)
    condition = np.tile(np.repeat(["control", "treated"], 10), 6)

    rates = np.full((n_cells, n_genes), 0.35, dtype=np.float64)
    for index, label in enumerate(("T", "B", "Myeloid")):
        mask = lineage == label
        rates[mask, 10 + index * 30 : 30 + index * 30] += 3.0
    rates[sample == "S2"] *= 1.15
    counts = rng.poisson(rates).astype(np.int32)

    obs = pd.DataFrame(
        {
            "sample": sample,
            "batch": sample,
            "condition": condition,
            "cell_type": lineage,
        },
        index=[f"cell_{index:03d}" for index in range(n_cells)],
    )
    var = pd.DataFrame(
        index=[
            *(f"MT-GENE{index}" for index in range(5)),
            *(f"GENE{index:03d}" for index in range(5, n_genes)),
        ]
    )
    adata = AnnData(X=sparse.csr_matrix(counts), obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    return adata


def _check(observed: Any, expected: Any, passed: bool) -> dict[str, Any]:
    return {"passed": bool(passed), "observed": observed, "expected": expected}


def _legacy_dense_pca_reference(adata: AnnData) -> np.ndarray:
    """Reproduce the pre-fix centered-PCA math without implicit densification."""
    mask = adata.var["discovery_feature"].to_numpy(bool)
    source = adata.layers["normalized_full"][:, mask]
    dense = source.toarray() if sparse.issparse(source) else np.asarray(source).copy()
    work = AnnData(
        X=dense,
        obs=adata.obs.copy(),
        var=adata.var.loc[mask].copy(),
    )
    sc.pp.scale(work, max_value=10.0, zero_center=True)
    n_comps = min(50, work.n_obs - 1, work.n_vars - 1)
    sc.tl.pca(work, n_comps=n_comps, svd_solver="arpack")
    return np.asarray(work.obsm["X_pca"])


def _geometry_metrics(observed: np.ndarray, reference: np.ndarray) -> dict[str, float]:
    observed_gram = observed @ observed.T
    reference_gram = reference @ reference.T
    denominator = max(float(np.linalg.norm(reference_gram)), np.finfo(float).eps)
    relative_gram_error = float(np.linalg.norm(observed_gram - reference_gram) / denominator)

    n_neighbors = min(15, observed.shape[0] - 1)
    observed_distances = np.sum(
        (observed[:, None, :] - observed[None, :, :]) ** 2,
        axis=2,
    )
    reference_distances = np.sum(
        (reference[:, None, :] - reference[None, :, :]) ** 2,
        axis=2,
    )
    observed_neighbors = np.argsort(observed_distances, axis=1)[:, 1 : n_neighbors + 1]
    reference_neighbors = np.argsort(reference_distances, axis=1)[:, 1 : n_neighbors + 1]
    overlap = np.mean(
        [
            len(set(left).intersection(right)) / n_neighbors
            for left, right in zip(observed_neighbors, reference_neighbors, strict=True)
        ]
    )
    return {
        "relative_pca_gram_error": relative_gram_error,
        "mean_neighbor_retention": float(overlap),
        "n_neighbors": int(n_neighbors),
    }


def run(output_path: Path = DEFAULT_OUTPUT, *, seed: int = 29) -> dict[str, Any]:
    """Run the controlled public-API contract benchmark and write JSON evidence."""
    adata = _controlled_umi_adata(seed=seed)
    input_snapshot = _snapshot(adata)
    input_counts = adata.layers["counts"].copy()
    context = scl.ProjectContext(
        dataset_type="cell_line",
        assay="scrna",
        input_provenance="filtered_counts",
        sample_key="sample",
        batch_key="batch",
        condition_key="condition",
        cell_type_key="cell_type",
        is_multi_sample=True,
    )

    first_card = scl.recommend_preprocess_policy(adata, context, consumer="exploration")
    second_card = scl.recommend_preprocess_policy(adata, context, consumer="exploration")
    integration_card = scl.recommend_preprocess_policy(adata, context, consumer="integration")
    inference_card = scl.recommend_preprocess_policy(
        adata,
        context,
        consumer="expression_inference",
    )
    after_review_snapshot = _snapshot(adata)

    with warnings.catch_warnings(record=True) as apply_warnings:
        warnings.simplefilter("always")
        first_evidence = scl.apply_preprocess_policy(adata, first_card.policy)
        second_evidence = scl.apply_preprocess_policy(adata, second_card.policy)
    first = first_evidence.adata
    second = second_evidence.adata
    after_apply_snapshot = _snapshot(adata)
    representation = first.uns["sclucid"]["preprocess"]["representation_contract"]
    resource_contract = first_evidence.result["discovery_temporary_contract"]
    implicit_densification_warnings = [
        str(item.message)
        for item in apply_warnings
        if "zero-centering a sparse array/matrix densifies it" in str(item.message)
    ]
    legacy_embedding = _legacy_dense_pca_reference(first)
    geometry = _geometry_metrics(
        np.asarray(first.obsm["X_pca"]),
        legacy_embedding,
    )

    mutated = adata.copy()
    mutated_counts = mutated.layers["counts"].tolil(copy=True)
    mutated_counts[0, 0] = int(mutated_counts[0, 0]) + 1
    mutated.layers["counts"] = mutated_counts.tocsr()
    fingerprint_rejected = False
    fingerprint_error = ""
    try:
        scl.apply_preprocess_policy(mutated, first_card.policy)
    except ValueError as exc:
        fingerprint_error = str(exc)
        fingerprint_rejected = "fingerprint" in fingerprint_error.lower()

    confounded = adata.copy()
    confounded.obs["condition"] = confounded.obs["batch"].astype(str)
    blocked_card = scl.recommend_preprocess_policy(
        confounded,
        context,
        consumer="integration",
    )

    first_pca_gram = np.asarray(first.obsm["X_pca"]) @ np.asarray(first.obsm["X_pca"]).T
    second_pca_gram = np.asarray(second.obsm["X_pca"]) @ np.asarray(second.obsm["X_pca"]).T
    integrated_keys = sorted(
        key for key in map(str, first.obsm.keys()) if key not in {"X_pca", "X_umap"}
    )

    checks = {
        "review_is_read_only": _check(
            after_review_snapshot,
            input_snapshot,
            after_review_snapshot == input_snapshot,
        ),
        "apply_does_not_mutate_input": _check(
            after_apply_snapshot,
            input_snapshot,
            after_apply_snapshot == input_snapshot,
        ),
        "fingerprint_rejects_changed_counts": _check(
            fingerprint_error,
            "ValueError mentioning fingerprint",
            fingerprint_rejected,
        ),
        "counts_are_permanent_and_exact": _check(
            {
                "layer_present": "counts" in first.layers,
                "source_digest": _matrix_digest(input_counts),
                "result_digest": _matrix_digest(first.layers["counts"]),
            },
            "layers[counts] exactly preserves input UMI counts",
            "counts" in first.layers and _matrix_equal(first.layers["counts"], input_counts),
        ),
        "normalized_full_is_full_gene_space": _check(
            {
                "present": "normalized_full" in first.layers,
                "shape": list(first.layers["normalized_full"].shape),
            },
            [int(adata.n_obs), int(adata.n_vars)],
            "normalized_full" in first.layers
            and first.layers["normalized_full"].shape == adata.shape,
        ),
        "persistent_expression_spaces_remain_sparse": _check(
            {
                "counts": sparse.issparse(first.layers["counts"]),
                "normalized_full": sparse.issparse(first.layers["normalized_full"]),
                "X": sparse.issparse(first.X),
                "raw": first.raw is not None and sparse.issparse(first.raw.X),
            },
            "counts, normalized_full, X, and raw remain sparse for sparse input",
            sparse.issparse(first.layers["counts"])
            and sparse.issparse(first.layers["normalized_full"])
            and sparse.issparse(first.X)
            and first.raw is not None
            and sparse.issparse(first.raw.X),
        ),
        "raw_is_normalized_full_snapshot": _check(
            {
                "raw_present": first.raw is not None,
                "raw_n_vars": int(first.raw.n_vars) if first.raw is not None else None,
            },
            "raw contains normalized_full for every input gene",
            first.raw is not None
            and first.raw.n_vars == adata.n_vars
            and _matrix_equal(first.raw.X, first.layers["normalized_full"], atol=1e-7),
        ),
        "X_is_interpretation_space": _check(
            representation.get("marker_program_source"),
            "layers[normalized_full]",
            _matrix_equal(first.X, first.layers["normalized_full"], atol=1e-7)
            and representation.get("marker_program_source") == "layers[normalized_full]",
        ),
        "discovery_rep_and_feature_mask_exist": _check(
            {
                "X_pca": "X_pca" in first.obsm,
                "discovery_feature": "discovery_feature" in first.var,
                "n_discovery_features": int(first.var["discovery_feature"].sum()),
            },
            "non-empty discovery_feature mask and obsm[X_pca]",
            "X_pca" in first.obsm
            and "discovery_feature" in first.var
            and int(first.var["discovery_feature"].sum()) >= 2,
        ),
        "discovery_densification_is_bounded_and_audited": _check(
            resource_contract,
            "temporary centered dense matrix is feature-bounded and PCA-only",
            resource_contract["densification_occurred"] is True
            and resource_contract["scope"] == "temporary_discovery_feature_matrix"
            and resource_contract["temporary_shape"]
            == [int(first.n_obs), int(first.var["discovery_feature"].sum())]
            and resource_contract["temporary_shape"][1]
            <= int(first_card.policy.execution["n_top_genes"])
            and resource_contract["estimated_peak_bytes"] >= resource_contract["dense_matrix_bytes"]
            and resource_contract["persistent"] is False
            and resource_contract["consumer"] == "PCA_and_neighbor_graph_only"
            and resource_contract["expression_inference_eligible"] is False,
        ),
        "no_implicit_sparse_zero_center_densification": _check(
            implicit_densification_warnings,
            [],
            not implicit_densification_warnings,
        ),
        "bounded_dense_pca_preserves_legacy_geometry": _check(
            geometry,
            {
                "relative_pca_gram_error_lte": 1e-5,
                "mean_neighbor_retention_gte": 0.99,
            },
            geometry["relative_pca_gram_error"] <= 1e-5
            and geometry["mean_neighbor_retention"] >= 0.99,
        ),
        "formal_models_point_to_counts": _check(
            representation.get("formal_count_model_source"),
            "layers[counts]",
            representation.get("formal_count_model_source") == "layers[counts]",
        ),
        "unintegrated_baseline_is_fail_safe": _check(
            {
                "normalization_method": first_card.policy.normalization_method,
                "feature_selection_method": first_card.policy.feature_selection_method,
                "run_integration": first_card.policy.run_integration,
                "integrated_rep": representation.get("integrated_rep"),
                "unexpected_integrated_keys": integrated_keys,
            },
            "standard + scanpy + unintegrated PCA until Pareto evidence exists",
            first_card.policy.normalization_method == "standard"
            and first_card.policy.feature_selection_method == "scanpy"
            and first_card.policy.run_integration is False
            and representation.get("integrated_rep") == "not_selected"
            and not integrated_keys,
        ),
        "integration_consumer_requires_pareto_evidence": _check(
            {
                "consumer": integration_card.policy.consumer,
                "status": integration_card.status,
                "integration_status": integration_card.details["integration_review"]["status"],
                "run_integration": integration_card.policy.run_integration,
            },
            "REVIEW and unintegrated until a Pareto comparison supports integration",
            integration_card.policy.consumer == "integration"
            and integration_card.status == "REVIEW"
            and integration_card.details["integration_review"]["status"] == "REVIEW"
            and integration_card.policy.run_integration is False,
        ),
        "expression_inference_consumer_uses_counts": _check(
            {
                "consumer": inference_card.policy.consumer,
                "formal_count_model_source": inference_card.policy.layer_contract.get(
                    "formal_count_model_source"
                ),
                "run_integration": inference_card.policy.run_integration,
            },
            "expression_inference points to layers[counts] and does not auto-integrate",
            inference_card.policy.consumer == "expression_inference"
            and inference_card.policy.layer_contract.get("formal_count_model_source")
            == "layers[counts]"
            and inference_card.policy.run_integration is False,
        ),
        "integration_consumer_respects_confounding_block": _check(
            {
                "status": blocked_card.status,
                "integration_status": blocked_card.details["integration_review"]["status"],
                "run_integration": blocked_card.policy.run_integration,
            },
            "BLOCKED and integration not selected when batch equals condition",
            blocked_card.status == "BLOCKED"
            and blocked_card.details["integration_review"]["status"] == "BLOCKED"
            and blocked_card.policy.run_integration is False,
        ),
        "repeat_policy_agreement": _check(
            {
                "first": [
                    first_card.policy.normalization_method,
                    first_card.policy.feature_selection_method,
                    first_card.policy.run_integration,
                ],
                "second": [
                    second_card.policy.normalization_method,
                    second_card.policy.feature_selection_method,
                    second_card.policy.run_integration,
                ],
            },
            "same method decisions for the same input and context",
            (
                first_card.policy.normalization_method,
                first_card.policy.feature_selection_method,
                first_card.policy.run_integration,
            )
            == (
                second_card.policy.normalization_method,
                second_card.policy.feature_selection_method,
                second_card.policy.run_integration,
            ),
        ),
        "repeat_execution_agreement": _check(
            {
                "normalized_full_equal": _matrix_equal(
                    first.layers["normalized_full"],
                    second.layers["normalized_full"],
                    atol=1e-7,
                ),
                "pca_geometry_equal": bool(
                    np.allclose(first_pca_gram, second_pca_gram, atol=1e-5, rtol=1e-5)
                ),
            },
            "same normalized expression and PCA geometry",
            _matrix_equal(
                first.layers["normalized_full"],
                second.layers["normalized_full"],
                atol=1e-7,
            )
            and bool(np.allclose(first_pca_gram, second_pca_gram, atol=1e-5, rtol=1e-5)),
        ),
    }

    all_passed = all(item["passed"] for item in checks.values())
    report = {
        "schema_name": "sclucid_controlled_preprocess_contract_benchmark",
        "schema_version": "1.1.0",
        "status": "CONTRACT_PASS_NOT_PERFORMANCE" if all_passed else "CONTRACT_FAIL",
        "seed": seed,
        "public_api": ["recommend_preprocess_policy", "apply_preprocess_policy"],
        "fixture": {
            "type": "deterministic_synthetic_umi",
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "samples": int(adata.obs["sample"].nunique()),
            "purpose": "execution contract only",
        },
        "checks": checks,
        "summary": {
            "n_checks": len(checks),
            "n_passed": sum(item["passed"] for item in checks.values()),
            "n_failed": sum(not item["passed"] for item in checks.values()),
        },
        "resource_contract": {
            **resource_contract,
            "legacy_geometry_compatibility": geometry,
        },
        "scalability": {
            "status": "BLOCKED_REAL_SCALE_NOT_RUN",
            "reason": (
                "The deterministic fixture validates an auditable bounded allocation, "
                "not peak RSS or runtime on representative large datasets and hardware."
            ),
            "next_evidence": (
                "Measure runtime and peak RSS on registered reference-scale sparse datasets."
            ),
        },
        "claim_boundary": {
            "supported": [
                "The public preprocessing review/apply API obeyed the tested execution and representation contracts on the controlled fixture."
            ],
            "unsupported": [
                "Scientific superiority over alternative preprocessing methods.",
                "Generalization to held-out mixology, primary tissue, or tumor projects.",
                "Biological fidelity of normalization, feature selection, or integration.",
                "Scalability on representative real datasets or reference hardware.",
            ],
            "next_evidence": [
                "Held-out mixology performance benchmark.",
                "Leave-one-project-out primary tissue and tumor validation.",
            ],
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--seed", type=int, default=29)
    args = parser.parse_args()
    report = run(args.output, seed=args.seed)
    print(json.dumps(report["summary"], indent=2))
    print(f"status={report['status']}")
    print(f"artifact={args.output}")
    return 0 if report["status"] == "CONTRACT_PASS_NOT_PERFORMANCE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
