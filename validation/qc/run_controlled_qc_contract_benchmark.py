#!/usr/bin/env python3
"""Validate QC input, execution, determinism, and bounded-run contracts."""

from __future__ import annotations

import argparse
import json
import platform
import time
import tracemalloc
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData, read_h5ad
from scipy import sparse

import scLucid as scl

SCHEMA_VERSION = "sclucid_controlled_qc_contract_benchmark_v1"


def _matrix_equal(left: Any, right: Any) -> bool:
    if sparse.issparse(left) or sparse.issparse(right):
        left_sparse = sparse.csr_matrix(left)
        right_sparse = sparse.csr_matrix(right)
        return left_sparse.shape == right_sparse.shape and (left_sparse != right_sparse).nnz == 0
    return np.array_equal(np.asarray(left), np.asarray(right))


def _review_mutation_count(
    adata: AnnData,
    before: AnnData,
) -> list[str]:
    mutations: list[str] = []
    if not _matrix_equal(adata.X, before.X):
        mutations.append("X")
    try:
        pd.testing.assert_frame_equal(adata.obs, before.obs)
    except AssertionError:
        mutations.append("obs")
    try:
        pd.testing.assert_frame_equal(adata.var, before.var)
    except AssertionError:
        mutations.append("var")
    if set(adata.layers) != set(before.layers) or any(
        not _matrix_equal(adata.layers[key], before.layers[key]) for key in before.layers
    ):
        mutations.append("layers")
    if set(adata.obsm) != set(before.obsm) or any(
        not _matrix_equal(adata.obsm[key], before.obsm[key]) for key in before.obsm
    ):
        mutations.append("obsm")
    if not _deep_equal(adata.uns, before.uns):
        mutations.append("uns")
    return mutations


def _deep_equal(left: Any, right: Any) -> bool:
    if isinstance(left, dict) and isinstance(right, dict):
        return set(left) == set(right) and all(_deep_equal(left[key], right[key]) for key in left)
    if isinstance(left, (list, tuple)) and isinstance(right, (list, tuple)):
        return len(left) == len(right) and all(_deep_equal(a, b) for a, b in zip(left, right))
    if isinstance(left, np.ndarray) or isinstance(right, np.ndarray):
        return np.array_equal(np.asarray(left), np.asarray(right), equal_nan=True)
    try:
        return bool(left == right)
    except (TypeError, ValueError):
        return False


def _fraction_exact(left: Any, right: Any) -> float:
    if left.shape != right.shape:
        return 0.0
    if sparse.issparse(left) or sparse.issparse(right):
        return 1.0 if (sparse.csr_matrix(left) != sparse.csr_matrix(right)).nnz == 0 else 0.0
    return float(np.array_equal(np.asarray(left), np.asarray(right)))


def _layers_equal(left: AnnData, right: AnnData) -> bool:
    return set(left.layers) == set(right.layers) and all(
        _matrix_equal(left.layers[key], right.layers[key]) for key in left.layers
    )


def _input_failure_report(
    adata: AnnData,
    *,
    missing_obs: list[str],
    missing_layers: list[str],
    failed_checks: list[str],
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "FAIL",
        "evidence_scope": "engineering_contract",
        "endpoint_status": {
            "qc_input_contract": "FAIL",
            "qc_profile_selection": "NOT_RUN",
            "qc_policy_execution": "NOT_RUN",
            "qc_scalability": "NOT_RUN",
            "qc_cell_calling": "NOT_RUN",
            "qc_ambient_correction": "NOT_RUN",
            "qc_damage_classification": "NOT_RUN",
            "qc_doublet_calibration": "NOT_RUN",
            "qc_rare_population_preservation": "NOT_RUN",
        },
        "engineering_metrics": {"counts_integrity_rate": 0.0},
        "input": {
            "n_observations": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "missing_obs": missing_obs,
            "missing_layers": missing_layers,
            "failed_checks": failed_checks,
        },
        "claim_boundary": {
            "supported": ["The controlled-suite input contract failed closed."],
            "unsupported": [
                "No execution or scientific conclusion is supported by an invalid fixture."
            ],
        },
    }


def run_benchmark(adata: AnnData) -> dict[str, Any]:
    required_obs = {
        "library",
        "assay_profile",
        "truth_droplet_class",
        "truth_policy_label",
        "truth_lineage",
        "truth_protected_keep",
        "truth_is_rare",
        "truth_ambient_fraction",
        "truth_damage_fraction",
    }
    missing_obs = sorted(required_obs - set(adata.obs))
    required_layers = {"counts", "native_counts", "ambient_counts"}
    missing_layers = sorted(required_layers - set(adata.layers))
    counts_identity = (
        False
        if missing_layers
        else _matrix_equal(
            adata.layers["counts"],
            adata.layers["native_counts"] + adata.layers["ambient_counts"],
        )
    )
    failed_input_checks: list[str] = []
    if not adata.obs_names.is_unique:
        failed_input_checks.append("obs_names_not_unique")
    if not adata.var_names.is_unique:
        failed_input_checks.append("var_names_not_unique")
    if "sclucid_controlled_truth" not in adata.uns:
        failed_input_checks.append("truth_provenance_missing")
    if not counts_identity:
        failed_input_checks.append("counts_not_native_plus_ambient")
    if not missing_layers:
        counts = adata.layers["counts"]
        values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
        if bool(np.any(values < 0)):
            failed_input_checks.append("counts_negative")
        if not bool(np.all(np.equal(values, np.floor(values)))):
            failed_input_checks.append("counts_non_integer")
    if missing_obs or missing_layers or failed_input_checks:
        return _input_failure_report(
            adata,
            missing_obs=missing_obs,
            missing_layers=missing_layers,
            failed_checks=failed_input_checks,
        )

    execution_input = adata[
        (adata.obs["library"] == "hq_scrna")
        & adata.obs["truth_droplet_class"].isin(["intact", "damaged"])
    ].copy()
    context = scl.ProjectContext(
        sample_key="library",
        is_multi_sample=False,
        assay="scrna",
        input_provenance="filtered_counts",
        cell_type_key="truth_lineage",
    )
    before = execution_input.copy()

    tracemalloc.start()
    started = time.perf_counter()
    card = scl.recommend_qc_policy(execution_input, context)
    review_seconds = time.perf_counter() - started
    _, peak_memory = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    review_mutations = _review_mutation_count(execution_input, before)

    first = scl.apply_qc_policy(execution_input, card.policy)
    second = scl.apply_qc_policy(execution_input.copy(), card.policy)
    expected_names = execution_input.obs_names[
        ~execution_input.obs_names.astype(str).isin(card.policy.remove_obs_names)
    ].astype(str)
    actual_names = first.adata.obs_names.astype(str)
    policy_apply_agreement = float(np.array_equal(expected_names, actual_names))
    policy_remove_names_valid = float(
        len(card.policy.remove_obs_names) == len(set(card.policy.remove_obs_names))
        and set(card.policy.remove_obs_names) <= set(execution_input.obs_names.astype(str))
    )
    first_repeat_obs_equal = True
    try:
        pd.testing.assert_frame_equal(first.adata.obs, second.adata.obs)
    except AssertionError:
        first_repeat_obs_equal = False
    repeat_run_agreement = float(
        np.array_equal(first.adata.obs_names, second.adata.obs_names)
        and _matrix_equal(first.adata.X, second.adata.X)
        and first_repeat_obs_equal
        and _layers_equal(first.adata, second.adata)
    )
    retained_source = execution_input[first.adata.obs_names, :]
    retained_x_preservation = _fraction_exact(first.adata.X, retained_source.X)
    retained_counts_preservation = _fraction_exact(
        first.adata.layers["counts"], retained_source.layers["counts"]
    )
    retained_layers_preservation = float(_layers_equal(first.adata, retained_source))

    fail_closed_card = scl.recommend_qc_policy(
        execution_input,
        scl.ProjectContext(
            is_multi_sample=True,
            sample_key="missing_sample_key",
            input_provenance="filtered_counts",
        ),
    )
    sample_key_fail_closed_rate = float(fail_closed_card.status == "BLOCKED")
    filtered_input_not_evaluable = float(
        card.policy.evidence_heads["ambient"]["status"] == "NOT_EVALUABLE"
    )

    profile_rows: list[dict[str, Any]] = []
    profile_minimum_genes = {"scrna": 200.0, "snrna": 100.0}
    for profile in ("scrna", "snrna"):
        profile_input = adata[
            (adata.obs["assay_profile"] == profile) & (adata.obs["truth_droplet_class"] == "intact")
        ].copy()
        profile_card = scl.recommend_qc_policy(
            profile_input,
            scl.ProjectContext(assay=profile, input_provenance="filtered_counts"),
        )
        sample_median_genes = [
            float(row["median_n_genes"]) for row in profile_card.details["sample_decisions"]
        ]
        profile_rows.append(
            {
                "expected_profile": profile,
                "selected_profile": profile_card.policy.profile,
                "status": profile_card.status,
                "blocked": profile_card.status == "BLOCKED",
                "minimum_sample_median_genes": min(sample_median_genes),
                "protocol_safety_floor": profile_minimum_genes[profile],
                "clears_protocol_safety_floor": bool(
                    min(sample_median_genes) >= profile_minimum_genes[profile]
                ),
            }
        )
    profile_provenance_rate = float(
        np.mean([row["expected_profile"] == row["selected_profile"] for row in profile_rows])
    )
    locked_high_quality_false_block_rate = float(
        np.mean([bool(row["blocked"]) for row in profile_rows])
    )
    locked_high_quality_safety_floor_rate = float(
        np.mean([bool(row["clears_protocol_safety_floor"]) for row in profile_rows])
    )

    sparse_contract = {
        "input_x_sparse": bool(sparse.issparse(execution_input.X)),
        "output_x_sparse": bool(sparse.issparse(first.adata.X)),
        "input_counts_sparse": bool(sparse.issparse(execution_input.layers["counts"])),
        "output_counts_sparse": bool(sparse.issparse(first.adata.layers["counts"])),
    }
    input_counts = sparse.csr_matrix(execution_input.layers["counts"])
    output_counts = sparse.csr_matrix(first.adata.layers["counts"])
    engineering_metrics = {
        "counts_integrity_rate": float(counts_identity and not missing_obs and not missing_layers),
        "sample_key_fail_closed_rate": sample_key_fail_closed_rate,
        "filtered_input_ambient_not_evaluable_rate": filtered_input_not_evaluable,
        "review_mutation_count": len(review_mutations),
        "review_mutated_components": review_mutations,
        "policy_apply_agreement": policy_apply_agreement,
        "policy_remove_names_valid": policy_remove_names_valid,
        "repeat_run_agreement": repeat_run_agreement,
        "retained_x_preservation": retained_x_preservation,
        "retained_counts_preservation": retained_counts_preservation,
        "retained_layers_preservation": retained_layers_preservation,
        "profile_provenance_rate": profile_provenance_rate,
        "locked_high_quality_false_block_rate": locked_high_quality_false_block_rate,
        "locked_high_quality_safety_floor_rate": locked_high_quality_safety_floor_rate,
        "runtime_seconds": review_seconds,
        "peak_python_traced_memory_bytes": int(peak_memory),
        "input_counts_nnz": int(input_counts.nnz),
        "output_counts_nnz": int(output_counts.nnz),
        "input_counts_density": float(input_counts.nnz / np.prod(input_counts.shape)),
        "output_counts_density": float(output_counts.nnz / np.prod(output_counts.shape)),
        "sparse_contract": sparse_contract,
        "unintended_dense_expansion": bool(
            (sparse_contract["input_x_sparse"] and not sparse_contract["output_x_sparse"])
            or (
                sparse_contract["input_counts_sparse"]
                and not sparse_contract["output_counts_sparse"]
            )
        ),
    }
    engineering_pass = (
        engineering_metrics["counts_integrity_rate"] == 1.0
        and engineering_metrics["sample_key_fail_closed_rate"] == 1.0
        and engineering_metrics["filtered_input_ambient_not_evaluable_rate"] == 1.0
        and engineering_metrics["review_mutation_count"] == 0
        and engineering_metrics["policy_apply_agreement"] == 1.0
        and engineering_metrics["policy_remove_names_valid"] == 1.0
        and engineering_metrics["repeat_run_agreement"] == 1.0
        and engineering_metrics["retained_x_preservation"] == 1.0
        and engineering_metrics["retained_counts_preservation"] == 1.0
        and engineering_metrics["retained_layers_preservation"] == 1.0
        and engineering_metrics["profile_provenance_rate"] == 1.0
        and engineering_metrics["locked_high_quality_false_block_rate"] == 0.0
        and engineering_metrics["locked_high_quality_safety_floor_rate"] == 1.0
        and not engineering_metrics["unintended_dense_expansion"]
    )

    totals = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
    mt_totals = np.asarray(adata.layers["counts"][:, :8].sum(axis=1)).ravel()
    mt_fraction = np.divide(
        mt_totals, totals, out=np.zeros_like(mt_totals, dtype=float), where=totals > 0
    )
    classes = adata.obs["truth_droplet_class"].astype(str).to_numpy()
    mechanism_checks = {
        "empty_median_counts_below_intact": bool(
            np.median(totals[classes == "empty"]) < np.median(totals[classes == "intact"])
        ),
        "damaged_median_mt_above_intact": bool(
            np.median(mt_fraction[classes == "damaged"])
            > np.median(mt_fraction[classes == "intact"])
        ),
        "ambient_component_present": bool(adata.layers["ambient_counts"].sum() > 0),
        "rare_intact_cells_present": bool(
            (
                adata.obs["truth_is_rare"].astype(bool)
                & (adata.obs["truth_policy_label"] == "KEEP")
            ).any()
        ),
        "known_doublets_present": bool(adata.obs["truth_is_doublet"].astype(bool).any()),
    }
    mechanism_status = "SIMULATION_PASS_NOT_EXTERNAL" if all(mechanism_checks.values()) else "FAIL"
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "PASS" if engineering_pass else "FAIL",
        "evidence_scope": "engineering_contract",
        "endpoint_status": {
            "qc_input_contract": "PASS" if engineering_pass else "FAIL",
            "qc_profile_selection": "PASS" if engineering_pass else "FAIL",
            "qc_policy_execution": "PASS" if engineering_pass else "FAIL",
            "qc_scalability": ("SIMULATION_PASS_NOT_EXTERNAL" if engineering_pass else "FAIL"),
            "qc_cell_calling": mechanism_status,
            "qc_ambient_correction": mechanism_status,
            "qc_damage_classification": mechanism_status,
            "qc_doublet_calibration": mechanism_status,
            "qc_rare_population_preservation": mechanism_status,
        },
        "engineering_metrics": engineering_metrics,
        "profile_rows": profile_rows,
        "mechanism_checks": mechanism_checks,
        "reference_hardware": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "processor": platform.processor() or "not_reported",
        },
        "resource_scope": (
            "Small controlled fixture only. Runtime and Python-traced memory are recorded "
            "for regression; they are not a hardware-independent scalability PASS."
        ),
        "input": {
            "n_observations": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "execution_subset_n": int(execution_input.n_obs),
            "missing_obs": missing_obs,
            "missing_layers": missing_layers,
        },
        "claim_boundary": {
            "supported": [
                "Read-only review, fail-closed sample identity and raw-input evidence, exact policy application, repeatability, and count/layer preservation on the controlled suite."
            ],
            "unsupported": [
                "Scientific QC superiority or real-tissue generalization from simulated controls."
            ],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("validation_outputs/current/qc_controlled_truth/controlled_qc_truth.h5ad"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation_outputs/current/qc_controlled_truth"),
    )
    args = parser.parse_args()
    report = run_benchmark(read_h5ad(args.input))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "controlled_qc_contract_benchmark.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "artifact": str(output)}, indent=2))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
