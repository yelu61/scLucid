from __future__ import annotations

import pandas as pd
from scipy import sparse

from validation.qc.generate_controlled_qc_truth_suite import build_controlled_truth
from validation.qc.run_controlled_qc_contract_benchmark import run_benchmark


def test_controlled_qc_truth_is_deterministic_and_has_explicit_mechanisms():
    first = build_controlled_truth(seed=17, n_genes=120)
    second = build_controlled_truth(seed=17, n_genes=120)

    assert first.shape == second.shape
    assert (sparse.csr_matrix(first.X) != sparse.csr_matrix(second.X)).nnz == 0
    pd.testing.assert_frame_equal(first.obs, second.obs)
    pd.testing.assert_frame_equal(first.var, second.var)
    for layer in ("counts", "native_counts", "ambient_counts"):
        assert (
            sparse.csr_matrix(first.layers[layer]) != sparse.csr_matrix(second.layers[layer])
        ).nnz == 0
    assert set(first.obs["truth_droplet_class"]) == {
        "intact",
        "damaged",
        "empty",
        "doublet",
    }
    assert {"KEEP", "REMOVE", "UNCERTAIN"} <= set(first.obs["truth_policy_label"])
    assert first.obs.loc[first.obs["truth_low_rna"], "truth_policy_label"].eq("KEEP").all()
    assert first.obs.loc[first.obs["truth_protected_keep"], "truth_policy_label"].eq("KEEP").all()
    assert (
        first.obs.loc[first.obs["truth_droplet_class"] == "empty", "truth_policy_label"]
        .eq("REMOVE")
        .all()
    )
    assert (
        first.obs.loc[first.obs["truth_droplet_class"] == "damaged", "truth_policy_label"]
        .eq("REMOVE")
        .all()
    )
    assert first.obs.loc[first.obs["truth_droplet_class"] == "doublet", "truth_is_doublet"].all()
    assert (
        sparse.csr_matrix(first.layers["counts"])
        != sparse.csr_matrix(first.layers["native_counts"] + first.layers["ambient_counts"])
    ).nnz == 0


def test_controlled_qc_contract_pass_is_scoped_to_engineering():
    report = run_benchmark(build_controlled_truth(seed=23, n_genes=240))

    assert report["status"] == "PASS"
    assert report["endpoint_status"]["qc_input_contract"] == "PASS"
    assert report["endpoint_status"]["qc_profile_selection"] == "PASS"
    assert report["endpoint_status"]["qc_policy_execution"] == "PASS"
    assert report["endpoint_status"]["qc_scalability"] == "SIMULATION_PASS_NOT_EXTERNAL"
    for endpoint in (
        "qc_cell_calling",
        "qc_ambient_correction",
        "qc_damage_classification",
        "qc_doublet_calibration",
        "qc_rare_population_preservation",
    ):
        assert report["endpoint_status"][endpoint] == "SIMULATION_PASS_NOT_EXTERNAL"
    assert report["engineering_metrics"]["review_mutation_count"] == 0
    assert report["engineering_metrics"]["sample_key_fail_closed_rate"] == 1.0
    assert report["engineering_metrics"]["filtered_input_ambient_not_evaluable_rate"] == 1.0
    assert report["engineering_metrics"]["policy_apply_agreement"] == 1.0
    assert report["engineering_metrics"]["repeat_run_agreement"] == 1.0
    assert report["engineering_metrics"]["retained_counts_preservation"] == 1.0
    assert report["engineering_metrics"]["retained_layers_preservation"] == 1.0
    assert report["engineering_metrics"]["unintended_dense_expansion"] is False
    assert report["engineering_metrics"]["locked_high_quality_safety_floor_rate"] == 1.0
    assert all(row["clears_protocol_safety_floor"] for row in report["profile_rows"])
    assert "Scientific QC superiority" in report["claim_boundary"]["unsupported"][0]


def test_controlled_qc_contract_fails_closed_on_invalid_counts_identity():
    adata = build_controlled_truth(seed=29, n_genes=120)
    adata.layers["counts"][0, 0] += 1

    report = run_benchmark(adata)

    assert report["status"] == "FAIL"
    assert report["endpoint_status"]["qc_input_contract"] == "FAIL"
    assert report["endpoint_status"]["qc_policy_execution"] == "NOT_RUN"
    assert "counts_not_native_plus_ambient" in report["input"]["failed_checks"]
