from __future__ import annotations

import json

from validation.preprocess.run_controlled_preprocess_contract_benchmark import run


def test_controlled_preprocess_contract_uses_public_api_and_calibrates_claim(tmp_path):
    output = tmp_path / "controlled_preprocess_contract_benchmark.json"

    report = run(output, seed=29)

    assert report["status"] == "CONTRACT_PASS_NOT_PERFORMANCE"
    assert report["public_api"] == [
        "recommend_preprocess_policy",
        "apply_preprocess_policy",
    ]
    assert report["summary"]["n_failed"] == 0
    assert all(item["passed"] for item in report["checks"].values())
    assert "Scientific superiority" in report["claim_boundary"]["unsupported"][0]
    assert json.loads(output.read_text())["status"] == "CONTRACT_PASS_NOT_PERFORMANCE"


def test_controlled_preprocess_contract_covers_fail_closed_and_four_spaces(tmp_path):
    report = run(tmp_path / "report.json", seed=29)
    checks = report["checks"]

    required = {
        "review_is_read_only",
        "apply_does_not_mutate_input",
        "fingerprint_rejects_changed_counts",
        "counts_are_permanent_and_exact",
        "normalized_full_is_full_gene_space",
        "persistent_expression_spaces_remain_sparse",
        "raw_is_normalized_full_snapshot",
        "discovery_rep_and_feature_mask_exist",
        "discovery_densification_is_bounded_and_audited",
        "no_implicit_sparse_zero_center_densification",
        "bounded_dense_pca_preserves_legacy_geometry",
        "formal_models_point_to_counts",
        "unintegrated_baseline_is_fail_safe",
        "integration_consumer_requires_pareto_evidence",
        "expression_inference_consumer_uses_counts",
        "integration_consumer_respects_confounding_block",
        "repeat_policy_agreement",
        "repeat_execution_agreement",
    }
    assert required <= set(checks)
    assert (
        checks["integration_consumer_respects_confounding_block"]["observed"]["status"] == "BLOCKED"
    )
    assert (
        checks["unintegrated_baseline_is_fail_safe"]["observed"]["integrated_rep"] == "not_selected"
    )
    assert report["resource_contract"]["densification_occurred"] is True
    assert report["resource_contract"]["persistent"] is False
    assert report["scalability"]["status"] == "BLOCKED_REAL_SCALE_NOT_RUN"
