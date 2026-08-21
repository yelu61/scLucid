from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from validation.qc_preprocess.locked_acceptance import (
    evaluate_preprocess_policy_acceptance,
    evaluate_qc_policy_acceptance,
    evaluate_real_project_ux_acceptance,
)


def test_qc_acceptance_excludes_uncertain_and_uses_grouped_bootstrap():
    names = [f"c{i}" for i in range(12)]
    truth = pd.Series(
        ["KEEP"] * 6 + ["REMOVE"] * 4 + ["UNCERTAIN"] * 2,
        index=names,
    )
    groups = pd.Series(["S1"] * 3 + ["S2"] * 3 + ["S1"] * 2 + ["S2"] * 4, index=names)
    policy = SimpleNamespace(
        remove_obs_names=names[6:10],
        candidate_policies=[
            {"name": "expert_global", "flagged_obs_names": ["c6", "c8"]},
            {"name": "per_sample_mad", "flagged_obs_names": ["c6", "c8"]},
        ],
    )

    result = evaluate_qc_policy_acceptance(
        policy,
        truth,
        groups,
        min_absolute_recall_gain=0.05,
        n_bootstrap=100,
    )

    assert result["status"] == "PASS"
    assert result["uncertain_excluded"] == 2
    assert all(row["grouped_bootstrap"]["status"] == "EVALUATED" for row in result["comparisons"])


def test_preprocess_acceptance_requires_baseline_or_pareto_gain():
    results = pd.DataFrame(
        [
            {
                "dataset": "heldout",
                "candidate": "standard_unintegrated",
                "selected": False,
                "preregistered_task_utility": 0.80,
                "biology_loss": 0.00,
            },
            {
                "dataset": "heldout",
                "candidate": "complex",
                "selected": True,
                "preregistered_task_utility": 0.81,
                "biology_loss": 0.01,
            },
        ]
    )
    result = evaluate_preprocess_policy_acceptance(results)
    assert result["status"] == "FAIL"
    assert result["datasets"][0]["simple_fallback_or_pareto_gain"] is False


def test_real_project_ux_acceptance_fails_closed_and_passes_complete_records():
    incomplete = pd.DataFrame({"project": ["P1"]})
    blocked = evaluate_real_project_ux_acceptance(
        incomplete,
        expected_projects=["P1", "P2"],
    )
    assert blocked["status"] == "BLOCKED"

    blank = pd.DataFrame(
        [
            {
                "project": project,
                "legacy_config_fields": "",
                "current_config_fields": "",
                "manual_predicted_doublet_deletion": "",
                "manual_review_summary_edit": "",
                "schema_bypass": "",
                "project_specific_patch_count": "",
                "run_evidence_status": "",
            }
            for project in ("P1", "P2")
        ]
    )
    blank_result = evaluate_real_project_ux_acceptance(
        blank,
        expected_projects=["P1", "P2"],
    )
    assert blank_result["status"] == "BLOCKED"

    records = pd.DataFrame(
        [
            {
                "project": project,
                "legacy_config_fields": 20,
                "current_config_fields": 5,
                "manual_predicted_doublet_deletion": False,
                "manual_review_summary_edit": False,
                "schema_bypass": False,
                "project_specific_patch_count": 0,
                "run_evidence_status": "REVIEW",
            }
            for project in ("P1", "P2")
        ]
    )
    passed = evaluate_real_project_ux_acceptance(
        records,
        expected_projects=["P1", "P2"],
    )
    assert passed["status"] == "PASS"
    assert all(row["config_field_reduction"] == 0.75 for row in passed["projects"])
