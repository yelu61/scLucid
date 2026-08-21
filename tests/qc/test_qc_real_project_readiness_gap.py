from __future__ import annotations

from validation.qc.build_qc_real_project_readiness_gap import summarize_bindings


def test_summarize_bindings_counts_exact_portfolio_without_scoring():
    heads = [
        {
            "endpoint_id": "a",
            "status": "PASS",
            "dataset_statuses": {"d1": "PASS", "d2": "PASS_BASELINE"},
        },
        {
            "endpoint_id": "b",
            "status": "BLOCKED",
            "dataset_statuses": {"d3": "FAIL", "d4": "NOT_RUN"},
        },
    ]

    result = summarize_bindings(heads)

    assert result["head_count"] == 2
    assert result["passed_head_count"] == 1
    assert result["binding_count"] == 4
    assert result["passing_binding_count"] == 2
    assert result["status_counts"] == {
        "FAIL": 1,
        "NOT_RUN": 1,
        "PASS": 1,
        "PASS_BASELINE": 1,
    }
    assert "not an aggregate quality score" in result["interpretation"]


def test_summarize_bindings_can_limit_to_a_use_specific_head_set():
    heads = [
        {"endpoint_id": "filtered", "status": "PASS", "dataset_statuses": {"d1": "PASS"}},
        {"endpoint_id": "raw_only", "status": "BLOCKED", "dataset_statuses": {"d2": "FAIL"}},
    ]

    result = summarize_bindings(heads, endpoint_ids={"filtered"})

    assert result["head_count"] == 1
    assert result["passed_head_count"] == 1
    assert result["binding_count"] == 1
