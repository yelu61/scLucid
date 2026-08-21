from __future__ import annotations

import pytest

from validation.qc.run_real_project_qc_concordance import compare_membership


def test_compare_membership_keeps_review_separate_from_remove():
    result = compare_membership(
        {"a", "b", "c", "d"},
        {"a", "b"},
        {"c"},
        {"a", "d"},
    )

    assert result["n_historical_removed"] == 2
    assert result["historical_removed_current_remove"] == 1
    assert result["historical_removed_current_review"] == 1
    assert result["historical_removed_flagged_fraction"] == 1.0
    assert result["historical_kept_current_review"] == 1


def test_compare_membership_fails_closed_on_invalid_membership():
    with pytest.raises(ValueError, match="exact subset"):
        compare_membership({"a"}, {"a", "missing"}, set(), set())

    with pytest.raises(ValueError, match="both REMOVE and REVIEW"):
        compare_membership({"a"}, {"a"}, {"a"}, {"a"})
