"""Tests for marker manager helpers and built-in resource composition."""

from scLucid.utils import get_marker_manager


def test_get_marker_manager_loads_scoped_state_resources():
    """Requested state markers should remain discoverable after resource reorganization."""
    mgr = get_marker_manager(
        species="human",
        states=["T cell exhaustion-like", "Stress-high"],
    )

    assert "T cell exhaustion-like" in mgr.CELLS
    assert "Stress-high" in mgr.CELLS
    assert mgr.CELLS["T cell exhaustion-like"].metadata["scope"] == "lineage_restricted"
    assert "T cells" in mgr.CELLS["T cell exhaustion-like"].metadata["applies_to"]
    assert mgr.CELLS["Stress-high"].metadata["scope"] == "all"
