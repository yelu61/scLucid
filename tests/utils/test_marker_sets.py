"""Tests for marker-set filtering utilities."""

import pytest

from scLucid.utils import flatten_marker_dict, filter_marker_dict


@pytest.mark.unit
def test_flatten_marker_dict_handles_nested_sets():
    markers = {
        "lymphoid": {
            "T": ["CD3D", "CD3E"],
            "NK": ("NKG7", "GNLY"),
        },
        "B": {"core": {"markers": ["MS4A1", "CD79A"]}},
    }

    flat = flatten_marker_dict(markers)

    assert flat["lymphoid.T"] == ["CD3D", "CD3E"]
    assert flat["lymphoid.NK"] == ["NKG7", "GNLY"]
    assert flat["B.core.markers"] == ["MS4A1", "CD79A"]


@pytest.mark.unit
def test_filter_marker_dict_returns_missing_and_drops_empty():
    markers = {
        "T": ["CD3D", "CD3E", "TRAC"],
        "NK": ["NKG7", "GNLY"],
        "mixed": {"cyto": ["PRF1", "BADGENE"]},
    }

    filtered, missing = filter_marker_dict(
        markers,
        ["cd3d", "TRAC", "NKG7", "PRF1"],
        return_missing=True,
    )

    assert filtered == {
        "T": ["CD3D", "TRAC"],
        "NK": ["NKG7"],
        "mixed": {"cyto": ["PRF1"]},
    }
    assert missing["T"] == ["CD3E"]
    assert missing["NK"] == ["GNLY"]
    assert missing["mixed.cyto"] == ["BADGENE"]


@pytest.mark.unit
def test_filter_marker_dict_keeps_empty_when_requested():
    markers = {"T": ["CD3D"], "B": ["MS4A1"]}
    filtered = filter_marker_dict(markers, ["CD3D"], drop_empty=False)
    assert filtered == {"T": ["CD3D"], "B": []}
