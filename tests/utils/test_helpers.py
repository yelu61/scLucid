"""Tests for scLucid.utils.helpers utility functions."""


import pytest

from scLucid.utils.helpers import build_metadata_dicts


class TestBuildMetadataDicts:
    def test_builds_group_and_batch_dicts(self):
        samples = ["S1", "S2", "S3"]
        group_dict = {"S1": "tumor", "S2": "normal", "S3": "tumor", "S4": "other"}
        batch_dict = {"S1": "A", "S2": "A", "S3": "B", "S4": "C"}

        result = build_metadata_dicts(
            samples, group_dict=group_dict, batch_dict=batch_dict
        )

        assert result == {
            "group": {"S1": "tumor", "S2": "normal", "S3": "tumor"},
            "batch": {"S1": "A", "S2": "A", "S3": "B"},
        }

    def test_empty_inputs_return_empty(self):
        samples = ["S1", "S2"]
        result = build_metadata_dicts(samples)
        assert result == {}

    def test_custom_keys(self):
        samples = ["S1", "S2"]
        group_dict = {"S1": "ctrl", "S2": "treat"}
        result = build_metadata_dicts(
            samples, group_dict=group_dict, group_key="condition"
        )
        assert result == {"condition": {"S1": "ctrl", "S2": "treat"}}

    def test_extra_dicts_and_default_value(self):
        samples = ["S1", "S2", "S3"]
        result = build_metadata_dicts(
            samples,
            extra_dicts={
                "patient": {"S1": "P1", "S2": "P2"},
                "timepoint": {"S1": "pre", "S3": "post"},
            },
            default_value="unknown",
        )

        assert result["patient"] == {"S1": "P1", "S2": "P2", "S3": "unknown"}
        assert result["timepoint"] == {"S1": "pre", "S2": "unknown", "S3": "post"}

    def test_strict_missing_sample_raises(self):
        with pytest.raises(KeyError, match="missing values"):
            build_metadata_dicts(
                ["S1", "S2"],
                group_dict={"S1": "ctrl"},
                strict=True,
            )
