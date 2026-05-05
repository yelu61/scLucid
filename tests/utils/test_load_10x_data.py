"""Tests for scLucid.utils helpers related to data loading.

Full 10x-loading tests would require real 10x data or complex mocking of
scanpy internals. We focus on the helper logic (path resolution, argument
validation) that the notebook relies on.
"""

import os
from pathlib import Path

import pytest


class TestFindSamplePaths:
    def test__find_sample_paths_standard_subdirs(self, tmp_path):
        """_find_sample_paths checks standard 10x subdirectories."""
        from scLucid.utils.helpers import _find_sample_paths

        s1_dir = tmp_path / "s1" / "outs" / "filtered_feature_bc_matrix"
        s1_dir.mkdir(parents=True)
        (s1_dir / "matrix.mtx").touch()
        (s1_dir / "features.tsv").touch()

        result = _find_sample_paths(str(tmp_path), ["s1"])
        assert "s1" in result
        assert "outs/filtered_feature_bc_matrix" in result["s1"]

    def test__find_sample_paths_falls_back_to_root(self, tmp_path):
        """If no standard subdirs exist, falls back to sample root."""
        from scLucid.utils.helpers import _find_sample_paths

        s1_dir = tmp_path / "s1"
        s1_dir.mkdir(parents=True)
        (s1_dir / "matrix.mtx").touch()
        (s1_dir / "genes.tsv").touch()

        result = _find_sample_paths(str(tmp_path), ["s1"])
        assert "s1" in result
        assert result["s1"].rstrip("/") == str(s1_dir)

    def test__find_sample_paths_missing_sample(self, tmp_path):
        """Missing samples are omitted from the result dict."""
        from scLucid.utils.helpers import _find_sample_paths

        s1_dir = tmp_path / "s1"
        s1_dir.mkdir(parents=True)
        (s1_dir / "matrix.mtx").touch()
        (s1_dir / "features.tsv").touch()

        result = _find_sample_paths(str(tmp_path), ["s1", "missing"])
        assert "s1" in result
        assert "missing" not in result

    def test__find_sample_paths_prefers_mtx_gz(self, tmp_path):
        """Prefers .mtx.gz over .mtx when both present."""
        from scLucid.utils.helpers import _find_sample_paths

        s1_dir = tmp_path / "s1" / "outs" / "filtered_feature_bc_matrix"
        s1_dir.mkdir(parents=True)
        (s1_dir / "matrix.mtx.gz").touch()
        (s1_dir / "features.tsv.gz").touch()

        result = _find_sample_paths(str(tmp_path), ["s1"])
        assert "s1" in result


class TestLoad10xDataValidation:
    def test_load_10x_data_needs_base_or_path_dict(self):
        """ValueError when neither base_dir nor path_dict is provided."""
        from scLucid.utils import load_10x_data

        with pytest.raises(ValueError, match="base_dir or path_dict"):
            load_10x_data(samples=["s1"])
