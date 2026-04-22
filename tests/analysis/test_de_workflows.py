"""Tests for high-level DE and characterization workflows."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.analysis.differential_expression.de_workflows import characterize_clusters


def _make_clustered_adata() -> AnnData:
    """Create minimal clustered AnnData for characterization workflow tests."""
    adata = AnnData(np.ones((6, 4), dtype=float))
    adata.obs_names = [f"cell_{i}" for i in range(6)]
    adata.var_names = [f"gene_{i}" for i in range(4)]
    adata.obs["leiden_clusters"] = pd.Categorical(["0", "0", "0", "1", "1", "1"])
    return adata


@pytest.mark.unit
def test_characterize_clusters_stores_review_tables_and_exports(tmp_path, monkeypatch):
    """characterize_clusters writes notebook-facing tables and markdown sidecars."""
    adata = _make_clustered_adata()

    fake_markers = pd.DataFrame(
        {
            "group": ["0", "0", "1", "1"],
            "names": ["CD3D", "LTB", "NKG7", "GNLY"],
            "logfoldchanges": [2.5, 1.8, 3.1, 2.2],
            "scores": [10.0, 8.0, 12.0, 9.0],
            "pvals_adj": [0.001, 0.005, 0.002, 0.01],
        }
    )
    fake_enrichment = {
        "0": {
            "ora": pd.DataFrame(
                {
                    "Term": ["T cell activation", "TCR signaling"],
                    "Adjusted P-value": [0.001, 0.02],
                }
            )
        },
        "1": {
            "ora": pd.DataFrame(
                {
                    "Term": ["NK cell mediated cytotoxicity", "Cytolysis"],
                    "Adjusted P-value": [0.003, 0.04],
                }
            )
        },
    }

    def fake_find_markers(adata_obj, config):
        adata_obj.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
            "de", {}
        )["rank_genes_groups_df"] = fake_markers.copy()
        return fake_markers.copy()

    monkeypatch.setattr(
        "scLucid.analysis.differential_expression.de_workflows.find_markers",
        fake_find_markers,
    )
    monkeypatch.setattr(
        "scLucid.analysis.differential_expression.de_workflows.run_enrichment",
        lambda adata_obj, groupby, config: fake_enrichment,
    )

    result = characterize_clusters(
        adata,
        groupby="leiden_clusters",
        save_path=tmp_path,
        n_top_markers=2,
        n_top_terms=2,
    )

    review = result.uns["cluster_characterization"]
    assert "summary_table" in review
    assert "top_markers" in review
    assert "enrichment_summary" in review
    assert "markdown_summary" in review
    assert {"cluster", "n_cells", "top_markers", "top_pathways"}.issubset(
        review["summary_table"].columns
    )
    assert review["top_markers"]["cluster"].astype(str).isin(["0", "1"]).all()
    assert review["enrichment_summary"]["term"].str.len().gt(0).all()

    export_paths = review["export_paths"]
    for path in export_paths.values():
        assert Path(path).exists(), f"Expected exported sidecar to exist: {path}"


@pytest.mark.unit
def test_characterize_clusters_accepts_base_csv_path(tmp_path, monkeypatch):
    """A CSV save path is used as the summary path and sidecars share its stem."""
    adata = _make_clustered_adata()

    fake_markers = pd.DataFrame(
        {
            "group": ["0"],
            "names": ["CD3D"],
            "logfoldchanges": [2.0],
            "scores": [9.0],
            "pvals_adj": [0.001],
        }
    )
    fake_enrichment = {"0": {"ora": pd.DataFrame({"Term": ["T cell activation"]})}}

    def fake_find_markers(adata_obj, config):
        adata_obj.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault(
            "de", {}
        )["rank_genes_groups_df"] = fake_markers.copy()
        return fake_markers.copy()

    monkeypatch.setattr(
        "scLucid.analysis.differential_expression.de_workflows.find_markers",
        fake_find_markers,
    )
    monkeypatch.setattr(
        "scLucid.analysis.differential_expression.de_workflows.run_enrichment",
        lambda adata_obj, groupby, config: fake_enrichment,
    )

    output = tmp_path / "review_bundle.csv"
    result = characterize_clusters(
        adata,
        groupby="leiden_clusters",
        save_path=output,
    )

    assert result.uns["cluster_characterization"]["export_paths"]["summary"] == str(output)
    assert output.exists()

