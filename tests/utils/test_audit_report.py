"""Tests for the HTML audit report exporter."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scLucid as scl
from scLucid.utils.audit_report import export_audit_report
from scLucid.utils.contracts import (
    REVIEW_SUMMARY_REQUIRED_KEYS,
    Modules,
    SCHEMA_VERSION,
    SCLUCID_ROOT,
    UnsKeys,
)


@pytest.fixture
def minimal_adata():
    """AnnData with no scLucid namespace at all."""
    return AnnData(np.random.poisson(5, size=(20, 50)).astype(float))


@pytest.fixture
def sclucid_adata():
    """AnnData with a populated scLucid namespace across QC / preprocess / analysis."""
    adata = AnnData(np.random.poisson(5, size=(40, 80)).astype(float))
    adata.layers["counts"] = adata.X.copy()
    adata.layers["normalized"] = adata.X.copy()
    adata.obsm["X_pca"] = np.random.randn(40, 5)
    adata.obs["leiden_clusters"] = pd.Categorical(["0"] * 20 + ["1"] * 20)
    adata.obs["cell_type_auto"] = pd.Categorical(["A"] * 25 + ["B"] * 15)

    def _review(module: str, steps: list[str], warnings_: list[str]) -> dict:
        return {
            "schema_version": SCHEMA_VERSION,
            "module": module,
            "workflow_name": f"standard_{module}",
            "steps_executed": steps,
            "data_shape": {"n_cells": adata.n_obs, "n_genes": adata.n_vars},
            "generated_at": "2026-05-13T10:00:00",
            "warnings": warnings_,
            "config": {"foo": "bar", "threshold": 0.05},
            "contract": {"valid": True, "schema_version": SCHEMA_VERSION},
            "config_lineage": {
                "global": {},
                "inherited": {"species": "human"},
                "stage": {},
                "effective": {"foo": "bar"},
            },
            "artifacts": {},
        }

    adata.uns[SCLUCID_ROOT] = {
        Modules.QC: {
            UnsKeys.WORKFLOW_CONFIG: {"min_genes": 200, "min_cells": 3},
            UnsKeys.STEPS_EXECUTED: ["metrics", "filter"],
            UnsKeys.REVIEW_SUMMARY: _review(
                "qc", ["metrics", "filter"], ["mito gene set inferred"]
            ),
            UnsKeys.CONTRACT: {"valid": True},
        },
        Modules.PREPROCESS: {
            UnsKeys.WORKFLOW_CONFIG: {"hvg": {"n_top_genes": 2000}},
            UnsKeys.STEPS_EXECUTED: ["normalize", "hvg", "scale", "pca"],
            UnsKeys.REVIEW_SUMMARY: _review(
                "preprocess", ["normalize", "hvg", "scale", "pca"], []
            ),
        },
        Modules.ANALYSIS: {
            UnsKeys.STEPS_EXECUTED: ["neighbors", "umap", "cluster", "annotate"],
            UnsKeys.REVIEW_SUMMARY: _review(
                "analysis",
                ["neighbors", "umap", "cluster", "annotate"],
                ["3 clusters had <10 cells; merged"],
            ),
            UnsKeys.ERRORS: [],
        },
        UnsKeys.PIPELINE_CONTEXT: {
            "species": "human",
            "tissue_type": "pbmc_or_blood",
            "batch_key": None,
        },
        UnsKeys.ANALYSIS_CONTEXT: {
            "species": "human",
            "tissue": "PBMC",
            "dataset_type": "pbmc_or_blood",
            "cancer_type": None,
        },
    }
    return adata


class TestExportAuditReportBasics:
    def test_writes_file(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "report.html"
        returned = export_audit_report(sclucid_adata, out_path)
        assert returned == out_path.resolve()
        assert out_path.exists()
        assert out_path.stat().st_size > 1000  # Non-trivial content

    def test_creates_parent_dirs(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "deep" / "nested" / "report.html"
        export_audit_report(sclucid_adata, out_path)
        assert out_path.exists()

    def test_returns_resolved_absolute_path(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "report.html"
        returned = export_audit_report(sclucid_adata, str(out_path))
        assert returned.is_absolute()


class TestRenderedContent:
    def test_includes_all_modules(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "report.html"
        export_audit_report(sclucid_adata, out_path)
        html = out_path.read_text()
        # Each module section's title is its capitalized name.
        for module in ("Qc", "Preprocess", "Analysis"):
            assert f">{module}<" in html, f"Missing section for {module}"

    def test_includes_dataset_shape_in_header(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "report.html"
        export_audit_report(sclucid_adata, out_path)
        html = out_path.read_text()
        assert "40" in html and "80" in html  # n_obs and n_vars

    def test_includes_pipeline_context_panel(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "report.html"
        export_audit_report(sclucid_adata, out_path)
        html = out_path.read_text()
        assert "Pipeline Context" in html
        assert "human" in html

    def test_includes_warnings_when_present(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "report.html"
        export_audit_report(sclucid_adata, out_path)
        html = out_path.read_text()
        assert "warning" in html.lower()
        assert "mito gene set inferred" in html

    def test_contract_badge_renders(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "report.html"
        export_audit_report(sclucid_adata, out_path)
        html = out_path.read_text()
        assert "contract" in html.lower()

    def test_escapes_html_special_chars(self, sclucid_adata, tmp_path):
        # Inject a value that would break HTML if not escaped.
        sclucid_adata.uns[SCLUCID_ROOT][Modules.QC][UnsKeys.REVIEW_SUMMARY][
            "warnings"
        ] = ["bad <script>alert(1)</script>"]
        out_path = tmp_path / "report.html"
        export_audit_report(sclucid_adata, out_path)
        html = out_path.read_text()
        assert "<script>alert(1)</script>" not in html
        assert "&lt;script&gt;" in html


class TestEmptyOrSparseInputs:
    def test_empty_namespace_does_not_crash(self, minimal_adata, tmp_path):
        out_path = tmp_path / "report.html"
        returned = export_audit_report(minimal_adata, out_path)
        assert returned.exists()
        html = out_path.read_text()
        assert "No module results were recorded" in html

    def test_non_dict_namespace_is_handled(self, minimal_adata, tmp_path):
        # Misconfigured namespace (e.g. dataframe instead of dict)
        minimal_adata.uns[SCLUCID_ROOT] = "not_a_dict"
        out_path = tmp_path / "report.html"
        returned = export_audit_report(minimal_adata, out_path)
        assert returned.exists()

    def test_partial_namespace_renders_available_modules(self, minimal_adata, tmp_path):
        minimal_adata.uns[SCLUCID_ROOT] = {
            Modules.QC: {
                UnsKeys.STEPS_EXECUTED: ["metrics"],
                UnsKeys.REVIEW_SUMMARY: {
                    key: f"value_for_{key}" for key in REVIEW_SUMMARY_REQUIRED_KEYS
                },
            }
        }
        out_path = tmp_path / "report.html"
        export_audit_report(minimal_adata, out_path)
        html = out_path.read_text()
        assert ">Qc<" in html
        assert "Preprocess" not in html or "Preprocess" not in html.split(">Qc<")[1]


class TestTopLevelAPI:
    def test_exposed_as_scl_attribute(self):
        assert hasattr(scl, "export_audit_report")
        assert callable(scl.export_audit_report)

    def test_exposed_under_utils(self):
        assert hasattr(scl.utils, "export_audit_report")

    def test_top_level_call_round_trip(self, sclucid_adata, tmp_path):
        out_path = tmp_path / "rt.html"
        result = scl.export_audit_report(sclucid_adata, out_path, title="Custom Title")
        assert result.exists()
        html = result.read_text()
        assert "Custom Title" in html


class TestOptions:
    def test_include_full_config_false_omits_workflow_config(
        self, sclucid_adata, tmp_path
    ):
        out_path = tmp_path / "compact.html"
        export_audit_report(sclucid_adata, out_path, include_full_config=False)
        html = out_path.read_text()
        # The collapsible "Effective workflow_config" detail should not appear
        assert "Effective workflow_config" not in html

    def test_include_full_config_true_includes_workflow_config(
        self, sclucid_adata, tmp_path
    ):
        out_path = tmp_path / "full.html"
        export_audit_report(sclucid_adata, out_path, include_full_config=True)
        html = out_path.read_text()
        assert "Effective workflow_config" in html


class TestAdditionalKeysRendered:
    """Module sections may carry non-canonical keys (e.g. execution_trace).

    Audit report must not silently drop them — that hides the recommendation
    engine's recommended_params snapshot, custom workflow metadata, and
    plugin-attached evidence.
    """

    def test_renders_execution_trace(self, sclucid_adata, tmp_path):
        from scLucid.utils.contracts import SCLUCID_ROOT, Modules

        sclucid_adata.uns[SCLUCID_ROOT][Modules.QC]["execution_trace"] = {
            "recommended_params": {"min_genes": 200},
            "user_overrides": {"min_genes": {"recommended": 200, "actual": 300}},
        }
        out_path = tmp_path / "trace.html"
        export_audit_report(sclucid_adata, out_path)
        html = out_path.read_text()
        assert "execution_trace" in html
        assert "user_overrides" in html
        # Recommended-vs-actual override key should appear
        assert "min_genes" in html
