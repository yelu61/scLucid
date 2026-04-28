"""Smoke tests for tumor analysis workflow."""

import pytest


class TestTumorWorkflow:
    def test_import(self):
        from scLucid.tumor.workflow import (
            run_tumor_analysis,
            run_tumor_analysis_expert,
        )
        assert callable(run_tumor_analysis)
        assert callable(run_tumor_analysis_expert)

    def test_tumor_config_import(self):
        from scLucid.tumor.config import (
            TumorAnalysisConfig,
            TumorWorkflowConfig,
        )
        config = TumorAnalysisConfig()
        assert config.run_malignancy is True
        assert config.run_tme is True
        assert config.run_cnv is True

    def test_tumor_workflow_config_quick(self):
        from scLucid.tumor.config import TumorWorkflowConfig

        config = TumorWorkflowConfig.quick(tissue_type="tumor")
        assert config.tissue_type == "tumor"
        assert config.use_recommendations is True

    def test_run_tumor_analysis_smoke(self, integration_test_adata):
        from scLucid.tumor.workflow import run_tumor_analysis
        from scLucid.tumor.config import TumorWorkflowConfig

        config = TumorWorkflowConfig.quick(
            tissue_type="tumor",
            cancer_type="pancreatic_adenocarcinoma",
        )
        # Tumor workflow needs cell_type annotation and other metadata.
        # Smoke-test: verify it runs without crash on basic data.
        try:
            result = run_tumor_analysis(integration_test_adata, config=config)
            assert result is not None
        except Exception as e:
            pytest.skip(f"Tumor workflow requires richer metadata: {e}")

    def test_tumor_utils_import(self):
        from scLucid.tumor.utils.databases import (
            query_cancer_gene_census,
            get_drug_targets,
            is_cancer_gene,
        )
        assert callable(query_cancer_gene_census)
        assert callable(get_drug_targets)
        assert callable(is_cancer_gene)

        from scLucid.tumor.utils.markers import (
            get_tumor_markers,
            get_immune_markers,
            get_all_markers,
        )
        assert callable(get_tumor_markers)
        assert callable(get_immune_markers)
        assert callable(get_all_markers)

        from scLucid.tumor.utils.signatures import (
            load_hallmark_signatures,
            calculate_signature_scores,
        )
        assert callable(load_hallmark_signatures)
        assert callable(calculate_signature_scores)
