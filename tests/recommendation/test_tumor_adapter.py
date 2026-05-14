"""Tests for tumor analysis recommendation adapter."""

import pytest
from anndata import AnnData

from scLucid.recommendation.tumor_adapter import adapt_tumor_recommendation
from scLucid.tumor.config import TumorAnalysisConfig


@pytest.fixture
def tumor_adata():
    """Synthetic AnnData with cell type annotations for tumor testing."""
    import numpy as np

    adata = AnnData(
        X=np.random.negative_binomial(5, 0.3, size=(800, 200)).astype(float)
    )
    adata.obs_names = [f"cell_{i}" for i in range(800)]
    adata.var_names = [f"gene_{i}" for i in range(200)]
    adata.obs["cell_type"] = ["T cell"] * 400 + ["Epithelial"] * 400
    return adata


@pytest.fixture
def small_tumor_adata():
    """Very small tumor AnnData (below malignancy minimum)."""
    import numpy as np

    adata = AnnData(
        X=np.random.negative_binomial(5, 0.3, size=(200, 100)).astype(float)
    )
    adata.obs_names = [f"cell_{i}" for i in range(200)]
    adata.var_names = [f"gene_{i}" for i in range(100)]
    return adata


class TestAdaptTumorRecommendation:
    def test_non_tumor_dataset_disabled(self, minimal_adata):
        section = adapt_tumor_recommendation(
            minimal_adata,
            config=TumorAnalysisConfig(),
        )
        assert section.name == "tumor"
        assert "disabled" in section.summary.lower()

        # All tumor modules should be False for non-tumor dataset
        for name in ["run_malignancy", "run_tme", "run_cnv", "run_therapy"]:
            param = section.get_parameter(name)
            assert param is not None
            assert param.value is False

    def test_tumor_dataset_enables_modules(self, tumor_adata):
        """Tumor module should enable when dataset type is tumor_tissue."""
        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="tumor_tissue", tissue="pancreas")

        section = adapt_tumor_recommendation(
            tumor_adata,
            config=TumorAnalysisConfig(),
            context=ctx,
            cancer_type="pancreatic_adenocarcinoma",
        )

        run_malignancy = section.get_parameter("run_malignancy")
        assert run_malignancy is not None
        assert run_malignancy.value is True

    def test_small_dataset_malignancy_disabled(self, small_tumor_adata):
        """Malignancy should be disabled for datasets < 500 cells."""
        section = adapt_tumor_recommendation(
            small_tumor_adata,
            config=TumorAnalysisConfig(run_malignancy=True),
        )
        # With n_cells < 500, malignancy should be disabled
        run_malignancy = section.get_parameter("run_malignancy")
        if run_malignancy is not None:
            assert run_malignancy.value is False or run_malignancy.confidence < 1.0

    def test_no_cell_type_tme_disabled(self, minimal_adata):
        """TME should be disabled without cell type annotations."""
        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="tumor_tissue", tissue="tumor")

        section = adapt_tumor_recommendation(
            minimal_adata,
            config=TumorAnalysisConfig(run_tme=True),
            context=ctx,
        )
        run_tme = section.get_parameter("run_tme")
        if run_tme is not None:
            # Without cell_type column, should be disabled or low confidence
            assert run_tme.value is False or run_tme.confidence <= 0.5

    def test_with_cnv_score_detected(self, tumor_adata):
        tumor_adata.obs["cnv_score"] = 1.0

        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="tumor_tissue", tissue="tumor")

        section = adapt_tumor_recommendation(
            tumor_adata,
            config=TumorAnalysisConfig(),
            context=ctx,
        )
        run_cnv = section.get_parameter("run_cnv")
        assert run_cnv is not None
        assert run_cnv.value is True
        assert run_cnv.confidence == 0.9

    def test_no_cancer_type_therapy_disabled(self, minimal_adata):
        """Therapy should be disabled without cancer_type."""
        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="tumor_tissue", tissue="tumor")

        section = adapt_tumor_recommendation(
            minimal_adata,
            config=TumorAnalysisConfig(),
            context=ctx,
            cancer_type=None,
        )
        run_therapy = section.get_parameter("run_therapy")
        assert run_therapy is not None
        assert run_therapy.value is False

    def test_section_structure(self, minimal_adata):
        section = adapt_tumor_recommendation(
            minimal_adata,
            config=TumorAnalysisConfig(),
        )
        assert section.name == "tumor"
        assert isinstance(section.summary, str)
        assert 0 <= section.confidence <= 1
        assert isinstance(section.parameters, list)
        assert isinstance(section.metadata, dict)
        assert "n_cells" in section.metadata

    def test_parameters_have_required_fields(self, minimal_adata):
        section = adapt_tumor_recommendation(minimal_adata)
        for param in section.parameters:
            assert param.name
            assert param.method
            assert param.rationale
            assert 0 <= param.confidence <= 1

    def test_uses_provided_context(self, minimal_adata):
        """When context is provided, it should be used directly."""
        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="tumor_tissue", tissue="pancreas")
        section = adapt_tumor_recommendation(
            minimal_adata,
            config=TumorAnalysisConfig(),
            context=ctx,
            cancer_type="pancreatic",
        )
        assert section.name == "tumor"
        assert section.metadata["cancer_type"] == "pancreatic"


class TestAlternativesPopulated:
    """Every recommended tumor parameter should expose ``alternatives``.

    The ``alternatives`` field is what tells users which other values the
    recommender considered before settling on the recommended one. If we
    leave it empty, downstream UIs (the audit report, plugin authors,
    notebook reviewers) have no way to show the user the option space.
    """

    def test_all_tumor_params_have_alternatives(self, tumor_adata):
        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="tumor_tissue", tissue="pancreas")
        section = adapt_tumor_recommendation(
            tumor_adata,
            config=TumorAnalysisConfig(),
            context=ctx,
            cancer_type="pancreatic",
        )
        # Each parameter should declare at least one alternative.
        for param in section.parameters:
            assert param.alternatives, (
                f"Parameter {param.name!r} on tumor adapter has empty "
                f"alternatives; users cannot see the option space."
            )

    def test_run_flags_alternatives_are_boolean(self, tumor_adata):
        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="tumor_tissue", tissue="pancreas")
        section = adapt_tumor_recommendation(
            tumor_adata,
            config=TumorAnalysisConfig(),
            context=ctx,
            cancer_type="pancreatic",
        )
        for name in ("run_malignancy", "run_tme", "run_cnv", "run_therapy"):
            param = section.get_parameter(name)
            assert param is not None
            assert set(param.alternatives) == {True, False}

    def test_malignancy_method_alternatives_cover_documented_methods(self, tumor_adata):
        from scLucid.utils.context import AnalysisContext

        ctx = AnalysisContext(dataset_type="tumor_tissue", tissue="pancreas")
        section = adapt_tumor_recommendation(
            tumor_adata,
            config=TumorAnalysisConfig(),
            context=ctx,
            cancer_type="pancreatic",
        )
        param = section.get_parameter("malignancy_method")
        assert param is not None
        assert set(param.alternatives) >= {"cnv", "threshold", "ml"}
