"""Public API checks for scLucid.preprocess."""

import pytest

import scLucid.preprocess as pp
from scLucid.preprocess.config import PreprocessingWorkflowConfig


@pytest.mark.unit
def test_preprocess_exports_resolve():
    for symbol in pp.__all__:
        assert hasattr(pp, symbol), f"scLucid.preprocess missing exported symbol: {symbol}"


@pytest.mark.unit
def test_normalization_config_success_and_reserved_layer_validation():
    cfg = pp.NormalizationConfig(target_sum=1e4, output_layer="normalized")
    assert cfg.target_sum == 1e4

    with pytest.raises(ValueError):
        pp.NormalizationConfig(output_layer="X")


@pytest.mark.unit
def test_preprocess_exports_gene_biotype_utilities():
    required = [
        "GeneBiotypeConfig",
        "annotate_gene_biotypes",
        "filter_genes_by_biotype",
        "get_biotype_statistics",
        "load_gene_biotypes",
        "recommend_biotype_strategy",
    ]
    for symbol in required:
        assert hasattr(pp, symbol), f"scLucid.preprocess missing gene biotype utility: {symbol}"


@pytest.mark.unit
def test_preprocess_legacy_aliases_are_not_top_level_api():
    removed = [
        "apply_gene_biotype_strategy",
        "get_gene_biotype_cache_dir",
        "list_gene_biotype_resources",
        "run_embedding_workflow",
        "adaptive_normalize",
        "quality_aware_normalize",
        "IntelligentPreprocessConfig",
        "IntelligentPreprocessRecommender",
        "PreprocessingStrategy",
        "run_intelligent_preprocessing",
        "recommend_intelligent_preprocessing",
        "PreprocessingBackend",
        "ScanpyBackend",
        "RapidsBackend",
        "get_backend",
        "set_backend",
        "list_available_backends",
    ]
    for symbol in removed:
        assert not hasattr(pp, symbol), f"legacy symbol should not be top-level API: {symbol}"
        assert symbol not in pp.__all__


@pytest.mark.unit
def test_iterative_preprocessing_is_public_default_advanced_entrypoint():
    assert hasattr(pp, "run_iterative_preprocessing")
    assert "run_iterative_preprocessing" in pp.__all__


@pytest.mark.unit
def test_from_simple_dict_does_not_mutate_input():
    simple = {
        "normalization_method": "standard",
        "hvg_n_top_genes": 1500,
        "results_dir": "./results",
        "run_regression": False,
    }
    original = dict(simple)

    config = PreprocessingWorkflowConfig.from_simple_dict(simple)

    assert simple == original
    assert config.normalization.method == "standard"
    assert config.hvg.n_top_genes == 1500
    assert config.save_dir == "./results"
    assert config.run_regression is False


@pytest.mark.unit
def test_from_simple_dict_accepts_gene_biotype_options():
    config = PreprocessingWorkflowConfig.from_simple_dict(
        {
            "gene_biotype_annotate": True,
            "gene_biotype_filter": True,
            "gene_biotype_method": "custom",
            "gene_biotype_custom_biotype_path": "gene_biotypes.csv",
            "gene_biotype_keep_biotypes": ["protein_coding"],
        }
    )

    assert config.gene_biotype.annotate is True
    assert config.gene_biotype.filter is True
    assert config.gene_biotype.method == "custom"
    assert config.gene_biotype.custom_biotype_path == "gene_biotypes.csv"
    assert config.gene_biotype.keep_biotypes == ["protein_coding"]
