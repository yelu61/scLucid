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
    if not hasattr(pp, "NormalizationConfig"):
        pytest.skip("NormalizationConfig not available in current environment")

    cfg = pp.NormalizationConfig(target_sum=1e4, output_layer="normalized")
    assert cfg.target_sum == 1e4

    with pytest.raises(ValueError):
        pp.NormalizationConfig(output_layer="X")


@pytest.mark.unit
def test_preprocess_exports_gene_biotype_utilities():
    required = [
        "apply_gene_biotype_strategy",
        "annotate_gene_biotypes",
        "filter_genes_by_biotype",
        "get_biotype_statistics",
        "get_gene_biotype_cache_dir",
        "list_gene_biotype_resources",
        "load_gene_biotypes",
        "recommend_biotype_strategy",
    ]
    for symbol in required:
        assert hasattr(pp, symbol), f"scLucid.preprocess missing gene biotype utility: {symbol}"


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
