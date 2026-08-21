"""Tests for Python-native ambient RNA diagnostics and correction."""

import subprocess

import numpy as np
import pytest
import scipy.sparse as sparse
from anndata import AnnData

import scLucid.qc as qc
from scLucid.qc import (
    AMBIENT_CORRECTED_COUNTS_LAYER,
    diagnose_ambient_rna,
    register_external_ambient_result,
)
from scLucid.qc.ambient import (
    build_ambient_layer_contract,
    correct_ambient_rna_linear,
    diagnose_empty_droplets,
    infer_ambient_input_context,
    record_ambient_correction_status,
)
from scLucid.qc.ambient_backends import (
    _probe_r_capability,
    _rpy2_available,
    correct_ambient_rna,
    decontx_available,
    list_ambient_backends,
    soupx_available,
)


def test_rpy2_probe_contains_a_crashed_embedded_r(monkeypatch):
    """A signal exit from embedded R must be reported without touching this process."""
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return subprocess.CompletedProcess(command, returncode=-11, stdout="", stderr="segfault")

    _probe_r_capability.cache_clear()
    monkeypatch.setattr("scLucid.qc.ambient_backends.subprocess.run", fake_run)
    try:
        assert _rpy2_available() is False
    finally:
        _probe_r_capability.cache_clear()

    assert len(calls) == 1
    assert "rpy2.robjects" in calls[0][0][2]
    assert calls[0][1]["check"] is False
    assert calls[0][1]["capture_output"] is True


def test_r_package_probe_timeout_is_unavailable(monkeypatch):
    """A wedged R initialization must time out and leave the test process alive."""
    def fake_run(command, **kwargs):
        raise subprocess.TimeoutExpired(command, timeout=kwargs["timeout"])

    _probe_r_capability.cache_clear()
    monkeypatch.setattr("scLucid.qc.ambient_backends.subprocess.run", fake_run)
    try:
        assert _probe_r_capability("SoupX") is False
    finally:
        _probe_r_capability.cache_clear()


def test_diagnose_ambient_rna_returns_risk_summary():
    X = np.ones((100, 30), dtype=float)
    X[:10, :] = 0
    X[:10, 0] = 50
    X[:10, 1] = 30
    adata = AnnData(X=X)
    adata.var_names = [f"Gene{i}" for i in range(30)]

    summary = diagnose_ambient_rna(adata, low_count_quantile=0.1, top_n_genes=5)

    assert summary["available"] is True
    assert summary["diagnostic_only"] is True
    assert summary["method"] == "python_heuristic_low_count_enrichment"
    assert "no expression correction" in summary["method_note"]
    assert summary["calibration_status"] == "heuristic_unvalidated"
    assert summary["risk_score_weights"] == {
        "top_gene_dominance": 0.4,
        "low_count_enrichment": 0.4,
        "enriched_gene_breadth": 0.2,
    }
    assert "not calibrated" in summary["risk_score_note"]
    assert summary["risk_level"] in {"low", "moderate", "high"}
    assert "top_genes" in summary
    assert summary["correction_status"]["corrected"] is False
    assert "input_context" in summary


def test_record_ambient_correction_status_is_used_by_diagnostic():
    adata = AnnData(X=np.ones((50, 10), dtype=float))
    record_ambient_correction_status(
        adata,
        corrected=True,
        backend="cellbender",
        output_layer="cellbender_corrected",
        details={"source": "external_cli"},
    )

    summary = diagnose_ambient_rna(adata, top_n_genes=3)

    assert summary["correction_status"]["corrected"] is True
    assert summary["correction_status"]["backend"] == "cellbender"


def test_register_external_ambient_result_copies_matching_matrix_to_layer():
    adata = AnnData(X=np.ones((5, 4), dtype=float))
    corrected = AnnData(X=np.full((5, 4), 2.0, dtype=float))
    adata.obs_names = corrected.obs_names = [f"cell{i}" for i in range(5)]
    adata.var_names = corrected.var_names = [f"gene{i}" for i in range(4)]

    status = register_external_ambient_result(
        adata,
        backend="cellbender",
        corrected_adata=corrected,
        output_layer="cellbender_corrected",
    )

    assert status["corrected"] is True
    assert status["backend"] == "cellbender"
    assert "cellbender_corrected" in adata.layers
    assert np.asarray(adata.layers["cellbender_corrected"]).mean() == 2.0


def test_register_external_ambient_result_defaults_to_contract_layer():
    adata = AnnData(X=np.ones((5, 4), dtype=float))
    corrected = AnnData(X=np.full((5, 4), 2.0, dtype=float))
    adata.obs_names = corrected.obs_names = [f"cell{i}" for i in range(5)]
    adata.var_names = corrected.var_names = [f"gene{i}" for i in range(4)]
    adata.layers["counts"] = adata.X.copy()

    status = register_external_ambient_result(
        adata,
        backend="cellbender",
        corrected_adata=corrected,
    )

    contract = adata.uns["sclucid"]["qc"]["ambient_layer_contract"]
    assert status["output_layer"] == AMBIENT_CORRECTED_COUNTS_LAYER
    assert AMBIENT_CORRECTED_COUNTS_LAYER in adata.layers
    assert contract["corrected_layer"] == AMBIENT_CORRECTED_COUNTS_LAYER
    assert contract["corrected_layer_present"] is True
    assert contract["recommended_preprocess_counts_layer"] == AMBIENT_CORRECTED_COUNTS_LAYER


def test_register_external_ambient_result_copies_canonical_obs_columns():
    adata = AnnData(X=np.ones((5, 4), dtype=float))
    corrected = AnnData(X=np.full((5, 4), 2.0, dtype=float))
    adata.obs_names = corrected.obs_names = [f"cell{i}" for i in range(5)]
    adata.var_names = corrected.var_names = [f"gene{i}" for i in range(4)]
    corrected.obs["cellbender_probability"] = np.linspace(0.9, 0.99, 5)
    corrected.obs["decontx_rho"] = np.linspace(0.01, 0.05, 5)

    status = register_external_ambient_result(
        adata,
        backend="cellbender",
        corrected_adata=corrected,
        obs_column_map={
            "cell_probability": "cellbender_probability",
            "ambient_fraction": "decontx_rho",
        },
    )

    assert "cell_probability" in adata.obs
    assert "ambient_fraction" in adata.obs
    assert status["details"]["obs_column_map"] == {
        "cell_probability": "cellbender_probability",
        "ambient_fraction": "decontx_rho",
    }
    assert np.isclose(float(adata.obs["cell_probability"].min()), 0.9)
    assert np.isclose(float(adata.obs["ambient_fraction"].max()), 0.05)


def test_infer_ambient_input_context_distinguishes_filtered_and_raw_like():
    filtered = AnnData(X=np.ones((5, 4), dtype=float))
    raw_like = AnnData(X=np.ones((5, 4), dtype=float))
    raw_like.obs["likely_empty_droplet"] = ["true", "false", "false", "false", "false"]

    assert infer_ambient_input_context(filtered)["matrix_type"] == "filtered_like"
    raw_context = infer_ambient_input_context(raw_like)
    assert raw_context["matrix_type"] == "raw_like"
    assert "cellbender" in raw_context["suitable_backends"]


def test_ambient_layer_contract_falls_back_to_counts_without_correction():
    adata = AnnData(X=np.ones((5, 4), dtype=float))
    adata.layers["counts"] = adata.X.copy()

    contract = build_ambient_layer_contract(
        adata,
        correction_summary={"corrected": False, "reason": "not_requested"},
    )

    assert contract["recommended_preprocess_counts_layer"] == "counts"
    assert "ambient_corrected_counts(optional)" in contract["canonical_flow"]


def test_diagnose_empty_droplets_is_pure_by_default_and_can_record_explicitly():
    X = np.ones((120, 20), dtype=float)
    X[:20, :] = 0
    X[:20, 0] = 10
    X[:20, 1] = 8
    adata = AnnData(X=X)
    adata.var_names = [f"Gene{i}" for i in range(20)]

    summary = diagnose_empty_droplets(adata, min_barcodes=50, top_n_genes=5)

    assert summary["available"] is True
    assert summary["diagnostic_only"] is True
    assert summary["method"] == "python_barcode_rank_background_profile"
    assert "not an EmptyDrops replacement" in summary["method_note"]
    assert "top_background_genes" in summary
    assert "sclucid" not in adata.uns or "qc" not in adata.uns.get("sclucid", {})

    recorded = diagnose_empty_droplets(adata, min_barcodes=50, top_n_genes=5, record=True)
    assert recorded["available"] is True
    assert "empty_droplet_summary" in adata.uns["sclucid"]["qc"]


def test_boolean_empty_droplet_indicator_uses_true_as_empty():
    empty = np.tile(np.array([[20, 5, 0, 0]], dtype=float), (20, 1))
    cells = np.tile(np.array([[1, 0, 5, 4]], dtype=float), (100, 1))
    adata = AnnData(X=np.vstack([empty, cells]))
    adata.obs["likely_empty_droplet"] = [True] * 20 + [False] * 100

    diagnostic = diagnose_empty_droplets(
        adata,
        cell_call_key="likely_empty_droplet",
        min_barcodes=50,
    )
    correction = correct_ambient_rna_linear(
        adata,
        empty_droplet_key="likely_empty_droplet",
        output_layer="ambient_corrected",
    )

    assert diagnostic["n_putative_empty_droplets"] == 20
    assert diagnostic["n_called_cells"] == 100
    assert correction["n_putative_empty_droplets"] == 20


def test_linear_correction_reduces_ambient_marker_expression():
    """Synthetic ambient profile should be partially removed by linear correction."""
    rng = np.random.default_rng(42)
    n_cells = 500
    n_genes = 50
    # True expression matrix with sparse counts.
    true_expr = rng.poisson(lam=2, size=(n_cells, n_genes)).astype(float)
    true_expr[true_expr < 1] = 0.0

    # Ambient profile concentrated in first 5 genes.
    ambient_profile = np.zeros(n_genes)
    ambient_profile[:5] = np.array([0.4, 0.3, 0.15, 0.1, 0.05])
    # Per-cell contamination fractions.
    rho = rng.uniform(0.05, 0.25, size=n_cells)
    ambient_counts = (rho[:, None] * true_expr.sum(axis=1)[:, None] * ambient_profile[None, :])
    observed = true_expr + ambient_counts
    observed = np.round(observed).astype(float)

    adata = AnnData(X=observed)
    adata.var_names = [f"Gene{i}" for i in range(n_genes)]

    summary = correct_ambient_rna_linear(adata, output_layer="ambient_corrected")

    assert summary["corrected"] is True
    assert "ambient_corrected" in adata.layers
    corrected = np.asarray(adata.layers["ambient_corrected"])
    # Ambient marker genes should lose counts on average.
    ambient_before = observed[:, :5].sum()
    ambient_after = corrected[:, :5].sum()
    assert ambient_after < ambient_before
    # Non-ambient genes should be nearly unchanged (allow some shrinkage noise
    # because the linear estimator uses all genes to scale rho, not just markers).
    non_ambient_before = observed[:, 5:].sum()
    non_ambient_after = corrected[:, 5:].sum()
    assert non_ambient_after > non_ambient_before * 0.75
    assert "residual_ambient_score" in summary
    assert 0.0 <= summary["residual_ambient_score"] <= 1.0


def test_linear_correction_preserves_sparse_output_without_dense_expected_matrix():
    empty = np.tile(np.array([[20, 5, 0, 0], [18, 4, 0, 1]], dtype=float), (5, 1))
    cells = np.tile(np.array([[1, 0, 5, 4], [1, 0, 4, 5]], dtype=float), (5, 1))
    X = sparse.csr_matrix(np.vstack([empty, cells]))
    adata = AnnData(X=X)
    adata.obs["empty_droplet"] = ["empty"] * empty.shape[0] + ["cell"] * cells.shape[0]

    summary = correct_ambient_rna_linear(
        adata,
        empty_droplet_key="empty_droplet",
        output_layer="ambient_corrected",
    )

    assert summary["corrected"] is True
    assert sparse.issparse(adata.layers["ambient_corrected"])
    assert adata.layers["ambient_corrected"].shape == adata.X.shape


def test_linear_correction_default_layer_contract():
    rng = np.random.default_rng(10)
    X = rng.poisson(lam=3, size=(200, 20)).astype(float)
    X[:30, :] = 0
    X[:30, :4] = rng.poisson(lam=5, size=(30, 4)).astype(float)
    adata = AnnData(X=X)
    adata.var_names = [f"Gene{i}" for i in range(20)]
    adata.layers["counts"] = adata.X.copy()

    summary = correct_ambient_rna_linear(adata)

    contract = adata.uns["sclucid"]["qc"]["ambient_layer_contract"]
    assert summary["output_layer"] == AMBIENT_CORRECTED_COUNTS_LAYER
    assert AMBIENT_CORRECTED_COUNTS_LAYER in adata.layers
    assert contract["recommended_preprocess_counts_layer"] == AMBIENT_CORRECTED_COUNTS_LAYER
    assert "artifact_contract" in adata.uns["sclucid"]["qc"]


def test_ambient_top_level_api_is_narrow_and_internal_api_stays_in_submodule():
    public = set(qc.__all__)
    for symbol in [
        "AMBIENT_CORRECTED_COUNTS_LAYER",
        "diagnose_ambient_rna",
        "register_external_ambient_result",
        "correct_ambient_rna",
    ]:
        assert symbol in public

    for hidden in [
        "build_ambient_layer_contract",
        "diagnose_empty_droplets",
        "infer_ambient_input_context",
        "record_ambient_correction_status",
        "record_ambient_layer_contract",
        "correct_ambient_rna_linear",
    ]:
        assert not hasattr(qc, hidden)
        assert hidden not in public


def test_correct_ambient_rna_auto_filtered_returns_diagnostic_without_backend():
    """Auto method on filtered matrices should not silently run linear correction."""
    rng = np.random.default_rng(7)
    adata = AnnData(X=rng.poisson(lam=5, size=(200, 20)).astype(float))
    adata.var_names = [f"Gene{i}" for i in range(20)]

    summary = correct_ambient_rna(
        adata, method="auto", backend="auto", output_layer="ambient_corrected"
    )
    if not (soupx_available() or decontx_available()):
        assert summary["corrected"] is False
        assert summary["method"] == "diagnostic_only"
        assert summary["matrix_type"] == "filtered_like"
        assert "ambient_corrected" not in adata.layers
        contract = adata.uns["sclucid"]["qc"]["ambient_layer_contract"]
        assert contract["recommended_preprocess_counts_layer"] in (None, "")


def test_correct_ambient_rna_linear_explicit():
    """Explicit linear method should always work without optional deps."""
    rng = np.random.default_rng(8)
    adata = AnnData(X=rng.poisson(lam=3, size=(200, 15)).astype(float))
    adata.var_names = [f"Gene{i}" for i in range(15)]

    summary = correct_ambient_rna(
        adata, method="linear", output_layer="ambient_corrected"
    )
    assert summary["corrected"] is True
    assert summary["method"] == "linear_background_subtraction"
    assert summary["calibration_status"] == "conservative_fallback_not_model_based"
    assert "not equivalent" in summary["correction_note"]
    assert "ambient_corrected" in adata.layers


def test_list_ambient_backends_reports_availability():
    backends = list_ambient_backends()
    for name in ("cellbender", "soupx", "decontx"):
        assert name in backends
        assert "available_now" in backends[name]
        assert "matrix_types" in backends[name]


def test_correct_ambient_rna_filtered_prefers_decontx_or_soupx():
    """Auto on a filtered matrix should choose a filtered-matrix backend if available."""
    rng = np.random.default_rng(9)
    adata = AnnData(X=rng.poisson(lam=5, size=(100, 20)).astype(float))
    adata.var_names = [f"Gene{i}" for i in range(20)]

    summary = correct_ambient_rna(
        adata,
        method="external",
        backend="auto",
        output_layer="ambient_corrected",
    )

    if not (soupx_available() or decontx_available()):
        assert summary["corrected"] is False
        assert summary["method"] == "diagnostic_only"
        assert summary["reason"] == "no_filtered_matrix_backend_available"
    else:
        assert summary["corrected"] is True


def test_correct_ambient_rna_unknown_backend_raises():
    rng = np.random.default_rng(10)
    adata = AnnData(X=rng.poisson(lam=2, size=(50, 10)).astype(float))
    adata.var_names = [f"Gene{i}" for i in range(10)]

    with pytest.raises(ValueError):
        correct_ambient_rna(adata, method="external", backend="not_a_backend")
