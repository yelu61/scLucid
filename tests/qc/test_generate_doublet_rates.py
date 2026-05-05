"""Tests for scLucid.qc.generate_doublet_rates."""

import numpy as np
import pytest
from anndata import AnnData


def _make_adata_with_samples(n_cells: int = 200, sample_key: str = "sampleID"):
    """Create a minimal AnnData with sample labels."""
    rng = np.random.default_rng(42)
    X = rng.integers(0, 10, size=(n_cells, 50)).astype(np.float32)
    adata = AnnData(X)
    adata.obs_names = [f"cell_{i:04d}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i:03d}" for i in range(50)]
    # Two samples
    adata.obs[sample_key] = np.where(np.arange(n_cells) < n_cells // 2, "s1", "s2")
    return adata


class TestGenerateDoubletRates:
    def test_generate_doublet_rates_v3(self):
        """v3 chemistry returns scale model with 0.008 per 1k cells."""
        from scLucid.qc.doublet.core import generate_doublet_rates

        adata = _make_adata_with_samples(n_cells=2000)
        rates = generate_doublet_rates(adata, chemistry="v3")

        assert "s1" in rates
        assert "s2" in rates
        # 1000 cells * 0.008 = 0.008 (capped at max 0.20)
        assert 0.0 < rates["s1"] <= 0.20
        assert rates["s1"] == pytest.approx(0.008, abs=0.001)

    def test_generate_doublet_rates_v2(self):
        """v2 chemistry uses 0.007 per 1k cells."""
        from scLucid.qc.doublet.core import generate_doublet_rates

        adata = _make_adata_with_samples(n_cells=2000)
        rates = generate_doublet_rates(adata, chemistry="v2")

        assert rates["s1"] == pytest.approx(0.007, abs=0.001)

    def test_generate_doublet_rates_bd_fixed(self):
        """BD chemistry returns a fixed 0.025 rate."""
        from scLucid.qc.doublet.core import generate_doublet_rates

        adata = _make_adata_with_samples(n_cells=5000)
        rates = generate_doublet_rates(adata, chemistry="BD")

        assert rates["s1"] == pytest.approx(0.025, abs=0.001)
        assert rates["s2"] == pytest.approx(0.025, abs=0.001)

    def test_generate_doublet_rates_custom_scale(self):
        """Custom scale model with user-provided rate factor."""
        from scLucid.qc.doublet.core import generate_doublet_rates

        adata = _make_adata_with_samples(n_cells=2000)
        rates = generate_doublet_rates(
            adata,
            chemistry="custom",
            custom_rate=0.010,
            custom_rate_model="scale",
        )

        assert rates["s1"] == pytest.approx(0.010, abs=0.001)

    def test_generate_doublet_rates_custom_fixed(self):
        """Custom fixed model applies the same rate to all samples."""
        from scLucid.qc.doublet.core import generate_doublet_rates

        adata = _make_adata_with_samples(n_cells=5000)
        rates = generate_doublet_rates(
            adata,
            chemistry="custom",
            custom_rate=0.05,
            custom_rate_model="fixed",
        )

        assert rates["s1"] == pytest.approx(0.05, abs=0.001)
        assert rates["s2"] == pytest.approx(0.05, abs=0.001)

    def test_generate_doublet_rates_max_rate_cap(self):
        """Rates are capped at max_rate."""
        from scLucid.qc.doublet.core import generate_doublet_rates

        # Very large sample should hit the max_rate cap
        adata = _make_adata_with_samples(n_cells=50_000)
        rates = generate_doublet_rates(adata, chemistry="v3", max_rate=0.10)

        assert rates["s1"] <= 0.10 + 1e-6
        assert rates["s2"] <= 0.10 + 1e-6

    def test_generate_doublet_rates_single_sample(self):
        """Works with a single sample."""
        from scLucid.qc.doublet.core import generate_doublet_rates

        rng = np.random.default_rng(42)
        X = rng.integers(0, 10, size=(100, 20)).astype(np.float32)
        adata = AnnData(X)
        adata.obs["sampleID"] = "only_one"

        rates = generate_doublet_rates(adata, chemistry="v3")
        assert "only_one" in rates
        assert 0.0 < rates["only_one"] <= 0.20

    def test_generate_doublet_rates_invalid_chemistry_fallback(self):
        """Invalid chemistry falls back to v3 with a warning."""
        from scLucid.qc.doublet.core import generate_doublet_rates

        adata = _make_adata_with_samples(n_cells=100)
        rates = generate_doublet_rates(adata, chemistry="not_real")
        # Falls back to v3 model
        assert "s1" in rates
        assert 0.0 < rates["s1"] <= 0.20
