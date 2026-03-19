"""
Tests for pyDWLS (R-free DWLS implementation)
"""

import pytest
import numpy as np
import pandas as pd
import sys

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

try:
    from scLucid.tools.pyDWLS import (
        DWLS,
        SignatureBuilder,
        DampenedWLS,
        MarkerSelector,
        CrossValidator,
        solve_nnls,
        normalize_data,
        filter_genes,
        create_pseudo_bulk,
    )
except Exception as exc:  # pragma: no cover - optional backend availability
    pytest.skip(f"Skipping pyDWLS tests: {exc}", allow_module_level=True)


@pytest.fixture
def sample_sc_data():
    """Generate sample single-cell data"""
    np.random.seed(42)
    n_genes = 100
    n_cells = 200

    sc_data = pd.DataFrame(
        np.random.poisson(5, (n_genes, n_cells)),
        index=[f"gene_{i}" for i in range(n_genes)],
        columns=[f"cell_{i}" for i in range(n_cells)],
    )

    cell_types = pd.Series(
        np.random.choice(['T_cell', 'B_cell', 'Macrophage'], n_cells),
        index=sc_data.columns,
    )

    return sc_data, cell_types


@pytest.fixture
def sample_bulk_data():
    """Generate sample bulk data"""
    np.random.seed(42)
    n_genes = 100
    n_samples = 10

    bulk_data = pd.DataFrame(
        np.random.poisson(500, (n_genes, n_samples)),
        index=[f"gene_{i}" for i in range(n_genes)],
        columns=[f"sample_{i}" for i in range(n_samples)],
    )

    return bulk_data


@pytest.mark.unit
class TestSignatureBuilder:
    """Test signature matrix building"""

    def test_build_signature(self, sample_sc_data):
        """Test signature building"""
        sc_data, cell_types = sample_sc_data

        builder = SignatureBuilder()
        signature = builder.build(sc_data, cell_types)

        assert isinstance(signature, pd.DataFrame)
        assert signature.shape[0] == sc_data.shape[0]
        assert signature.shape[1] == len(cell_types.unique())

    def test_trimmed_mean(self, sample_sc_data):
        """Test trimmed mean aggregation"""
        sc_data, cell_types = sample_sc_data

        builder = SignatureBuilder(trim_percent=0.1)
        signature = builder.build(sc_data, cell_types, method="trimmed_mean")

        assert signature.shape[0] == sc_data.shape[0]

    def test_min_cells_filter(self, sample_sc_data):
        """Test minimum cells filtering"""
        sc_data, cell_types = sample_sc_data

        builder = SignatureBuilder(min_cells=100)
        signature = builder.build(sc_data, cell_types)

        # Should only include cell types with >= 100 cells
        for ct in signature.columns:
            assert (cell_types == ct).sum() >= 100


@pytest.mark.unit
class TestDampenedWLS:
    """Test Dampened WLS solver"""

    def test_solve(self):
        """Test basic solving"""
        np.random.seed(42)

        n_genes = 50
        n_cell_types = 3

        # Create simple signature and bulk
        S = np.random.rand(n_genes, n_cell_types)
        b = S @ np.array([0.3, 0.4, 0.3]) + np.random.randn(n_genes) * 0.01

        solver = DampenedWLS(dampen_factor=1.0)
        proportions = solver.solve(S, b)

        assert proportions.shape == (n_cell_types,)
        assert np.all(proportions >= 0)
        assert np.abs(proportions.sum() - 1.0) < 0.01

    def test_dampening_effect(self):
        """Test that dampening affects results"""
        np.random.seed(42)

        n_genes = 50
        n_cell_types = 3

        S = np.random.rand(n_genes, n_cell_types)
        S[0, :] *= 10  # Make first gene highly expressed
        b = S @ np.array([0.3, 0.4, 0.3])

        solver_no_dampen = DampenedWLS(dampen_factor=0.0)
        solver_dampen = DampenedWLS(dampen_factor=2.0)

        props_no_dampen = solver_no_dampen.solve(S, b)
        props_dampen = solver_dampen.solve(S, b)

        # Results should be different
        assert not np.allclose(props_no_dampen, props_dampen)


@pytest.mark.unit
class TestMarkerSelector:
    """Test marker gene selection"""

    def test_select_markers(self, sample_sc_data):
        """Test marker selection"""
        sc_data, cell_types = sample_sc_data

        selector = MarkerSelector()
        markers = selector.select(sc_data, cell_types, n_markers=10)

        assert len(markers) <= 30  # Max 10 per cell type * 3 types
        assert len(markers) > 0

    def test_different_methods(self, sample_sc_data):
        """Test different selection methods"""
        sc_data, cell_types = sample_sc_data

        selector = MarkerSelector()

        markers_ratio = selector.select(sc_data, cell_types, method="ratio")
        markers_diff = selector.select(sc_data, cell_types, method="difference")
        markers_fc = selector.select(sc_data, cell_types, method="fold_change")

        assert len(markers_ratio) > 0
        assert len(markers_diff) > 0
        assert len(markers_fc) > 0


@pytest.mark.unit
class TestDWLS:
    """Test main DWLS class"""

    def test_initialization(self):
        """Test DWLS initialization"""
        dwls = DWLS(dampen_factor=0.5, use_nonneg=True)

        assert dwls.dampen_factor == 0.5
        assert dwls.use_nonneg is True

    def test_build_signature(self, sample_sc_data):
        """Test signature building via DWLS"""
        sc_data, cell_types = sample_sc_data

        dwls = DWLS()
        signature = dwls.build_signature_matrix(sc_data, cell_types)

        assert dwls.signature_matrix is not None
        assert signature.shape[1] == len(cell_types.unique())

    def test_deconvolve(self, sample_sc_data, sample_bulk_data):
        """Test full deconvolution"""
        sc_data, cell_types = sample_sc_data

        dwls = DWLS()
        dwls.build_signature_matrix(sc_data, cell_types)

        # Use subset of genes that exist in both
        common_genes = sc_data.index.intersection(sample_bulk_data.index)
        bulk_subset = sample_bulk_data.loc[common_genes]

        proportions = dwls.deconvolve(bulk_subset, verbose=False)

        assert proportions.shape == (bulk_subset.shape[1], len(cell_types.unique()))
        assert np.all(proportions >= 0)

        # Check rows sum to ~1
        row_sums = proportions.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=0.01)

    def test_select_and_build(self, sample_sc_data):
        """Test marker selection and signature building"""
        sc_data, cell_types = sample_sc_data

        dwls = DWLS()
        markers = dwls.select_marker_genes(sc_data, cell_types, n_markers=10)
        signature = dwls.build_signature_matrix(
            sc_data, cell_types, genes_to_use=markers
        )

        assert signature.shape[0] == len(markers)


@pytest.mark.unit
class TestUtils:
    """Test utility functions"""

    def test_normalize_data(self):
        """Test data normalization"""
        data = pd.DataFrame(np.random.poisson(100, (50, 10)))

        normalized = normalize_data(data, method="cpm")

        # CPM should sum to 1e6
        col_sums = normalized.sum(axis=0)
        assert np.allclose(col_sums, 1e6)

    def test_filter_genes(self):
        """Test gene filtering"""
        data = pd.DataFrame(np.random.poisson(5, (100, 20)))

        filtered = filter_genes(data, min_cells=5, min_expression=1)

        assert filtered.shape[0] <= data.shape[0]

    def test_create_pseudo_bulk(self, sample_sc_data):
        """Test pseudo-bulk creation"""
        sc_data, cell_types = sample_sc_data

        pseudo_bulk, true_props = create_pseudo_bulk(
            sc_data, cell_types, n_cells=50, random_state=42
        )

        assert len(pseudo_bulk) == sc_data.shape[0]
        assert len(true_props) <= len(cell_types.unique())
        assert np.abs(true_props.sum() - 1.0) < 0.01


@pytest.mark.unit
class TestCrossValidator:
    """Test cross-validation"""

    def test_cv_initialization(self):
        """Test CV initialization"""
        cv = CrossValidator(n_folds=5, n_cells_per_bulk=100)

        assert cv.n_folds == 5
        assert cv.n_cells_per_bulk == 100

    def test_cross_validate(self, sample_sc_data):
        """Test cross-validation"""
        sc_data, cell_types = sample_sc_data

        cv = CrossValidator(n_folds=3, n_cells_per_bulk=50)
        results = cv.run(sc_data, cell_types, verbose=False)

        assert 'mean_correlation' in results
        assert 'mean_rmse' in results
        assert len(results['fold_correlations']) == 3


@pytest.mark.integration
class TestFullWorkflow:
    """Integration test for complete workflow"""

    def test_complete_pipeline(self):
        """Test complete DWLS pipeline"""
        np.random.seed(42)

        # Generate synthetic data
        n_genes = 200
        n_cells = 300
        n_cell_types = 3

        cell_types = pd.Series(
            np.random.choice(['A', 'B', 'C'], n_cells),
            index=[f"cell_{i}" for i in range(n_cells)],
        )

        sc_data = pd.DataFrame(
            np.random.poisson(5, (n_genes, n_cells)),
            index=[f"gene_{i}" for i in range(n_genes)],
            columns=cell_types.index,
        )

        # Create bulk data
        n_samples = 5
        bulk_data = pd.DataFrame(
            np.random.poisson(500, (n_genes, n_samples)),
            index=sc_data.index,
            columns=[f"sample_{i}" for i in range(n_samples)],
        )

        # Full workflow
        dwls = DWLS()

        # Select markers
        markers = dwls.select_marker_genes(sc_data, cell_types, n_markers=20)

        # Build signature with markers
        signature = dwls.build_signature_matrix(
            sc_data, cell_types, genes_to_use=markers
        )

        # Deconvolve
        proportions = dwls.deconvolve(bulk_data, verbose=False)

        assert proportions.shape == (n_samples, n_cell_types)
        assert np.all(proportions >= 0)
        assert np.allclose(proportions.sum(axis=1), 1.0, atol=0.01)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
