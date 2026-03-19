"""
Tests for pyBayesPrism (R-free BayesPrism implementation)
"""

import pytest
import numpy as np
import pandas as pd
from scipy import sparse
import sys

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

from scLucid.tools.pyBayesPrism import (
    PrismConfig,
    ReferenceConfig,
    BayesPrismReference,
    BayesPrism,
    BayesPrismEmbedding,
    GibbsSampler,
    cleanup_genes,
    compute_correlation,
    compute_rmse,
    validate_inputs,
)


@pytest.fixture
def sample_reference():
    """Generate sample reference data"""
    np.random.seed(42)
    n_genes = 100
    n_cells = 200

    # Generate sparse reference data
    reference = sparse.random(n_genes, n_cells, density=0.1, format='csr')
    reference.data = np.random.poisson(5, size=reference.data.shape)

    # Cell type labels
    cell_types = pd.Series(
        np.random.choice(['T_cell', 'B_cell', 'Macrophage', 'Tumor'], size=n_cells)
    )

    gene_names = [f'Gene_{i}' for i in range(n_genes)]

    return reference, cell_types, gene_names


@pytest.fixture
def sample_mixture():
    """Generate sample bulk mixture data"""
    np.random.seed(42)
    n_genes = 100
    n_samples = 20

    mixture = pd.DataFrame(
        np.random.poisson(1000, (n_genes, n_samples)),
        index=[f'Gene_{i}' for i in range(n_genes)],
        columns=[f'Sample_{i}' for i in range(n_samples)],
    )

    return mixture


@pytest.fixture
def bayes_prism_obj(sample_reference, sample_mixture):
    """Create BayesPrism object"""
    reference, cell_types, gene_names = sample_reference

    ref = BayesPrismReference(
        reference=reference,
        cell_type_labels=cell_types,
    )

    bp = BayesPrism(
        reference=ref,
        mixture=sample_mixture,
    )

    return bp


@pytest.mark.unit
class TestPrismConfig:
    """Test PrismConfig"""

    def test_default_config(self):
        """Test default configuration"""
        config = PrismConfig()

        assert config.n_iter == 100
        assert config.burnin == 50
        assert config.n_chains == 4

    def test_config_validation(self):
        """Test config validation"""
        config = PrismConfig(n_iter=100, burnin=50)
        validated = config.validate()
        assert isinstance(validated, PrismConfig)
        assert validated.n_iter == 100

        with pytest.raises(ValueError):
            bad_config = PrismConfig(n_iter=10, burnin=20)
            bad_config.validate()

    def test_config_to_dict(self):
        """Test config serialization"""
        config = PrismConfig(n_iter=200)
        config_dict = config.to_dict()

        assert config_dict['n_iter'] == 200
        assert 'gibbs_control' in config_dict


@pytest.mark.unit
class TestBayesPrismReference:
    """Test BayesPrismReference"""

    def test_initialization(self, sample_reference):
        """Test reference initialization"""
        reference, cell_types, gene_names = sample_reference

        ref = BayesPrismReference(
            reference=reference,
            cell_type_labels=cell_types,
        )

        assert ref.reference_matrix.shape == (100, 200)
        assert len(ref.cell_types) == 4
        assert ref.phi.shape == (100, 4)

    def test_get_marker_genes(self, sample_reference):
        """Test marker gene selection"""
        reference, cell_types, gene_names = sample_reference

        ref = BayesPrismReference(
            reference=reference,
            cell_type_labels=cell_types,
        )

        markers = ref.get_marker_genes(n_markers=10)

        assert len(markers) == 4  # One per cell type
        for cell_type, genes in markers.items():
            assert len(genes) <= 10

    def test_filter_genes(self, sample_reference):
        """Test gene filtering"""
        reference, cell_types, gene_names = sample_reference

        ref = BayesPrismReference(
            reference=reference,
            cell_type_labels=cell_types,
        )

        genes_to_keep = gene_names[:50]
        ref.filter_genes(genes_to_keep)

        assert ref.reference_matrix.shape[0] == 50
        assert ref.phi.shape[0] == 50


@pytest.mark.unit
class TestGibbsSampler:
    """Test GibbsSampler"""

    def test_sampler_initialization(self):
        """Test sampler initialization"""
        sampler = GibbsSampler(n_iter=50, burnin=20)

        assert sampler.n_iter == 50
        assert sampler.burnin == 20

    def test_sample(self):
        """Test Gibbs sampling"""
        np.random.seed(42)

        n_genes = 50
        n_cell_types = 3

        mixture = np.random.poisson(100, n_genes).astype(float)
        theta_init = np.array([0.3, 0.4, 0.3])
        phi = np.random.dirichlet(np.ones(n_cell_types), n_genes)

        sampler = GibbsSampler(n_iter=20, burnin=10, use_numba=False)
        theta_samples, Z_samples = sampler.sample(
            mixture=mixture,
            theta_init=theta_init,
            phi=phi,
            verbose=False,
        )

        assert theta_samples.shape[1] == n_cell_types
        assert Z_samples.shape[1] == n_genes
        assert Z_samples.shape[2] == n_cell_types

    def test_posterior_mean(self):
        """Test posterior mean computation"""
        np.random.seed(42)

        sampler = GibbsSampler(n_iter=20, burnin=10, use_numba=False)

        mixture = np.random.poisson(100, 50).astype(float)
        theta_init = np.array([0.3, 0.4, 0.3])
        phi = np.random.dirichlet(np.ones(3), 50)

        sampler.sample(mixture, theta_init, phi, verbose=False)
        theta_mean, Z_mean = sampler.get_posterior_mean()

        assert theta_mean.shape == (3,)
        assert Z_mean.shape == (50, 3)
        assert np.allclose(theta_mean.sum(), 1.0, atol=0.01)


@pytest.mark.unit
class TestBayesPrism:
    """Test BayesPrism main class"""

    def test_initialization(self, sample_reference, sample_mixture):
        """Test BayesPrism initialization"""
        reference, cell_types, gene_names = sample_reference

        ref = BayesPrismReference(
            reference=reference,
            cell_type_labels=cell_types,
        )

        bp = BayesPrism(reference=ref, mixture=sample_mixture)

        assert bp.reference is ref
        assert bp.mixture.shape == (100, 20)
        assert len(bp.aligned_genes_) == 100

    def test_cleanup_genes(self, bayes_prism_obj):
        """Test gene cleanup"""
        bp = bayes_prism_obj

        # Add some ribosomal genes
        bp.mixture.index = [f'RPS{i}' if i < 10 else f'Gene_{i}' for i in range(100)]
        bp.aligned_genes_ = bp.mixture.index.tolist()

        initial_genes = len(bp.aligned_genes_)
        bp.cleanup_genes(remove_ribo=True, remove_mito=True)

        assert len(bp.aligned_genes_) < initial_genes

    def test_run_deconvolution(self, bayes_prism_obj):
        """Test deconvolution workflow"""
        bp = bayes_prism_obj

        # Use smaller config for testing
        bp.config.burnin = 5
        bp.config.n_iter = 10

        bp.run_deconvolution(n_cores=1, verbose=False, use_numba=False)

        assert bp.theta_initial_ is not None
        assert bp.theta_updated_ is not None
        assert bp.Z_ is not None

    def test_get_fraction(self, bayes_prism_obj):
        """Test getting cell type fractions"""
        bp = bayes_prism_obj
        bp.config.burnin = 5
        bp.config.n_iter = 10

        bp.run_deconvolution(n_cores=1, verbose=False, use_numba=False)
        fractions = bp.get_fraction(updated=True)

        assert fractions.shape == (20, 4)  # samples x cell_types
        assert list(fractions.columns) == bp.reference.cell_types

        # Check proportions sum to ~1
        row_sums = fractions.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=0.1)


@pytest.mark.unit
class TestBayesPrismEmbedding:
    """Test BayesPrismEmbedding"""

    def test_initialization(self, bayes_prism_obj):
        """Test embedding initialization"""
        bp = bayes_prism_obj

        # Need to run deconvolution first
        bp.config.burnin = 5
        bp.config.n_iter = 10
        bp.run_deconvolution(n_cores=1, verbose=False, use_numba=False)

        embedding = BayesPrismEmbedding(
            prism=bp,
            tumor_key='Tumor',
            n_programs=3,
        )

        assert embedding.n_programs == 3
        assert embedding.tumor_key == 'Tumor'

    def test_run_nmf(self, bayes_prism_obj):
        """Test NMF gene program learning"""
        bp = bayes_prism_obj

        bp.config.burnin = 5
        bp.config.n_iter = 10
        bp.run_deconvolution(n_cores=1, verbose=False, use_numba=False)

        embedding = BayesPrismEmbedding(
            prism=bp,
            tumor_key='Tumor',
            n_programs=2,
        )

        embedding.run_nmf(max_iter=50, verbose=False)

        assert embedding.W_ is not None
        assert embedding.H_ is not None
        assert embedding.W_df_ is not None
        assert embedding.H_df_ is not None

    def test_get_top_genes(self, bayes_prism_obj):
        """Test getting top genes for programs"""
        bp = bayes_prism_obj

        bp.config.burnin = 5
        bp.config.n_iter = 10
        bp.run_deconvolution(n_cores=1, verbose=False, use_numba=False)

        embedding = BayesPrismEmbedding(
            prism=bp,
            tumor_key='Tumor',
            n_programs=2,
        )

        embedding.run_nmf(max_iter=50, verbose=False)
        top_genes = embedding.get_top_genes('Program_1', n=10)

        assert len(top_genes) == 10


@pytest.mark.unit
class TestUtils:
    """Test utility functions"""

    def test_cleanup_genes(self):
        """Test gene cleanup"""
        genes = ['RPS14', 'MT-CO1', 'TP53', 'XIST', 'BRCA1', 'RPL5']

        cleaned = cleanup_genes(
            genes,
            remove_ribo=True,
            remove_mito=True,
            remove_sex=True,
        )

        assert 'TP53' in cleaned
        assert 'BRCA1' in cleaned
        assert 'RPS14' not in cleaned
        assert 'MT-CO1' not in cleaned
        assert 'XIST' not in cleaned

    def test_compute_correlation(self):
        """Test correlation computation"""
        predicted = pd.DataFrame({
            'Type_A': [0.2, 0.3, 0.4],
            'Type_B': [0.8, 0.7, 0.6],
        })

        actual = pd.DataFrame({
            'Type_A': [0.25, 0.35, 0.45],
            'Type_B': [0.75, 0.65, 0.55],
        })

        corr = compute_correlation(predicted, actual, method='pearson')

        assert len(corr) == 2
        assert 'pearson_r' in corr.columns

    def test_compute_rmse(self):
        """Test RMSE computation"""
        predicted = pd.DataFrame({
            'Type_A': [0.2, 0.3, 0.4],
            'Type_B': [0.8, 0.7, 0.6],
        })

        actual = pd.DataFrame({
            'Type_A': [0.3, 0.3, 0.3],
            'Type_B': [0.7, 0.7, 0.7],
        })

        rmse = compute_rmse(predicted, actual)

        assert len(rmse) == 2
        assert 'rmse' in rmse.columns


@pytest.mark.integration
class TestFullWorkflow:
    """Integration test for complete workflow"""

    def test_complete_pipeline(self):
        """Test complete BayesPrism pipeline"""
        np.random.seed(42)

        # Generate data
        n_genes = 50
        n_cells = 100
        n_samples = 10

        reference = sparse.random(n_genes, n_cells, density=0.1, format='csr')
        reference.data = np.random.poisson(5, size=reference.data.shape)

        cell_types = pd.Series(
            np.random.choice(['T_cell', 'B_cell', 'Tumor'], size=n_cells)
        )

        mixture = pd.DataFrame(
            np.random.poisson(500, (n_genes, n_samples)),
            index=[f'Gene_{i}' for i in range(n_genes)],
            columns=[f'Sample_{i}' for i in range(n_samples)],
        )

        # Create reference
        ref = BayesPrismReference(
            reference=reference,
            cell_type_labels=cell_types,
        )

        # Run deconvolution
        config = PrismConfig(n_iter=10, burnin=5, n_chains=2)
        bp = BayesPrism(reference=ref, mixture=mixture, config=config)
        bp.run_deconvolution(n_cores=1, verbose=False, use_numba=False)

        # Get results
        fractions = bp.get_fraction()
        expression = bp.get_expression()

        assert fractions.shape == (n_samples, 3)
        assert isinstance(expression, dict)
        assert len(expression) == 3


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
