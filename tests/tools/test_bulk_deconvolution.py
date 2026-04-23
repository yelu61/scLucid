"""
Tests for R-free bulk deconvolution tools.

Tests BayesPrism, DWLS, Bisque, and NNLS implementations.
"""

import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

try:
    from scLucid.tools.bulk import deconvolve_bulk, differential_abundance
    from scLucid.tools.pyBayesPrism import BayesPrismReference, PrismConfig
    from scLucid.tools.pyDWLS import DWLS
except Exception as exc:  # pragma: no cover - optional backend availability
    pytest.skip(f"Skipping bulk deconvolution tests: {exc}", allow_module_level=True)

from tests.fixtures.synthetic_data import generate_minimal_adata


@pytest.fixture
def sc_reference():
    """Generate single-cell reference data for testing."""
    adata = generate_minimal_adata(n_cells=200, n_genes=500)

    # Add cell types
    cell_types = ["T_cell", "B_cell", "Monocyte"]
    adata.obs["cell_type"] = np.random.choice(cell_types, size=adata.n_obs)
    adata.obs["sampleID"] = "sample_1"

    return adata


@pytest.fixture
def bulk_data(sc_reference):
    """Generate synthetic bulk data from sc reference."""
    # Create bulk by summing expression
    cell_types = sc_reference.obs["cell_type"].unique()

    bulk_dict = {}
    for sample in ["bulk_1", "bulk_2", "bulk_3"]:
        # Random proportions
        props = np.random.dirichlet(np.ones(len(cell_types)))

        # Create bulk expression
        bulk_expr = np.zeros(sc_reference.n_vars)
        for ct, prop in zip(cell_types, props):
            mask = sc_reference.obs["cell_type"] == ct
            ct_mean = sc_reference[mask].X.mean(axis=0)
            if hasattr(ct_mean, "flatten"):
                ct_mean = ct_mean.flatten()
            bulk_expr += prop * ct_mean

        bulk_dict[sample] = bulk_expr

    bulk_df = pd.DataFrame(bulk_dict, index=sc_reference.var_names)
    return bulk_df


@pytest.mark.unit
class TestBulkDeconvolution:
    """Test bulk deconvolution methods."""

    def test_deconvolve_bulk_nnls(self, sc_reference, bulk_data):
        """Test NNLS deconvolution."""
        result = deconvolve_bulk(sc_reference, bulk_data, cell_type_key="cell_type", method="NNLS")

        # Check results stored
        assert "bulk_deconvolution" in result.uns["sclucid"]["tools"]

        props = result.uns["sclucid"]["tools"]["bulk_deconvolution"]["proportions"]
        assert isinstance(props, pd.DataFrame)
        assert props.shape == (3, 3)  # 3 samples x 3 cell types

    def test_deconvolve_bulk_dwls(self, sc_reference, bulk_data):
        """Test DWLS deconvolution."""
        result = deconvolve_bulk(sc_reference, bulk_data, cell_type_key="cell_type", method="DWLS")

        props = result.uns["sclucid"]["tools"]["bulk_deconvolution"]["proportions"]

        # Check proportions sum to ~1
        row_sums = props.sum(axis=1)
        np.testing.assert_allclose(row_sums, 1.0, atol=0.1)

    def test_deconvolve_bulk_bisque(self, sc_reference, bulk_data):
        """Test Bisque-like deconvolution."""
        result = deconvolve_bulk(
            sc_reference, bulk_data, cell_type_key="cell_type", method="Bisque"
        )

        props = result.uns["sclucid"]["tools"]["bulk_deconvolution"]["proportions"]
        assert props.shape[0] == 3  # 3 samples

    def test_deconvolve_bulk_bayesprism(self, sc_reference, bulk_data):
        """Test BayesPrism deconvolution."""
        result = deconvolve_bulk(
            sc_reference,
            bulk_data,
            cell_type_key="cell_type",
            method="BayesPrism",
            n_iter=10,  # Small for testing
            n_chains=2,
            burnin=5,
        )

        props = result.uns["sclucid"]["tools"]["bulk_deconvolution"]["proportions"]
        assert props.shape == (3, 3)

    def test_proportions_sum_to_one(self, sc_reference, bulk_data):
        """Test that proportions sum to approximately 1."""
        result = deconvolve_bulk(sc_reference, bulk_data, cell_type_key="cell_type", method="NNLS")

        props = result.uns["sclucid"]["tools"]["bulk_deconvolution"]["proportions"]

        for sample in props.index:
            total = props.loc[sample].sum()
            assert 0.99 <= total <= 1.01, f"Proportions sum to {total}"

    def test_invalid_method_raises_error(self, sc_reference, bulk_data):
        """Test that invalid method raises ValueError."""
        with pytest.raises(ValueError):
            deconvolve_bulk(
                sc_reference, bulk_data, cell_type_key="cell_type", method="invalid_method"
            )

    def test_insufficient_common_genes(self, sc_reference):
        """Test error when too few common genes."""
        # Create bulk data with different genes
        bulk_data_mismatch = pd.DataFrame(
            np.random.rand(10, 2),  # Only 10 genes
            index=[f"gene_{i}" for i in range(1000, 1010)],  # Different gene names
        )

        with pytest.raises(ValueError, match="common genes"):
            deconvolve_bulk(sc_reference, bulk_data_mismatch, cell_type_key="cell_type")


@pytest.mark.unit
class TestDifferentialAbundance:
    """Test differential abundance analysis."""

    def test_differential_abundance_ttest(self):
        """Test t-test differential abundance."""
        # Create synthetic proportions
        np.random.seed(42)

        # Group 1: Higher T_cell proportion
        group1_props = pd.DataFrame(
            {
                "T_cell": np.random.beta(8, 2, 10),
                "B_cell": np.random.beta(2, 8, 10),
                "Monocyte": np.random.beta(2, 2, 10),
            },
            index=[f"sample_{i}" for i in range(10)],
        )

        # Group 2: Lower T_cell proportion
        group2_props = pd.DataFrame(
            {
                "T_cell": np.random.beta(2, 8, 10),
                "B_cell": np.random.beta(8, 2, 10),
                "Monocyte": np.random.beta(2, 2, 10),
            },
            index=[f"sample_{i}" for i in range(10, 20)],
        )

        proportions = pd.concat([group1_props, group2_props])

        metadata = pd.DataFrame(
            {"group": ["control"] * 10 + ["treatment"] * 10}, index=proportions.index
        )

        result = differential_abundance(
            proportions,
            metadata,
            group_col="group",
            group1="control",
            group2="treatment",
            method="ttest",
        )

        # Check results
        assert "T_cell" in result["cell_type"].values
        assert "pvalue" in result.columns

        # T_cell should be significantly different
        t_cell_result = result[result["cell_type"] == "T_cell"].iloc[0]
        assert t_cell_result["pvalue"] < 0.05

    def test_differential_abundance_wilcoxon(self):
        """Test Wilcoxon differential abundance."""
        np.random.seed(42)

        props = pd.DataFrame(
            {"Type_A": np.random.rand(20), "Type_B": np.random.rand(20)},
            index=[f"sample_{i}" for i in range(20)],
        )

        metadata = pd.DataFrame({"condition": ["A"] * 10 + ["B"] * 10}, index=props.index)

        result = differential_abundance(
            props, metadata, group_col="condition", group1="A", group2="B", method="wilcoxon"
        )

        assert len(result) > 0
        assert "pvalue" in result.columns


@pytest.mark.unit
class TestDWLSClass:
    """Test DWLS class directly."""

    def test_dwls_initialization(self):
        """Test DWLS initialization."""
        signature = pd.DataFrame(np.random.rand(100, 3), columns=["Type_A", "Type_B", "Type_C"])

        bulk = pd.DataFrame(np.random.rand(100, 5), columns=[f"sample_{i}" for i in range(5)])

        dwls = DWLS(signature_matrix=signature, bulk_data=bulk)
        assert dwls.signature_matrix is not None
        assert dwls.bulk_data is not None


@pytest.mark.unit
class TestBayesPrismClass:
    """Test BayesPrism class directly."""

    def test_prism_reference_initialization(self):
        """Test BayesPrismReference initialization."""
        reference = pd.DataFrame(
            np.random.poisson(5, (100, 50)),
            index=[f"gene_{i}" for i in range(100)],
            columns=[f"cell_{i}" for i in range(50)],
        )

        cell_types = pd.Series(np.random.choice(["A", "B", "C"], 50), index=reference.columns)

        prism_ref = BayesPrismReference(reference=reference, cell_type_labels=cell_types)

        assert prism_ref.cell_types == ["A", "B", "C"]
        assert prism_ref.phi is not None

    def test_prism_config(self):
        """Test PrismConfig."""
        config = PrismConfig(n_iter=200, n_chains=8)
        assert config.n_iter == 200
        assert config.n_chains == 8
        assert "chain_length" in config.gibbs_control


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
