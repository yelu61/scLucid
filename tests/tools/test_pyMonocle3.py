"""
Tests for pyMonocle3 (R-free Monocle3 implementation)
"""

import pytest
import numpy as np
import pandas as pd
import scipy.sparse as sp
import sys

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

from scLucid.tools.pyMonocle3 import (
    CellDataSet,
    new_cell_data_set,
    create_cds_from_scanpy,
    export_to_scanpy,
    detect_genes,
    estimate_size_factors,
    preprocess_cds,
    reduce_dimension,
    run_pca,
    run_umap,
    cluster_cells,
    partition_cells,
    learn_graph,
    order_cells,
    graph_test,
    top_markers,
    compare_genes,
    plot_cells,
    plot_trajectory,
    validate_cds,
)


@pytest.fixture
def sample_cds():
    """Generate sample CellDataSet for testing"""
    np.random.seed(42)
    n_genes = 100
    n_cells = 50

    # Generate expression data
    expr = np.random.poisson(5, (n_genes, n_cells))
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    cell_names = [f"cell_{i}" for i in range(n_cells)]

    # Create metadata
    cell_meta = pd.DataFrame({
        'cell_type': np.random.choice(['A', 'B', 'C'], n_cells),
        'batch': np.random.choice(['batch1', 'batch2'], n_cells),
    }, index=cell_names)

    gene_meta = pd.DataFrame({
        'gene_short_name': gene_names,
    }, index=gene_names)

    return CellDataSet(
        expression_data=expr,
        cell_metadata=cell_meta,
        gene_metadata=gene_meta,
    )


@pytest.fixture
def preprocessed_cds(sample_cds):
    """Preprocessed CellDataSet"""
    return preprocess_cds(sample_cds, num_dim=20)


@pytest.mark.unit
class TestCellDataSet:
    """Test CellDataSet core functionality"""

    def test_initialization(self):
        """Test CellDataSet creation"""
        expr = np.random.poisson(5, (50, 20))
        cell_meta = pd.DataFrame(index=[f"cell_{i}" for i in range(20)])
        gene_meta = pd.DataFrame(index=[f"gene_{i}" for i in range(50)])

        cds = CellDataSet(
            expression_data=expr,
            cell_metadata=cell_meta,
            gene_metadata=gene_meta,
        )

        assert cds.n_genes == 50
        assert cds.n_cells == 20
        assert len(cds.reducedDims) == 0

    def test_validation_error(self):
        """Test validation catches mismatched dimensions"""
        expr = np.random.poisson(5, (50, 20))
        # Wrong number of cells in metadata
        cell_meta = pd.DataFrame(index=[f"cell_{i}" for i in range(10)])
        gene_meta = pd.DataFrame(index=[f"gene_{i}" for i in range(50)])

        with pytest.raises(ValueError):
            CellDataSet(
                expression_data=expr,
                cell_metadata=cell_meta,
                gene_metadata=gene_meta,
            )

    def test_copy(self, sample_cds):
        """Test CellDataSet copy"""
        cds_copy = sample_cds.copy()

        assert cds_copy.n_genes == sample_cds.n_genes
        assert cds_copy.n_cells == sample_cds.n_cells
        assert cds_copy.expression_data is not sample_cds.expression_data

    def test_save_load(self, sample_cds, tmp_path):
        """Test save and load"""
        filepath = tmp_path / "test_cds.pkl"
        sample_cds.save(str(filepath))

        assert filepath.exists()

        loaded = CellDataSet.load(str(filepath))
        assert loaded.n_genes == sample_cds.n_genes
        assert loaded.n_cells == sample_cds.n_cells


@pytest.mark.unit
class TestNewCellDataSet:
    """Test new_cell_data_set factory function"""

    def test_from_matrix(self):
        """Test creating from matrix"""
        expr = np.random.poisson(5, (50, 20))
        cds = new_cell_data_set(expr)

        assert cds.n_genes == 50
        assert cds.n_cells == 20

    def test_from_dataframe(self):
        """Test creating from DataFrame"""
        expr = pd.DataFrame(
            np.random.poisson(5, (50, 20)),
            index=[f"gene_{i}" for i in range(50)],
            columns=[f"cell_{i}" for i in range(20)],
        )
        cds = new_cell_data_set(expr)

        assert cds.n_genes == 50
        assert cds.n_cells == 20


@pytest.mark.unit
class TestPreprocessing:
    """Test preprocessing functions"""

    def test_detect_genes(self, sample_cds):
        """Test gene detection"""
        cds = detect_genes(sample_cds, min_expr=0.1, min_cells=5)

        assert 'mean_expression' in cds.gene_metadata.columns
        assert 'num_cells_expressed' in cds.gene_metadata.columns
        assert 'use_for_ordering' in cds.gene_metadata.columns

    def test_estimate_size_factors(self, sample_cds):
        """Test size factor estimation"""
        cds = estimate_size_factors(sample_cds)

        assert 'Size_Factor' in cds.cell_metadata.columns
        assert len(cds.cell_metadata['Size_Factor']) == cds.n_cells
        assert all(cds.cell_metadata['Size_Factor'] > 0)

    def test_preprocess_cds(self, sample_cds):
        """Test full preprocessing pipeline"""
        cds = preprocess_cds(sample_cds, num_dim=10)

        assert 'PCA' in cds.reducedDims
        assert cds.reducedDims['PCA'].shape[0] == cds.n_cells
        assert cds.reducedDims['PCA'].shape[1] <= 10
        assert len(cds.preprocessing_params) > 0


@pytest.mark.unit
class TestDimensionality:
    """Test dimensionality reduction"""

    def test_reduce_dimension_umap(self, preprocessed_cds):
        """Test UMAP reduction"""
        pytest.importorskip("umap")

        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')

        assert 'UMAP' in cds.reducedDims
        assert cds.reducedDims['UMAP'].shape == (cds.n_cells, 2)

    def test_reduce_dimension_tsne(self, preprocessed_cds):
        """Test tSNE reduction"""
        cds = reduce_dimension(preprocessed_cds, reduction_method='tSNE')

        assert 'tSNE' in cds.reducedDims
        assert cds.reducedDims['tSNE'].shape == (cds.n_cells, 2)

    def test_run_pca(self, sample_cds):
        """Test PCA"""
        cds = run_pca(sample_cds, n_components=10)

        assert 'PCA' in cds.reducedDims
        assert cds.reducedDims['PCA'].shape[1] <= 10

    def test_run_umap(self, preprocessed_cds):
        """Test standalone UMAP"""
        pytest.importorskip("umap")

        cds = run_umap(preprocessed_cds, reduction_key='PCA')

        assert 'UMAP_2D' in cds.reducedDims
        assert cds.reducedDims['UMAP_2D'].shape == (cds.n_cells, 2)


@pytest.mark.unit
class TestClustering:
    """Test clustering functions"""

    def test_cluster_cells_leiden(self, preprocessed_cds):
        """Test Leiden clustering"""
        pytest.importorskip("leidenalg")
        pytest.importorskip("igraph")

        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')
        cds = cluster_cells(cds, cluster_method='leiden')

        assert cds.clusters is not None
        assert 'cluster' in cds.cell_metadata.columns

    def test_cluster_cells_louvain(self, preprocessed_cds):
        """Test Louvain clustering"""
        pytest.importorskip("louvain")
        pytest.importorskip("igraph")

        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')
        cds = cluster_cells(cds, cluster_method='louvain')

        assert cds.clusters is not None
        assert 'cluster' in cds.cell_metadata.columns

    def test_partition_cells(self, preprocessed_cds):
        """Test cell partitioning"""
        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')
        cds = cluster_cells(cds)
        cds = partition_cells(cds)

        assert cds.partitions is not None
        assert 'partition' in cds.cell_metadata.columns


@pytest.mark.unit
class TestTrajectory:
    """Test trajectory inference"""

    def test_learn_graph(self, preprocessed_cds):
        """Test graph learning"""
        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')
        cds = cluster_cells(cds)
        cds = partition_cells(cds)
        cds = learn_graph(cds)

        assert cds.principal_graph is not None
        assert 'adj_matrix' in cds.principal_graph

    def test_order_cells(self, preprocessed_cds):
        """Test pseudotime calculation"""
        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')
        cds = cluster_cells(cds)
        cds = partition_cells(cds)
        cds = learn_graph(cds)
        cds = order_cells(cds)

        assert 'pseudotime' in cds.cell_metadata.columns

    def test_graph_test(self, preprocessed_cds):
        """Test graph-based differential expression"""
        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')
        cds = cluster_cells(cds)
        cds = partition_cells(cds)
        cds = learn_graph(cds)

        results = graph_test(cds)

        assert isinstance(results, pd.DataFrame)
        assert 'gene' in results.columns
        assert 'morans_i' in results.columns


@pytest.mark.unit
class TestDifferential:
    """Test differential expression"""

    def test_top_markers(self, preprocessed_cds):
        """Test marker gene finding"""
        cds = reduce_dimension(preprocessed_cds, reduction_method='PCA')
        cds = cluster_cells(cds, reduction_method='PCA')

        markers = top_markers(cds)

        assert isinstance(markers, pd.DataFrame)
        assert len(markers) > 0

    def test_compare_genes(self, preprocessed_cds):
        """Test gene comparison"""
        cell_list1 = preprocessed_cds.cell_metadata.index[:10].tolist()
        cell_list2 = preprocessed_cds.cell_metadata.index[10:20].tolist()

        results = compare_genes(preprocessed_cds, cell_list1, cell_list2)

        assert isinstance(results, pd.DataFrame)
        assert 'gene' in results.columns
        assert 'log2fc' in results.columns


@pytest.mark.unit
class TestVisualization:
    """Test visualization functions"""

    def test_plot_cells(self, preprocessed_cds):
        """Test cell plotting"""
        pytest.importorskip("matplotlib")

        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')
        cds = cluster_cells(cds)

        fig, ax = plot_cells(cds, color_cells_by='cluster')

        assert fig is not None
        assert ax is not None

    def test_plot_trajectory(self, preprocessed_cds):
        """Test trajectory plotting"""
        pytest.importorskip("matplotlib")

        cds = reduce_dimension(preprocessed_cds, reduction_method='UMAP')
        cds = cluster_cells(cds)
        cds = partition_cells(cds)
        cds = learn_graph(cds)
        cds = order_cells(cds)

        fig, ax = plot_trajectory(cds)

        assert fig is not None
        assert ax is not None


@pytest.mark.unit
class TestUtils:
    """Test utility functions"""

    def test_validate_cds(self, sample_cds):
        """Test CDS validation"""
        is_valid, message = validate_cds(sample_cds)

        assert is_valid
        assert message == "Valid"


@pytest.mark.integration
class TestFullWorkflow:
    """Integration test for full workflow"""

    def test_complete_pipeline(self, sample_cds):
        """Test complete Monocle3 workflow"""
        pytest.importorskip("umap")
        pytest.importorskip("leidenalg")

        cds = sample_cds

        # Preprocessing
        cds = preprocess_cds(cds, num_dim=10)
        assert 'PCA' in cds.reducedDims

        # Dimensionality reduction
        cds = reduce_dimension(cds, reduction_method='UMAP')
        assert 'UMAP' in cds.reducedDims

        # Clustering
        cds = cluster_cells(cds)
        assert cds.clusters is not None

        # Trajectory
        cds = learn_graph(cds)
        assert cds.principal_graph is not None

        # Pseudotime
        cds = order_cells(cds)
        assert 'pseudotime' in cds.cell_metadata.columns

        # Validation
        is_valid, _ = validate_cds(cds)
        assert is_valid


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
