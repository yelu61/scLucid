"""
Tests for pyCellChat (R-free CellChat implementation)
"""

import sys

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

sys.path.insert(0, "/Users/luye/Scripts/scLucid/src")

from scLucid.tools.pyCellChat import (
    CellChat,
    CellChatDB,
    create_cellchat_from_scanpy,
    get_default_database,
    plot_heatmap,
)


@pytest.fixture
def sample_data():
    """Generate sample data for testing"""
    np.random.seed(42)
    n_genes = 200
    n_cells = 100

    # Generate expression data
    expr = np.random.poisson(5, (n_genes, n_cells))
    db_genes = sorted(list(CellChatDB(species="human").get_all_genes()))
    if len(db_genes) >= n_genes:
        gene_names = db_genes[:n_genes]
    else:
        extra = [f"Gene_{i}" for i in range(n_genes - len(db_genes))]
        gene_names = db_genes + extra
    cell_names = [f"Cell_{i}" for i in range(n_cells)]

    expr_df = pd.DataFrame(expr, index=gene_names, columns=cell_names)

    # Create metadata
    cell_types = np.random.choice(["T_cell", "B_cell", "Monocyte"], n_cells)
    meta = pd.DataFrame({"cell_type": cell_types}, index=cell_names)

    return expr_df, meta


@pytest.fixture
def sample_anndata():
    """Generate sample AnnData for testing"""
    np.random.seed(42)
    n_genes = 200
    n_cells = 100

    X = np.random.poisson(5, (n_cells, n_genes))
    db_genes = sorted(list(CellChatDB(species="human").get_all_genes()))
    if len(db_genes) >= n_genes:
        var_names = db_genes[:n_genes]
    else:
        extra = [f"Gene_{i}" for i in range(n_genes - len(db_genes))]
        var_names = db_genes + extra
    obs = pd.DataFrame({"cell_type": np.random.choice(["T_cell", "B_cell", "Monocyte"], n_cells)})
    var = pd.DataFrame(index=var_names)

    adata = AnnData(X=X, obs=obs, var=var)
    return adata


@pytest.mark.unit
class TestCellChatDB:
    """Test CellChatDB database"""

    def test_db_initialization(self):
        """Test database initialization"""
        db = CellChatDB(species="human")
        assert db.species == "human"
        assert len(db.interaction) > 0

    def test_db_get_genes(self):
        """Test getting genes from database"""
        db = CellChatDB(species="human")
        genes = db.get_all_genes()
        assert len(genes) > 0
        assert "TNF" in genes or "TGFB1" in genes

    def test_db_get_pathways(self):
        """Test getting pathways"""
        db = CellChatDB(species="human")
        pathways = db.get_pathways()
        assert len(pathways) > 0
        assert "TGFB" in pathways or "TNF" in pathways

    def test_db_subset(self):
        """Test subsetting database"""
        db = CellChatDB(species="human")
        subset = db.subset_db(pathways=["TGFB", "TNF"])
        assert len(subset.interaction) > 0
        assert all(p in ["TGFB", "TNF"] for p in subset.interaction["pathway_name"])

    def test_db_mouse(self):
        """Test mouse database"""
        db = CellChatDB(species="mouse")
        assert db.species == "mouse"
        assert len(db.interaction) > 0


@pytest.mark.unit
class TestCellChatCore:
    """Test CellChat core functionality"""

    def test_cellchat_initialization(self, sample_data):
        """Test CellChat initialization"""
        expr_df, meta = sample_data
        cellchat = CellChat(data=expr_df, meta=meta, group_by="cell_type")

        assert cellchat.n_cells == 100
        assert cellchat.n_genes == 200
        assert len(cellchat.unique_groups) == 3

    def test_set_database(self, sample_data):
        """Test setting database"""
        expr_df, meta = sample_data
        cellchat = CellChat(data=expr_df, meta=meta, group_by="cell_type")

        db = CellChatDB(species="human")
        cellchat.set_database(db)

        assert cellchat.db is not None
        assert cellchat.LR is not None

    def test_preprocess(self, sample_data):
        """Test data preprocessing"""
        expr_df, meta = sample_data
        cellchat = CellChat(data=expr_df, meta=meta, group_by="cell_type")

        db = CellChatDB(species="human")
        cellchat.set_database(db)

        # Should work without error
        cellchat.preprocess_data()

    def test_identify_overexpressed_genes(self, sample_data):
        """Test identifying overexpressed genes"""
        expr_df, meta = sample_data
        cellchat = CellChat(data=expr_df, meta=meta, group_by="cell_type")

        db = CellChatDB(species="human")
        cellchat.set_database(db)
        cellchat.preprocess_data()

        group_means = cellchat.identify_overexpressed_genes()
        assert group_means.shape[1] == 3  # 3 cell types

    def test_compute_communication_prob(self, sample_data):
        """Test computing communication probability"""
        expr_df, meta = sample_data
        cellchat = CellChat(data=expr_df, meta=meta, group_by="cell_type")

        db = CellChatDB(species="human")
        cellchat.set_database(db)
        cellchat.preprocess_data()

        prob = cellchat.compute_communication_prob()
        assert prob.shape[0] == len(db.interaction)
        assert prob.shape[1] == 3  # 3 cell types
        assert prob.shape[2] == 3  # 3 cell types

    def test_save_load(self, sample_data, tmp_path):
        """Test saving and loading"""
        expr_df, meta = sample_data
        cellchat = CellChat(data=expr_df, meta=meta, group_by="cell_type")

        # Save
        filepath = tmp_path / "test_cellchat.pkl"
        cellchat.save(str(filepath))
        assert filepath.exists()

        # Load
        loaded = CellChat.load(str(filepath))
        assert loaded.n_cells == cellchat.n_cells
        assert loaded.n_genes == cellchat.n_genes


@pytest.mark.unit
class TestCellChatUtils:
    """Test utility functions"""

    def test_create_from_scanpy(self, sample_anndata):
        """Test creating CellChat from AnnData"""
        cellchat = create_cellchat_from_scanpy(sample_anndata, group_by="cell_type")

        assert cellchat.n_cells == sample_anndata.n_obs
        assert cellchat.n_genes == sample_anndata.n_vars

    def test_get_default_database(self):
        """Test getting default database"""
        db = get_default_database(species="human")
        assert isinstance(db, CellChatDB)
        assert len(db.interaction) > 0


@pytest.mark.unit
class TestCellChatVisualization:
    """Test visualization functions"""

    def test_plot_heatmap(self, sample_data):
        """Test heatmap plotting"""
        expr_df, meta = sample_data
        cellchat = CellChat(data=expr_df, meta=meta, group_by="cell_type")

        db = CellChatDB(species="human")
        cellchat.set_database(db)
        cellchat.preprocess_data()
        cellchat.compute_communication_prob()

        fig, ax = plot_heatmap(cellchat)
        assert fig is not None
        assert ax is not None


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
