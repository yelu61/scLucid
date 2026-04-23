"""Tests for differential expression core functions."""

import pandas as pd
import pytest
import scanpy as sc
from anndata import AnnData

from scLucid.analysis.config import (
    CompareConditionsConfig,
    CompareGroupsConfig,
    DifferentialConfig,
    FilterMarkersConfig,
)
from scLucid.analysis.differential_expression.de_core import (
    compare_conditions,
    compare_groups,
    filter_markers,
    find_markers,
)


def _preprocess_for_de(adata: AnnData) -> AnnData:
    """Minimal preprocessing for DE tests."""
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata


@pytest.fixture
def de_adata(minimal_adata):
    """Preprocessed AnnData with cell types for DE testing."""
    adata = _preprocess_for_de(minimal_adata)
    return adata


class TestFindMarkers:
    """Tests for find_markers."""

    def test_basic_find_markers(self, de_adata):
        """find_markers runs and returns a DataFrame."""
        config = DifferentialConfig(groupby="cell_type", method="wilcoxon")
        df = find_markers(de_adata, config)
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        assert "group" in df.columns
        assert "names" in df.columns
        assert "logfoldchanges" in df.columns

    def test_kwargs_override_config(self, de_adata):
        """Kwargs override config values."""
        config = DifferentialConfig(groupby="cell_type", method="wilcoxon")
        df = find_markers(de_adata, config, method="t-test")
        # Verify the override was applied by checking stored params
        stored = de_adata.uns["sclucid"]["analysis"]["de"]["rank_genes_groups_params"]
        assert stored["method"] == "t-test"

    def test_stores_in_uns(self, de_adata):
        """Results are stored in adata.uns under sclucid path."""
        config = DifferentialConfig(groupby="cell_type", key_added="my_markers")
        find_markers(de_adata, config)
        root = de_adata.uns["sclucid"]["analysis"]["de"]
        assert "my_markers_df" in root
        assert "my_markers" in root
        assert "my_markers_params" in root

    def test_pval_cutoff_filtering(self, de_adata):
        """pval_cutoff removes rows above threshold."""
        config = DifferentialConfig(groupby="cell_type", pval_cutoff=1e-20)
        df = find_markers(de_adata, config)
        assert (df["pvals_adj"] <= 1e-20).all()

    def test_missing_groupby_raises(self, de_adata):
        """KeyError when groupby column doesn't exist."""
        config = DifferentialConfig(groupby="nonexistent")
        with pytest.raises(KeyError, match="nonexistent"):
            find_markers(de_adata, config)

    def test_config_none_uses_defaults(self, de_adata):
        """Passing config=None uses default DifferentialConfig."""
        # Need to set groupby on adata if using defaults
        de_adata.obs["leiden_clusters"] = de_adata.obs["cell_type"]
        df = find_markers(de_adata, config=None)
        assert isinstance(df, pd.DataFrame)


class TestFilterMarkers:
    """Tests for filter_markers."""

    @pytest.fixture
    def marked_adata(self, de_adata):
        """AnnData with find_markers already run."""
        config = DifferentialConfig(groupby="cell_type", method="wilcoxon")
        find_markers(de_adata, config)
        return de_adata

    def test_basic_filter(self, marked_adata):
        """filter_markers returns a DataFrame."""
        config = FilterMarkersConfig(key="rank_genes_groups", min_log2fc=0.5)
        df = filter_markers(marked_adata, config)
        assert isinstance(df, pd.DataFrame)

    def test_min_log2fc_filter(self, marked_adata):
        """min_log2fc removes low fold-change genes."""
        config = FilterMarkersConfig(key="rank_genes_groups", min_log2fc=2.0)
        df = filter_markers(marked_adata, config)
        assert (df["logfoldchanges"] >= 2.0).all()

    def test_max_padj_filter(self, marked_adata):
        """max_padj removes high p-value genes."""
        config = FilterMarkersConfig(key="rank_genes_groups", max_padj=0.01, min_log2fc=0.0)
        df = filter_markers(marked_adata, config)
        assert (df["pvals_adj"] <= 0.01).all()

    def test_keep_top_n(self, marked_adata):
        """keep_top_n limits rows per group."""
        config = FilterMarkersConfig(
            key="rank_genes_groups", keep_top_n=5, min_log2fc=0.0, max_padj=1.0
        )
        df = filter_markers(marked_adata, config)
        for group in df["group"].unique():
            assert len(df[df["group"] == group]) <= 5

    def test_missing_source_raises(self, de_adata):
        """KeyError when source DE results not found."""
        config = FilterMarkersConfig(key="missing_key")
        with pytest.raises(KeyError, match="missing_key_df"):
            filter_markers(de_adata, config)

    def test_empty_df_returns_empty(self, marked_adata):
        """Empty source returns empty DataFrame without error."""
        # Replace stored df with empty one
        root = marked_adata.uns["sclucid"]["analysis"]["de"]
        root["rank_genes_groups_df"] = pd.DataFrame()
        config = FilterMarkersConfig(key="rank_genes_groups")
        df = filter_markers(marked_adata, config)
        assert df.empty

    def test_abs_log2fc(self, marked_adata):
        """use_abs_log2fc keeps both up and down regulated genes."""
        config = FilterMarkersConfig(
            key="rank_genes_groups",
            use_abs_log2fc=True,
            min_log2fc=1.0,
            max_padj=1.0,
        )
        df = filter_markers(marked_adata, config)
        assert (df["logfoldchanges"].abs() >= 1.0).all()


class TestCompareGroups:
    """Tests for compare_groups."""

    def test_basic_comparison(self, de_adata):
        """compare_groups returns filtered DE genes."""
        cell_types = de_adata.obs["cell_type"].unique()
        config = CompareGroupsConfig(
            groupby="cell_type",
            group1=cell_types[0],
            group2=cell_types[1],
            min_log2fc=0.0,
            max_padj=1.0,
            n_top_genes=10,
        )
        df = compare_groups(de_adata, config)
        assert isinstance(df, pd.DataFrame)

    def test_missing_groupby_raises(self, de_adata):
        """KeyError when groupby column doesn't exist."""
        config = CompareGroupsConfig(groupby="missing", group1="A", group2="B")
        with pytest.raises(KeyError, match="missing"):
            compare_groups(de_adata, config)

    def test_missing_groups_raises(self, de_adata):
        """ValueError when neither group exists."""
        config = CompareGroupsConfig(
            groupby="cell_type", group1="NonExistent", group2="AlsoMissing"
        )
        with pytest.raises(ValueError, match="No cells found"):
            compare_groups(de_adata, config)

    def test_stores_results(self, de_adata):
        """Results stored in adata.uns."""
        cell_types = de_adata.obs["cell_type"].unique()
        config = CompareGroupsConfig(
            groupby="cell_type",
            group1=cell_types[0],
            group2=cell_types[1],
            min_log2fc=0.0,
            max_padj=1.0,
            key_added="my_comparison",
        )
        compare_groups(de_adata, config)
        root = de_adata.uns["sclucid"]["analysis"]["de"]
        assert "my_comparison" in root.keys()


class TestCompareConditions:
    """Tests for compare_conditions."""

    def test_basic_condition_comparison(self, de_adata):
        """compare_conditions within a cell type."""
        # Add fake condition
        de_adata.obs["condition"] = ["ctrl", "treat"] * (de_adata.n_obs // 2)
        cell_type = de_adata.obs["cell_type"].unique()[0]
        config = CompareConditionsConfig(
            groupby="cell_type",
            group_name=cell_type,
            condition_key="condition",
            condition1="ctrl",
            condition2="treat",
            min_log2fc=0.0,
            max_padj=1.0,
        )
        df = compare_conditions(de_adata, config)
        assert isinstance(df, pd.DataFrame)

    def test_missing_group_raises(self, de_adata):
        """ValueError when group not found."""
        de_adata.obs["condition"] = ["ctrl", "treat"] * (de_adata.n_obs // 2)
        config = CompareConditionsConfig(
            groupby="cell_type",
            group_name="Missing",
            condition_key="condition",
            condition1="ctrl",
            condition2="treat",
        )
        with pytest.raises(ValueError, match="Missing"):
            compare_conditions(de_adata, config)
