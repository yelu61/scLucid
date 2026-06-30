"""Tests for differential expression core functions."""

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from anndata import AnnData

from scLucid.analysis.config import (
    CompareConditionsConfig,
    CompareGroupsConfig,
    DifferentialConfig,
    FilterMarkersConfig,
    PseudobulkDEConfig,
)
from scLucid.analysis.differential_expression import de_core
from scLucid.analysis.differential_expression.de_core import (
    compare_conditions,
    compare_groups,
    filter_markers,
    find_markers,
    run_mixedlm_de,
    run_pseudobulk_de,
)
from scLucid.analysis.differential_expression.de_utils import _hierarchical_fdr_correction


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
        assert set(df["inference_level"]) == {"cell_level_marker_discovery"}
        assert set(df["claim_level"]) == {"exploratory_marker_screen"}
        assert set(df["recommended_formal_inference_api"]) == {"run_pseudobulk_de"}
        assert not df["valid_for_publication_inference"].any()
        assert df["pseudoreplication_warning"].all()

    def test_kwargs_override_config(self, de_adata):
        """Kwargs override config values."""
        config = DifferentialConfig(groupby="cell_type", method="wilcoxon")
        find_markers(de_adata, config, method="t-test")
        # Verify the override was applied by checking stored params
        stored = de_adata.uns["sclucid"]["analysis"]["de"]["rank_genes_groups_params"]
        assert stored["method"] == "t-test"
        assert stored["claim_level"] == "exploratory_marker_screen"
        assert stored["valid_for_publication_inference"] is False

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
        if not df.empty:
            assert "claim_level" in df.columns

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
            acknowledge_exploratory=True,
        )
        df = compare_groups(de_adata, config)
        assert isinstance(df, pd.DataFrame)
        if not df.empty:
            assert set(df["inference_level"]) == {"exploratory_cell_level"}
            assert set(df["claim_level"]) == {"exploratory_hypothesis_generation"}
            assert not df["valid_for_publication_inference"].any()
            assert df["pseudoreplication_warning"].all()

    def test_acknowledge_exploratory_required(self, de_adata):
        """compare_groups raises if exploratory use is not acknowledged."""
        cell_types = de_adata.obs["cell_type"].unique()
        config = CompareGroupsConfig(
            groupby="cell_type",
            group1=cell_types[0],
            group2=cell_types[1],
            min_log2fc=0.0,
            max_padj=1.0,
            n_top_genes=10,
        )
        with pytest.raises(ValueError, match="acknowledge_exploratory"):
            compare_groups(de_adata, config)

    def test_missing_groupby_raises(self, de_adata):
        """KeyError when groupby column doesn't exist."""
        config = CompareGroupsConfig(
            groupby="missing", group1="A", group2="B", acknowledge_exploratory=True
        )
        with pytest.raises(KeyError, match="missing"):
            compare_groups(de_adata, config)

    def test_missing_groups_raises(self, de_adata):
        """ValueError when neither group exists."""
        config = CompareGroupsConfig(
            groupby="cell_type",
            group1="NonExistent",
            group2="AlsoMissing",
            acknowledge_exploratory=True,
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
            acknowledge_exploratory=True,
        )
        compare_groups(de_adata, config)
        root = de_adata.uns["sclucid"]["analysis"]["de"]
        assert "my_comparison" in root.keys()
        assert root["my_comparison_params"]["inference_level"] == "exploratory_cell_level"
        assert root["my_comparison_params"]["claim_level"] == "exploratory_hypothesis_generation"
        assert root["my_comparison_params"]["valid_for_publication_inference"] is False


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
            acknowledge_exploratory=True,
        )
        df = compare_conditions(de_adata, config)
        assert isinstance(df, pd.DataFrame)
        if not df.empty:
            assert set(df["inference_level"]) == {"exploratory_cell_level"}
            assert set(df["claim_level"]) == {"exploratory_hypothesis_generation"}
            assert not df["valid_for_publication_inference"].any()
            assert df["pseudoreplication_warning"].all()

    def test_acknowledge_exploratory_required(self, de_adata):
        """compare_conditions raises if exploratory use is not acknowledged."""
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
        with pytest.raises(ValueError, match="acknowledge_exploratory"):
            compare_conditions(de_adata, config)

    def test_missing_group_raises(self, de_adata):
        """ValueError when group not found."""
        de_adata.obs["condition"] = ["ctrl", "treat"] * (de_adata.n_obs // 2)
        config = CompareConditionsConfig(
            groupby="cell_type",
            group_name="Missing",
            condition_key="condition",
            condition1="ctrl",
            condition2="treat",
            acknowledge_exploratory=True,
        )
        with pytest.raises(ValueError, match="Missing"):
            compare_conditions(de_adata, config)


class TestRunPseudobulkDE:
    """Tests for sample-level pseudobulk DE."""

    def _make_pseudobulk_adata(self, n_reps=3):
        import numpy as np

        rows = []
        obs = []
        var_names = ["GeneA", "GeneB", "GeneC"]
        for cell_type in ["T", "B"]:
            for condition in ["ctrl", "treat"]:
                for rep in range(n_reps):
                    sample = f"{cell_type}_{condition}_{rep}"
                    for _ in range(4):
                        if cell_type == "T" and condition == "treat":
                            counts = [80, 10, 20]
                        elif cell_type == "T":
                            counts = [10, 80, 20]
                        else:
                            counts = [20, 20, 40]
                        rows.append(counts)
                        obs.append(
                            {
                                "sample": sample,
                                "condition": condition,
                                "cell_type": cell_type,
                            }
                        )
        adata = AnnData(X=np.asarray(rows, dtype=float), obs=pd.DataFrame(obs))
        adata.var_names = var_names
        return adata

    def test_pseudobulk_de_multiple_groups_and_contrasts(self):
        adata = self._make_pseudobulk_adata(n_reps=3)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T", "B"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            key_added="pb_de",
        )

        df = run_pseudobulk_de(adata, config)

        assert isinstance(df, pd.DataFrame)
        assert {"names", "logfoldchanges", "pvals_adj", "group", "contrast"}.issubset(df.columns)
        assert set(df["group"]) == {"T", "B"}
        t_gene_a = df[(df["group"] == "T") & (df["names"] == "GeneA")].iloc[0]
        assert t_gene_a["logfoldchanges"] > 0
        assert t_gene_a["direction"] == "treat - ctrl"
        assert t_gene_a["n_samples_condition1"] == 3
        assert t_gene_a["inference_level"] == "sample_level"
        assert t_gene_a["claim_level"] == "replicate_aware_sample_level_condition_inference"
        assert bool(t_gene_a["valid_for_publication_inference"]) is True
        assert t_gene_a["replicate_status"] == "replicated"
        assert adata.uns["sclucid"]["analysis"]["de"]["pb_de"].equals(df)

    def test_pseudobulk_de_n2_uses_welch_logcpm(self):
        adata = self._make_pseudobulk_adata(n_reps=2)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="welch_logcpm",
        )

        df = run_pseudobulk_de(adata, config)

        assert not df.empty
        assert set(df["method"]) == {"welch_logcpm_n2"}

    def test_pseudobulk_de_welch_excludes_zero_library_samples(self):
        import numpy as np

        adata = self._make_pseudobulk_adata(n_reps=2)
        zero_sample = "T_ctrl_0"
        adata.X[adata.obs["sample"].to_numpy() == zero_sample, :] = 0
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="welch_logcpm",
        )

        df = run_pseudobulk_de(adata, config)

        assert not df.empty
        assert np.isfinite(df["logfoldchanges"]).all()
        assert set(df["n_samples_condition1"]) == {1}
        assert set(df["n_samples_condition2"]) == {2}

    def test_pseudobulk_de_auto_prefers_deseq2(self, monkeypatch):
        adata = self._make_pseudobulk_adata(n_reps=3)

        def fake_deseq2(*args, **kwargs):
            return pd.DataFrame(
                {
                    "names": ["GeneA"],
                    "gene": ["GeneA"],
                    "logfoldchanges": [1.0],
                    "log2fc": [1.0],
                    "scores": [3.0],
                    "statistic": [3.0],
                    "pvals": [0.01],
                    "pval": [0.01],
                    "pvals_adj": [0.02],
                    "padj": [0.02],
                    "method": ["deseq2"],
                }
            )

        monkeypatch.setattr(de_core, "_run_pydeseq2_de", fake_deseq2)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
        )

        df = run_pseudobulk_de(adata, config)

        assert not df.empty
        assert set(df["method"]) == {"deseq2"}

    def test_pseudobulk_de_single_replicate_returns_descriptive_pseudobulk(self):
        adata = self._make_pseudobulk_adata(n_reps=1)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
        )

        df = run_pseudobulk_de(adata, config)

        assert not df.empty
        assert set(df["method"]) == {"descriptive_pseudobulk"}
        assert set(df["inference_level"]) == {"descriptive_single_sample"}
        assert set(df["claim_level"]) == {"descriptive_effect_size_only"}
        assert not df["valid_for_publication_inference"].any()
        assert df["pvals"].isna().all()
        assert df["pvals_adj"].isna().all()
        assert "pseudobulk_warning" in df.columns

    def test_pseudobulk_de_forced_welch_skips_single_replicate(self):
        adata = self._make_pseudobulk_adata(n_reps=1)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="welch_logcpm",
        )

        df = run_pseudobulk_de(adata, config)

        assert df.empty

    def test_pseudobulk_de_can_force_cell_level_fallback(self):
        adata = self._make_pseudobulk_adata(n_reps=3)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="cell_level_fallback",
        )

        df = run_pseudobulk_de(adata, config)

        assert not df.empty
        assert set(df["method"]) == {"cell_level_fallback"}
        assert set(df["inference_level"]) == {"exploratory_cell_level"}
        assert set(df["claim_level"]) == {"exploratory_hypothesis_generation"}
        assert df["pseudoreplication_warning"].all()
        assert not df["valid_for_publication_inference"].any()
        assert df["pseudobulk_warning"].str.contains("forced").all()

    def test_pseudobulk_de_parallel_multiple_contrasts(self):
        adata = self._make_pseudobulk_adata(n_reps=3)
        adata.obs["condition3"] = adata.obs["condition"].replace({"ctrl": "A", "treat": "B"})
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition3",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("A", "B"), ("B", "A")],
            min_cells_per_sample=1,
            min_counts=0,
            n_jobs=2,
        )

        df = run_pseudobulk_de(adata, config)

        assert set(df["contrast"]) == {"B_vs_A", "A_vs_B"}

    def test_pseudobulk_de_linear_model_includes_batch_covariate(self):
        adata = self._make_pseudobulk_adata(n_reps=3)
        sample_meta = adata.obs[["sample", "condition"]].drop_duplicates().copy()
        sample_meta["batch"] = ["b1", "b2", "b1", "b2", "b1", "b2"] * 2
        batch_map = dict(zip(sample_meta["sample"], sample_meta["batch"]))
        adata.obs["batch"] = adata.obs["sample"].map(batch_map)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="linear_model_logcpm",
            design_covariates=["batch"],
        )

        df = run_pseudobulk_de(adata, config)

        assert not df.empty
        assert set(df["method"]) == {"linear_model_logcpm"}
        assert set(df["inference_level"]) == {"sample_level"}
        assert set(df["claim_level"]) == {"replicate_aware_sample_level_condition_inference"}
        assert df["valid_for_publication_inference"].all()
        assert df["design_covariates"].str.contains("batch").all()
        assert df["design_formula"].str.contains("C\\(__condition\\)").all()
        assert set(df["covariance_type"]) == {"HC3"}
        assert df["model_warning"].str.contains("HC3 covariance").all()

    def test_pseudobulk_de_linear_model_can_use_nonrobust_covariance(self):
        adata = self._make_pseudobulk_adata(n_reps=3)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="linear_model_logcpm",
            robust_cov_type="nonrobust",
        )

        df = run_pseudobulk_de(adata, config)

        assert not df.empty
        assert set(df["covariance_type"]) == {"nonrobust"}

    def test_pseudobulk_de_config_serializes_robust_covariance_type(self):
        cfg = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            contrasts=[("ctrl", "treat")],
            robust_cov_type="HC1",
        )

        assert cfg.model_dump()["robust_cov_type"] == "HC1"

    def test_pseudobulk_de_block_col_enters_linear_model(self):
        adata = self._make_pseudobulk_adata(n_reps=3)
        sample_meta = adata.obs[["sample", "condition"]].drop_duplicates().copy()
        sample_meta["patient"] = [f"p{idx % 3}" for idx in range(len(sample_meta))]
        patient_map = dict(zip(sample_meta["sample"], sample_meta["patient"]))
        adata.obs["patient"] = adata.obs["sample"].map(patient_map)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="linear_model_logcpm",
            block_col="patient",
        )

        df = run_pseudobulk_de(adata, config)

        assert not df.empty
        assert set(df["method"]) == {"linear_model_logcpm"}
        assert df["design_covariates"].str.contains("patient").all()
        assert set(df["block_col"]) == {"patient"}

    def test_pseudobulk_de_rejects_nonunique_sample_covariate(self):
        adata = self._make_pseudobulk_adata(n_reps=2)
        adata.obs["batch"] = "b1"
        first_sample = adata.obs["sample"].iloc[0]
        mask = adata.obs["sample"] == first_sample
        adata.obs.loc[adata.obs.index[mask.to_numpy().nonzero()[0][0]], "batch"] = "b2"
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="linear_model_logcpm",
            design_covariates=["batch"],
        )

        with pytest.raises(ValueError, match="multiple values in covariate"):
            run_pseudobulk_de(adata, config)


class TestStatsmodelsBackends:
    """Tests for Python-native statsmodels GLM/GEE pseudobulk backends."""

    def _make_pseudobulk_adata(self, n_reps=3):
        import numpy as np

        rows = []
        obs = []
        var_names = ["GeneA", "GeneB", "GeneC"]
        for cell_type in ["T", "B"]:
            for condition in ["ctrl", "treat"]:
                for rep in range(n_reps):
                    sample = f"{cell_type}_{condition}_{rep}"
                    for _ in range(4):
                        if cell_type == "T" and condition == "treat":
                            counts = [80, 10, 20]
                        elif cell_type == "T":
                            counts = [10, 80, 20]
                        else:
                            counts = [20, 20, 40]
                        rows.append(counts)
                        obs.append(
                            {
                                "sample": sample,
                                "condition": condition,
                                "cell_type": cell_type,
                            }
                        )
        adata = AnnData(X=np.asarray(rows, dtype=float), obs=pd.DataFrame(obs))
        adata.var_names = var_names
        return adata

    def test_statsmodels_glm_negative_binomial_runs(self):
        pytest.importorskip("statsmodels")
        adata = self._make_pseudobulk_adata(n_reps=3)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="statsmodels_glm",
            statsmodels_glm_family="NegativeBinomial",
        )
        df = run_pseudobulk_de(adata, config)
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        assert {"names", "logFC", "log2fc", "pvals", "padj", "group", "contrast"}.issubset(
            df.columns
        )
        t_gene_a = df[(df["group"] == "T") & (df["names"] == "GeneA")].iloc[0]
        assert t_gene_a["logFC"] > 0
        assert t_gene_a["method"] == "statsmodels_glm_negativebinomial"
        assert t_gene_a["inference_level"] == "sample_level"

    def test_statsmodels_glm_gaussian_runs(self):
        pytest.importorskip("statsmodels")
        adata = self._make_pseudobulk_adata(n_reps=3)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="statsmodels_glm",
            statsmodels_glm_family="Gaussian",
        )
        df = run_pseudobulk_de(adata, config)
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        t_gene_a = df[(df["group"] == "T") & (df["names"] == "GeneA")].iloc[0]
        assert t_gene_a["logFC"] > 0
        assert t_gene_a["method"] == "statsmodels_glm_gaussian"

    def test_statsmodels_glm_gamma_runs(self):
        pytest.importorskip("statsmodels")
        adata = self._make_pseudobulk_adata(n_reps=3)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="statsmodels_glm",
            statsmodels_glm_family="Gamma",
        )
        df = run_pseudobulk_de(adata, config)
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        assert set(df["method"]) == {"statsmodels_glm_gamma"}

    def test_statsmodels_gee_runs(self):
        pytest.importorskip("statsmodels")
        adata = self._make_pseudobulk_adata(n_reps=3)
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            groupby="cell_type",
            group_names=["T"],
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="statsmodels_gee",
            statsmodels_glm_family="NegativeBinomial",
        )
        df = run_pseudobulk_de(adata, config)
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        assert {"names", "logFC", "pvals", "padj", "group", "contrast"}.issubset(df.columns)
        assert df["method"].iloc[0].startswith("statsmodels_gee_")


class TestHierarchicalFDR:
    """Tests for hierarchical multiple-testing correction."""

    def test_hierarchical_fdr_produces_expected_columns(self):
        df = pd.DataFrame(
            {
                "names": ["g1", "g2", "g3", "g4"],
                "cell_type": ["A", "A", "B", "B"],
                "contrast": ["c1", "c1", "c1", "c1"],
                "pvals": [0.01, 0.04, 0.001, 0.10],
            }
        )
        corrected = _hierarchical_fdr_correction(df, group_cols=["cell_type", "contrast"])
        assert "padj_within_group" in corrected.columns
        assert "padj_global_by" in corrected.columns
        assert corrected["padj_within_group"].notna().all()
        assert corrected["padj_global_by"].notna().all()

    def test_hierarchical_fdr_through_pseudobulk_de(self):
        pytest.importorskip("statsmodels")
        rows = []
        obs = []
        var_names = ["GeneA", "GeneB", "GeneC"]
        for condition in ["ctrl", "treat"]:
            for rep in range(3):
                sample = f"{condition}_{rep}"
                for _ in range(4):
                    counts = [80, 10, 20] if condition == "treat" else [10, 80, 20]
                    rows.append(counts)
                    obs.append({"sample": sample, "condition": condition})
        adata = AnnData(X=np.asarray(rows, dtype=float), obs=pd.DataFrame(obs))
        adata.var_names = var_names
        config = PseudobulkDEConfig(
            sample_col="sample",
            condition_key="condition",
            contrasts=[("ctrl", "treat")],
            min_cells_per_sample=1,
            min_counts=0,
            method="statsmodels_glm",
            statsmodels_glm_family="Gaussian",
            hierarchical_correction=True,
        )
        df = run_pseudobulk_de(adata, config)
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        assert "padj_within_group" in df.columns
        assert "padj_global_by" in df.columns


class TestMixedLM:
    """Tests for cell-level sample-aware MixedLM DE."""

    def _make_mixedlm_adata(self):
        import numpy as np

        rng = np.random.default_rng(42)
        n_cells = 60
        n_genes = 10
        samples = [f"sample_{i % 6}" for i in range(n_cells)]
        conditions = ["ctrl" if i % 6 < 3 else "treat" for i in range(n_cells)]
        X = rng.poisson(2, size=(n_cells, n_genes)).astype(float)
        # Upregulate first gene in treat
        for i, cond in enumerate(conditions):
            if cond == "treat":
                X[i, 0] += 5
        obs = pd.DataFrame({"sample": samples, "condition": conditions})
        adata = AnnData(X=X, obs=obs)
        adata.var_names = [f"gene_{i}" for i in range(n_genes)]
        return adata

    def test_run_mixedlm_de_returns_sample_aware_tags(self):
        pytest.importorskip("statsmodels")
        adata = self._make_mixedlm_adata()
        df = run_mixedlm_de(
            adata,
            sample_col="sample",
            condition_col="condition",
            condition1="ctrl",
            condition2="treat",
            max_genes=10,
        )
        assert isinstance(df, pd.DataFrame)
        assert not df.empty
        assert {"names", "logFC", "pvals", "padj"}.issubset(df.columns)
        assert set(df["inference_level"]) == {"cell_level_sample_aware"}
        assert set(df["model"]) == {"statsmodels_mixedlm"}
        assert not df["valid_for_publication_inference"].any()
        assert "mixedlm_de" in adata.uns["sclucid"]["analysis"]["de"]

    def test_run_mixedlm_de_warns_and_truncates_over_max_genes(self, caplog):
        pytest.importorskip("statsmodels")
        import logging

        adata = self._make_mixedlm_adata()
        # Add a bunch of fake genes to exceed max_genes
        extra = np.zeros((adata.n_obs, 100))
        full_X = np.concatenate([adata.X, extra], axis=1)
        adata = AnnData(X=full_X, obs=adata.obs)
        adata.var_names = [f"gene_{i}" for i in range(full_X.shape[1])]
        with caplog.at_level(logging.WARNING):
            df = run_mixedlm_de(
                adata,
                sample_col="sample",
                condition_col="condition",
                condition1="ctrl",
                condition2="treat",
                max_genes=50,
            )
        assert len(df) <= 50
        assert "exceeding max_genes" in caplog.text
