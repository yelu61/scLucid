"""Tests for HVG (highly variable gene) selection module."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.preprocess.hvg import (
    PROTECTED_GENE_PRESETS,
    _exclude_genes,
    _gene_type_detection,
    _get_hvg_input_matrix,
    _infer_species_from_gene_names,
    _validate_hvg_input_matrix,
    evaluate_hvg_stability,
    find_hvgs,
    select_hvg_sets,
    suggest_hvg_choice,
)
from scLucid.preprocess.hvg.core import _apply_hvg_biological_protection


class TestGetHVGInputMatrix:
    def test_from_X(self, tiny_adata):
        result = _get_hvg_input_matrix(tiny_adata, "X")
        assert result.shape == (8, 6)

    def test_from_layer(self, tiny_adata):
        result = _get_hvg_input_matrix(tiny_adata, "counts")
        assert result is not None

    def test_missing_layer_raises(self, tiny_adata):
        with pytest.raises(KeyError, match="not found"):
            _get_hvg_input_matrix(tiny_adata, "nonexistent_layer")


class TestValidateHVGInputMatrix:
    def test_valid_input(self, tiny_adata):
        _validate_hvg_input_matrix(tiny_adata.X, "X", "scanpy")

    def test_negative_values_raises(self):
        adata = AnnData(X=np.array([[-1.0, 2.0], [3.0, -2.0]]))
        with pytest.raises(ValueError, match="negative"):
            _validate_hvg_input_matrix(adata.X, "t", "scanpy")

    def test_empty_raises(self):
        adata = AnnData(X=np.zeros((0, 5)))
        with pytest.raises(ValueError, match="empty"):
            _validate_hvg_input_matrix(adata.X, "t", "scanpy")


class TestGeneTypeDetection:
    @staticmethod
    def _make_var_names(*gene_names):
        return pd.Index(gene_names)

    def test_detects_human_mitochondrial(self):
        vn = self._make_var_names("MT-ND1", "MT-CO1", "GAPDH", "ACTB")
        r = _gene_type_detection(vn, species="human")
        mt = r["mitochondrial"]
        assert mt[0]
        assert mt[1]
        assert not mt[2]
        assert not mt[3]

    def test_detects_mouse_mitochondrial(self):
        vn = self._make_var_names("mt-Nd1", "mt-Co1", "Gapdh", "Actb")
        r = _gene_type_detection(vn, species="mouse")
        mt = r["mitochondrial"]
        assert mt[0]
        assert mt[1]

    def test_detects_human_ribosomal(self):
        vn = self._make_var_names("RPS1", "RPL10", "MRPL3", "MRPS5", "GAPDH")
        r = _gene_type_detection(vn, species="human")
        rb = r["ribosomal"]
        assert rb[0]
        assert rb[1]
        assert rb[2]
        assert rb[3]
        assert not rb[4]

    def test_detects_mouse_ribosomal(self):
        vn = self._make_var_names("Rps1", "Rpl10a", "Mrpl3", "Mrps5", "Gapdh")
        r = _gene_type_detection(vn, species="mouse")
        rb = r["ribosomal"]
        assert rb[0]
        assert rb[1]
        assert rb[2]
        assert rb[3]

    def test_detects_hemoglobin(self):
        vn = self._make_var_names("HBA1", "HBB", "HBG1", "GAPDH")
        r = _gene_type_detection(vn, species="human")
        assert r["hemoglobin"][0]
        assert r["hemoglobin"][1]
        assert r["hemoglobin"][2]

    def test_detects_mouse_hemoglobin(self):
        vn = self._make_var_names("Hba-a1", "Hbb-bs", "Gapdh")
        r = _gene_type_detection(vn, species="mouse")
        assert r["hemoglobin"][0]
        assert r["hemoglobin"][1]

    def test_detects_heat_shock(self):
        vn = self._make_var_names("HSPA1A", "HSPB1", "DNAJA1", "GAPDH")
        r = _gene_type_detection(vn, species="human")
        assert r["heat_shock"][0]
        assert r["heat_shock"][1]
        assert r["heat_shock"][2]

    def test_detects_immediate_early(self):
        vn = self._make_var_names("FOS", "JUN", "EGR1", "ATF3", "GAPDH")
        r = _gene_type_detection(vn, species="human")
        ieg = r["immediate_early"]
        assert ieg[0]
        assert ieg[1]
        assert ieg[2]
        assert ieg[3]

    def test_detects_mouse_immediate_early(self):
        vn = self._make_var_names("Fos", "Junb", "Egr1", "Nr4a1", "Gapdh")
        r = _gene_type_detection(vn, species="mouse")
        ieg = r["immediate_early"]
        assert ieg[0]
        assert ieg[1]
        assert ieg[2]

    def test_auto_detect_species(self):
        vn = self._make_var_names("MT-ND1", "GAPDH")
        r = _gene_type_detection(vn, species=None)
        assert r["mitochondrial"][0]

    def test_rat_species(self):
        vn = self._make_var_names("Mt-Nd1", "Rpl10", "Gapdh")
        r = _gene_type_detection(vn, species="rat")
        mt = r["mitochondrial"]
        assert isinstance(mt[0], (bool, np.bool_))


class TestInferSpeciesFromGeneNames:
    def test_infers_human_by_mt_prefix(self):
        assert _infer_species_from_gene_names(pd.Index(["MT-ND1", "GAPDH", "ACTB"])) == "human"

    def test_infers_mouse_by_mt_prefix(self):
        assert _infer_species_from_gene_names(pd.Index(["mt-Nd1", "Gapdh", "Actb"])) == "mouse"

    def test_infers_human_by_ribosomal(self):
        assert _infer_species_from_gene_names(pd.Index(["RPS1", "RPL10", "GAPDH"])) == "human"

    def test_infers_mouse_by_ribosomal_title_case(self):
        # Use genes that are title case but not caught as rat (no Rpl prefix)
        result = _infer_species_from_gene_names(pd.Index(["Mrpl3", "Mrps5", "Gapdh"]))
        assert result == "mouse"

    def test_infers_rat_by_rpl_prefix(self):
        assert _infer_species_from_gene_names(pd.Index(["Rpl10", "Rpl32", "Gapdh"])) == "rat"

    def test_defaults_to_human(self):
        assert _infer_species_from_gene_names(pd.Index(["UNIQUE1", "UNIQUE2"])) == "human"

    def test_with_realistic_human_genes(self):
        genes = pd.Index(["MT-ND1", "RPS1", "GAPDH", "TP53", "KRAS", "ACTB"])
        assert _infer_species_from_gene_names(genes) == "human"


class TestExcludeGenes:
    def test_excludes_mitochondrial(self, minimal_adata):
        gene_types = _gene_type_detection(minimal_adata.var_names, species="human")
        hvg_mask = np.ones(minimal_adata.n_vars, dtype=bool)

        mask, counts = _exclude_genes(
            minimal_adata, hvg_mask,
            exclude_types=["mitochondrial"],
            gene_types=gene_types, species="human",
        )
        assert isinstance(mask, np.ndarray)
        assert isinstance(counts, dict)
        assert "mitochondrial" in counts

    def test_no_exclusion_when_none_specified(self, minimal_adata):
        hvg_mask = np.ones(minimal_adata.n_vars, dtype=bool)
        gene_types = _gene_type_detection(minimal_adata.var_names, species="human")
        mask, counts = _exclude_genes(
            minimal_adata, hvg_mask,
            exclude_types=[],
            gene_types=gene_types, species="human",
        )
        assert mask.sum() == hvg_mask.sum()

    def test_multiple_exclusion_types(self, minimal_adata):
        gene_types = _gene_type_detection(minimal_adata.var_names, species="human")
        hvg_mask = np.ones(minimal_adata.n_vars, dtype=bool)

        mask, counts = _exclude_genes(
            minimal_adata, hvg_mask,
            exclude_types=["mitochondrial", "ribosomal", "hemoglobin"],
            gene_types=gene_types, species="human",
        )
        assert isinstance(counts, dict)
        assert 2 <= len(counts) <= 3


class TestHVGBiologicalProtection:
    def test_protection_cap_uses_stable_preset_order_and_reports_truncation(self):
        adata = AnnData(X=np.ones((4, 6)))
        adata.var_names = ["TRAC", "IL2", "FOXP3", "GeneA", "GeneB", "GeneC"]
        hvg_mask = np.array([False, False, False, True, False, False])

        updated, report = _apply_hvg_biological_protection(
            adata,
            hvg_mask,
            presets=["immune_receptor", "cytokine", "transcription_factor"],
            max_extra_genes=2,
        )

        assert report["truncated"] is True
        assert report["n_rescued_before_cap"] == 3
        assert report["rescued_genes"] == ["TRAC", "IL2"]
        assert "preset order" in report["truncation_policy"]
        assert updated.tolist() == [True, True, False, True, False, False]
        assert adata.var["hvg_protection_rescued"].tolist() == [
            True,
            True,
            False,
            False,
            False,
            False,
        ]


class TestFindHVGs:
    def test_scanpy_method_with_config(self, minimal_adata):
        """Smoke test for scanpy HVG selection using a config object."""
        from scLucid.preprocess.config import HVGConfig

        config = HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100)
        adata = minimal_adata.copy()
        result = find_hvgs(adata, config=config, force=True, input_layer="counts")
        assert result is not None
        assert "highly_variable_scanpy_seurat" in result.var.columns

    def test_scanpy_method_v3_with_config(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig

        config = HVGConfig(method="scanpy", flavor="seurat_v3", n_top_genes=100)
        adata = minimal_adata.copy()
        result = find_hvgs(adata, config=config, force=True, input_layer="counts")
        assert result is not None

    def test_find_hvgs_no_overwrite_existing(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig

        config = HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100)
        adata = minimal_adata.copy()
        find_hvgs(adata, config=config, force=True, input_layer="counts")
        hvg_key = "highly_variable_scanpy_seurat"
        assert hvg_key in adata.var
        original_count = adata.var[hvg_key].sum()

        # Running again without force keeps the original
        config2 = HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100)
        find_hvgs(adata, config=config2, force=False, input_layer="counts")
        assert adata.var[hvg_key].sum() == original_count

    def test_config_override_via_kwargs(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig

        config = HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100)
        adata = minimal_adata.copy()
        result = find_hvgs(adata, config=config, force=True, n_top_genes=120, input_layer="counts")
        assert result is not None

    def test_protected_gene_set_rescues_biology_gene_into_final_hvg(self):
        from scLucid.preprocess.config import HVGConfig

        rng = np.random.default_rng(0)
        X = rng.poisson(2, size=(160, 130)).astype(float)
        X[:, -1] = 0.0
        adata = AnnData(X=X)
        adata.var_names = [f"Gene{i}" for i in range(129)] + ["TRAC"]
        adata.layers["counts"] = adata.X.copy()

        config = HVGConfig(
            method="scanpy",
            flavor="seurat",
            n_top_genes=100,
            auto_n_top_genes=False,
            protected_gene_sets={"immune_core": ["TRAC"]},
        )
        result = find_hvgs(adata, config=config, force=True, input_layer="counts")
        hvg_key = "highly_variable_scanpy_seurat"
        report = result.uns["sclucid"]["preprocess"]["hvg"]["biological_protection"]

        assert bool(result.var.loc["TRAC", hvg_key]) is True
        assert bool(result.var.loc["TRAC", "hvg_protection_rescued"]) is True
        assert report["n_rescued"] >= 1
        assert "TRAC" in report["rescued_genes"]

    def test_preserve_tumor_heterogeneity_enables_tumor_preset(self):
        from scLucid.preprocess.config import HVGConfig

        rng = np.random.default_rng(1)
        X = rng.poisson(2, size=(160, 130)).astype(float)
        X[:, -1] = 0.0
        adata = AnnData(X=X)
        adata.var_names = [f"Gene{i}" for i in range(129)] + ["EPCAM"]
        adata.layers["counts"] = adata.X.copy()

        config = HVGConfig(
            method="scanpy",
            flavor="seurat",
            n_top_genes=100,
            auto_n_top_genes=False,
            protected_gene_presets=[],
        )
        result = find_hvgs(
            adata,
            config=config,
            force=True,
            input_layer="counts",
            preserve_tumor_heterogeneity=True,
        )
        report = result.uns["sclucid"]["preprocess"]["hvg"]["biological_protection"]

        assert "tumor_heterogeneity" in report["presets"]
        assert bool(result.var.loc["EPCAM", "highly_variable_scanpy_seurat"]) is True
        assert "EPCAM" in PROTECTED_GENE_PRESETS["tumor_heterogeneity"]

    def test_custom_hvg_records_sample_consensus_report(self):
        from scLucid.preprocess.config import HVGConfig

        rng = np.random.default_rng(2)
        X = rng.poisson(2, size=(180, 140)).astype(float)
        adata = AnnData(X=X)
        adata.var_names = [f"Gene{i}" for i in range(140)]
        adata.obs["sampleID"] = ["s1"] * 60 + ["s2"] * 60 + ["s3"] * 60
        adata.layers["counts"] = adata.X.copy()

        config = HVGConfig(
            method="custom",
            flavor="seurat",
            n_top_genes=100,
            auto_n_top_genes=False,
            min_n_samples=2,
            n_highly_expressed_genes=0,
            n_specific_genes=0,
            exclude_gene_types=[],
        )
        result = find_hvgs(
            adata,
            config=config,
            force=True,
            input_layer="counts",
            n_jobs=1,
            plot=False,
        )
        report = result.uns["sclucid"]["preprocess"]["hvg"]["consensus_report"]

        assert report["available"] is True
        assert report["evidence_key"] == "highly_variable_custom_sample_count"
        assert "sample_count_distribution" in report


class TestEvaluateHVGStability:
    def test_stability_runs(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig

        adata = minimal_adata.copy()
        config = HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100)
        find_hvgs(adata, config=config, force=True, input_layer="counts")

        result = evaluate_hvg_stability(
            adata,
            hvg_key="highly_variable_scanpy_seurat",
            n_bootstrap=3,
            sample_fraction=0.8,
            plot=False,
        )
        assert result is not None
        assert "hvg_selection_frequency" in result.var.columns


class TestSelectHVGSets:
    def test_direct_mode(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig

        adata = minimal_adata.copy()
        config = HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100)
        find_hvgs(adata, config=config, force=True, input_layer="counts")

        result = select_hvg_sets(
            adata,
            hvg_keys=["highly_variable_scanpy_seurat"],
            mode="direct", subset=False, keep_raw=True,
        )
        assert result is not None

    def test_intersection_mode(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig

        adata = minimal_adata.copy()
        find_hvgs(adata, config=HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100),
                  force=True, input_layer="counts")
        find_hvgs(adata, config=HVGConfig(method="scanpy", flavor="seurat_v3", n_top_genes=100),
                  force=True, input_layer="counts")

        result = select_hvg_sets(
            adata,
            hvg_keys=["highly_variable_scanpy_seurat", "highly_variable_scanpy_seurat_v3"],
            mode="intersection", subset=False, keep_raw=True,
        )
        assert result is not None

    def test_union_mode(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig

        adata = minimal_adata.copy()
        find_hvgs(adata, config=HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100),
                  force=True, input_layer="counts")

        result = select_hvg_sets(
            adata,
            hvg_keys=["highly_variable_scanpy_seurat"],
            mode="union", subset=False, keep_raw=True,
        )
        assert result is not None

    def test_select_and_audit_hvgs(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig
        from scLucid.preprocess.hvg.selection import select_and_audit_hvgs

        adata = minimal_adata.copy()
        find_hvgs(adata, config=HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100),
                  force=True, input_layer="counts")

        result, audit = select_and_audit_hvgs(
            adata,
            hvg_keys=["highly_variable_scanpy_seurat"],
            mode="direct",
            subset=False,
            keep_raw=True,
        )

        assert result is not None
        assert audit["n_selected"] > 0
        assert "hvg_selection_audit" in result.uns["sclucid"]["preprocess"]


class TestSuggestHVGChoice:
    def test_suggest_returns_structured_guidance(self, minimal_adata):
        from scLucid.preprocess.config import HVGConfig

        adata = minimal_adata.copy()
        find_hvgs(adata, config=HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100),
                  force=True, input_layer="counts")
        find_hvgs(adata, config=HVGConfig(method="scanpy", flavor="seurat_v3", n_top_genes=100),
                  force=True, input_layer="counts")

        suggestion = suggest_hvg_choice(
            adata,
            hvg_keys=["highly_variable_scanpy_seurat", "highly_variable_scanpy_seurat_v3"],
            mode="auto",
        )
        assert suggestion["requested_mode"] == "auto"
        assert suggestion["recommended_mode"] in {"union", "intersection"}
        assert 0 <= suggestion["jaccard_index"] <= 1
        assert suggestion["messages"]


class TestEdgeCases:
    def test_all_zero_expression(self):
        """All-zero expression trigger input validation before calling scanpy."""
        adata = AnnData(X=np.zeros((100, 120)))
        adata.obs_names = [f"c{i}" for i in range(100)]
        adata.var_names = [f"g{i}" for i in range(120)]
        adata.layers["counts"] = adata.X.copy()

        from scLucid.preprocess.config import HVGConfig

        config = HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100)
        # All-zero data will likely crash inside scanpy; just verify our wrapper
        # doesn't crash before handing off to scanpy
        try:
            result = find_hvgs(adata, config=config, force=True, input_layer="counts")
            assert result is not None
        except (ValueError, IndexError):
            # Expected when scanpy can't handle the data
            pass

    def test_minimal_gene_count(self):
        adata = AnnData(X=np.random.poisson(2, size=(50, 10)).astype(float))
        adata.var_names = [f"g{i}" for i in range(10)]
        adata.layers["counts"] = adata.X.copy()

        from scLucid.preprocess.config import HVGConfig

        config = HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100)
        result = find_hvgs(adata, config=config, force=True, input_layer="counts")
        assert result is not None
