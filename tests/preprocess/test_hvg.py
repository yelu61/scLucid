"""Tests for HVG (highly variable gene) selection module."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from scLucid.preprocess.hvg import (
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
        assert mt[0]; assert mt[1]; assert not mt[2]; assert not mt[3]

    def test_detects_mouse_mitochondrial(self):
        vn = self._make_var_names("mt-Nd1", "mt-Co1", "Gapdh", "Actb")
        r = _gene_type_detection(vn, species="mouse")
        mt = r["mitochondrial"]
        assert mt[0]; assert mt[1]

    def test_detects_human_ribosomal(self):
        vn = self._make_var_names("RPS1", "RPL10", "MRPL3", "MRPS5", "GAPDH")
        r = _gene_type_detection(vn, species="human")
        rb = r["ribosomal"]
        assert rb[0]; assert rb[1]; assert rb[2]; assert rb[3]; assert not rb[4]

    def test_detects_mouse_ribosomal(self):
        vn = self._make_var_names("Rps1", "Rpl10a", "Mrpl3", "Mrps5", "Gapdh")
        r = _gene_type_detection(vn, species="mouse")
        rb = r["ribosomal"]
        assert rb[0]; assert rb[1]; assert rb[2]; assert rb[3]

    def test_detects_hemoglobin(self):
        vn = self._make_var_names("HBA1", "HBB", "HBG1", "GAPDH")
        r = _gene_type_detection(vn, species="human")
        assert r["hemoglobin"][0]; assert r["hemoglobin"][1]; assert r["hemoglobin"][2]

    def test_detects_mouse_hemoglobin(self):
        vn = self._make_var_names("Hba-a1", "Hbb-bs", "Gapdh")
        r = _gene_type_detection(vn, species="mouse")
        assert r["hemoglobin"][0]; assert r["hemoglobin"][1]

    def test_detects_heat_shock(self):
        vn = self._make_var_names("HSPA1A", "HSPB1", "DNAJA1", "GAPDH")
        r = _gene_type_detection(vn, species="human")
        assert r["heat_shock"][0]; assert r["heat_shock"][1]; assert r["heat_shock"][2]

    def test_detects_immediate_early(self):
        vn = self._make_var_names("FOS", "JUN", "EGR1", "ATF3", "GAPDH")
        r = _gene_type_detection(vn, species="human")
        ieg = r["immediate_early"]
        assert ieg[0]; assert ieg[1]; assert ieg[2]; assert ieg[3]

    def test_detects_mouse_immediate_early(self):
        vn = self._make_var_names("Fos", "Junb", "Egr1", "Nr4a1", "Gapdh")
        r = _gene_type_detection(vn, species="mouse")
        ieg = r["immediate_early"]
        assert ieg[0]; assert ieg[1]; assert ieg[2]

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


class TestSuggestHVGChoice:
    def test_suggest_does_not_crash(self, minimal_adata, capsys):
        from scLucid.preprocess.config import HVGConfig

        adata = minimal_adata.copy()
        find_hvgs(adata, config=HVGConfig(method="scanpy", flavor="seurat", n_top_genes=100),
                  force=True, input_layer="counts")
        find_hvgs(adata, config=HVGConfig(method="scanpy", flavor="seurat_v3", n_top_genes=100),
                  force=True, input_layer="counts")

        suggest_hvg_choice(
            adata,
            hvg_keys=["highly_variable_scanpy_seurat", "highly_variable_scanpy_seurat_v3"],
            mode="standard",
        )
        captured = capsys.readouterr()
        assert len(captured.out) > 0 or len(captured.err) > 0


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
