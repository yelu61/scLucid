"""Smoke tests for tumor microenvironment analysis."""

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


class TestMicroenvironment:
    def test_import(self):
        from scLucid.tumor.microenvironment.deconvolution import (
            TMEProfiler,
            analyze_immune_infiltration,
            deconvolve_tme,
            estimate_stromal_content,
        )
        assert callable(deconvolve_tme)
        assert callable(estimate_stromal_content)
        assert callable(analyze_immune_infiltration)

    def test_tme_profiler_init(self):
        from scLucid.tumor.microenvironment.deconvolution import TMEProfiler

        profiler = TMEProfiler()
        assert profiler is not None

    def test_tme_profiler_normalizes_label_variants(self):
        from scLucid.tumor.microenvironment.deconvolution import TMEProfiler

        adata = AnnData(np.zeros((8, 2)))
        adata.obs["cell_type"] = pd.Categorical(
            ["T_cell", "T cell", "CD8_T", "B_cell", "B cell", "Fibroblast", "CAF", "Tumor"]
        )
        profiler = TMEProfiler()
        profiler.fit(adata)

        assert profiler.immune_score_ == pytest.approx(5 / 8, abs=1e-6)
        assert profiler.stromal_score_ == pytest.approx(2 / 8, abs=1e-6)
        assert profiler.malignant_score_ == pytest.approx(1 / 8, abs=1e-6)
        assert set(profiler.proportions_.index.tolist()) == {"immune", "stromal", "malignant"}

    def test_tme_profiler_tracks_unmapped_labels(self):
        from scLucid.tumor.microenvironment.deconvolution import TMEProfiler

        adata = AnnData(np.zeros((3, 2)))
        adata.obs["cell_type"] = pd.Categorical([" WeirdType", "weird_type", "Tumor"])
        profiler = TMEProfiler()
        profiler.fit(adata)

        assert "weird type" in profiler.unmapped_labels_
        assert profiler.malignant_score_ == pytest.approx(1 / 3, abs=1e-6)
        assert profiler.proportions_["other"] == pytest.approx(2 / 3, abs=1e-6)

    def test_tme_profiler_does_not_call_normal_epithelial_malignant(self):
        from scLucid.tumor.microenvironment.deconvolution import TMEProfiler

        adata = AnnData(np.zeros((3, 2)))
        adata.obs["cell_type"] = pd.Categorical(["normal_epithelial", "Tumor", "Cancer"])
        profiler = TMEProfiler()
        profiler.fit(adata)

        assert profiler.malignant_score_ == pytest.approx(2 / 3, abs=1e-6)
        assert profiler.proportions_["other"] == pytest.approx(1 / 3, abs=1e-6)
        assert "normal epithelial" not in profiler.unmapped_labels_

    def test_deconvolve_tme_adds_compartment_claim(self):
        from scLucid.tumor.microenvironment.deconvolution import deconvolve_tme

        adata = AnnData(np.zeros((4, 2)))
        adata.obs["cell_type"] = pd.Categorical(["T_cell", "Fibroblast", "Tumor", "Tumor"])
        result = deconvolve_tme(adata, copy=True)

        assert "tme_compartment_claim" in result.uns
        assert "annotation-derived" in result.uns["tme_compartment_claim"]
        assert "tme_compartment" in result.obs.columns
        assert set(result.obs["tme_compartment"].unique().tolist()) == {"immune", "stromal", "malignant"}

    def test_deconvolve_tme_smoke(self, qc_test_adata):
        from scLucid.tumor.microenvironment.deconvolution import deconvolve_tme

        result = deconvolve_tme(qc_test_adata)
        assert result is not None

    def test_immune_infiltration(self, qc_test_adata):
        from scLucid.tumor.microenvironment.deconvolution import analyze_immune_infiltration

        result = analyze_immune_infiltration(qc_test_adata)
        assert result is not None

    def test_stromal_content(self, qc_test_adata):
        from scLucid.tumor.microenvironment.deconvolution import estimate_stromal_content

        result = estimate_stromal_content(qc_test_adata)
        assert result is not None

    def test_ecosystem_import(self):
        from scLucid.tumor.microenvironment.ecosystem import (
            EcosystemAnalyzer,
        )
        assert EcosystemAnalyzer

    def test_interaction_import(self):
        from scLucid.tumor.microenvironment.interaction import (
            InteractionAnalyzer,
        )
        assert InteractionAnalyzer
