"""Smoke tests for tumor microenvironment analysis."""

import pytest


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
