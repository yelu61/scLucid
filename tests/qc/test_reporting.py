"""Smoke tests for QC reporting module."""

import os

from scLucid.qc.reporting import EnhancedQCReport, generate_qc_html_report


class TestEnhancedQCReport:
    def test_init(self, qc_test_adata):
        reporter = EnhancedQCReport(qc_test_adata)
        assert reporter is not None
        assert reporter.adata is qc_test_adata

    def test_init_with_before(self, qc_test_adata):
        reporter = EnhancedQCReport(qc_test_adata, adata_before=qc_test_adata)
        assert reporter.adata_before is qc_test_adata

    def test_generate_html_report_smoke(self, qc_test_adata, temp_output_dir):
        output = os.path.join(temp_output_dir, "qc_report.html")
        reporter = EnhancedQCReport(qc_test_adata)
        reporter.generate_html_report(
            output_path=output,
            include_plots=False,
            include_recommendations=True,
        )
        assert os.path.exists(output)

class TestConvenienceFunctions:
    def test_generate_qc_html_report(self, qc_test_adata, temp_output_dir):
        output = os.path.join(temp_output_dir, "qc_conv.html")
        result = generate_qc_html_report(
            qc_test_adata,
            output_path=output,
            title="Test QC Report",
        )
        assert result == output
        assert os.path.exists(output)
