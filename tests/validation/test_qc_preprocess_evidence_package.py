"""Tests for combined QC+preprocess evidence extraction."""

from pathlib import Path

from validation.qc_preprocess.build_review_summary_evidence_tables import _decision_rows


def test_combined_evidence_rows_include_qc_and_preprocess_handoffs():
    sclucid = {
        "qc": {
            "review_summary": {
                "data": {
                    "filtering_summary": {"removed_fraction": 0.1},
                    "qc_handoff_readiness": {
                        "status": "review_required",
                        "recommended_preprocess_counts_layer": "counts",
                        "handoff_cell_counts": {"review_required": 3, "sensitivity_only": 1},
                        "handoff_cell_fractions": {"review_required": 0.2},
                        "warnings": ["carry review cells into sensitivity analysis"],
                    },
                    "doublet_evidence_summary": {"review_required": False},
                }
            }
        },
        "preprocess": {
            "review_summary": {
                "data": {
                    "analysis_handoff_readiness": {
                        "status": "ready",
                        "graph_representation": "X_pca",
                        "expression_for_markers": "adata.raw.X",
                        "expression_for_de": "adata.layers['normalized']",
                    }
                }
            }
        },
    }

    rows = _decision_rows("demo", Path("demo.h5ad"), sclucid)
    decisions = {row["decision"]: row for row in rows}

    assert decisions["qc_to_preprocess_handoff"]["value"] == "counts"
    assert decisions["review_sensitivity_cell_tracking"]["review_required"] == "yes"
    assert decisions["preprocess_to_analysis_handoff"]["value"] == "X_pca"
    assert decisions["marker_de_expression_source"]["value"] == "adata.layers['normalized']"
