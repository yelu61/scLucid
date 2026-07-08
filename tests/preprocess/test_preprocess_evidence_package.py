"""Tests for preprocess validation evidence package helpers."""

import pandas as pd

from validation.preprocess.build_preprocess_evidence_package import (
    _source_rows_from_module_table,
)


def test_preprocess_source_rows_accept_legacy_figure_panel_column():
    table = pd.DataFrame(
        [
            {
                "figure_panel": "3B",
                "dataset": "lin2020.pdac",
                "strategy": "semantic_auto_protected_union",
                "inclusion_rate": 0.9,
                "genes_present": 20,
                "context": '{"auto_risk": "low"}',
            }
        ]
    )

    rows = _source_rows_from_module_table(
        table,
        evidence_domain="hvg_marker_program_preservation",
        source_file="figure3_hvg_data.tsv",
        default_panel="preprocess_hvg",
        metric_columns=("inclusion_rate", "genes_present"),
    )

    assert len(rows) == 2
    assert rows[0]["module_panel"] == "3B"
    assert rows[0]["strategy_or_method"] == "semantic_auto_protected_union"
    assert "auto_risk" in rows[0]["context"]
