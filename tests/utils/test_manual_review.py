"""Tests for manual review finalization helpers."""

import numpy as np
from anndata import AnnData

from scLucid.utils import finalize_manual_review_summary


def test_finalize_manual_review_summary_writes_contract(tmp_path):
    adata = AnnData(np.ones((3, 2)))
    summary = finalize_manual_review_summary(
        adata,
        module="qc",
        workflow_name="manual_test",
        steps=["qc_metrics", "filtering"],
        config={"sample_key": "sampleID"},
        summary={"custom_note": "reviewed"},
        save_dir=tmp_path,
        warnings=["manual warning"],
    )

    qc_ns = adata.uns["sclucid"]["qc"]
    assert qc_ns["workflow_config"] == {"sample_key": "sampleID"}
    assert qc_ns["steps_executed"] == ["qc_metrics", "filtering"]
    assert qc_ns["review_summary"] is summary
    assert summary["module"] == "qc"
    assert summary["workflow_name"] == "manual_test"
    assert summary["steps_executed"] == ["qc_metrics", "filtering"]
    assert summary["data_shape"] == {"n_cells": 3, "n_genes": 2}
    assert summary["warnings"] == ["manual warning"]
    assert summary["data"]["custom_note"] == "reviewed"
    assert (tmp_path / "qc_review_summary.json").exists()
    assert (tmp_path / "qc_review_summary.md").exists()
