"""Tests that legacy fixtures degrade gracefully when cell_type is missing."""

import numpy as np
import pandas as pd
from anndata import AnnData

from validation.qc.run_threshold_benchmark import _retention_rows
from validation.qc.run_tumor_biological_fidelity_benchmark import _retention_bias_rows
from validation.preprocess.run_graph_stability_benchmark import _rare_population_rows


def test_threshold_retention_records_missing_cell_type():
    rng = np.random.default_rng(0)
    X = rng.poisson(2, (50, 20)).astype(float)
    adata = AnnData(X)
    adata.obs["sample"] = ["A"] * 25 + ["B"] * 25
    keep = pd.Series([True] * 50, index=adata.obs_names)
    rows = _retention_rows(adata, "demo", "scanpy_fixed_threshold", keep)
    statuses = {r["group_type"]: r["annotation_status"] for r in rows}
    assert statuses["cell_type"] == "missing_cell_type"
    assert statuses["sample"] == "ok"


def test_tumor_fidelity_records_missing_cell_type():
    rng = np.random.default_rng(0)
    X = rng.poisson(2, (50, 20)).astype(float)
    adata = AnnData(X)
    adata.obs["sample"] = ["A"] * 25 + ["B"] * 25
    keep = pd.Series([True] * 50, index=adata.obs_names)
    rows = _retention_bias_rows(adata, "demo", "sclucid_tumor_aware", keep)
    cell_type_rows = [r for r in rows if r["group_type"] == "cell_type"]
    assert len(cell_type_rows) == 1
    assert cell_type_rows[0]["review_reason"] == "missing_cell_type"


def test_graph_stability_records_missing_labels():
    labels = pd.Series(name="cell_type")
    clusters = np.array([0] * 30)
    rows = _rare_population_rows("demo", labels, clusters, n_pcs=20)
    assert len(rows) == 1
    assert rows[0]["rare_population"] == "__missing__"
    assert rows[0]["preservation_proxy"] == "missing_annotation"
