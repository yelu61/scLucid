from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from validation.preprocess.run_mixology_preprocess_benchmark import run


def _synthetic_mixology(path) -> None:
    rng = np.random.default_rng(41)
    protocols = np.repeat(["P1", "P2", "P3"], 45)
    identities = np.tile(np.repeat(["A", "B", "C"], 15), 3)
    matrix = rng.poisson(0.4, size=(len(protocols), 90)).astype(np.int32)
    for index, identity in enumerate(("A", "B", "C")):
        matrix[identities == identity, index * 20 : (index + 1) * 20] += rng.poisson(
            3.0, size=((identities == identity).sum(), 20)
        )
    adata = AnnData(
        X=sparse.csr_matrix(matrix),
        obs=pd.DataFrame(
            {"protocol": protocols, "mixology_identity": identities},
            index=[f"cell_{idx}" for idx in range(len(protocols))],
        ),
    )
    adata.layers["counts"] = adata.X.copy()
    adata.write_h5ad(path)


def test_mixology_benchmark_emits_heldout_and_locked_evidence(tmp_path):
    source = tmp_path / "mixology.h5ad"
    _synthetic_mixology(source)

    report = run(
        source,
        tmp_path / "output",
        n_top_genes=45,
        n_pcs=8,
        seed=7,
        run_harmony=False,
    )

    assert report["status"] == "REVIEW"
    assert report["experimental_unit"] == "protocol"
    assert report["release_gate"]["status"] == "BLOCKED"
    assert report["candidate_acceptance"]["status"] in {"PASS", "FAIL"}
    assert report["product_policy"]["selected_candidate"] == "standard_unintegrated"
    assert report["product_policy"]["representation_contract"]["counts_preserved"] is True
    assert report["product_policy"]["representation_contract"]["normalized_full_present"] is True
    assert report["product_policy"]["representation_contract"]["discovery_feature_present"] is True
    assert report["product_policy"]["representation_contract"]["discovery_rep_present"] is True
    assert report["product_policy"]["representation_contract"]["integrated_rep_selected"] is False
    fold = pd.read_csv(report["artifacts"]["fold_metrics"], sep="\t")
    candidates = pd.read_csv(report["artifacts"]["candidate_metrics"], sep="\t")
    assert len(fold) == 12
    assert candidates["selected"].sum() == 1
    assert set(candidates["candidate"]) == {
        "standard_unintegrated",
        "pearson_residuals",
        "multinomial_deviance",
        "pearson_residuals_deviance",
    }
