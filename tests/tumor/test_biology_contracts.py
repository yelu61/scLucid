"""Biology-oriented contract tests for tumor outputs.

These tests assert directional biological behavior on synthetic data. They are
not benchmarks; they guard against regressions where a workflow output remains
present but no longer carries the expected biological signal.
"""

import numpy as np
import pandas as pd
from anndata import AnnData


def test_cnv_inference_detects_synthetic_aneuploid_shift():
    """Tumor cells with chromosome-wide amplification should look more aneuploid."""
    from scLucid.tumor.cnv.infercnv import infer_cnv

    rng = np.random.default_rng(42)
    n_ref = 60
    n_tumor = 60
    n_genes = 120
    X = rng.poisson(8, size=(n_ref + n_tumor, n_genes)).astype(np.float32)
    X[n_ref:, :60] += rng.poisson(12, size=(n_tumor, 60)).astype(np.float32)

    adata = AnnData(X)
    adata.obs_names = [f"cell_{i:03d}" for i in range(adata.n_obs)]
    adata.var_names = [f"gene_{i:03d}" for i in range(n_genes)]
    adata.layers["counts"] = X.copy()
    adata.obs["cell_type"] = pd.Categorical(["Normal"] * n_ref + ["Tumor"] * n_tumor)
    adata.var["chromosome"] = ["1"] * 60 + ["2"] * 60
    adata.var["start"] = list(range(n_genes))
    adata.var["end"] = list(range(1, n_genes + 1))

    result = infer_cnv(
        adata,
        reference_cells="Normal",
        reference_key="cell_type",
        threshold_mad=2.0,
        use_gmm=False,
    )

    obs = result.obs
    tumor = obs["cell_type"].astype(str) == "Tumor"
    normal = obs["cell_type"].astype(str) == "Normal"
    assert obs.loc[tumor, "cnv_extreme_frac"].mean() > obs.loc[normal, "cnv_extreme_frac"].mean()
    tumor_aneuploid = (obs.loc[tumor, "cnv_predicted_class"].astype(str) == "aneuploid").mean()
    normal_aneuploid = (obs.loc[normal, "cnv_predicted_class"].astype(str) == "aneuploid").mean()
    assert tumor_aneuploid > normal_aneuploid


def test_malignancy_interpretation_combines_annotation_cnv_and_signature_evidence():
    """High tumor annotation/CNV/signature evidence should call malignant cells."""
    from scLucid.tumor.malignancy.interpretation import run_malignancy_interpretation

    adata = AnnData(np.ones((6, 3), dtype=np.float32))
    adata.obs_names = [f"cell_{i}" for i in range(6)]
    adata.var_names = ["MKI67", "KRAS", "PTPRC"]
    adata.obs["cell_type_auto"] = pd.Categorical(
        ["T_cell", "T_cell", "immune", "tumor epithelial", "carcinoma", "malignant"]
    )
    adata.obs["cnv_score"] = [0.02, 0.05, 0.03, 0.85, 0.90, 0.95]
    adata.obs["analysis_malignancy"] = [0.05, 0.10, 0.08, 0.80, 0.88, 0.92]

    table = run_malignancy_interpretation(
        adata,
        annotation_key="cell_type_auto",
        cnv_score_key="cnv_score",
        run_cnv=False,
        run_malignancy_score=False,
        malignancy_score_key="analysis_malignancy",
        threshold=0.55,
        suspect_threshold=0.35,
    )

    calls = adata.obs["malignancy_call"].astype(str)
    assert (calls.iloc[:3] == "non_malignant").all()
    assert (calls.iloc[3:] == "malignant").all()
    assert {"annotation_prior", "cnv_score", "malignancy_signature_score"}.issubset(
        set(
            adata.uns["sclucid"]["analysis"]["malignancy"][
                "malignancy_interpretation_summary"
            ]["evidence_sources"]
        )
    )
    assert not table.empty


def test_tme_compartment_contract_keeps_unknown_labels_reviewable():
    """Unknown labels should become other/reviewable, not silently malignant."""
    from scLucid.tumor.microenvironment.deconvolution import deconvolve_tme

    adata = AnnData(np.zeros((5, 2), dtype=np.float32))
    adata.obs["cell_type"] = pd.Categorical(
        ["T_cell", "CAF", "normal_epithelial", "Tumor", "RareUnknown"]
    )

    result = deconvolve_tme(adata, copy=True)
    compartments = result.obs["tme_compartment"].astype(str).tolist()
    assert compartments == ["immune", "stromal", "other", "malignant", "other"]
    assert result.uns["tme_malignant_score"] == 1 / 5
    assert "rareunknown" in result.uns["tme_unmapped_labels"]
