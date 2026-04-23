"""Adapter for generating tumor analysis recommendations."""

from __future__ import annotations

import importlib.util
import logging
from typing import List, Optional

import numpy as np
from anndata import AnnData

from ..tumor.config import TumorAnalysisConfig
from .schema import ParameterRecommendation, RecommendationSection

log = logging.getLogger(__name__)


def adapt_tumor_recommendation(
    adata: AnnData,
    config: Optional[TumorAnalysisConfig] = None,
    cancer_type: Optional[str] = None,
) -> RecommendationSection:
    """Build a RecommendationSection for tumor-specific analysis."""
    if config is None:
        config = TumorAnalysisConfig()

    parameters: List[ParameterRecommendation] = []
    concerns: List[str] = []
    notes: List[str] = []

    n_cells = int(adata.n_obs)
    has_cnv_score = "cnv_score" in adata.obs.columns or "cnv" in adata.obs.columns
    infercnvpy_available = importlib.util.find_spec("infercnvpy") is not None

    # Malignancy recommendation
    run_malignancy = config.run_malignancy
    malignancy_confidence = 0.8
    if n_cells < 500:
        run_malignancy = False
        concerns.append("Dataset too small for reliable malignancy classification (<500 cells).")
        malignancy_confidence = 0.4
    elif n_cells < 1000:
        notes.append("Small dataset; malignancy classification may have reduced power.")
        malignancy_confidence = 0.65

    parameters.append(
        ParameterRecommendation(
            name="run_malignancy",
            value=run_malignancy,
            method="dataset_size_heuristic",
            confidence=malignancy_confidence,
            rationale="Malignancy scoring is recommended for tumor datasets with sufficient cells.",
            evidence={"n_cells": n_cells},
        )
    )

    # TME recommendation
    run_tme = config.run_tme
    tme_confidence = 0.75
    if n_cells < 300:
        run_tme = False
        concerns.append("Dataset too small for meaningful TME deconvolution (<300 cells).")
        tme_confidence = 0.4
    cell_type_key = config.tme_cell_type_key
    if cell_type_key not in adata.obs.columns:
        if "cell_type" in adata.obs.columns:
            cell_type_key = "cell_type"
            notes.append("Using 'cell_type' as fallback for TME deconvolution.")
        else:
            run_tme = False
            concerns.append("No cell type annotations found for TME deconvolution.")
            tme_confidence = 0.3

    parameters.append(
        ParameterRecommendation(
            name="run_tme",
            value=run_tme,
            method="dataset_size_and_annotation_check",
            confidence=tme_confidence,
            rationale="TME deconvolution is recommended when cell type annotations are available.",
            evidence={"n_cells": n_cells, "cell_type_key": cell_type_key},
        )
    )
    if cell_type_key != config.tme_cell_type_key:
        parameters.append(
            ParameterRecommendation(
                name="tme_cell_type_key",
                value=cell_type_key,
                method="fallback_annotation_key",
                confidence=0.7,
                rationale="Fallback cell type key selected for TME analysis.",
            )
        )

    # CNV recommendation
    run_cnv = config.run_cnv
    cnv_confidence = 0.7
    if has_cnv_score:
        run_cnv = True
        cnv_confidence = 0.9
        notes.append("Existing CNV scores detected; enabling CNV analysis.")
    elif infercnvpy_available:
        if config.cnv_reference_key is None or (config.cnv_reference_key not in adata.obs.columns):
            concerns.append(
                "CNV inference recommended but no reference cell key specified. "
                "Consider setting cnv_reference_key."
            )
            cnv_confidence = 0.55
    else:
        run_cnv = False
        cnv_confidence = 0.5
        notes.append("infercnvpy is not available; CNV inference disabled.")

    parameters.append(
        ParameterRecommendation(
            name="run_cnv",
            value=run_cnv,
            method="dependency_and_data_check",
            confidence=cnv_confidence,
            rationale="CNV inference is enabled when infercnvpy is available or CNV scores exist.",
            evidence={
                "has_cnv_score": has_cnv_score,
                "infercnvpy_available": infercnvpy_available,
            },
        )
    )

    # Therapy recommendation
    run_therapy = config.run_therapy
    therapy_confidence = 0.6
    if cancer_type is None:
        run_therapy = False
        therapy_confidence = 0.4
        concerns.append("Therapy response prediction requires cancer_type to be specified.")
    else:
        notes.append(f"Cancer type '{cancer_type}' provided for therapy prediction.")

    parameters.append(
        ParameterRecommendation(
            name="run_therapy",
            value=run_therapy,
            method="cancer_type_availability",
            confidence=therapy_confidence,
            rationale="Therapy prediction is recommended when a specific cancer type is provided.",
            evidence={"cancer_type": cancer_type},
        )
    )

    # Malignancy method recommendation
    recommended_method = config.malignancy_method
    if recommended_method == "cnv" and not run_cnv:
        recommended_method = "threshold"
        notes.append(
            "CNV-based malignancy method unavailable without CNV inference; fallback to threshold."
        )

    parameters.append(
        ParameterRecommendation(
            name="malignancy_method",
            value=recommended_method,
            method="cnv_availability_fallback",
            confidence=0.75 if recommended_method == config.malignancy_method else 0.6,
            rationale="Malignancy classification method selected based on available data.",
            alternatives=["cnv", "threshold", "ml"],
        )
    )

    overall_confidence = float(
        np.mean([malignancy_confidence, tme_confidence, cnv_confidence, therapy_confidence])
    )

    recommended_config = config.model_copy(
        update={
            "run_malignancy": run_malignancy,
            "run_tme": run_tme,
            "run_cnv": run_cnv,
            "run_therapy": run_therapy,
            "malignancy_method": recommended_method,
            "tme_cell_type_key": cell_type_key,
        }
    )

    return RecommendationSection(
        name="tumor",
        summary=(
            f"Tumor analysis: malignancy={run_malignancy}, tme={run_tme}, "
            f"cnv={run_cnv}, therapy={run_therapy}."
        ),
        confidence=overall_confidence,
        parameters=parameters,
        concerns=concerns,
        notes=notes,
        raw_result=recommended_config,
        metadata={
            "cancer_type": cancer_type,
            "n_cells": n_cells,
            "infercnvpy_available": infercnvpy_available,
            "has_cnv_score": has_cnv_score,
        },
    )
