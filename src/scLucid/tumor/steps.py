"""Lightweight AnalysisStep adapters for the tumor stage.

These adapters wrap existing tumor functions so the tumor workflow can use the
shared ``AnalysisStep`` interface without a full refactor of the underlying
implementations. They provide:

- ``validate_input``: pre-flight checks with structured failure reasons
- ``run``: execute the wrapped function and capture outputs
- ``get_summary``: return a compact, audit-friendly summary

The adapters are registered with ``AnalysisStepFactory`` so callers can
instantiate them by name if desired.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from anndata import AnnData

from ..base_interfaces import AnalysisStep, AnalysisStepFactory
from ..utils import StepResult, step_results_to_storage
from .config import TumorAnalysisConfig

log = logging.getLogger(__name__)

__all__ = [
    "CNVInferenceStep",
    "MalignancyInterpretationStep",
    "MalignancyScoringStep",
    "TMEDeconvolutionStep",
    "TherapyPredictionStep",
]


def _get_cfg_attr(config: Any, name: str, default: Any = None) -> Any:
    """Safely get an attribute from a config object or dict."""
    if config is None:
        return default
    if isinstance(config, dict):
        return config.get(name, default)
    return getattr(config, name, default)


class _TumorAdapterBase(AnalysisStep):
    """Base adapter providing common validation and summary patterns."""

    step_name: str = "tumor_step"
    default_evidence_level: str = "heuristic"

    def __init__(self, config: Optional[TumorAnalysisConfig] = None):
        super().__init__(config=config)
        self._validation_warnings: List[str] = []
        self._validation_errors: List[str] = []
        self._outputs: Dict[str, Any] = {}

    def validate_input(self, adata: AnnData) -> bool:
        """Default validation; subclasses should extend."""
        self._validation_warnings = []
        self._validation_errors = []
        return True

    def get_summary(self) -> Dict[str, Any]:
        """Return an audit-friendly summary of the last run."""
        return {
            "step": self.step_name,
            "status": "completed" if self._results is not None else "not_run",
            "outputs": dict(self._outputs),
            "warnings": list(self._validation_warnings),
            "errors": list(self._validation_errors),
        }

    def make_step_result(self, status: str = "completed") -> StepResult:
        return StepResult(
            name=self.step_name,
            status=status,  # type: ignore[arg-type]
            evidence_level=self.default_evidence_level,  # type: ignore[arg-type]
            outputs=dict(self._outputs),
            warnings=list(self._validation_warnings),
        )


class MalignancyScoringStep(_TumorAdapterBase):
    """Adapter for malignancy scoring + classification."""

    step_name = "malignancy_scoring"
    default_evidence_level = "heuristic"

    def validate_input(self, adata: AnnData) -> bool:
        valid = super().validate_input(adata)
        # score_malignancy is robust to missing markers; no hard pre-checks needed.
        return valid

    def run(self, adata: AnnData, **kwargs: Any) -> AnnData:
        from .malignancy.classification import classify_malignant_cells
        from .malignancy.scoring import score_malignancy

        self.validate_input(adata)
        cancer_type = kwargs.get(
            "cancer_type", _get_cfg_attr(self.config, "cancer_type")
        )
        method = _get_cfg_attr(self.config, "malignancy_method", "combined")
        ref_key = _get_cfg_attr(self.config, "malignancy_reference_key")

        log.info("Tumor stage: scoring malignancy")
        adata = score_malignancy(adata, key_added="malignancy", cancer_type=cancer_type)

        log.info("Tumor stage: classifying malignant cells")
        ref_adata = None
        if method in ("threshold", "ml") and ref_key and ref_key in adata.obs.columns:
            ref_mask = (
                adata.obs[ref_key]
                .astype(str)
                .str.lower()
                .isin({"normal", "healthy", "reference", "immune", "stromal"})
            )
            if ref_mask.any():
                ref_adata = adata[ref_mask].copy()
            else:
                self._validation_warnings.append(
                    f"No reference cells found via '{ref_key}'. Falling back to unsupervised."
                )
        elif method in ("threshold", "ml"):
            self._validation_warnings.append(
                f"malignancy_method='{method}' may require reference cells."
            )

        classify_malignant_cells(
            adata,
            method=method,
            reference_adata=ref_adata,
            key_added="is_malignant",
        )

        self._outputs["method"] = method
        self._outputs["score_key"] = "malignancy"
        self._outputs["classification_key"] = "is_malignant"
        self._outputs["reference_used"] = ref_adata is not None
        self._results = self.get_summary()
        return adata


class CNVInferenceStep(_TumorAdapterBase):
    """Adapter for expression-based CNV inference."""

    step_name = "cnv_inference"
    default_evidence_level = "heuristic"

    def validate_input(self, adata: AnnData) -> bool:
        valid = super().validate_input(adata)

        has_counts = False
        count_layer: Optional[str] = None
        for layer in ("counts", "raw_counts", "raw"):
            if layer in adata.layers:
                has_counts = True
                count_layer = layer
                break
        if not has_counts:
            import numpy as np
            from scipy.sparse import issparse

            sample = adata.X[: min(100, adata.n_obs), : min(100, adata.n_vars)]
            dense = sample.toarray() if issparse(sample) else np.asarray(sample)
            has_counts = bool(np.max(dense) >= 20 and not np.any(dense % 1 != 0))

        if not has_counts:
            self._validation_warnings.append(
                "No clear raw count layer found; CNV inference requires raw counts."
            )

        has_coords = all(c in adata.var.columns for c in ("chromosome", "start", "end"))
        if not has_coords:
            self._validation_warnings.append(
                "Genomic coordinates missing; CNV inference will run without "
                "chromosome-aware smoothing."
            )

        ref_key = _get_cfg_attr(self.config, "cnv_reference_key")
        if ref_key and ref_key not in adata.obs.columns:
            self._validation_warnings.append(
                f"Reference key '{ref_key}' not found in adata.obs."
            )

        self._outputs["has_count_input"] = has_counts
        self._outputs["count_layer_used"] = count_layer
        self._outputs["has_genomic_coordinates"] = has_coords
        return valid

    def run(self, adata: AnnData, **kwargs: Any) -> AnnData:
        from .cnv.infercnv import infer_cnv

        self.validate_input(adata)
        ref_key = _get_cfg_attr(self.config, "cnv_reference_key")
        key_added = kwargs.get("key_added", "cnv")

        ref_cells: Optional[Any] = None
        reference_key_to_use = "cell_type"
        has_ref = False
        if ref_key and ref_key in adata.obs.columns:
            reference_key_to_use = ref_key
            ref_mask = (
                adata.obs[ref_key]
                .astype(str)
                .str.lower()
                .isin({"normal", "healthy", "reference", "immune", "stromal"})
            )
            has_ref = bool(ref_mask.any())
            if has_ref:
                ref_cells = adata.obs.loc[ref_mask, ref_key].unique().tolist()

        adata = infer_cnv(
            adata,
            reference_cells=ref_cells,
            reference_key=reference_key_to_use,
            key_added=key_added,
        )
        cnv_summary = adata.uns.get(f"{key_added}_summary", {})
        input_quality = cnv_summary.get("input_quality", {})
        has_coords = input_quality.get(
            "has_genomic_coordinates",
            all(c in adata.var.columns for c in ("chromosome", "start", "end")),
        )
        evidence_level = "validated_core" if (has_ref and has_coords) else "heuristic"
        self.default_evidence_level = evidence_level  # type: ignore[misc]

        self._outputs["key_added"] = key_added
        self._outputs["score_key"] = f"{key_added}_score"
        self._outputs["predicted_class_key"] = f"{key_added}_predicted_class"
        self._outputs["n_aneuploid"] = cnv_summary.get("n_aneuploid")
        self._outputs["n_diploid"] = cnv_summary.get("n_diploid")
        self._outputs["threshold"] = cnv_summary.get("threshold")
        self._outputs["mean_cnv_score"] = cnv_summary.get("mean_cnv_score")
        self._outputs["has_reference_cells"] = has_ref
        self._outputs["has_genomic_coordinates"] = has_coords
        if not has_ref:
            self._validation_warnings.append(
                "No reference cells provided; threshold-based calls less reliable."
            )
        self._results = self.get_summary()
        return adata


class MalignancyInterpretationStep(_TumorAdapterBase):
    """Adapter for malignancy interpretation."""

    step_name = "malignancy_interpretation"
    default_evidence_level = "heuristic"

    def validate_input(self, adata: AnnData) -> bool:
        valid = super().validate_input(adata)
        annotation_key = _get_cfg_attr(
            self.config, "malignancy_annotation_key", None
        )
        if annotation_key is None:
            annotation_key = (
                "cell_type_auto"
                if "cell_type_auto" in adata.obs.columns
                else "cell_type"
            )
        if annotation_key not in adata.obs.columns:
            for fallback in ("cell_type_auto", "cell_type"):
                if fallback in adata.obs.columns:
                    self._validation_warnings.append(
                        f"Annotation key '{annotation_key}' not found; using '{fallback}'."
                    )
                    annotation_key = fallback
                    break
            else:
                self._validation_errors.append(
                    f"Annotation key '{annotation_key}' not found in adata.obs."
                )
                valid = False
        self._outputs["annotation_key"] = annotation_key
        return valid

    def run(self, adata: AnnData, **kwargs: Any) -> AnnData:
        from .malignancy.interpretation import run_malignancy_interpretation

        self.validate_input(adata)
        annotation_key = self._outputs.get(
            "annotation_key",
            kwargs.get("annotation_key", "cell_type_auto"),
        )
        cluster_key = kwargs.get("cluster_key", None)
        cancer_type = kwargs.get(
            "cancer_type", _get_cfg_attr(self.config, "cancer_type")
        )

        adata = run_malignancy_interpretation(
            adata,
            annotation_key=annotation_key,
            cluster_key=cluster_key,
            cancer_type=cancer_type,
            run_cnv=bool(_get_cfg_attr(self.config, "run_cnv", True)),
            run_malignancy_score=True,
        )
        self._outputs["call_key"] = "malignancy_call"
        self._outputs["score_key"] = "malignancy_interpretation_score"
        malignancy_ns = (
            adata.uns.get("sclucid", {}).get("analysis", {}).get("malignancy", {})
        )
        summary = malignancy_ns.get("malignancy_interpretation_summary", {})
        evidence_sources = summary.get("evidence_sources", [])
        if len(evidence_sources) >= 3:
            self.default_evidence_level = "validated_core"  # type: ignore[misc]

        if summary:
            self._outputs["summary"] = {
                k: v
                for k, v in summary.items()
                if k
                in (
                    "n_malignant",
                    "n_suspect_malignant",
                    "n_non_malignant",
                    "n_unresolved",
                    "tumor_purity_estimate",
                    "evidence_sources",
                    "review_required",
                    "mean_score",
                    "threshold",
                    "suspect_threshold",
                )
            }
        self._results = self.get_summary()
        return adata


class TMEDeconvolutionStep(_TumorAdapterBase):
    """Adapter for TME composition profiling."""

    step_name = "tme_deconvolution"
    default_evidence_level = "heuristic"

    def validate_input(self, adata: AnnData) -> bool:
        valid = super().validate_input(adata)
        cell_type_key = _get_cfg_attr(
            self.config, "tme_cell_type_key", "cell_type_auto"
        )
        if cell_type_key not in adata.obs.columns:
            for fallback in ("cell_type_auto", "cell_type"):
                if fallback in adata.obs.columns:
                    self._validation_warnings.append(
                        f"TME cell type key '{cell_type_key}' not found; using '{fallback}'."
                    )
                    cell_type_key = fallback
                    break
            else:
                self._validation_errors.append(
                    "No cell type annotation found for TME deconvolution."
                )
                valid = False
        self._outputs["cell_type_key"] = cell_type_key
        return valid

    def run(self, adata: AnnData, **kwargs: Any) -> AnnData:
        from .microenvironment.deconvolution import deconvolve_tme

        self.validate_input(adata)
        cell_type_key = self._outputs.get("cell_type_key", "cell_type_auto")
        key_added = kwargs.get("key_added", "tme")
        adata = deconvolve_tme(adata, cell_type_key=cell_type_key, key_added=key_added)
        self._outputs["key_added"] = key_added
        proportions = adata.uns.get(f"{key_added}_proportions", {})
        if hasattr(proportions, "to_dict"):
            proportions = proportions.to_dict()
        self._outputs["proportions"] = proportions
        self._outputs["immune_score"] = float(
            adata.uns.get(f"{key_added}_immune_score", 0.0)
        )
        self._outputs["stromal_score"] = float(
            adata.uns.get(f"{key_added}_stromal_score", 0.0)
        )
        self._outputs["malignant_score"] = float(
            adata.uns.get(f"{key_added}_malignant_score", 0.0)
        )
        self._outputs["claim"] = (
            "annotation-derived TME composition; not a bulk deconvolution"
        )
        self._validation_warnings.append(
            "TME composition is annotation-derived, not bulk deconvolution."
        )
        self._results = self.get_summary()
        return adata


class TherapyPredictionStep(_TumorAdapterBase):
    """Adapter for signature-based therapy response prediction."""

    step_name = "therapy_prediction"
    default_evidence_level = "exploratory"

    def validate_input(self, adata: AnnData) -> bool:
        valid = super().validate_input(adata)
        drugs = _get_cfg_attr(self.config, "therapy_drugs") or ["chemotherapy"]
        if not drugs:
            self._validation_warnings.append("No therapy drugs configured.")
        self._outputs["drugs_requested"] = list(drugs)
        return valid

    def run(self, adata: AnnData, **kwargs: Any) -> AnnData:
        from .therapy.prediction import predict_therapy_response

        self.validate_input(adata)
        drugs = _get_cfg_attr(self.config, "therapy_drugs") or ["chemotherapy"]
        drug_results: List[StepResult] = []

        for drug in drugs:
            drug_step_name = f"therapy_prediction_{drug}"
            try:
                predict_therapy_response(
                    adata,
                    therapy_type=drug,
                    method="signature",
                    key_added=f"therapy_response_{drug}",
                )
                drug_results.append(
                    StepResult(
                        name=drug_step_name,
                        status="completed",
                        evidence_level="exploratory",
                        outputs={"drug": drug, "key_added": f"therapy_response_{drug}"},
                        warnings=["Signature-based therapy prediction is exploratory."],
                    )
                )
            except Exception as drug_exc:
                msg = f"Therapy prediction failed for {drug}: {drug_exc}"
                log.warning(f"{msg}. Skipping drug.")
                drug_results.append(
                    StepResult.from_exception(
                        name=drug_step_name,
                        exc=drug_exc,
                        degraded=True,
                        evidence_level="unavailable",
                    )
                )

        n_completed = sum(1 for r in drug_results if r.status == "completed")
        status = "completed" if n_completed == len(drug_results) else "degraded"
        self._outputs["drugs_completed"] = n_completed
        self._outputs["drug_results"] = [r.name for r in drug_results]
        self._outputs["per_drug_step_results"] = step_results_to_storage(drug_results)
        self._validation_warnings.append(
            "Therapy prediction uses built-in signatures only."
        )
        self._results = self.get_summary()
        # Expose per-drug results as extra step results via a side channel
        self._drug_step_results = drug_results
        return adata


# Register adapters with the factory so they can be instantiated by name.
AnalysisStepFactory.register("malignancy_scoring", MalignancyScoringStep)
AnalysisStepFactory.register("cnv_inference", CNVInferenceStep)
AnalysisStepFactory.register("malignancy_interpretation", MalignancyInterpretationStep)
AnalysisStepFactory.register("tme_deconvolution", TMEDeconvolutionStep)
AnalysisStepFactory.register("therapy_prediction", TherapyPredictionStep)
