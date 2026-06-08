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
from ..utils import StepResult
from .config import TumorAnalysisConfig

log = logging.getLogger(__name__)

__all__ = [
    "CNVInferenceStep",
    "MalignancyInterpretationStep",
    "TMEDeconvolutionStep",
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
        if ref_key and ref_key in adata.obs.columns:
            reference_key_to_use = ref_key
            ref_mask = (
                adata.obs[ref_key]
                .astype(str)
                .str.lower()
                .isin({"normal", "healthy", "reference", "immune", "stromal"})
            )
            if ref_mask.any():
                ref_cells = adata.obs.loc[ref_mask, ref_key].unique().tolist()

        adata = infer_cnv(
            adata,
            reference_cells=ref_cells,
            reference_key=reference_key_to_use,
            key_added=key_added,
        )
        self._outputs["key_added"] = key_added
        self._outputs["score_key"] = f"{key_added}_score"
        self._outputs["predicted_class_key"] = f"{key_added}_predicted_class"
        if f"{key_added}_summary" in adata.uns:
            self._outputs["summary"] = adata.uns[f"{key_added}_summary"]
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
            annotation_key = "cell_type_auto" if "cell_type_auto" in adata.obs.columns else "cell_type"
        if annotation_key not in adata.obs.columns:
            # Try fallback keys
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
            # Try fallback to common keys
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
        self._outputs["claim"] = (
            "annotation-derived TME composition; not a bulk deconvolution"
        )
        self._results = self.get_summary()
        return adata


# Register adapters with the factory so they can be instantiated by name.
AnalysisStepFactory.register("cnv_inference", CNVInferenceStep)
AnalysisStepFactory.register("malignancy_interpretation", MalignancyInterpretationStep)
AnalysisStepFactory.register("tme_deconvolution", TMEDeconvolutionStep)
