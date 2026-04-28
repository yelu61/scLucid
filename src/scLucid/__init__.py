"""
scLucid: A Comprehensive System for Single-Cell Analysis
=========================================================

scLucid is a powerful and flexible Python toolkit for the analysis of
single-cell RNA-sequencing data.
"""

import os
import warnings
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from typing import Any, Dict, Optional

from .runtime import setup_runtime_environment

setup_runtime_environment()

# Default to non-interactive plotting to prevent pop-up windows in scripts/CI.
# Users can enable interactivity via scl.set_interactive_mode(True).
if not os.environ.get("MPLBACKEND"):
    import matplotlib

    matplotlib.use("Agg", force=False)
    warnings.filterwarnings("ignore", message="FigureCanvasAgg is non-interactive")

try:
    __version__ = version("sclucid")
except PackageNotFoundError:
    __version__ = "0.1"


def _import_optional(module: str, *, hint: Optional[str] = None):
    """Import module defensively so package import remains usable."""
    try:
        return import_module(f".{module}", __name__)
    except Exception as exc:
        msg = f"Could not import module '{module}': {exc}"
        if hint:
            msg = f"{msg}. {hint}"
        warnings.warn(msg, ImportWarning)
        return None


# --- Core Modules ---
qc = _import_optional("qc")
preprocess = _import_optional("preprocess")
analysis = _import_optional("analysis")
plotting = _import_optional("plotting")
utils = _import_optional("utils")
recommendation = _import_optional("recommendation")
tools = _import_optional(
    "tools",
    hint="Install with 'pip install sclucid[tools]' to use advanced features",
)

# --- Optional Web Module ---
_web = _import_optional(
    "web",
    hint="Install with 'pip install sclucid[web]' to use web features",
)
launch_web_app = getattr(_web, "launch_web_app", None)

# --- Convenient Aliases ---
pp = preprocess
al = analysis
tl = tools
ut = utils
pl = plotting
rc = recommendation

# --- Configuration and Settings ---
try:
    from .settings import (
        is_interactive_mode,
        reset_figure_params,
        set_figure_params,
        set_interactive_mode,
        setup_logging,
    )
except Exception as exc:
    warnings.warn(f"Could not import settings module: {exc}", ImportWarning)

    def _settings_unavailable(*args, **kwargs):
        raise RuntimeError(
            "scLucid settings are unavailable because optional plotting dependencies "
            "failed to import."
        )

    setup_logging = _settings_unavailable
    set_figure_params = _settings_unavailable
    reset_figure_params = _settings_unavailable
    set_interactive_mode = _settings_unavailable

    def is_interactive_mode() -> bool:
        return False

from .config import get_config, reset_config, set_config
from .utils.context import (
    AnalysisContext,
    DatasetProfile,
    infer_analysis_context,
    infer_dataset_profile,
    normalize_dataset_type,
)
from .utils.contracts import (
    Modules,
    UnsKeys,
    build_config_lineage,
    record_contract_result,
    validate_stage_contract,
)
from .utils.validation import (
    ValidationError,
    assert_analysis_ready,
    assert_preprocessing_ready,
    assert_qc_ready,
)


def _stage_kwargs(prefix: str, kwargs: Dict[str, Any]) -> Dict[str, Any]:
    """Return kwargs for one pipeline stage with the stage prefix stripped."""
    token = f"{prefix}_"
    return {key[len(token) :]: value for key, value in kwargs.items() if key.startswith(token)}


# --- Abstract Base Classes for Extensibility ---
AnalysisStep = None
QCFilter = None
CellAnnotator = None
ScoringMethod = None
PlottingBackend = None
ProportionAnalysisMethod = None
AnalysisStepFactory = None

_base_interfaces = _import_optional("base_interfaces")
if _base_interfaces is not None:
    AnalysisStep = getattr(_base_interfaces, "AnalysisStep", None)
    QCFilter = getattr(_base_interfaces, "QCFilter", None)
    CellAnnotator = getattr(_base_interfaces, "CellAnnotator", None)
    ScoringMethod = getattr(_base_interfaces, "ScoringMethod", None)
    PlottingBackend = getattr(_base_interfaces, "PlottingBackend", None)
    ProportionAnalysisMethod = getattr(_base_interfaces, "ProportionMethod", None)
    AnalysisStepFactory = getattr(_base_interfaces, "AnalysisStepFactory", None)

# --- Academic Font Style Constants ---
FONT_NATURE = "nature"
FONT_CELL = "cell"
FONT_TRADITIONAL = "traditional"

# --- High-Level Workflows ---
run_standard_qc = None
run_advanced_qc = None
run_preprocessing = None
run_annotation = None
characterize_clusters = None
recommend_analysis_parameters = None
run_tumor_analysis = None

_qc_workflow = _import_optional("qc.workflow")
if _qc_workflow is not None:
    run_standard_qc = getattr(_qc_workflow, "run_standard_qc", None)
    run_advanced_qc = getattr(_qc_workflow, "run_advanced_qc", None)

_preprocess_workflow = _import_optional("preprocess.workflow")
if _preprocess_workflow is not None:
    run_preprocessing = getattr(_preprocess_workflow, "run_preprocessing", None)

_analysis_workflow = _import_optional("analysis.workflow")
if _analysis_workflow is not None:
    run_standard_analysis = getattr(_analysis_workflow, "run_standard_analysis", None)
    run_custom_analysis = getattr(_analysis_workflow, "run_custom_analysis", None)

if analysis is not None:
    run_annotation = getattr(analysis, "run_annotation", None)
    characterize_clusters = getattr(analysis, "characterize_clusters", None)

if recommendation is not None:
    recommend_analysis_parameters = getattr(recommendation, "recommend_analysis_parameters", None)

_tumor_workflow = _import_optional("tumor.workflow")
if _tumor_workflow is not None:
    run_tumor_analysis = getattr(_tumor_workflow, "run_tumor_analysis", None)


# --- Unified Pipeline Entry Point ---


def run_pipeline(
    adata,
    stages: list[str] = None,
    *,
    context=None,
    dataset_type: str = None,
    qc_config=None,
    preprocess_config=None,
    analysis_config=None,
    show_progress: bool = True,
    tissue_type: str = "unknown",
    tissue: str = None,
    species: str = "human",
    cancer_type: str = None,
    sample_key: str = None,
    batch_key: str = None,
    condition_key: str = None,
    **kwargs,
):
    """
    Run a multi-stage analysis pipeline in sequence.

    This is the unified framework entry point that chains QC, preprocessing,
    and analysis stages with automatic state bridging between modules.

    Parameters
    ----------
    adata : AnnData
        Input data matrix
    stages : list of str, optional
        Stages to run in order. Valid stages: ``"qc"``, ``"preprocess"``, ``"analysis"``.
        Default: all three.
    qc_config : QCWorkflowConfig, optional
        Configuration for QC stage
    preprocess_config : PreprocessingWorkflowConfig, optional
        Configuration for preprocessing stage
    analysis_config : AnalysisWorkflowConfig, optional
        Configuration for analysis stage
    show_progress : bool, default=True
        Show progress bars
    context : AnalysisContext or dict, optional
        Shared dataset context used to tune defaults and document assumptions.
    dataset_type : str, optional
        Canonical or alias dataset type, e.g. ``"pbmc_or_blood"``,
        ``"normal_tissue"``, ``"tumor_tissue"``, ``"cell_line"``,
        ``"organoid"``, or ``"spatial"``. Multi-sample status is tracked
        separately as ``AnalysisContext.is_multi_sample``.
    tissue_type : str, default="unknown"
        Backward-compatible tissue type hint passed to QC when no context is provided.
    **kwargs
        Additional parameters passed to individual stages

    Returns:
    -------
    AnnData
        Processed data with all stage results stored under
        ``adata.uns["sclucid"][stage]``.

    Examples:
    --------
    >>> # Run full pipeline
    >>> adata = scl.run_pipeline(adata, stages=["qc", "preprocess", "analysis"])

    >>> # Run only QC and preprocessing
    >>> adata = scl.run_pipeline(adata, stages=["qc", "preprocess"])

    >>> # Skip QC, start from preprocessed data
    >>> adata = scl.run_pipeline(adata, stages=["analysis"])
    """
    if stages is None:
        stages = ["qc", "preprocess", "analysis"]

    valid_stages = {"qc", "preprocess", "analysis"}
    invalid = set(stages) - valid_stages
    if invalid:
        raise ValueError(f"Invalid stages: {invalid}. Valid stages are: {valid_stages}")

    global_config = get_config()
    effective_dataset_type = (
        dataset_type
        if dataset_type is not None
        else getattr(global_config, "default_dataset_type", "unknown")
    )
    effective_species = species or getattr(global_config, "default_species", "human")

    analysis_context = infer_analysis_context(
        adata,
        context=context,
        dataset_type=effective_dataset_type,
        species=effective_species,
        tissue=tissue,
        tissue_type=tissue_type,
        cancer_type=cancer_type,
        sample_key=sample_key,
        batch_key=batch_key,
        condition_key=condition_key,
    )

    # Cross-module context to propagate settings across stages
    pipeline_context: Dict[str, Any] = analysis_context.to_dict()
    pipeline_context["config_lineage"] = build_config_lineage(
        global_config=global_config.to_dict() if hasattr(global_config, "to_dict") else {},
        inherited=analysis_context.to_dict(),
    )

    # Run QC
    if "qc" in stages:
        if run_standard_qc is None:
            raise RuntimeError("QC module not available")
        assert_qc_ready(adata)
        adata = run_standard_qc(
            adata,
            config=qc_config,
            tissue_type=analysis_context.qc_tissue_type,
            show_progress=show_progress,
            **_stage_kwargs("qc", kwargs),
        )
        qc_contract = validate_stage_contract(adata, "qc", when="output")
        record_contract_result(adata, Modules.QC, qc_contract)
        # Capture upstream context from QC results
        qc_uns = adata.uns.get("sclucid", {}).get(Modules.QC, {})
        qc_config_dict = qc_uns.get(UnsKeys.WORKFLOW_CONFIG, {})
        pipeline_context["species"] = qc_config_dict.get(
            "species", getattr(qc_config, "species", analysis_context.species)
        )
        pipeline_context["tissue_type"] = analysis_context.qc_tissue_type

    # Run Preprocessing
    if "preprocess" in stages:
        if run_preprocessing is None:
            raise RuntimeError("Preprocess module not available")
        if "qc" not in stages:
            try:
                assert_preprocessing_ready(adata)
            except ValidationError as e:
                raise RuntimeError(
                    f"Preprocessing stage requires QC results, but QC was skipped. "
                    f"Include 'qc' in stages or run QC manually first. Details: {e}"
                )
        adata = run_preprocessing(
            adata,
            config=preprocess_config,
            tissue_type=analysis_context.qc_tissue_type,
            show_progress=show_progress,
            **_stage_kwargs("preprocess", kwargs),
        )
        preprocess_contract = validate_stage_contract(adata, "preprocess", when="output")
        record_contract_result(adata, Modules.PREPROCESS, preprocess_contract)
        # Capture upstream context from preprocess results
        pp_uns = adata.uns.get("sclucid", {}).get(Modules.PREPROCESS, {})
        pp_config_dict = pp_uns.get(UnsKeys.WORKFLOW_CONFIG, {})
        integration_cfg = pp_config_dict.get("integration", {})
        if integration_cfg.get("batch_key"):
            pipeline_context["batch_key"] = integration_cfg["batch_key"]
        if integration_cfg.get("method"):
            pipeline_context["integration_method"] = integration_cfg["method"]

    # Run Analysis
    if "analysis" in stages:
        if run_standard_analysis is None:
            raise RuntimeError("Analysis module not available")
        if "preprocess" not in stages:
            try:
                assert_analysis_ready(adata)
            except ValidationError as e:
                raise RuntimeError(
                    f"Analysis stage requires preprocessing results, but preprocessing was skipped. "
                    f"Include 'preprocess' in stages or run preprocessing manually first. Details: {e}"
                )
        # Inherit species from QC context if analysis_config is not explicitly provided
        effective_analysis_config = analysis_config
        if effective_analysis_config is None and "species" in pipeline_context:
            from .analysis.config import AnalysisWorkflowConfig

            effective_analysis_config = AnalysisWorkflowConfig()
            if effective_analysis_config.annotation is not None:
                new_annotation = effective_analysis_config.annotation.model_copy(
                    update={"marker_species": pipeline_context["species"]}
                )
                effective_analysis_config = effective_analysis_config.model_copy(
                    update={"annotation": new_annotation}
                )
        adata = run_standard_analysis(
            adata,
            config=effective_analysis_config,
            show_progress=show_progress,
            **_stage_kwargs("analysis", kwargs),
        )
        analysis_contract = validate_stage_contract(adata, "analysis", when="output")
        record_contract_result(adata, Modules.ANALYSIS, analysis_contract)

    # Store pipeline context for downstream inspection
    adata.uns.setdefault("sclucid", {})[UnsKeys.PIPELINE_CONTEXT] = pipeline_context
    adata.uns.setdefault("sclucid", {})[UnsKeys.ANALYSIS_CONTEXT] = analysis_context.to_dict()

    return adata


__all__ = [
    "pp",
    "al",
    "tl",
    "ut",
    "pl",
    "qc",
    "preprocess",
    "analysis",
    "tools",
    "plotting",
    "utils",
    "recommendation",
    "AnalysisContext",
    "DatasetProfile",
    "infer_analysis_context",
    "infer_dataset_profile",
    "normalize_dataset_type",
    "setup_logging",
    "set_figure_params",
    "reset_figure_params",
    "set_interactive_mode",
    "is_interactive_mode",
    "FONT_NATURE",
    "FONT_CELL",
    "FONT_TRADITIONAL",
    "run_standard_qc",
    "run_advanced_qc",
    "run_preprocessing",
    "run_standard_analysis",
    "run_custom_analysis",
    "run_pipeline",
    "run_annotation",
    "characterize_clusters",
    "recommend_analysis_parameters",
    "get_config",
    "set_config",
    "reset_config",
    "AnalysisStep",
    "QCFilter",
    "CellAnnotator",
    "ScoringMethod",
    "PlottingBackend",
    "ProportionAnalysisMethod",
    "AnalysisStepFactory",
    "launch_web_app",
    "rc",
]
