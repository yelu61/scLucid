"""
scLucid: A Comprehensive System for Single-Cell Analysis
=========================================================

scLucid is a powerful and flexible Python toolkit for the analysis of
single-cell RNA-sequencing data.
"""

from importlib import import_module
from importlib.metadata import version, PackageNotFoundError
import warnings
from typing import Optional

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

# --- Configuration and Settings ---
try:
    from .settings import setup_logging, set_figure_params, reset_figure_params
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

from .config import get_config, set_config, reset_config

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

_qc_workflow = _import_optional("qc.workflow")
if _qc_workflow is not None:
    run_standard_qc = getattr(_qc_workflow, "run_standard_qc", None)
    run_advanced_qc = getattr(_qc_workflow, "run_advanced_qc", None)

_preprocess_workflow = _import_optional("preprocess.workflow")
if _preprocess_workflow is not None:
    run_preprocessing = getattr(_preprocess_workflow, "run_preprocessing", None)

if analysis is not None:
    run_annotation = getattr(analysis, "run_annotation", None)
    characterize_clusters = getattr(analysis, "characterize_clusters", None)

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
    "setup_logging",
    "set_figure_params",
    "reset_figure_params",
    "FONT_NATURE",
    "FONT_CELL",
    "FONT_TRADITIONAL",
    "run_standard_qc",
    "run_advanced_qc",
    "run_preprocessing",
    "run_annotation",
    "characterize_clusters",
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
]
