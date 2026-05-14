"""Canonical AnnData contracts for scLucid workflows.

This module centralizes the names and structural expectations that define the
main scLucid workflow boundary. It is intentionally conservative: contracts
describe the stable core while modules can still store richer module-specific
details under ``adata.uns["sclucid"][module]``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime
from typing import Any, Dict, Iterable, Literal, Optional

from anndata import AnnData

SCHEMA_VERSION = "1.0"
SCLUCID_ROOT = "sclucid"
STAGE_ORDER = ("qc", "preprocess", "analysis")
API_LAYER_ORDER = ("workflow", "simple_api", "advanced")
REVIEW_SUMMARY_REQUIRED_KEYS = (
    "schema_version",
    "module",
    "workflow_name",
    "steps_executed",
    "data_shape",
)
REVIEW_SUMMARY_RECOMMENDED_KEYS = (
    "generated_at",
    "warnings",
    "config",
    "contract",
    "config_lineage",
    "artifacts",
)


class LayerKeys:
    """Canonical layer keys."""

    COUNTS = "counts"
    NORMALIZED = "normalized"
    SCALED = "scaled"


class ObsKeys:
    """Canonical obs column keys used across workflows."""

    SAMPLE = "sampleID"
    QC_N_GENES = "n_genes_by_counts"
    QC_TOTAL_COUNTS = "total_counts"
    QC_PCT_MT = "pct_counts_mt"
    LOW_QUALITY = "low_quality"
    CLUSTER = "leiden_clusters"
    CELL_TYPE = "cell_type_auto"


class VarKeys:
    """Canonical var column keys."""

    HIGHLY_VARIABLE = "highly_variable"


class ObsmKeys:
    """Canonical obsm keys."""

    PCA = "X_pca"
    UMAP = "X_umap"
    SPATIAL = "spatial"


class UnsKeys:
    """Canonical uns keys under adata.uns['sclucid'][module]."""

    NAMESPACE_METADATA = "_metadata"
    WORKFLOW_CONFIG = "workflow_config"
    STEPS_EXECUTED = "steps_executed"
    REVIEW_SUMMARY = "review_summary"
    PIPELINE_CONTEXT = "pipeline_context"
    ANALYSIS_CONTEXT = "analysis_context"
    CONFIG_LINEAGE = "config_lineage"
    CONTRACT = "contract"
    ARTIFACTS = "artifacts"
    ERRORS = "errors"


class Modules:
    """Canonical scLucid result namespaces."""

    QC = "qc"
    PREPROCESS = "preprocess"
    ANALYSIS = "analysis"
    TUMOR = "tumor"
    TOOLS = "tools"


StageName = Literal["qc", "preprocess", "analysis"]
APILayerName = Literal["workflow", "simple_api", "advanced"]


@dataclass(frozen=True)
class APILayerContract:
    """Stable public API layer contract."""

    name: APILayerName
    purpose: str
    primary_entrypoints: tuple[str, ...]
    example_artifacts: tuple[str, ...] = ()
    expected_outputs: tuple[str, ...] = ()


API_LAYER_CONTRACTS: dict[APILayerName, APILayerContract] = {
    "workflow": APILayerContract(
        name="workflow",
        purpose="Run the supported baseline workflow with minimal user code.",
        primary_entrypoints=(
            "scLucid.run_pipeline",
            "scLucid.qc.run_standard_qc",
            "scLucid.preprocess.run_preprocessing",
            "scLucid.analysis.run_standard_analysis",
        ),
        example_artifacts=("examples/01_workflow/basic_pipeline.py",),
        expected_outputs=(
            'adata.uns["sclucid"]["qc"]["review_summary"]',
            'adata.uns["sclucid"]["preprocess"]["review_summary"]',
            'adata.uns["sclucid"]["analysis"]["review_summary"]',
        ),
    ),
    "simple_api": APILayerContract(
        name="simple_api",
        purpose="Expose composable stage-level functions for inspection and overrides.",
        primary_entrypoints=(
            "scLucid.qc.calculate_qc_metric",
            "scLucid.qc.recommend_intelligent_qc",
            "scLucid.qc.mark_low_quality_cell",
            "scLucid.qc.filter_cells",
            "scLucid.preprocess.normalize_data",
            "scLucid.preprocess.find_hvgs",
            "scLucid.preprocess.scale_data",
            "scLucid.preprocess.batch_correction",
        ),
        example_artifacts=(
            "examples/02_simple_api/qc_step_by_step.py",
            "examples/02_simple_api/preprocess_step_by_step.py",
            "examples/02_simple_api/qc_preprocess_review.py",
        ),
        expected_outputs=(
            "inspectable intermediate AnnData state",
            "stage review summaries when workflow entrypoints are used",
            "reviewable tables and reports for manual decisions",
        ),
    ),
    "advanced": APILayerContract(
        name="advanced",
        purpose="Support full audit trails for real exploratory analysis projects.",
        primary_entrypoints=(
            "examples/03_advanced_notebooks/Step1A-QC_Audit.ipynb",
            "examples/03_advanced_notebooks/Step1B-Preprocessing_Audit.ipynb",
            "examples/03_advanced_notebooks/Step2-Annotation_and_Malignancy.ipynb",
            "scripts/run_pbmc_golden_path.py",
        ),
        example_artifacts=(
            "examples/03_advanced_notebooks/",
            "scripts/run_pbmc_golden_path.py",
        ),
        expected_outputs=(
            "final .h5ad",
            "stage review summaries",
            "figures and sidecar artifacts",
            "machine-readable manifest for golden paths",
        ),
    ),
}

MINIMAL_WORKFLOW_CONTRACT: dict[str, Any] = {
    "schema_version": SCHEMA_VERSION,
    "storage_root": SCLUCID_ROOT,
    "stage_order": STAGE_ORDER,
    "stages": STAGE_ORDER,
    "required_stage_namespace_keys": (
        UnsKeys.WORKFLOW_CONFIG,
        UnsKeys.STEPS_EXECUTED,
        UnsKeys.REVIEW_SUMMARY,
    ),
    "pipeline_records": (
        UnsKeys.PIPELINE_CONTEXT,
        UnsKeys.ANALYSIS_CONTEXT,
    ),
    "contract_validation": "run_pipeline records input and output validation under each stage namespace",
}


@dataclass(frozen=True)
class StageContract:
    """Input/output contract for one workflow stage."""

    name: StageName
    input_layers: tuple[str, ...] = ()
    input_obs: tuple[str, ...] = ()
    input_obsm: tuple[str, ...] = ()
    input_uns: tuple[tuple[str, ...], ...] = ()
    output_layers: tuple[str, ...] = ()
    output_obs: tuple[str, ...] = ()
    output_obsm: tuple[str, ...] = ()
    output_uns: tuple[tuple[str, ...], ...] = ()
    notes: tuple[str, ...] = ()


STAGE_CONTRACTS: dict[StageName, StageContract] = {
    "qc": StageContract(
        name="qc",
        output_obs=(ObsKeys.QC_N_GENES, ObsKeys.QC_TOTAL_COUNTS),
        output_uns=((SCLUCID_ROOT, Modules.QC, UnsKeys.REVIEW_SUMMARY),),
        notes=("QC input requires non-empty AnnData and raw counts in X or layers['counts'].",),
    ),
    "preprocess": StageContract(
        name="preprocess",
        input_layers=(LayerKeys.COUNTS,),
        input_uns=((SCLUCID_ROOT, Modules.QC),),
        output_layers=(LayerKeys.NORMALIZED,),
        output_obsm=(ObsmKeys.PCA,),
        output_uns=((SCLUCID_ROOT, Modules.PREPROCESS, UnsKeys.REVIEW_SUMMARY),),
    ),
    "analysis": StageContract(
        name="analysis",
        input_obsm=(ObsmKeys.PCA,),
        input_uns=((SCLUCID_ROOT, Modules.PREPROCESS),),
        output_uns=((SCLUCID_ROOT, Modules.ANALYSIS, UnsKeys.REVIEW_SUMMARY),),
    ),
}


@dataclass
class ContractValidationResult:
    """Result returned by contract validation."""

    valid: bool
    stage: str
    when: Literal["input", "output"]
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    checked: dict[str, list[str]] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serializable representation."""
        return {
            "valid": self.valid,
            "stage": self.stage,
            "when": self.when,
            "errors": self.errors,
            "warnings": self.warnings,
            "checked": self.checked,
        }


class ContractError(ValueError):
    """Raised when an AnnData object violates a workflow contract."""

    def __init__(self, result: ContractValidationResult):
        self.result = result
        super().__init__(format_contract_error(result))


def _has_uns_path(adata: AnnData, path: tuple[str, ...]) -> bool:
    current: Any = adata.uns
    for part in path:
        if not isinstance(current, dict) or part not in current:
            return False
        current = current[part]
    return True


def _path_label(path: Iterable[str]) -> str:
    parts = list(path)
    if not parts:
        return "adata.uns"
    if parts[0] == SCLUCID_ROOT:
        return 'adata.uns["sclucid"]' + "".join(f'["{part}"]' for part in parts[1:])
    return "adata.uns" + "".join(f'["{part}"]' for part in parts)


def _check_keys(
    available: Iterable[str],
    required: Iterable[str],
    *,
    container_label: str,
    errors: list[str],
) -> list[str]:
    available_set = set(available)
    required_list = list(required)
    for key in required_list:
        if key not in available_set:
            errors.append(f"Missing required {container_label} key: {key!r}")
    return required_list


def validate_stage_contract(
    adata: AnnData,
    stage: StageName,
    *,
    when: Literal["input", "output"],
    raise_on_error: bool = False,
) -> ContractValidationResult:
    """Validate an AnnData object against a stage input or output contract."""
    if stage not in STAGE_CONTRACTS:
        valid = ", ".join(STAGE_ORDER)
        raise ValueError(f"Unknown workflow stage {stage!r}. Valid stages are: {valid}.")

    contract = STAGE_CONTRACTS[stage]
    errors: list[str] = []
    warnings: list[str] = []
    checked: dict[str, list[str]] = {}

    if adata is None:
        errors.append("AnnData object is None.")
    elif adata.n_obs == 0 or adata.n_vars == 0:
        errors.append(f"AnnData must be non-empty, got shape={adata.shape}.")

    layer_keys = contract.input_layers if when == "input" else contract.output_layers
    obs_keys = contract.input_obs if when == "input" else contract.output_obs
    obsm_keys = contract.input_obsm if when == "input" else contract.output_obsm
    uns_paths = contract.input_uns if when == "input" else contract.output_uns

    checked["layers"] = _check_keys(
        adata.layers.keys(), layer_keys, container_label="layer", errors=errors
    )
    checked["obs"] = _check_keys(adata.obs.columns, obs_keys, container_label="obs", errors=errors)
    checked["obsm"] = _check_keys(
        adata.obsm.keys(), obsm_keys, container_label="obsm", errors=errors
    )

    checked_uns = []
    for path in uns_paths:
        label = _path_label(path)
        checked_uns.append(label)
        if not _has_uns_path(adata, path):
            errors.append(f"Missing required uns path: {label}")
    checked["uns"] = checked_uns

    if stage == "qc" and when == "input" and LayerKeys.COUNTS not in adata.layers:
        warnings.append("layers['counts'] is absent; QC will treat adata.X as the count source.")

    result = ContractValidationResult(
        valid=len(errors) == 0,
        stage=stage,
        when=when,
        errors=errors,
        warnings=warnings,
        checked=checked,
    )
    if raise_on_error and not result.valid:
        raise ContractError(result)
    return result


def stage_contract_to_dict(contract: StageContract) -> Dict[str, Any]:
    """Return a JSON-serializable representation of a stage contract."""
    data = asdict(contract)
    for key, value in list(data.items()):
        if isinstance(value, tuple):
            data[key] = list(value)
        if key.endswith("_uns"):
            data[key] = [list(path) for path in value]
    return data


def api_layer_contract_to_dict(contract: APILayerContract) -> Dict[str, Any]:
    """Return a JSON-serializable representation of one API layer contract."""
    data = asdict(contract)
    for key, value in list(data.items()):
        if isinstance(value, tuple):
            data[key] = list(value)
    return data


def get_stage_contract(stage: StageName) -> Dict[str, Any]:
    """Return the canonical contract for one workflow stage."""
    if stage not in STAGE_CONTRACTS:
        valid = ", ".join(STAGE_ORDER)
        raise ValueError(f"Unknown workflow stage {stage!r}. Valid stages are: {valid}.")
    return stage_contract_to_dict(STAGE_CONTRACTS[stage])


def get_api_layer_spec(layer: Optional[APILayerName] = None) -> Dict[str, Any]:
    """Return the frozen public API layer specification."""
    if layer is not None:
        if layer not in API_LAYER_CONTRACTS:
            valid = ", ".join(API_LAYER_ORDER)
            raise ValueError(f"Unknown API layer {layer!r}. Valid layers are: {valid}.")
        return api_layer_contract_to_dict(API_LAYER_CONTRACTS[layer])

    return {
        "schema_version": SCHEMA_VERSION,
        "layer_order": list(API_LAYER_ORDER),
        "layers": {
            name: api_layer_contract_to_dict(contract)
            for name, contract in API_LAYER_CONTRACTS.items()
        },
    }


def get_minimal_workflow_contract() -> Dict[str, Any]:
    """Return the frozen minimal workflow contract."""
    return {
        key: list(value) if isinstance(value, tuple) else value
        for key, value in MINIMAL_WORKFLOW_CONTRACT.items()
    }


def get_contract_spec() -> Dict[str, Any]:
    """Return the public scLucid data-contract specification."""
    return {
        "schema_version": SCHEMA_VERSION,
        "storage_root": SCLUCID_ROOT,
        "stage_order": list(STAGE_ORDER),
        "api_layers": get_api_layer_spec(),
        "minimal_workflow": get_minimal_workflow_contract(),
        "canonical_keys": {
            "layers": {
                "counts": LayerKeys.COUNTS,
                "normalized": LayerKeys.NORMALIZED,
                "scaled": LayerKeys.SCALED,
            },
            "obs": {
                "sample": ObsKeys.SAMPLE,
                "qc_n_genes": ObsKeys.QC_N_GENES,
                "qc_total_counts": ObsKeys.QC_TOTAL_COUNTS,
                "qc_pct_mt": ObsKeys.QC_PCT_MT,
                "low_quality": ObsKeys.LOW_QUALITY,
                "cluster": ObsKeys.CLUSTER,
                "cell_type": ObsKeys.CELL_TYPE,
            },
            "var": {
                "highly_variable": VarKeys.HIGHLY_VARIABLE,
            },
            "obsm": {
                "pca": ObsmKeys.PCA,
                "umap": ObsmKeys.UMAP,
                "spatial": ObsmKeys.SPATIAL,
            },
            "uns": {
                "namespace_metadata": UnsKeys.NAMESPACE_METADATA,
                "workflow_config": UnsKeys.WORKFLOW_CONFIG,
                "steps_executed": UnsKeys.STEPS_EXECUTED,
                "review_summary": UnsKeys.REVIEW_SUMMARY,
                "pipeline_context": UnsKeys.PIPELINE_CONTEXT,
                "analysis_context": UnsKeys.ANALYSIS_CONTEXT,
                "config_lineage": UnsKeys.CONFIG_LINEAGE,
                "contract": UnsKeys.CONTRACT,
                "artifacts": UnsKeys.ARTIFACTS,
                "errors": UnsKeys.ERRORS,
            },
            "modules": {
                "qc": Modules.QC,
                "preprocess": Modules.PREPROCESS,
                "analysis": Modules.ANALYSIS,
                "tumor": Modules.TUMOR,
                "tools": Modules.TOOLS,
            },
        },
        "review_summary": {
            "required_keys": list(REVIEW_SUMMARY_REQUIRED_KEYS),
            "recommended_keys": list(REVIEW_SUMMARY_RECOMMENDED_KEYS),
        },
        "stages": {
            stage: stage_contract_to_dict(contract)
            for stage, contract in STAGE_CONTRACTS.items()
        },
    }


def validate_all_stage_contracts(
    adata: AnnData,
    *,
    stages: Optional[Iterable[StageName]] = None,
    when: Literal["input", "output"] = "output",
    raise_on_error: bool = False,
) -> Dict[str, ContractValidationResult]:
    """Validate multiple workflow stage contracts against one AnnData object."""
    selected_stages = tuple(stages) if stages is not None else STAGE_ORDER
    results: Dict[str, ContractValidationResult] = {}
    failures: list[str] = []

    for stage in selected_stages:
        result = validate_stage_contract(adata, stage, when=when, raise_on_error=False)
        results[stage] = result
        if not result.valid:
            failures.append(format_contract_error(result))

    if failures and raise_on_error:
        first_invalid = next(result for result in results.values() if not result.valid)
        raise ContractError(first_invalid)

    return results


def format_contract_error(result: ContractValidationResult) -> str:
    """Format a contract validation result as an actionable error message."""
    errors = "; ".join(result.errors) if result.errors else "unknown contract violation"
    return f"[{result.stage}:{result.when}] AnnData contract failed: {errors}"


def ensure_sclucid_namespace(adata: AnnData) -> Dict[str, Any]:
    """Ensure and return ``adata.uns['sclucid']``."""
    existing = adata.uns.get(SCLUCID_ROOT)
    if existing is None:
        root: Dict[str, Any] = {}
        adata.uns[SCLUCID_ROOT] = root
    elif not isinstance(existing, dict):
        raise TypeError(
            f"adata.uns[{SCLUCID_ROOT!r}] must be a dictionary, "
            f"got {type(existing).__name__}."
        )
    else:
        root = existing

    now = datetime.now().isoformat()
    metadata = root.setdefault(UnsKeys.NAMESPACE_METADATA, {})
    metadata.setdefault("schema_version", SCHEMA_VERSION)
    metadata.setdefault("created_at", now)
    metadata["updated_at"] = now
    return root


def module_namespace(adata: AnnData, module: str, *, create: bool = True) -> Dict[str, Any]:
    """Return ``adata.uns['sclucid'][module]`` with optional creation."""
    root = ensure_sclucid_namespace(adata) if create else adata.uns.get(SCLUCID_ROOT, {})
    if module not in root:
        if not create:
            return {}
        root[module] = {}
    namespace = root[module]
    if not isinstance(namespace, dict):
        raise TypeError(
            f"adata.uns[{SCLUCID_ROOT!r}][{module!r}] must be a dictionary, "
            f"got {type(namespace).__name__}."
        )
    if create:
        now = datetime.now().isoformat()
        metadata = namespace.setdefault(UnsKeys.NAMESPACE_METADATA, {})
        metadata.setdefault("schema_version", SCHEMA_VERSION)
        metadata.setdefault("module", module)
        metadata.setdefault("created_at", now)
        metadata["updated_at"] = now
    return namespace


def build_config_lineage(
    *,
    global_config: Optional[Dict[str, Any]] = None,
    inherited: Optional[Dict[str, Any]] = None,
    stage_config: Optional[Dict[str, Any]] = None,
    effective_config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Build a standard record explaining configuration inheritance."""
    return {
        "schema_version": SCHEMA_VERSION,
        "global": global_config or {},
        "inherited": inherited or {},
        "stage": stage_config or {},
        "effective": effective_config or {},
        "precedence": ["explicit stage config", "inherited pipeline context", "global defaults"],
    }


def normalize_review_summary(
    summary: Dict[str, Any],
    *,
    module: str,
    workflow_name: str,
    adata: Optional[AnnData] = None,
    steps_executed: Optional[list[str]] = None,
    config: Optional[Dict[str, Any]] = None,
    warnings: Optional[list[str]] = None,
    config_lineage: Optional[Dict[str, Any]] = None,
    artifacts: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Add the standard review summary envelope while preserving existing fields."""
    normalized = dict(summary or {})
    normalized.setdefault("schema_version", SCHEMA_VERSION)
    normalized.setdefault("module", module)
    normalized.setdefault("workflow_name", workflow_name)
    normalized.setdefault("generated_at", datetime.now().isoformat())
    if steps_executed is not None:
        normalized.setdefault("steps_executed", list(steps_executed))
    else:
        normalized.setdefault("steps_executed", [])
    if adata is not None:
        normalized.setdefault(
            "data_shape", {"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)}
        )
    else:
        normalized.setdefault("data_shape", {})
    normalized.setdefault("warnings", warnings or normalized.get("warnings", []))
    if config is not None:
        normalized.setdefault("config", config)
    normalized.setdefault("config_lineage", config_lineage or {})
    normalized.setdefault("artifacts", artifacts or {})
    normalized.setdefault(
        "contract",
        {
            "schema_version": SCHEMA_VERSION,
            "required_keys": list(REVIEW_SUMMARY_REQUIRED_KEYS),
        },
    )
    # Backward-compatible read view for older tests/notebooks that accessed
    # ``adata.uns["sclucid"][module]["review_summary"]["data"]``. The canonical
    # contract remains the flat envelope above; ``data`` is a shallow mirror so
    # nested artifacts/contract updates stay shared without creating recursion.
    if "data" not in normalized or not isinstance(normalized.get("data"), dict):
        normalized["data"] = {key: value for key, value in normalized.items() if key != "data"}
    return normalized


def validate_review_summary_schema(
    summary: Dict[str, Any],
    *,
    module: Optional[str] = None,
    raise_on_error: bool = False,
) -> ContractValidationResult:
    """Validate the standard review summary envelope."""
    errors: list[str] = []
    warnings: list[str] = []

    if not isinstance(summary, dict):
        errors.append("Review summary must be a dictionary.")
        summary = {}

    for key in REVIEW_SUMMARY_REQUIRED_KEYS:
        if key not in summary:
            errors.append(f"Review summary missing required key: {key!r}")

    for key in REVIEW_SUMMARY_RECOMMENDED_KEYS:
        if key not in summary:
            warnings.append(f"Review summary missing recommended key: {key!r}")

    if module is not None and summary.get("module") != module:
        errors.append(
            f"Review summary module mismatch: expected {module!r}, "
            f"got {summary.get('module')!r}"
        )

    if "schema_version" in summary and str(summary["schema_version"]) != SCHEMA_VERSION:
        warnings.append(
            f"Review summary schema_version is {summary['schema_version']!r}; "
            f"current contract schema is {SCHEMA_VERSION!r}."
        )
    if "steps_executed" in summary and not isinstance(summary["steps_executed"], list):
        errors.append("Review summary 'steps_executed' must be a list.")
    if "data_shape" in summary:
        data_shape = summary["data_shape"]
        if not isinstance(data_shape, dict):
            errors.append("Review summary 'data_shape' must be a dictionary.")
        else:
            for key in ("n_cells", "n_genes"):
                if key in data_shape and not isinstance(data_shape[key], int):
                    errors.append(f"Review summary data_shape[{key!r}] must be an integer.")
    if "warnings" in summary and not isinstance(summary["warnings"], list):
        errors.append("Review summary 'warnings' must be a list.")
    if "config" in summary and not isinstance(summary["config"], dict):
        errors.append("Review summary 'config' must be a dictionary.")
    if "contract" in summary and not isinstance(summary["contract"], dict):
        errors.append("Review summary 'contract' must be a dictionary.")
    if "config_lineage" in summary and not isinstance(summary["config_lineage"], dict):
        errors.append("Review summary 'config_lineage' must be a dictionary.")
    if "artifacts" in summary and not isinstance(summary["artifacts"], dict):
        errors.append("Review summary 'artifacts' must be a dictionary.")

    result = ContractValidationResult(
        valid=len(errors) == 0,
        stage=module or str(summary.get("module", "unknown")),
        when="output",
        errors=errors,
        warnings=warnings,
        checked={"summary": list(REVIEW_SUMMARY_REQUIRED_KEYS)},
    )
    if raise_on_error and not result.valid:
        raise ContractError(result)
    return result


def record_contract_result(
    adata: AnnData,
    module: str,
    result: ContractValidationResult,
) -> None:
    """Store a contract validation result in the module namespace."""
    namespace = module_namespace(adata, module, create=True)
    contract_ns = namespace.setdefault(UnsKeys.CONTRACT, {})
    contract_ns[f"{result.when}_validation"] = result.to_dict()
    summary = namespace.get(UnsKeys.REVIEW_SUMMARY)
    if isinstance(summary, dict):
        summary.setdefault(UnsKeys.CONTRACT, {})
        summary[UnsKeys.CONTRACT][f"{result.when}_validation"] = result.to_dict()


def record_config_lineage(
    adata: AnnData,
    module: str,
    lineage: Dict[str, Any],
) -> None:
    """Store config lineage in the module namespace and review summary."""
    namespace = module_namespace(adata, module, create=True)
    namespace[UnsKeys.CONFIG_LINEAGE] = lineage
    summary = namespace.get(UnsKeys.REVIEW_SUMMARY)
    if isinstance(summary, dict):
        summary[UnsKeys.CONFIG_LINEAGE] = lineage
        if isinstance(summary.get("data"), dict):
            summary["data"][UnsKeys.CONFIG_LINEAGE] = lineage


def record_artifact(
    adata: AnnData,
    module: str,
    key: str,
    path: str,
    *,
    kind: str = "file",
    description: Optional[str] = None,
) -> None:
    """Record a saved artifact path under the module namespace."""
    namespace = module_namespace(adata, module, create=True)
    artifacts = namespace.setdefault(UnsKeys.ARTIFACTS, {})
    artifact = {
        "path": str(path),
        "kind": kind,
        "description": description,
        "recorded_at": datetime.now().isoformat(),
    }
    artifacts[key] = artifact
    summary = namespace.get(UnsKeys.REVIEW_SUMMARY)
    if isinstance(summary, dict):
        summary.setdefault(UnsKeys.ARTIFACTS, {})
        summary[UnsKeys.ARTIFACTS][key] = artifact


def record_error(
    adata: AnnData,
    module: str,
    error: BaseException,
    *,
    step_name: str = "unknown",
    context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Record a structured workflow error under the module namespace."""
    namespace = module_namespace(adata, module, create=True)
    errors = namespace.setdefault(UnsKeys.ERRORS, [])
    record = {
        "schema_version": SCHEMA_VERSION,
        "module": module,
        "step_name": step_name,
        "error_type": type(error).__name__,
        "message": str(error),
        "context": context or {},
        "recorded_at": datetime.now().isoformat(),
    }
    errors.append(record)
    return record


__all__ = [
    "API_LAYER_CONTRACTS",
    "API_LAYER_ORDER",
    "APILayerContract",
    "ContractError",
    "ContractValidationResult",
    "LayerKeys",
    "MINIMAL_WORKFLOW_CONTRACT",
    "Modules",
    "ObsmKeys",
    "ObsKeys",
    "REVIEW_SUMMARY_RECOMMENDED_KEYS",
    "REVIEW_SUMMARY_REQUIRED_KEYS",
    "SCHEMA_VERSION",
    "SCLUCID_ROOT",
    "STAGE_CONTRACTS",
    "STAGE_ORDER",
    "StageContract",
    "UnsKeys",
    "VarKeys",
    "api_layer_contract_to_dict",
    "build_config_lineage",
    "ensure_sclucid_namespace",
    "format_contract_error",
    "get_contract_spec",
    "get_api_layer_spec",
    "get_minimal_workflow_contract",
    "get_stage_contract",
    "module_namespace",
    "normalize_review_summary",
    "record_artifact",
    "record_contract_result",
    "record_config_lineage",
    "record_error",
    "stage_contract_to_dict",
    "validate_all_stage_contracts",
    "validate_review_summary_schema",
    "validate_stage_contract",
]
