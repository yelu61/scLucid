"""
Unified storage utilities for scLucid.

Provides standardized access to adata.uns['sclucid'] for storing analysis results,
configurations, and metadata across all modules (qc, preprocess, analysis, tools).

Usage:
    >>> from scLucid.utils import get_storage, save_result, load_result
    >>>
    >>> # Save results in standardized location
    >>> save_result(adata, module='qc', key='metrics', data=metrics_dict)
    >>>
    >>> # Retrieve results
    >>> metrics = load_result(adata, module='qc', key='metrics')
    >>>
    >>> # Get entire module storage
    >>> qc_storage = get_storage(adata, module='qc')
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union

from anndata import AnnData

from .contracts import (
    SCHEMA_VERSION,
    SCLUCID_ROOT,
    Modules,
    UnsKeys,
    _restore_review_sequences,
    module_namespace,
    record_artifact,
)
from .sanitize import sanitize_for_hdf5

log = logging.getLogger(__name__)

# Storage hierarchy: adata.uns['sclucid'][module][key]
STORAGE_ROOT = SCLUCID_ROOT

# Valid modules for storage organization
VALID_MODULES = {
    Modules.QC,
    Modules.PREPROCESS,
    Modules.ANALYSIS,
    "clustering",
    "annotation",
    "de",
    "enrichment",
    "proportion",
    "scenic",
    Modules.TOOLS,
    "checkpoint",
    Modules.TUMOR,
}


def get_storage(adata: AnnData, module: str, create: bool = True) -> Dict[str, Any]:
    """
    Get storage dictionary for a specific module.

    Args:
        adata: AnnData object
        module: Module name (e.g., 'qc', 'preprocess', 'analysis')
        create: If True, create storage structure if it doesn't exist

    Returns:
        Dictionary for the module's storage space

    Example:
        >>> qc_storage = get_storage(adata, 'qc')
        >>> qc_storage['my_key'] = my_data
    """
    if module not in VALID_MODULES:
        log.warning(f"Unknown module '{module}'. Valid modules: {VALID_MODULES}")

    return module_namespace(adata, module, create=create)


def save_result(
    adata: AnnData,
    module: str,
    key: str,
    data: Any,
    config: Optional[Dict[str, Any]] = None,
    overwrite: bool = True,
) -> None:
    """
    Save analysis result to standardized storage location.

    Args:
        adata: AnnData object
        module: Module name (e.g., 'qc', 'preprocess', 'analysis')
        key: Unique key for this result within the module
        data: Data to store (must be HDF5-serializable)
        config: Optional configuration dict to store alongside result
        overwrite: If False, raise error if key already exists

    Raises:
        KeyError: If key exists and overwrite=False

    Example:
        >>> save_result(adata, 'qc', 'metrics', {'n_cells': 1000})
        >>> save_result(adata, 'qc', 'metrics_config', config_dict)
    """
    storage = get_storage(adata, module, create=True)

    if key in storage and not overwrite:
        raise KeyError(
            f"Key '{key}' already exists in {module} storage. " f"Use overwrite=True to replace."
        )

    canonical_direct_keys = {
        UnsKeys.WORKFLOW_CONFIG,
        UnsKeys.STEPS_EXECUTED,
        UnsKeys.REVIEW_SUMMARY,
        UnsKeys.CONFIG_LINEAGE,
        UnsKeys.CONTRACT,
        UnsKeys.ARTIFACTS,
        UnsKeys.ERRORS,
    }
    if key in canonical_direct_keys:
        sanitized = sanitize_for_hdf5(data)
        if key == UnsKeys.REVIEW_SUMMARY:
            sanitized = _restore_review_sequences(sanitized)
        storage[key] = sanitized
        log.debug(f"Saved canonical result '{key}' to {module} storage")
        return

    # Store result with metadata
    storage[key] = {
        "schema_version": SCHEMA_VERSION,
        "module": module,
        "key": key,
        "data": sanitize_for_hdf5(data),
        "timestamp": datetime.now().isoformat(),
    }

    if config is not None:
        storage[f"{key}_config"] = {
            "schema_version": SCHEMA_VERSION,
            "module": module,
            "key": f"{key}_config",
            "data": sanitize_for_hdf5(config),
            "timestamp": datetime.now().isoformat(),
        }

    log.debug(f"Saved result '{key}' to {module} storage")


def load_result(adata: AnnData, module: str, key: str, default: Any = None) -> Any:
    """
    Load analysis result from standardized storage.

    Args:
        adata: AnnData object
        module: Module name
        key: Key for the stored result
        default: Default value if key not found

    Returns:
        Stored data, or default if not found

    Example:
        >>> metrics = load_result(adata, 'qc', 'metrics', default={})
    """
    storage = get_storage(adata, module, create=False)

    if key not in storage:
        return default

    result = storage[key]

    # Handle new format with metadata wrapper
    if isinstance(result, dict) and "data" in result:
        return result["data"]

    # Handle legacy format (direct storage)
    return result


def load_config(adata: AnnData, module: str, key: str) -> Optional[Dict[str, Any]]:
    """
    Load configuration stored alongside a result.

    Args:
        adata: AnnData object
        module: Module name
        key: Base key for the result

    Returns:
        Configuration dict, or None if not found
    """
    config_key = f"{key}_config"
    result = load_result(adata, module, config_key)

    # Handle legacy format where config was stored directly
    if result is None:
        storage = get_storage(adata, module, create=False)
        if config_key in storage:
            return storage[config_key]

    return result


def has_result(adata: AnnData, module: str, key: str) -> bool:
    """Check if a result exists in storage."""
    storage = get_storage(adata, module, create=False)
    return key in storage


def list_results(adata: AnnData, module: Optional[str] = None) -> Dict[str, list]:
    """
    List all stored results.

    Args:
        adata: AnnData object
        module: If specified, only list results for this module

    Returns:
        Dict mapping module names to lists of result keys
    """
    if STORAGE_ROOT not in adata.uns:
        return {}

    root = adata.uns[STORAGE_ROOT]

    if module:
        if module in root:
            return {module: list(root[module].keys())}
        return {}

    return {
        mod: list(keys.keys())
        for mod, keys in root.items()
        if mod != UnsKeys.NAMESPACE_METADATA and isinstance(keys, dict)
    }


def clear_storage(
    adata: AnnData, module: Optional[str] = None, keys: Optional[list] = None, dry_run: bool = False
) -> Dict[str, Any]:
    """
    Clear stored results.

    Args:
        adata: AnnData object
        module: If specified, only clear this module
        keys: If specified, only clear these keys
        dry_run: If True, only report what would be deleted

    Returns:
        Dict with information about cleared items
    """
    if STORAGE_ROOT not in adata.uns:
        return {"cleared": [], "modules_cleared": []}

    root = adata.uns[STORAGE_ROOT]
    cleared = []
    modules_cleared = []

    if module:
        if module in root:
            if keys:
                for key in keys:
                    if key in root[module]:
                        if not dry_run:
                            del root[module][key]
                        cleared.append(f"{module}.{key}")
            else:
                if not dry_run:
                    del root[module]
                modules_cleared.append(module)
    else:
        # Clear everything
        all_modules = list(root.keys())
        if not dry_run:
            adata.uns[STORAGE_ROOT] = {}
        modules_cleared = all_modules

    return {"cleared": cleared, "modules_cleared": modules_cleared}


def migrate_legacy_storage(adata: AnnData, dry_run: bool = False) -> Dict[str, list]:
    """
    Migrate legacy storage formats to standardized format.

    Handles:
    - Top-level 'qc' key -> 'sclucid.qc'
    - Direct result storage -> wrapped with metadata

    Args:
        adata: AnnData object
        dry_run: If True, only report what would be migrated

    Returns:
        Dict with migration summary
    """
    migrated = []

    # Migrate top-level 'qc' key
    if "qc" in adata.uns and isinstance(adata.uns["qc"], dict):
        if not dry_run:
            storage = get_storage(adata, "qc", create=True)
            for key, value in adata.uns["qc"].items():
                if key not in storage:
                    storage[key] = value
            del adata.uns["qc"]
        migrated.append("qc (top-level -> sclucid.qc)")

    # Migrate top-level analysis keys
    for old_key in ["clustering", "annotation"]:
        if old_key in adata.uns:
            if not dry_run:
                storage = get_storage(adata, "analysis", create=True)
                if old_key not in storage:
                    storage[old_key] = adata.uns[old_key]
                del adata.uns[old_key]
            migrated.append(f"{old_key} (top-level -> sclucid.analysis)")

    return {"migrated": migrated}


# Convenience functions for common patterns


def save_workflow_result(
    adata: AnnData,
    module: str,
    workflow_name: str,
    steps: list,
    config: Dict[str, Any],
    *,
    artifacts: Optional[Dict[str, Union[str, Path]]] = None,
) -> None:
    """
    Save workflow completion metadata.

    Standardized format used by all workflow modules.
    """
    storage = get_storage(adata, module, create=True)
    storage[UnsKeys.WORKFLOW_CONFIG] = config
    storage[UnsKeys.STEPS_EXECUTED] = list(steps)

    save_result(
        adata,
        module,
        f"{workflow_name}_workflow",
        {
            "name": workflow_name,
            "steps_executed": steps,
            "completed_at": datetime.now().isoformat(),
        },
        config=config,
    )
    if artifacts:
        for key, path in artifacts.items():
            record_artifact(adata, module, key, str(path), description=f"{workflow_name} artifact")


def load_workflow_result(
    adata: AnnData, module: str, workflow_name: str
) -> Optional[Dict[str, Any]]:
    """Load workflow completion metadata."""
    return load_result(adata, module, f"{workflow_name}_workflow")


def export_review_summary(
    review_summary: Dict[str, Any],
    save_dir: Union[str, Path],
    module: str,
    title: Optional[str] = None,
    adata: Optional[AnnData] = None,
) -> Dict[str, Path]:
    """
    Export a review summary as JSON and Markdown sidecars.

    This is the unified entry point for writing review-facing artifacts
    from any scLucid module (qc, preprocess, analysis).

    Args:
        review_summary: Structured review summary dict
        save_dir: Directory to write files
        module: Module name (used in filenames, e.g. 'qc', 'preprocess', 'analysis')
        title: Optional title override for Markdown header

    Returns:
        Dict mapping 'json' and 'md' to written Path objects

    Example:
        >>> export_review_summary(summary, "./results", "qc")
        {'json': Path('.../qc_review_summary.json'), 'md': Path('.../qc_review_summary.md')}
    """
    import json as _json

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    # JSON sidecar
    json_path = save_dir / f"{module}_review_summary.json"
    json_path.write_text(
        _json.dumps(review_summary, indent=2, default=str),
        encoding="utf-8",
    )

    # Markdown sidecar
    md_title = title or f"{module.capitalize()} Review Summary"
    md_lines = [f"# {md_title}", ""]

    def _render_dict(d: Dict[str, Any], indent: int = 0) -> list[str]:
        lines = []
        prefix = "  " * indent
        for key, val in d.items():
            if isinstance(val, dict):
                lines.append(f"{prefix}- **{key}**:")
                lines.extend(_render_dict(val, indent + 1))
            elif isinstance(val, list):
                lines.append(f"{prefix}- **{key}**:")
                for item in val:
                    if isinstance(item, dict):
                        lines.append(f"{prefix}  -")
                        lines.extend(_render_dict(item, indent + 2))
                    else:
                        lines.append(f"{prefix}  - {item}")
            else:
                lines.append(f"{prefix}- **{key}**: {val}")
        return lines

    md_lines.extend(_render_dict(review_summary))
    md_lines.append("")

    md_path = save_dir / f"{module}_review_summary.md"
    md_path.write_text("\n".join(md_lines), encoding="utf-8")

    log.info(f"Review summary exported to {json_path} and {md_path}")
    artifacts = {"json": json_path, "md": md_path}
    if adata is not None:
        for key, path in artifacts.items():
            record_artifact(
                adata,
                module,
                f"{module}_review_summary_{key}",
                str(path),
                kind=key,
                description=f"{module} review summary {key.upper()} sidecar",
            )
    return artifacts


def write_h5ad_safe(
    adata: AnnData,
    path: Union[str, Path],
    *,
    compression: Optional[str] = "gzip",
    sanitize_uns: bool = True,
    atomic: bool = True,
    lightweight: bool = False,
    drop_layers: Iterable[str] = ("regressed", "scaled"),
    x_layer: Optional[str] = None,
) -> Path:
    """Write an AnnData object with scLucid-safe provenance handling.

    Parameters
    ----------
    adata
        Object to save.
    path
        Output ``.h5ad`` path.
    compression
        Compression passed to ``AnnData.write_h5ad``.
    sanitize_uns
        Sanitize ``adata.uns['sclucid']`` before writing so nested provenance
        records are HDF5-compatible.
    atomic
        Write to a temporary file in the same directory and rename on success.
    lightweight
        Save a copy with selected dense working layers removed.
    drop_layers
        Layers to remove when ``lightweight=True``.
    x_layer
        Optional layer to copy into ``.X`` before writing the lightweight copy.

    Returns
    -------
    pathlib.Path
        Final written path.
    """
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata_to_write = adata.copy() if lightweight else adata

    if sanitize_uns and "sclucid" in adata_to_write.uns:
        adata_to_write.uns["sclucid"] = sanitize_for_hdf5(adata_to_write.uns["sclucid"])

    if lightweight:
        for layer in drop_layers:
            if layer in adata_to_write.layers:
                del adata_to_write.layers[layer]
        if x_layer is not None:
            if x_layer not in adata_to_write.layers:
                raise KeyError(f"x_layer '{x_layer}' not found in adata.layers")
            adata_to_write.X = adata_to_write.layers[x_layer].copy()

    write_path = output_path
    temp_path: Optional[Path] = None
    if atomic:
        temp_path = output_path.with_name(f".{output_path.name}.tmp")
        write_path = temp_path
        if temp_path.exists():
            temp_path.unlink()

    try:
        adata_to_write.write_h5ad(write_path, compression=compression)
        if temp_path is not None:
            temp_path.replace(output_path)
    except Exception:
        if temp_path is not None and temp_path.exists():
            temp_path.unlink()
        raise

    return output_path


__all__ = [
    # Core storage functions
    "get_storage",
    "save_result",
    "load_result",
    "load_config",
    "has_result",
    "list_results",
    "clear_storage",
    "migrate_legacy_storage",
    # Convenience functions
    "save_workflow_result",
    "load_workflow_result",
    "export_review_summary",
    "write_h5ad_safe",
    # Constants
    "STORAGE_ROOT",
    "VALID_MODULES",
]
