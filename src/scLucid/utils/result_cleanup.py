"""Utilities for cleaning up sclucid results stored in AnnData objects.

Long-running workflows can leave large nested structures in
``adata.uns['sclucid']``.  This module provides helpers to prune them
and reclaim memory / disk space.
"""

import logging
from typing import Iterable, Optional

import numpy as np
from anndata import AnnData

log = logging.getLogger(__name__)

# Modules that scLucid stores under adata.uns["sclucid"]
KNOWN_MODULES = {"qc", "preprocess", "analysis", "proportion", "recommendation"}
COMPACTABLE_MODULES = ("qc", "preprocess", "analysis")
COMPACT_KEEP_KEYS = {
    "_metadata",
    "workflow_config",
    "steps_executed",
    "review_summary",
    "config_lineage",
    "contract",
    "artifacts",
    "errors",
    "context",
    "warnings",
    "original_config",
    "applied_config",
    "user_overrides",
    "pipeline_context",
    "analysis_context",
    "tumor_aware_flags",
    "tumor_aware_notes",
}


def _estimate_size(obj) -> int:
    """Rough byte-size estimate for nested dicts/DataFrames/arrays."""
    try:
        import pandas as pd

        if isinstance(obj, pd.DataFrame):
            return int(obj.memory_usage(deep=True).sum())
    except Exception:
        pass

    if isinstance(obj, np.ndarray):
        return obj.nbytes

    if isinstance(obj, dict):
        return sum(_estimate_size(v) for v in obj.values())

    if isinstance(obj, (list, tuple)):
        return sum(_estimate_size(v) for v in obj)

    # Fallback — very rough
    return 64


def clear_sclucid_results(
    adata: AnnData,
    module: Optional[str] = None,
    dry_run: bool = False,
) -> dict:
    """Remove sclucid result data from ``adata.uns``.

    Parameters
    ----------
    adata : AnnData
        The object whose ``uns["sclucid"]`` entries should be cleared.
    module : {None, 'qc', 'preprocess', 'analysis', 'proportion', 'recommendation'}, optional
        If *None*, all sclucid modules are removed.
        If a module name is given, only that sub-key is removed.
    dry_run : bool, default=False
        If *True*, report what would be removed without actually deleting.

    Returns:
    -------
    dict
        Summary with ``removed_modules`` (list) and ``estimated_bytes`` (int).

    Examples:
    --------
    >>> clear_sclucid_results(adata, module='qc')
    {'removed_modules': ['qc'], 'estimated_bytes': 123456}

    >>> clear_sclucid_results(adata)  # clear everything
    """
    if "sclucid" not in adata.uns:
        log.info("No sclucid results found in adata.uns.")
        return {"removed_modules": [], "estimated_bytes": 0}

    sclucid = adata.uns["sclucid"]
    if not isinstance(sclucid, dict):
        log.warning("adata.uns['sclucid'] is not a dict; skipping.")
        return {"removed_modules": [], "estimated_bytes": 0}

    if module is not None and module not in sclucid:
        log.info(f"Module '{module}' not found in adata.uns['sclucid'].")
        return {"removed_modules": [], "estimated_bytes": 0}

    to_remove = [module] if module is not None else list(sclucid.keys())
    estimated_bytes = 0
    removed = []

    for key in to_remove:
        if key not in sclucid:
            continue
        size = _estimate_size(sclucid[key])
        estimated_bytes += size
        removed.append(key)

        if not dry_run:
            del sclucid[key]
            log.info(f"Removed adata.uns['sclucid']['{key}'] (~{size:,} bytes)")

    # If the sclucid dict is now empty, remove the top-level key entirely
    if not dry_run and not sclucid:
        del adata.uns["sclucid"]
        log.info("Removed empty adata.uns['sclucid']")

    action = "Would remove" if dry_run else "Removed"
    log.info(
        f"{action} {len(removed)} module(s) from adata.uns['sclucid']: {removed} "
        f"(~{estimated_bytes:,} bytes)"
    )

    return {"removed_modules": removed, "estimated_bytes": estimated_bytes}


def list_sclucid_modules(adata: AnnData) -> dict:
    """Inspect what sclucid data is stored in *adata*.

    Returns a dict mapping module names to estimated byte sizes.
    """
    if "sclucid" not in adata.uns:
        return {}

    sclucid = adata.uns["sclucid"]
    if not isinstance(sclucid, dict):
        return {}

    return {k: _estimate_size(v) for k, v in sclucid.items()}


def compact_sclucid_uns(
    adata: AnnData,
    modules: Optional[Iterable[str]] = None,
    *,
    keep_keys: Optional[Iterable[str]] = None,
    dry_run: bool = False,
) -> dict:
    """Prune heavy intermediate ``adata.uns['sclucid']`` entries stage by stage.

    The compact form is intended for handoff to the next workflow stage or for
    smaller ``.h5ad`` files after sidecar reports have been exported. It keeps
    the audit contract and compact review material, while removing large result
    tables, model outputs, diagnostic payloads, and duplicated summaries.

    Parameters
    ----------
    adata
        AnnData object to compact.
    modules
        Module namespaces to compact. Defaults to ``("qc", "preprocess", "analysis")``.
    keep_keys
        Additional per-module top-level keys to keep.
    dry_run
        If ``True``, report what would be removed without modifying ``adata``.

    Returns
    -------
    dict
        Summary with removed paths, preserved paths, and rough byte estimates.
    """
    if "sclucid" not in adata.uns or not isinstance(adata.uns["sclucid"], dict):
        return {
            "modules_compacted": [],
            "removed": [],
            "preserved": [],
            "estimated_bytes_removed": 0,
            "dry_run": dry_run,
        }

    sclucid = adata.uns["sclucid"]
    modules_to_compact = tuple(modules or COMPACTABLE_MODULES)
    keys_to_keep = set(COMPACT_KEEP_KEYS)
    if keep_keys is not None:
        keys_to_keep.update(str(key) for key in keep_keys)

    removed: list[str] = []
    preserved: list[str] = []
    modules_compacted: list[str] = []
    estimated_bytes = 0

    for module in modules_to_compact:
        namespace = sclucid.get(module)
        if not isinstance(namespace, dict):
            continue

        modules_compacted.append(module)
        compact_namespace = {}
        for key, value in namespace.items():
            path = f"{module}.{key}"
            if key in keys_to_keep:
                compact_namespace[key] = value
                preserved.append(path)
            else:
                estimated_bytes += _estimate_size(value)
                removed.append(path)

        if not dry_run:
            sclucid[module] = compact_namespace

    log.info(
        "%s sclucid uns compact: modules=%s, removed=%d key(s), preserved=%d key(s), ~%d bytes",
        "Dry-run" if dry_run else "Applied",
        list(modules_to_compact),
        len(removed),
        len(preserved),
        estimated_bytes,
    )

    return {
        "modules_compacted": modules_compacted,
        "removed": removed,
        "preserved": preserved,
        "estimated_bytes_removed": estimated_bytes,
        "dry_run": dry_run,
    }


__all__ = [
    "clear_sclucid_results",
    "compact_sclucid_uns",
    "list_sclucid_modules",
    "KNOWN_MODULES",
    "COMPACTABLE_MODULES",
    "COMPACT_KEEP_KEYS",
]
