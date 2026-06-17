"""Manual workflow review finalization helpers.

These helpers bridge low-level API usage and the workflow-layer contract. They
are intended for advanced notebooks or project scripts that execute individual
steps manually but still need the same review-summary handoff as standard
workflows.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Sequence, Union

from anndata import AnnData

from .contracts import normalize_review_summary, validate_review_summary_schema
from .sanitize import sanitize_for_hdf5
from .storage import export_review_summary, save_result


def _config_to_dict(config: Any) -> Optional[Dict[str, Any]]:
    if config is None:
        return None
    if isinstance(config, dict):
        return dict(config)
    if hasattr(config, "to_dict"):
        return config.to_dict()
    if hasattr(config, "model_dump"):
        return config.model_dump(mode="json")
    return {"repr": repr(config)}


def finalize_manual_review_summary(
    adata: AnnData,
    *,
    module: str,
    workflow_name: str,
    steps: Sequence[str],
    config: Any,
    summary: Optional[Dict[str, Any]] = None,
    save_dir: Optional[Union[str, Path]] = None,
    warnings: Optional[Sequence[str]] = None,
    title: Optional[str] = None,
) -> Dict[str, Any]:
    """Finalize a manual workflow path with the standard scLucid review contract.

    Parameters
    ----------
    adata
        AnnData object to annotate with ``adata.uns["sclucid"][module]``.
    module
        Module namespace, such as ``"qc"`` or ``"preprocess"``.
    workflow_name
        Stable workflow name for this manual path.
    steps
        Ordered step identifiers executed manually.
    config
        Config object or dict used for the manual path.
    summary
        Module-specific review payload to wrap in the standard review envelope.
    save_dir
        Optional directory for JSON/Markdown sidecars.
    warnings
        Optional review warnings to place in the standard envelope.
    title
        Optional Markdown title for sidecar export.

    Returns
    -------
    dict
        Normalized review summary stored under
        ``adata.uns["sclucid"][module]["review_summary"]``.
    """

    config_dict = _config_to_dict(config)
    steps_executed = list(steps)
    review_summary = normalize_review_summary(
        summary or {},
        module=module,
        workflow_name=workflow_name,
        adata=adata,
        steps_executed=steps_executed,
        config=config_dict,
        warnings=list(warnings or []),
    )
    validate_review_summary_schema(review_summary, module=module, raise_on_error=True)

    if config_dict is not None:
        save_result(adata, module, "workflow_config", config_dict)
    save_result(adata, module, "steps_executed", steps_executed)
    save_result(adata, module, "review_summary", review_summary)
    module_ns = adata.uns["sclucid"][module]
    review_summary = module_ns["review_summary"]

    if save_dir is not None:
        export_review_summary(
            review_summary,
            save_dir=save_dir,
            module=module,
            title=title,
            adata=adata,
        )
        review_summary = module_ns["review_summary"]

    return review_summary


__all__ = ["finalize_manual_review_summary"]
