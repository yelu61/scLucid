"""Final QC cell filtering and subsetting utilities."""

import ast
import logging
from typing import Dict, Optional

import pandas as pd
from anndata import AnnData

from .artifacts import record_filter_result
from .config import FilterConfig

log = logging.getLogger(__name__)

__all__ = [
    "filter_cells",
]


def _evaluate_filter_logic_expr(
    expr: str,
    namespace: Dict[str, pd.Series],
    *,
    index: pd.Index,
) -> pd.Series:
    """Evaluate a restricted boolean expression over QC criterion masks.

    Supported syntax is intentionally small: criterion names, parentheses,
    ``&`` / ``|`` and ``~`` plus the word forms ``and`` / ``or`` / ``not``.
    Function calls, comparisons, attributes, subscripts and constants are
    rejected so custom filtering cannot execute arbitrary Python.
    """

    def _eval(node: ast.AST) -> pd.Series:
        if isinstance(node, ast.Expression):
            return _eval(node.body)
        if isinstance(node, ast.Name):
            if node.id not in namespace:
                raise ValueError(f"Unknown criterion in custom logic: {node.id!r}")
            return namespace[node.id].reindex(index).fillna(False).astype(bool)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.Invert, ast.Not)):
            return ~_eval(node.operand)
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.BitAnd, ast.BitOr)):
            left = _eval(node.left)
            right = _eval(node.right)
            return left & right if isinstance(node.op, ast.BitAnd) else left | right
        if isinstance(node, ast.BoolOp) and isinstance(node.op, (ast.And, ast.Or)):
            values = [_eval(value) for value in node.values]
            result = values[0]
            for value in values[1:]:
                result = result & value if isinstance(node.op, ast.And) else result | value
            return result
        raise ValueError(
            "Custom filter logic only supports criterion names combined with &, |, ~, and/or/not."
        )

    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as exc:
        raise ValueError(f"Invalid custom filter logic syntax: {expr!r}") from exc
    return _eval(tree).reindex(index).fillna(False).astype(bool)





def filter_cells(
    adata: AnnData,
    config: Optional[FilterConfig] = None,
    copy: bool = False,
    **kwargs,
) -> Optional[AnnData]:
    """
    Enhanced cell filtering with flexible logical combinations and detailed reporting.

    This function filters cells based on previously calculated QC and doublet boolean flags.
    It supports multiple combination strategies and provides comprehensive statistics.

    Args:
        adata: AnnData object with QC metrics calculated
        config: FilterConfig object with filtering parameters
        copy: Whether to return a new filtered AnnData object

    Returns:
        Filtered AnnData object if copy=True, otherwise filters in place and returns None

    Example:
        # Basic usage
        adata_filtered = filter_cells(adata, copy=True)

        # Custom logic - remove cells with both MT and QC issues
        config = FilterConfig(
            combination_logic="custom",
            custom_logic_expr="outlier_mt & outlier_qc_metrics"
        )
        filter_cells(adata, config=config)

        # Threshold-based - remove cells with at least 2 issues
        config = FilterConfig(
            combination_logic="threshold",
            min_criteria_for_removal=2
        )
        filter_cells(adata, config=config)
    """
    # === 1. CONFIGURATION SETUP ===
    cfg = FilterConfig()
    if config is not None:
        cfg = config.model_copy(deep=True)
        if kwargs:
            cfg = cfg.model_copy(update=kwargs, deep=True)
    elif kwargs:
        cfg = cfg.model_copy(update=kwargs, deep=True)

    # --- Use cfg.criteria_to_filter instead of building a new list ---
    criteria = cfg.criteria_to_filter

    # === 2. FILTERING ===
    # Build criteria list if not provided
    if criteria is None:
        criteria = []

        # Map config attributes to column names
        criteria_mapping = {
            "filter_by_outlier_min_genes": "outlier_min_genes",
            "filter_by_outlier_max_genes": "outlier_max_genes",
            "filter_by_outlier_min_counts": "outlier_min_counts",
            "filter_by_outlier_max_counts": "outlier_max_counts",
            "filter_by_outlier_mt": "outlier_mt",
            "filter_by_outlier_hb": "outlier_hb",
            "filter_by_outlier_qc_metrics": "outlier_qc_metrics",
            "filter_by_scrublet_predicted": "scrublet_predicted",
            "filter_by_heuristic_predicted": "heuristic_predicted",
            "filter_by_predicted_doublet": "predicted_doublet",
        }

        for config_attr, col_name in criteria_mapping.items():
            if getattr(cfg, config_attr, False) and col_name in adata.obs.columns:
                criteria.append(col_name)

    # Filter out criteria that don't exist
    valid_criteria = [c for c in criteria if c in adata.obs.columns]
    missing_criteria = set(criteria) - set(valid_criteria)

    if missing_criteria:
        log.warning(f"Criteria not found in adata.obs and will be ignored: {missing_criteria}")

    if not valid_criteria:
        log.warning("No valid filtering criteria selected. Returning original object.")
        return adata.copy() if copy else None

    initial_cells = adata.n_obs
    log.info(f"Starting cell filtering with {initial_cells} cells")
    log.info(f"Using criteria: {', '.join(valid_criteria)}")
    log.info(f"Combination logic: {cfg.combination_logic}")

    # Apply metadata filters first if specified
    metadata_mask = pd.Series(True, index=adata.obs_names)
    if cfg.metadata_filters:
        for key, value in cfg.metadata_filters.items():
            if key in adata.obs.columns:
                if isinstance(value, list):
                    metadata_mask &= adata.obs[key].isin(value)
                else:
                    metadata_mask &= adata.obs[key] == value
                log.info(f"Applied metadata filter {key}={value}")
            else:
                log.warning(f"Metadata key '{key}' not found in adata.obs")

    # Calculate individual criteria masks
    criteria_masks = {}
    criteria_counts = {}

    for col in valid_criteria:
        col_mask = adata.obs[col].fillna(False).astype(bool)
        criteria_masks[col] = col_mask
        criteria_counts[col] = col_mask.sum()

    # Apply combination logic
    if cfg.combination_logic == "any":
        # Remove if ANY criterion is true (default behavior)
        combined_removal_mask = pd.Series(False, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            combined_removal_mask |= col_mask

    elif cfg.combination_logic == "all":
        # Remove only if ALL criteria are true
        combined_removal_mask = pd.Series(True, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            combined_removal_mask &= col_mask

    elif cfg.combination_logic == "custom":
        # Use custom expression
        if not cfg.custom_logic_expr:
            raise ValueError("custom_logic_expr must be provided when combination_logic='custom'")
        try:
            namespace = {col: criteria_masks[col] for col in valid_criteria}
            combined_removal_mask = _evaluate_filter_logic_expr(
                cfg.custom_logic_expr,
                namespace,
                index=adata.obs_names,
            )

        except Exception as e:
            log.error(f"Error evaluating custom logic expression: {e}")
            raise ValueError(f"Invalid custom logic expression: {cfg.custom_logic_expr}") from e

    elif cfg.combination_logic == "threshold":
        # Remove if at least min_criteria_for_removal criteria are true
        criteria_sum = pd.Series(0, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            criteria_sum += col_mask.astype(int)
        combined_removal_mask = criteria_sum >= cfg.min_criteria_for_removal

    else:
        raise ValueError(f"Unknown combination logic: {cfg.combination_logic}")

    # Apply metadata filter
    combined_removal_mask = combined_removal_mask & metadata_mask

    # Calculate final keep mask
    keep_mask = ~combined_removal_mask

    # Report statistics
    log.info("\n" + "=" * 40)
    log.info("CELL FILTERING STATISTICS")
    log.info("=" * 40)

    # Individual criteria counts
    total_cells = len(adata.obs_names)
    for col, count in criteria_counts.items():
        percentage = count / total_cells * 100
        log.info(f"{col}: {count} cells ({percentage:.2f}%)")

    # Final filtering results
    removed_count = combined_removal_mask.sum()
    kept_count = keep_mask.sum()

    log.info("\nFiltering results:")
    log.info(f"  Initial cells: {initial_cells}")
    log.info(f"  Cells removed: {removed_count} ({removed_count / initial_cells:.2%})")
    log.info(f"  Cells retained: {kept_count} ({kept_count / initial_cells:.2%})")

    # Analyze overlap between criteria
    if len(valid_criteria) > 1:
        log.info("\nOverlap analysis:")

        # Count cells with multiple issues
        criteria_sum = pd.Series(0, index=adata.obs_names)
        for col_mask in criteria_masks.values():
            criteria_sum += col_mask.astype(int)

        for i in range(1, len(valid_criteria) + 1):
            count = (criteria_sum == i).sum()
            if count > 0:
                percentage = count / total_cells * 100
                log.info(f"  Cells with exactly {i} issues: {count} ({percentage:.2f}%)")

    log.info("=" * 40)

    # === Store filtering results in the unified namespace ===
    if "sclucid" not in adata.uns:
        adata.uns["sclucid"] = {}
    if "qc" not in adata.uns["sclucid"]:
        adata.uns["sclucid"]["qc"] = {}

    removal_reasons = {}
    for col, mask in criteria_masks.items():
        removal_reasons[col] = int((mask & combined_removal_mask).sum())

    # Store stats in a dictionary
    stats = {
        "schema_version": "qc_filter_result_v1",
        "initial_cells": initial_cells,
        "final_cells": kept_count,
        "removed_cells": removed_count,
        "removed_fraction": removed_count / initial_cells if initial_cells > 0 else 0,
        "criteria_used": valid_criteria,
        "criteria_requested": list(criteria),
        "criteria_missing": sorted(missing_criteria),
        "combination_logic": cfg.combination_logic,
        "min_criteria_for_removal": cfg.min_criteria_for_removal,
        "criteria_counts": criteria_counts,
        "removal_reason_counts": removal_reasons,
        "config": cfg.to_dict(),
    }

    # Perform filtering
    if copy:
        adata_filtered = adata[keep_mask, :].copy()
        record_filter_result(adata_filtered, stats=stats)
        return adata_filtered
    else:
        record_filter_result(adata, stats=stats)
        adata._inplace_subset_obs(keep_mask)
        return None
