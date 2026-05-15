"""Shared utilities for annotation modules."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Dict, List, Optional

import pandas as pd
from anndata import AnnData

from ...utils import Manager, get_marker_manager
from ..config import AnnotationConfig

_ANNOTATION_NOISE_EXACT_GENES = {
    "MALAT1": "housekeeping",
    "NEAT1": "housekeeping",
    "XIST": "housekeeping",
}

_ANNOTATION_NOISE_PREFIXES = {
    "ribosomal": ("RPL", "RPS", "MRPL", "MRPS"),
    "mitochondrial": ("MT-",),
    "stress": (
        "HSP",
        "HSPA",
        "HSPB",
        "HSPC",
        "HSPD",
        "HSPE",
        "HSPH",
        "DNAJ",
        "FOS",
        "JUN",
        "ATF3",
        "IER",
        "DDIT3",
        "PPP1R15A",
    ),
}


def _get_default_compartment_map() -> Dict[str, str]:
    """Return a default compartment mapper for common cell type labels."""
    return {
        # Tumor / malignant
        "Malignant": "tumor",
        "Tumor": "tumor",
        "Cancer": "tumor",
        "Carcinoma": "tumor",
        "Epithelial": "tumor",
        "Pan-Cancer": "tumor",
        "LUAD": "tumor",
        "LUSC": "tumor",
        "SCLC": "tumor",
        # Immune
        "T cells": "immune",
        "B cells": "immune",
        "NK cells": "immune",
        "Macrophage": "immune",
        "Monocyte": "immune",
        "DC": "immune",
        "Neutrophil": "immune",
        "Mast": "immune",
        "Plasma": "immune",
        # Stromal
        "Fibroblast": "stromal",
        "Endothelial": "stromal",
        "Pericyte": "stromal",
        "Smooth muscle": "stromal",
        "Stromal": "stromal",
        # Uncertain
        "Unknown": "uncertain",
        "Unassigned": "uncertain",
    }


def _map_compartments(
    labels: pd.Series,
    compartment_map: Optional[Dict[str, str]] = None,
) -> pd.Series:
    """Map cell type labels to broad compartments (tumor/immune/stromal/uncertain/mixed)."""
    if compartment_map is None:
        compartment_map = _get_default_compartment_map()

    def _map_label(label: str) -> str:
        if pd.isna(label):
            return "uncertain"
        label_str = str(label)
        # Exact match first
        if label_str in compartment_map:
            return compartment_map[label_str]
        # Partial match fallback
        for key, compartment in compartment_map.items():
            if key.lower() in label_str.lower():
                return compartment
        return "uncertain"

    mapped = labels.apply(_map_label)
    return mapped


def _read_table_file(path: str) -> pd.DataFrame:
    """
    Read a tabular mapping file with robust handling:
    - Supports .xlsx/.xls (preferred) and .csv
    - For CSV, tries common encodings, then chardet probing as fallback
    Returns a DataFrame.
    """
    import os

    ext = os.path.splitext(path)[1].lower()
    if ext in [".xlsx", ".xls"]:
        # Excel is robust against encoding issues
        return pd.read_excel(path)
    elif ext == ".csv":
        # Try multiple encodings
        last_err = None
        encodings = ["utf-8", "utf-8-sig", "gbk", "gb2312", "big5", "cp1252", "latin1"]
        for enc in encodings:
            try:
                return pd.read_csv(path, encoding=enc)
            except Exception as e:
                last_err = e
        # chardet fallback
        try:
            import chardet

            with open(path, "rb") as f:
                raw = f.read(200000)
            guess = chardet.detect(raw).get("encoding") or "utf-8"
            return pd.read_csv(path, encoding=guess, errors="replace")
        except Exception:
            pass
        # final fallback
        try:
            return pd.read_csv(path, encoding="latin1", errors="replace")
        except Exception:
            raise (
                last_err
                if last_err
                else RuntimeError(f"Failed to read CSV with multiple encodings: {path}")
            )
    else:
        raise ValueError("Unsupported table format. Use .xlsx, .xls, or .csv")


def _read_json_file(path: str) -> dict:
    """
    Read JSON with robust encoding handling.
    """
    import json

    last_err = None
    for enc in ["utf-8", "utf-8-sig", "gbk", "cp1252", "latin1"]:
        try:
            with open(path, encoding=enc) as f:
                return json.load(f)
        except Exception as e:
            last_err = e
    try:
        import chardet

        with open(path, "rb") as fb:
            raw = fb.read(200000)
        guess = chardet.detect(raw).get("encoding") or "utf-8"
        with open(path, encoding=guess, errors="replace") as f:
            return json.load(f)
    except Exception:
        pass
    raise (last_err if last_err else RuntimeError("Failed to read JSON with multiple encodings."))


def _classify_annotation_marker(gene: Any) -> Optional[str]:
    """Classify common non-informative marker genes seen during manual annotation."""
    if pd.isna(gene):
        return None

    gene_upper = str(gene).upper()
    if gene_upper in _ANNOTATION_NOISE_EXACT_GENES:
        return _ANNOTATION_NOISE_EXACT_GENES[gene_upper]

    for category, prefixes in _ANNOTATION_NOISE_PREFIXES.items():
        if any(gene_upper.startswith(prefix) for prefix in prefixes):
            return category
    return None


def _resolve_score_columns(
    adata: AnnData,
    score_cols: Optional[Sequence[str]] = None,
) -> List[str]:
    """Resolve module score columns from obs."""
    if score_cols is not None:
        return [col for col in score_cols if col in adata.obs.columns]
    return [
        col
        for col in adata.obs.columns
        if col.endswith("_score") and pd.api.types.is_numeric_dtype(adata.obs[col])
    ]


def _resolve_annotation_manager(
    *,
    species: str,
    tissue: Optional[str],
    states: Optional[List[str]] = None,
    marker_config: Optional[str] = None,
    view: Optional[str] = None,
):
    """Resolve either a custom marker file or the built-in combined marker manager."""
    if marker_config:
        return Manager(marker_config, case_sensitive=True)
    return get_marker_manager(species=species, tissue=tissue, states=states, view=view)


def _label_matches_target_lineage(labels: pd.Series, target_lineage: Optional[str]) -> pd.Series:
    """Return a boolean mask for cells matching the requested target lineage."""
    if not target_lineage:
        return pd.Series(True, index=labels.index)

    target = str(target_lineage).strip().lower()
    values = labels.astype(str).str.lower()
    return values.eq(target) | values.str.contains(target, regex=False)


def _build_modular_annotation_label(
    lineage: pd.Series,
    subtype: pd.Series,
    state: pd.Series,
) -> pd.Series:
    """Construct a compact modular display label from lineage/subtype/state columns."""
    result = []
    for lin, sub, st in zip(lineage.astype(str), subtype.astype(str), state.astype(str)):
        label = lin
        if sub not in {"Unknown", "Not_applicable", "nan", ""} and sub != lin:
            label = sub
        if st not in {"Unknown", "Not_applicable", "nan", ""}:
            label = f"{label} | {st}"
        result.append(label)
    return pd.Series(result, index=lineage.index, dtype="object")


def _rename_score_columns_for_manager(
    adata: AnnData,
    manager: Manager,
    suffix: str,
) -> None:
    """Rename manager-derived *_score columns to a scoped suffix."""
    rename_map = {}
    for cell_type in manager.CELLS:
        source = f"{cell_type}_score"
        if source in adata.obs.columns:
            rename_map[source] = f"{cell_type}{suffix}"
    if rename_map:
        adata.obs.rename(columns=rename_map, inplace=True)


def _collect_state_signatures(
    config: AnnotationConfig,
) -> tuple[Dict[str, List[str]], Dict[str, Dict[str, object]]]:
    """Collect state/program signatures plus optional scope metadata."""
    signatures: Dict[str, List[str]] = {}
    metadata: Dict[str, Dict[str, object]] = {}

    if config.marker_states:
        state_mgr = _resolve_annotation_manager(
            species=config.marker_species,
            tissue=None,
            marker_config=config.state_marker_config or f"cell_state_{config.marker_species}",
        )
        selected = state_mgr.select_cells(config.marker_states, include_children=True)
        for name, cell in selected.CELLS.items():
            if cell.markers:
                signatures[name] = list(cell.markers)
                metadata[name] = dict(cell.metadata)

    if config.custom_state_signatures:
        signatures.update(config.custom_state_signatures)
        if config.custom_state_metadata:
            metadata.update(config.custom_state_metadata)

    if config.state_signature_names or config.state_signature_categories:
        manager = get_marker_manager(
            species=config.marker_species,
            view="state_annotation",
            include_functional=True,
        )
        requested_categories = set(config.state_signature_categories)
        for category in config.state_signature_categories:
            matched = {
                name: list(cell.markers)
                for name, cell in manager.CELLS.items()
                if cell.metadata.get("category") == category and cell.markers
            }
            if not matched:
                raise KeyError(f"State signature category '{category}' not found")
            signatures.update(matched)
            for name in matched:
                metadata[name] = dict(manager.CELLS[name].metadata)
        for name in config.state_signature_names:
            if name not in manager.CELLS:
                raise KeyError(f"State signature '{name}' not found")
            cell = manager.CELLS[name]
            signatures[name] = list(cell.markers)
            metadata[name] = {
                **dict(cell.metadata),
                "requested_categories": sorted(requested_categories),
            }

    return signatures, metadata


def _state_applies_to_cell(
    state_meta: Optional[Dict[str, object]],
    lineage_label: str,
    subtype_label: str,
) -> bool:
    """Check whether a state program is valid for the current lineage/subtype context."""
    if not state_meta:
        return True

    scope = str(state_meta.get("scope", "all")).lower()
    applies_to = state_meta.get("applies_to", ["all"])
    if isinstance(applies_to, str):
        applies = [applies_to]
    else:
        applies = [str(x) for x in applies_to]

    applies_lower = [x.lower() for x in applies]
    if "all" in applies_lower or scope == "all":
        return True

    lineage_lower = str(lineage_label).lower()
    subtype_lower = str(subtype_label).lower()

    for allowed in applies_lower:
        if allowed and (allowed == lineage_lower or allowed in lineage_lower):
            return True
        if allowed and subtype_lower not in {"not_applicable", "unknown", "nan"}:
            if allowed == subtype_lower or allowed in subtype_lower:
                return True
    return False
