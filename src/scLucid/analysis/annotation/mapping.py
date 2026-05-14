"""Manual annotation mapping and label remapping."""

from __future__ import annotations

from typing import Dict, Optional, Union

import logging
import pandas as pd
from anndata import AnnData

from ...utils import sanitize_for_hdf5
from .utils import _read_json_file, _read_table_file

log = logging.getLogger(__name__)


def apply_annotation_mapping(
    adata: AnnData,
    cluster_key: str,
    mapping: Union[Dict[str, str], str],
    key_added: str = "cell_type",
) -> AnnData:
    """
    Apply AI/manual cluster-to-celltype mapping from dict, Excel (.xlsx/.xls), CSV, or JSON.
    Robust to type mismatches and common header conventions.

    Supported table schemas:
    - Columns named ['cluster','cell_type'] (preferred)
    - Or the first two columns are treated as [cluster, cell_type]
    """
    import datetime
    import os

    # 0) Validate cluster_key
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"'{cluster_key}' not found in adata.obs.")

    # 1) Load mapping
    if isinstance(mapping, str):
        ext = os.path.splitext(mapping)[1].lower()
        if ext in [".xlsx", ".xls", ".csv"]:
            df = _read_table_file(mapping)
            if {"cluster", "cell_type"}.issubset(df.columns):
                mapping_dict = pd.Series(df["cell_type"].values, index=df["cluster"]).to_dict()
            elif len(df.columns) >= 2:
                mapping_dict = pd.Series(df.iloc[:, 1].values, index=df.iloc[:, 0]).to_dict()
                log.info(
                    "Mapping file has no standard headers; used first two columns as [cluster, cell_type]."
                )
            else:
                raise ValueError(
                    "Mapping table must have either ['cluster','cell_type'] columns or at least 2 columns."
                )
            mapping_source = {"type": "file", "path": mapping}
        elif ext == ".json":
            mapping_dict = _read_json_file(mapping)
            mapping_source = {"type": "file", "path": mapping}
        else:
            raise ValueError("Unsupported mapping file format. Use .xlsx, .xls, .csv or .json.")
    elif isinstance(mapping, dict):
        mapping_dict = mapping
        mapping_source = {"type": "dict"}
    else:
        raise TypeError("mapping must be a dictionary or a file path string.")

    # 2) Ensure robust type matching (convert keys to string)
    source_clusters = adata.obs[cluster_key].astype(str)
    mapping_dict_str_keys = {str(k): v for k, v in mapping_dict.items()}

    # 3) Apply the mapping
    new_series = source_clusters.map(mapping_dict_str_keys)
    adata.obs[key_added] = pd.Categorical(new_series)

    # 4) Unmapped handling
    unmapped_mask = adata.obs[key_added].isnull()
    if unmapped_mask.any():
        unmapped_ids = source_clusters[unmapped_mask].unique().tolist()
        log.warning(
            f"{unmapped_mask.sum()} cells could not be mapped. "
            f"Missing cluster IDs in mapping: {unmapped_ids}"
        )
        adata.obs[key_added] = (
            adata.obs[key_added].cat.add_categories("Unmapped").fillna("Unmapped")
        )

    # 5) Store metadata snapshot
    annot_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    annot_ns[f"{key_added}_mapping"] = sanitize_for_hdf5(mapping_dict)
    annot_ns[f"{key_added}_mapping_meta"] = sanitize_for_hdf5(
        {
            "cluster_key": cluster_key,
            "created_at": datetime.datetime.now().isoformat(timespec="seconds"),
            **mapping_source,
        }
    )

    log.info(
        f"Applied mapping to column '{key_added}'. Categories: {list(adata.obs[key_added].cat.categories)}"
    )
    return adata


def remap_labels(
    adata: AnnData,
    column: str,
    mapping: Optional[Dict[str, str]] = None,
    where: Optional[pd.Series] = None,
    to: Optional[str] = None,
    in_place: bool = True,
    key_added: Optional[str] = None,
    tidy_categories: bool = True,
) -> AnnData:
    """
    Partially remap or correct labels in an existing categorical column of adata.obs.

    Modes:
    1) Dictionary-based rename:
       - mapping={'OldA':'NewA', 'OldB':'NewB'}
    2) Condition-based assign:
       - where: boolean Series aligned to adata.obs (True = replace), to='NewLabel'

    Args:
        column: Existing obs column to modify (e.g., 'celltype').
        mapping: Old -> New label mapping.
        where: Boolean mask of cells to set to `to`.
        to: Target label for cells where mask is True.
        in_place: If False, write to `key_added` and keep original intact.
        key_added: Required if in_place=False.
        tidy_categories: Drop unused categories after remap (clean legend).
    """
    if column not in adata.obs.columns:
        raise KeyError(f"Column '{column}' not found in adata.obs.")
    if not in_place and not key_added:
        raise ValueError("When in_place=False, you must provide key_added for the new column name.")

    src = adata.obs[column].astype(str)
    s = src.copy()

    changed_by_mapping = 0
    changed_by_where = 0

    if mapping is not None:
        before = s.copy()
        s = s.replace(mapping)
        changed_by_mapping = int((s != before).sum())

    if where is not None:
        if to is None:
            raise ValueError("When providing 'where', you must also provide 'to' label.")
        if not isinstance(where, pd.Series) or not where.index.equals(adata.obs.index):
            raise ValueError("'where' must be a boolean Series aligned to adata.obs index.")
        before = s.copy()
        s.loc[where] = to
        changed_by_where = int((s != before).sum())

    # Cast back to categorical
    cat = pd.Categorical(s)
    target_col = column if in_place else key_added
    adata.obs[target_col] = cat

    # Optionally tidy categories (drop unused)
    if tidy_categories:
        adata.obs[target_col] = adata.obs[target_col].cat.remove_unused_categories()

    # Audit
    annot_ns = (
        adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    )
    audit = annot_ns.setdefault("remap_audit", [])
    audit_entry = sanitize_for_hdf5(
        {
            "source_column": column,
            "target_column": target_col,
            "mode": (
                "mapping+where"
                if (mapping is not None and where is not None)
                else ("mapping" if mapping is not None else "where")
            ),
            "mapping": mapping if mapping is not None else None,
            "to": to if where is not None else None,
            "changed_by_mapping": changed_by_mapping,
            "changed_by_where": changed_by_where,
            "final_categories": list(adata.obs[target_col].cat.categories),
        }
    )
    audit.append(audit_entry)

    log.info(
        f"Remapped labels in '{target_col}': +{changed_by_mapping} (mapping), +{changed_by_where} (where). Categories: {list(adata.obs[target_col].cat.categories)}"
    )
    return adata
