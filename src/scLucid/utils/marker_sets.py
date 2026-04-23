"""Utilities for sanitizing and filtering marker or gene-set dictionaries."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any, Dict, List, Tuple, Union

import numpy as np
import pandas as pd

GeneSetTree = Dict[str, Any]


def flatten_marker_dict(markers: Mapping[str, Any], prefix: str = "") -> Dict[str, List[str]]:
    """
    Flatten a nested marker dictionary into ``{path: [genes...]}``.

    Parameters
    ----------
    markers : mapping
        Nested marker dictionary.
    prefix : str
        Optional prefix used during recursive traversal.

    Returns:
    -------
    dict
        Flattened dictionary where nested keys are joined by ``.``.
    """
    flat: Dict[str, List[str]] = {}
    for key, value in markers.items():
        name = str(key) if not prefix else f"{prefix}.{key}"
        if isinstance(value, Mapping):
            flat.update(flatten_marker_dict(value, prefix=name))
        elif _is_gene_collection(value):
            flat[name] = _normalize_gene_collection(value)
    return flat


def filter_marker_dict(
    markers: Mapping[str, Any],
    var_names: Union[pd.Index, List[str], np.ndarray],
    *,
    uppercase: bool = True,
    drop_empty: bool = True,
    return_missing: bool = False,
) -> Union[GeneSetTree, Tuple[GeneSetTree, Dict[str, List[str]]]]:
    """
    Filter marker dictionaries against available feature names.

    Parameters
    ----------
    markers : mapping
        Marker dictionary. Leaf values can be lists, tuples, sets, numpy arrays,
        or pandas Index objects.
    var_names : sequence-like
        Available gene names to match against.
    uppercase : bool
        If True, perform case-insensitive matching by uppercasing both sides.
    drop_empty : bool
        If True, remove branches with no retained genes.
    return_missing : bool
        If True, also return a ``{path: missing_genes}`` dictionary.

    Returns:
    -------
    dict or tuple(dict, dict)
        Filtered marker dictionary, optionally with missing-gene report.
    """
    present = {str(g).upper() if uppercase else str(g) for g in pd.Index(var_names).astype(str)}
    missing: Dict[str, List[str]] = {}

    def _recurse(node: Any, path: str = "") -> Any:
        if isinstance(node, Mapping):
            out: Dict[str, Any] = {}
            for key, value in node.items():
                child_path = f"{path}.{key}" if path else str(key)
                child = _recurse(value, path=child_path)
                if not (drop_empty and _is_empty_container(child)):
                    out[str(key)] = child
            return out

        if _is_gene_collection(node):
            kept: List[str] = []
            missed: List[str] = []
            for gene in _normalize_gene_collection(node):
                probe = gene.upper() if uppercase else gene
                if probe in present:
                    kept.append(gene)
                else:
                    missed.append(gene)
            if return_missing and missed:
                missing[path] = missed
            return kept

        return node

    filtered = _recurse(markers)
    if return_missing:
        return filtered, missing
    return filtered


def _is_gene_collection(value: Any) -> bool:
    """Return True when a node should be treated as a gene collection leaf."""
    return isinstance(value, (list, tuple, set, np.ndarray, pd.Index))


def _normalize_gene_collection(value: Iterable[Any]) -> List[str]:
    """Normalize an iterable of genes into a clean string list."""
    genes: List[str] = []
    for item in value:
        if item is None:
            continue
        gene = str(item).strip()
        if gene:
            genes.append(gene)
    return genes


def _is_empty_container(value: Any) -> bool:
    """Return True for empty dict/list leaves produced during filtering."""
    return value == {} or value == []
