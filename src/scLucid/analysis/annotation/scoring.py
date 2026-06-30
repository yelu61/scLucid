"""Cell-type marker scoring helpers."""

from __future__ import annotations

import logging
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from scipy.stats import rankdata

from ...utils import Manager, sanitize_for_hdf5, use_layer_as_X

log = logging.getLogger(__name__)


def _score_genes_aucell(adata: AnnData, gene_list: list[str], **kwargs) -> pd.Series:
    """
    Compute a classic AUCell rank-based AUC score per cell for ``gene_list``.

    Genes are ranked per cell by expression (rank 1 = highest). The AUCell
    score is the area under the recovery curve obtained by walking the ranked
    list. It is equivalent to the normalized Mann-Whitney U statistic: a
    random signature scores ~0.5, a perfectly top-ranked signature scores ~1.0,
    and a bottom-ranked signature scores ~0.0.

    Returns a pandas Series indexed by ``adata.obs_names``.
    """
    X = _extract_expression_matrix(adata)
    gene_idx = np.array(
        [adata.var_names.get_loc(g) for g in gene_list if g in adata.var_names]
    )
    if len(gene_idx) == 0:
        return pd.Series(np.nan, index=adata.obs_names, name="aucell_score")

    n_genes = X.shape[1]
    set_size = len(gene_idx)
    # rank 1 = lowest expression (ascending); highly expressed genes get large ranks
    ranks = np.apply_along_axis(lambda row: rankdata(row, method="average"), 1, X)
    set_ranks = ranks[:, gene_idx]
    # Classic AUCell AUC (normalized Mann-Whitney U). Random ~0.5, top ~1.0, bottom ~0.0.
    auc = (
        set_ranks.sum(axis=1) - set_size * (set_size + 1.0) / 2.0
    ) / (set_size * (n_genes - set_size))
    return pd.Series(auc, index=adata.obs_names, name="aucell_score")


def _score_genes_ucell(
    adata: AnnData, gene_list: list[str], max_rank: int = 1000, **kwargs
) -> pd.Series:
    """
    Compute a UCell rank-sum score per cell for ``gene_list``.

    Genes are ranked per cell by ascending expression (rank 1 = lowest), so
    highly expressed genes receive large ranks. The UCell score is the sum of
    the ranks of the signature genes divided by ``max_rank * |S|`` and capped
    at 1. No control genes are required, making it suitable for small or
    sparse signatures.

    Returns a pandas Series indexed by ``adata.obs_names``.
    """
    X = _extract_expression_matrix(adata)
    gene_idx = np.array(
        [adata.var_names.get_loc(g) for g in gene_list if g in adata.var_names]
    )
    if len(gene_idx) == 0:
        return pd.Series(np.nan, index=adata.obs_names, name="ucell_score")

    set_size = len(gene_idx)
    # rank 1 = lowest expression (ascending)
    ranks = np.apply_along_axis(lambda row: rankdata(row, method="average"), 1, X)
    set_ranks = ranks[:, gene_idx]
    ucell = np.minimum(1.0, set_ranks.sum(axis=1) / (max_rank * set_size))
    return pd.Series(ucell, index=adata.obs_names, name="ucell_score")


def _extract_expression_matrix(adata: AnnData) -> np.ndarray:
    """Return a dense numpy expression matrix from ``adata.X``."""
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def score_cell_types(
    adata: AnnData,
    marker_config: Union[str, Manager],
    layer: Optional[str] = "normalized",
    use_raw: bool = True,
    min_genes: int = 3,
    ctrl_size: int = 50,
    score_name_suffix: str = "_score",
    scoring_backend: Literal["scanpy", "aucell", "ucell"] = "scanpy",
    copy: bool = False,
) -> AnnData:
    """
    Score cells for cell type marker gene sets.

    Adds score columns to ``adata.obs``. Three backends are supported:

    - ``"scanpy"`` (default): ``scanpy.tl.score_genes`` with control-gene
      sampling. Best for larger signatures.
    - ``"aucell"``: rank-based AUC without control genes.
    - ``"ucell"``: Mann-Whitney rank-sum statistic without control genes.
      Preferred for small or sparse signatures because it avoids control-gene
      sampling noise.

    Parameters
    ----------
    scoring_backend
        One of ``"scanpy"``, ``"aucell"``, ``"ucell"``.
    """
    if scoring_backend not in {"scanpy", "aucell", "ucell"}:
        raise ValueError(f"Unknown scoring_backend: {scoring_backend}")
    if copy:
        adata = adata.copy()
    if isinstance(marker_config, str):
        mgr = Manager(marker_config)
    elif isinstance(marker_config, Manager):
        mgr = marker_config
    else:
        raise TypeError("marker_config must be a file path or Manager instance.")
    mgr.intersect_with(adata.raw if use_raw and adata.raw is not None else adata)
    n_scored, n_skipped = 0, 0

    expression_adata = (
        adata.raw.to_adata() if use_raw and adata.raw is not None else adata
    )

    if scoring_backend == "scanpy":
        if use_raw:
            for cell_type, cell in mgr.CELLS.items():
                if len(cell.markers) >= min_genes:
                    sc.tl.score_genes(
                        adata,
                        cell.markers,
                        score_name=f"{cell_type}{score_name_suffix}",
                        use_raw=True,
                        ctrl_size=ctrl_size,
                    )
                    n_scored += 1
                else:
                    n_skipped += 1
        else:
            with use_layer_as_X(adata, layer):
                for cell_type, cell in mgr.CELLS.items():
                    if len(cell.markers) >= min_genes:
                        sc.tl.score_genes(
                            adata,
                            cell.markers,
                            score_name=f"{cell_type}{score_name_suffix}",
                            ctrl_size=ctrl_size,
                        )
                        n_scored += 1
                    else:
                        n_skipped += 1
    else:
        scorer = _score_genes_aucell if scoring_backend == "aucell" else _score_genes_ucell
        if use_raw:
            for cell_type, cell in mgr.CELLS.items():
                present_markers = [
                    g for g in cell.markers if g in expression_adata.var_names
                ]
                if len(present_markers) >= min_genes:
                    scores = scorer(expression_adata, present_markers)
                    adata.obs[f"{cell_type}{score_name_suffix}"] = scores.reindex(
                        adata.obs_names
                    ).values
                    n_scored += 1
                else:
                    n_skipped += 1
        else:
            with use_layer_as_X(expression_adata, layer):
                for cell_type, cell in mgr.CELLS.items():
                    present_markers = [
                        g for g in cell.markers if g in expression_adata.var_names
                    ]
                    if len(present_markers) >= min_genes:
                        scores = scorer(expression_adata, present_markers)
                        adata.obs[f"{cell_type}{score_name_suffix}"] = scores.reindex(
                            adata.obs_names
                        ).values
                        n_scored += 1
                    else:
                        n_skipped += 1

    log.info(f"Scored {n_scored} cell types ({n_skipped} skipped).")
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    scoring_params = sanitize_for_hdf5(
        {
            "backend": scoring_backend,
            "use_raw": use_raw,
            "layer": layer,
            "min_genes": min_genes,
            "ctrl_size": ctrl_size,
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"]["scoring_params"] = scoring_params
    return adata
