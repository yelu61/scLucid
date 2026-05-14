"""Cell-type marker scoring helpers."""

from __future__ import annotations

from typing import Optional, Union

import logging
import scanpy as sc
from anndata import AnnData

from ...utils import Manager, sanitize_for_hdf5, use_layer_as_X

log = logging.getLogger(__name__)


def score_cell_types(
    adata: AnnData,
    marker_config: Union[str, Manager],
    layer: Optional[str] = "normalized",
    use_raw: bool = True,
    min_genes: int = 3,
    ctrl_size: int = 50,
    score_name_suffix: str = "_score",
    copy: bool = False,
) -> AnnData:
    """
    Score cells for cell type marker gene sets.
    Adds score columns to adata.obs.
    """
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
    log.info(f"Scored {n_scored} cell types ({n_skipped} skipped).")
    adata.uns.setdefault("sclucid", {}).setdefault("analysis", {}).setdefault("annotation", {})
    scoring_params = sanitize_for_hdf5(
        {
            "use_raw": use_raw,
            "layer": layer,
            "min_genes": min_genes,
            "ctrl_size": ctrl_size,
        }
    )
    adata.uns["sclucid"]["analysis"]["annotation"]["scoring_params"] = scoring_params
    return adata
