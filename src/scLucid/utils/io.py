"""
User-facing ``.h5ad`` reader for scLucid.

A thin wrapper around :func:`anndata.read_h5ad` that guarantees the
resulting AnnData satisfies the scLucid input contract: a populated
``layers["counts"]`` slot (when ``X`` looks integer-valued), uniqued gene
names, and any optional sample-level metadata lifted into
``adata.uns["sclucid"]["analysis_context"]``.

The 10x Genomics readers (single-sample and multi-sample) live in
:mod:`scLucid.utils.helpers` to share the multi-sample concatenation /
fallback parser with the legacy ``load_10x_data`` API.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union

import scanpy as sc
from anndata import AnnData

from .helpers import _attach_counts_layer, _attach_sample_metadata

log = logging.getLogger(__name__)


def read_h5ad(
    path: Union[str, Path],
    *,
    ensure_counts_layer: bool = True,
    sample_id: Optional[str] = None,
    species: Optional[str] = None,
    tissue: Optional[str] = None,
    tissue_type: Optional[str] = None,
    cancer_type: Optional[str] = None,
) -> AnnData:
    """
    Load a ``.h5ad`` file and ensure it is ready for ``run_pipeline()``.

    Wrapper around :func:`anndata.read_h5ad` that:

    - copies ``X`` into ``layers["counts"]`` if no counts layer exists and the
      matrix looks integer-valued (controlled by ``ensure_counts_layer``);
    - attaches optional sample-level metadata to ``.obs`` and the
      ``analysis_context``;
    - leaves the rest of the AnnData untouched.

    Parameters
    ----------
    path : str or Path
        Path to a ``.h5ad`` file.
    ensure_counts_layer : bool, default=True
        If ``True`` and ``adata.layers["counts"]`` is missing, copy ``X`` to it
        when ``X`` appears to contain raw counts. If ``X`` is already
        normalized/log-transformed, no copy happens and a warning is logged.
    sample_id, species, tissue, tissue_type, cancer_type : str, optional
        Sample-level annotations recorded on every cell of ``adata.obs`` and
        lifted into ``adata.uns["sclucid"]["analysis_context"]``.

    Returns:
    -------
    AnnData

    Examples:
    --------
    >>> import scLucid as scl
    >>> adata = scl.read_h5ad("data/pbmc3k.h5ad", species="human")
    >>> adata_processed = scl.run_pipeline(adata)
    """
    path_obj = Path(path).expanduser().resolve()
    if not path_obj.exists():
        raise FileNotFoundError(f"h5ad path does not exist: {path_obj}")

    log.info("Reading h5ad file: %s", path_obj)
    adata = sc.read_h5ad(path_obj)

    if ensure_counts_layer:
        _attach_counts_layer(adata)

    _attach_sample_metadata(
        adata,
        sample_id=sample_id,
        species=species,
        tissue=tissue,
        tissue_type=tissue_type,
        cancer_type=cancer_type,
    )
    return adata


__all__ = ["read_h5ad"]
