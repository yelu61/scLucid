"""Spatial neighbor graph construction."""

from __future__ import annotations

import logging
from typing import Optional

import numpy as np
import scipy.sparse as sp
from anndata import AnnData
from sklearn.neighbors import NearestNeighbors, radius_neighbors_graph

from .config import SpatialDiagnosticsConfig, SpatialNeighborsConfig
from .diagnostics import diagnose_spatial_data_quality

log = logging.getLogger(__name__)


def build_spatial_neighbors(
    adata: AnnData,
    config: Optional[SpatialNeighborsConfig] = None,
) -> AnnData:
    """Build a spatial neighbor graph from ``adata.obsm[spatial_key]``.

    Parameters
    ----------
    adata
        AnnData with spatial coordinates.
    config
        Neighbor configuration.

    Returns
    -------
    AnnData
        ``adata`` with ``obsp["spatial_connectivities"]`` and
        ``obsp["spatial_distances"]`` populated.
    """
    if config is None:
        config = SpatialNeighborsConfig()

    diag = diagnose_spatial_data_quality(adata, config=SpatialDiagnosticsConfig(spatial_key=config.spatial_key))
    if not diag["passed"]:
        log.warning("Spatial data quality warnings: %s", diag["warnings"])

    coords = np.asarray(adata.obsm[config.spatial_key])[:, :2]
    n = coords.shape[0]

    if config.method == "knn":
        nbrs = NearestNeighbors(n_neighbors=min(config.n_neigh, n - 1), metric="euclidean")
        nbrs.fit(coords)
        distances, indices = nbrs.kneighbors(coords)
        row_idx = np.repeat(np.arange(n), distances.shape[1])
        col_idx = indices.ravel()
        data = distances.ravel()
        dist_mat = sp.csr_matrix((data, (row_idx, col_idx)), shape=(n, n))
        conn_mat = dist_mat.copy()
        conn_mat.data = np.ones_like(conn_mat.data)
    elif config.method == "radius":
        if config.radius is None:
            raise ValueError("radius method requires config.radius > 0")
        dist_mat = radius_neighbors_graph(coords, radius=config.radius, mode="distance")
        conn_mat = radius_neighbors_graph(coords, radius=config.radius, mode="connectivity")
    else:
        raise ValueError(f"Unknown spatial neighbor method: {config.method}")

    # Symmetrize
    conn_mat = conn_mat.maximum(conn_mat.T)
    dist_mat = dist_mat.maximum(dist_mat.T)

    adata.obsp["spatial_connectivities"] = conn_mat
    adata.obsp["spatial_distances"] = dist_mat
    adata.uns[config.key_added] = {
        "params": config.to_dict(),
        "n_neighbors_total": int(conn_mat.sum()),
    }

    return adata
