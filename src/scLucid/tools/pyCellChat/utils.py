"""
Utility functions for CellChat (R-free)
"""

import pandas as pd
import numpy as np
from typing import Optional, List
import logging

log = logging.getLogger(__name__)


def create_cellchat_from_scanpy(
    adata,
    group_by: str,
    use_raw: bool = False,
    spatial_key: Optional[str] = None
):
    """
    Create CellChat object from Scanpy AnnData object

    Parameters
    ----------
    adata : anndata.AnnData
        Scanpy AnnData object
    group_by : str
        Column in adata.obs for cell grouping
    use_raw : bool
        Use raw counts
    spatial_key : Optional[str]
        Key for spatial coordinates in adata.obsm
    """
    # Get expression data
    if use_raw and adata.raw is not None:
        expr = adata.raw.X.T
        gene_names = adata.raw.var_names.tolist()
    else:
        expr = adata.X.T
        gene_names = adata.var_names.tolist()

    # Create expression DataFrame
    expr_df = pd.DataFrame(
        expr,
        index=gene_names,
        columns=adata.obs_names
    )

    # Get metadata
    if group_by not in adata.obs.columns:
        raise ValueError(f"Column '{group_by}' not found in adata.obs")

    meta = adata.obs[[group_by]].copy()

    # Get spatial coordinates if available
    spatial_coords = None
    if spatial_key is not None and spatial_key in adata.obsm:
        spatial_coords = pd.DataFrame(
            adata.obsm[spatial_key],
            index=adata.obs_names
        )

    # Create CellChat object
    from .core import CellChat

    cellchat = CellChat(
        data=expr_df,
        meta=meta,
        group_by=group_by,
        spatial_coords=spatial_coords
    )

    log.info(f"Created CellChat from AnnData: {adata.n_obs} cells, {adata.n_vars} genes")
    return cellchat


def merge_cellchat_objects(
    cellchat_list: list,
    add_names: Optional[List] = None
):
    """
    Merge multiple CellChat objects for comparison

    Parameters
    ----------
    cellchat_list : list
        List of CellChat objects
    add_names : Optional[List]
        Names for each condition
    """
    if add_names is None:
        add_names = [f"Condition_{i+1}" for i in range(len(cellchat_list))]

    if len(add_names) != len(cellchat_list):
        raise ValueError("Length of add_names must match cellchat_list")

    merged = {
        'objects': cellchat_list,
        'names': add_names
    }

    return merged


def export_to_cytoscape(
    cellchat_obj,
    pathway: str,
    filename: str,
    thresh: float = 0.05
):
    """
    Export network to Cytoscape format

    Parameters
    ----------
    cellchat_obj : CellChat
        CellChat object
    pathway : str
        Pathway to export
    filename : str
        Output filename
    thresh : float
        Threshold for edges
    """
    if pathway not in cellchat_obj.netP['prob']:
        raise ValueError(f"Pathway {pathway} not found")

    prob = cellchat_obj.netP['prob'][pathway]
    groups = cellchat_obj.unique_groups

    # Create edge list
    edges = []
    for i, source in enumerate(groups):
        for j, target in enumerate(groups):
            if prob[i, j] > thresh:
                edges.append({
                    'source': source,
                    'target': target,
                    'weight': prob[i, j],
                    'interaction': 'pp'
                })

    edge_df = pd.DataFrame(edges)
    edge_df.to_csv(filename, index=False)

    log.info(f"Exported {len(edges)} edges to {filename}")


def save_cellchat(cellchat_obj, filename: str):
    """Save CellChat object to file"""
    import pickle
    with open(filename, 'wb') as f:
        pickle.dump(cellchat_obj, f)
    log.info(f"Saved CellChat object to {filename}")


def load_cellchat(filename: str):
    """Load CellChat object from file"""
    import pickle
    with open(filename, 'rb') as f:
        cellchat_obj = pickle.load(f)
    log.info(f"Loaded CellChat object from {filename}")
    return cellchat_obj
