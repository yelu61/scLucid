"""
Tumor progression trajectory analysis.

This module provides tools for analyzing tumor progression
trajectories and identifying transition states.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from anndata import AnnData
import logging

log = logging.getLogger(__name__)


class ProgressionAnalyzer:
    """
    Analyze tumor progression trajectories.

    Parameters
    ----------
    method : str
        Trajectory inference method ("pseudotime", "velocity")

    Attributes
    ----------
    pseudotime_ : pd.Series
        Pseudotime assignments
    """

    def __init__(
        self,
        method: str = "pseudotime",
    ):
        self.method = method
        self.pseudotime_: Optional[pd.Series] = None

    def analyze_tumor_progression(
        self,
        adata: AnnData,
        start_cell_type: Optional[str] = None,
        end_cell_type: Optional[str] = None,
        cell_type_key: str = "cell_type",
        embedding_key: str = "X_umap",
    ) -> pd.Series:
        """
        Analyze tumor progression using pseudotime.

        Parameters
        ----------
        adata : AnnData
            Expression data
        start_cell_type : str, optional
            Cell type at start of trajectory
        end_cell_type : str, optional
            Cell type at end of trajectory
        cell_type_key : str
            Column containing cell type annotations
        embedding_key : str
            Key for embedding coordinates

        Returns
        -------
        pd.Series
            Pseudotime assignments
        """
        # Use diffusion pseudotime or similar
        if embedding_key not in adata.obsm:
            raise ValueError(f"Embedding {embedding_key} not found. Run dimensionality reduction first.")

        X = adata.obsm[embedding_key]

        # Identify start cells
        if start_cell_type and cell_type_key in adata.obs.columns:
            root_cells = adata.obs[cell_type_key] == start_cell_type
            root_idx = np.where(root_cells)[0][0] if root_cells.any() else 0
        else:
            root_idx = 0

        # Calculate diffusion pseudotime
        from scipy.spatial.distance import cdist

        # Simple pseudotime based on distance from root
        distances = cdist(X, X[[root_idx]])
        pseudotime = distances.flatten()

        # Normalize to [0, 1]
        pseudotime = (pseudotime - pseudotime.min()) / (pseudotime.max() - pseudotime.min())

        self.pseudotime_ = pd.Series(pseudotime, index=adata.obs_names)

        return self.pseudotime_

    def identify_transition_states(
        self,
        adata: AnnData,
        pseudotime_key: str = "pseudotime",
        cell_type_key: str = "cell_type",
        window_size: float = 0.1,
    ) -> pd.DataFrame:
        """
        Identify transition states along progression.

        Parameters
        ----------
        adata : AnnData
            Expression data with pseudotime
        pseudotime_key : str
            Column containing pseudotime
        cell_type_key : str
            Column containing cell type annotations
        window_size : float
            Window size for sliding analysis

        Returns
        -------
        pd.DataFrame
            Transition state information
        """
        if pseudotime_key not in adata.obs.columns:
            raise ValueError(f"Pseudotime column {pseudotime_key} not found")

        pseudotime = adata.obs[pseudotime_key]

        transitions = []
        n_windows = int(1 / window_size)

        for i in range(n_windows):
            start_pt = i * window_size
            end_pt = (i + 1) * window_size

            mask = (pseudotime >= start_pt) & (pseudotime < end_pt)

            if mask.sum() < 5:
                continue

            window_adata = adata[mask]

            # Cell type composition
            if cell_type_key in window_adata.obs.columns:
                composition = window_adata.obs[cell_type_key].value_counts(normalize=True)

                transition = {
                    "pseudotime_start": start_pt,
                    "pseudotime_end": end_pt,
                    "n_cells": mask.sum(),
                    "dominant_type": composition.index[0],
                    "dominant_fraction": composition.iloc[0],
                    "n_types": len(composition),
                    "entropy": -np.sum(composition * np.log(composition + 1e-10)),
                }

                transitions.append(transition)

        return pd.DataFrame(transitions)

    def find_branch_points(
        self,
        adata: AnnData,
        embedding_key: str = "X_umap",
        min_branch_size: int = 50,
    ) -> pd.DataFrame:
        """
        Find branch points in tumor progression.

        Parameters
        ----------
        adata : AnnData
            Expression data
        embedding_key : str
            Key for embedding coordinates
        min_branch_size : int
            Minimum size for a branch

        Returns
        -------
        pd.DataFrame
            Branch point information
        """
        from sklearn.cluster import DBSCAN

        X = adata.obsm[embedding_key]

        # Use DBSCAN to identify distinct populations
        clustering = DBSCAN(eps=0.5, min_samples=min_branch_size).fit(X)
        labels = clustering.labels_

        # Find transition regions (high density between clusters)
        branches = []
        for label in np.unique(labels):
            if label == -1:
                continue

            mask = labels == label
            center = X[mask].mean(axis=0)

            branches.append({
                "branch_id": label,
                "n_cells": mask.sum(),
                "center_x": center[0],
                "center_y": center[1],
            })

        return pd.DataFrame(branches)


def analyze_tumor_progression(
    adata: AnnData,
    start_cell_type: Optional[str] = None,
    end_cell_type: Optional[str] = None,
    cell_type_key: str = "cell_type",
    embedding_key: str = "X_umap",
    key_added: str = "pseudotime",
) -> pd.Series:
    """
    Analyze tumor progression using pseudotime.

    Parameters
    ----------
    adata : AnnData
        Expression data
    start_cell_type : str, optional
        Cell type at start
    end_cell_type : str, optional
        Cell type at end
    cell_type_key : str
        Column with cell types
    embedding_key : str
        Key for embedding
    key_added : str
        Key for storing pseudotime

    Returns
    -------
    pd.Series
        Pseudotime assignments
    """
    analyzer = ProgressionAnalyzer()
    pseudotime = analyzer.analyze_tumor_progression(
        adata, start_cell_type, end_cell_type, cell_type_key, embedding_key
    )

    adata.obs[key_added] = pseudotime

    log.info(f"Calculated pseudotime for {len(pseudotime)} cells")

    return pseudotime


def identify_transition_states(
    adata: AnnData,
    pseudotime_key: str = "pseudotime",
    cell_type_key: str = "cell_type",
) -> pd.DataFrame:
    """
    Identify transition states along tumor progression.

    Parameters
    ----------
    adata : AnnData
        Expression data
    pseudotime_key : str
        Column with pseudotime
    cell_type_key : str
        Column with cell types

    Returns
    -------
    pd.DataFrame
        Transition states
    """
    analyzer = ProgressionAnalyzer()
    return analyzer.identify_transition_states(adata, pseudotime_key, cell_type_key)


def align_progression_trajectories(
    adata_list: List[AnnData],
    pseudotime_key: str = "pseudotime",
    gene_subset: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Align progression trajectories across multiple samples.

    Parameters
    ----------
    adata_list : list
        List of AnnData objects
    pseudotime_key : str
        Column with pseudotime
    gene_subset : list, optional
        Genes to use for alignment

    Returns
    -------
    pd.DataFrame
        Aligned trajectory data
    """
    aligned_data = []

    for i, adata in enumerate(adata_list):
        if pseudotime_key not in adata.obs.columns:
            continue

        pseudotime = adata.obs[pseudotime_key]

        # Bin pseudotime
        bins = np.linspace(0, 1, 11)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        for j, (start, end) in enumerate(zip(bins[:-1], bins[1:])):
            mask = (pseudotime >= start) & (pseudotime < end)

            if mask.sum() < 5:
                continue

            bin_adata = adata[mask]

            if gene_subset:
                genes = [g for g in gene_subset if g in bin_adata.var_names]
                expr = bin_adata[:, genes].X.mean(axis=0)
            else:
                expr = bin_adata.X.mean(axis=0)

            if hasattr(expr, 'toarray'):
                expr = expr.toarray().flatten()

            aligned_data.append({
                "sample": i,
                "pseudotime_bin": bin_centers[j],
                "mean_expression": expr.mean(),
            })

    return pd.DataFrame(aligned_data)
