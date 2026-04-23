"""
Clustering functions for pyMonocle3 (R-free)
"""

import logging
from typing import List, Optional

import numpy as np
import pandas as pd

from .core import CellDataSet

log = logging.getLogger(__name__)


def cluster_cells(
    cds: CellDataSet,
    reduction_method: str = "UMAP",
    k: int = 20,
    cluster_method: str = "leiden",
    num_iter: int = 2,
    partition_qval: float = 0.05,
    weight: bool = False,
    resolution: Optional[float] = None,
    random_seed: int = 42,
) -> CellDataSet:
    """
    Cluster cells using graph-based clustering

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    reduction_method : str
        Reduction method to use for clustering
    k : int
        Number of nearest neighbors
    cluster_method : str
        Clustering method ("leiden" or "louvain")
    num_iter : int
        Number of iterations for Leiden algorithm
    partition_qval : float
        Q-value threshold for partition
    weight : bool
        Weight edges by distance
    resolution : float, optional
        Resolution parameter for clustering
    random_seed : int
        Random seed

    Returns:
    -------
    CellDataSet
        CellDataSet with cluster assignments
    """
    if reduction_method not in cds.reducedDims:
        raise ValueError(f"{reduction_method} not found. " f"Run reduce_dimension first.")

    reduced_data = cds.reducedDims[reduction_method]

    # Build nearest neighbor graph
    from sklearn.neighbors import kneighbors_graph

    adj_matrix = kneighbors_graph(
        reduced_data,
        n_neighbors=k,
        mode="connectivity",
        include_self=False,
    )

    # Convert to symmetric graph
    adj_matrix = adj_matrix.maximum(adj_matrix.T)

    # Run clustering
    if cluster_method == "leiden":
        try:
            import igraph as ig
            import leidenalg as la
        except ImportError:
            raise ImportError(
                "leidenalg and igraph are required. " "Install with: pip install leidenalg igraph"
            )

        # Convert to igraph
        g = ig.Graph.Adjacency((adj_matrix > 0).toarray().tolist(), mode="UNDIRECTED")

        # Run Leiden algorithm
        resolution = resolution or 1.0

        partition = la.find_partition(
            g,
            la.RBConfigurationVertexPartition,
            weights=adj_matrix[adj_matrix.nonzero()].tolist() if weight else None,
            n_iterations=num_iter,
            resolution_parameter=resolution,
            seed=random_seed,
        )

        clusters = pd.Series(partition.membership, index=cds.cell_metadata.index, name="cluster")

    elif cluster_method == "louvain":
        try:
            import igraph as ig
            import louvain
        except ImportError:
            raise ImportError(
                "louvain and igraph are required. " "Install with: pip install louvain igraph"
            )

        # Convert to igraph
        g = ig.Graph.Adjacency((adj_matrix > 0).toarray().tolist(), mode="UNDIRECTED")

        # Run Louvain algorithm
        partition = louvain.find_partition(
            g,
            louvain.ModularityVertexPartition,
            weights=adj_matrix[adj_matrix.nonzero()].tolist() if weight else None,
        )

        clusters = pd.Series(partition.membership, index=cds.cell_metadata.index, name="cluster")

    else:
        raise ValueError(f"Unknown cluster_method: {cluster_method}")

    cds.clusters = clusters
    cds.cell_metadata["cluster"] = clusters

    n_clusters = len(np.unique(clusters))
    log.info(f"Clustering completed: {n_clusters} clusters found using {cluster_method}")

    return cds


def partition_cells(
    cds: CellDataSet,
    reduction_method: str = "UMAP",
    partition_list: Optional[List] = None,
) -> CellDataSet:
    """
    Partition cells into trajectory segments

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    reduction_method : str
        Reduction method to use
    partition_list : list, optional
        Pre-defined partitions

    Returns:
    -------
    CellDataSet
        CellDataSet with partition assignments
    """
    if partition_list is not None:
        partitions = pd.Series(partition_list, index=cds.cell_metadata.index, name="partition")
    else:
        # Use clusters as partitions by default
        if cds.clusters is None:
            raise ValueError("No clusters found. Run cluster_cells first.")

        partitions = cds.clusters.copy()
        partitions.name = "partition"

    cds.partitions = partitions
    cds.cell_metadata["partition"] = partitions

    n_partitions = len(np.unique(partitions))
    log.info(f"Partitioning completed: {n_partitions} partitions")

    return cds


def group_cells(
    cds: CellDataSet,
    group_cells_by: str = "cluster",
) -> pd.DataFrame:
    """
    Group cells by a metadata column and compute statistics

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    group_cells_by : str
        Column to group by

    Returns:
    -------
    pd.DataFrame
        Group statistics
    """
    if group_cells_by not in cds.cell_metadata.columns:
        raise ValueError(f"Column '{group_cells_by}' not found in cell metadata")

    groups = cds.cell_metadata.groupby(group_cells_by)

    stats = pd.DataFrame(
        {
            "n_cells": groups.size(),
            "mean_size_factor": (
                groups["Size_Factor"].mean() if "Size_Factor" in cds.cell_metadata.columns else None
            ),
        }
    )

    return stats


def find_cluster_markers(
    cds: CellDataSet,
    group_cells_by: str = "cluster",
    test_type: str = "wilcox",
    verbose: bool = False,
) -> pd.DataFrame:
    """
    Find marker genes for each cluster

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    group_cells_by : str
        Column containing cluster assignments
    test_type : str
        Statistical test type
    verbose : bool
        Verbose output

    Returns:
    -------
    pd.DataFrame
        Marker genes with statistics
    """
    from scipy.stats import mannwhitneyu, ranksums

    if group_cells_by not in cds.cell_metadata.columns:
        raise ValueError(f"Column '{group_cells_by}' not found")

    clusters = cds.cell_metadata[group_cells_by].unique()
    expr = (
        cds.expression_data.toarray()
        if hasattr(cds.expression_data, "toarray")
        else cds.expression_data
    )

    markers = []

    for cluster in clusters:
        cluster_cells = cds.cell_metadata[group_cells_by] == cluster
        other_cells = ~cluster_cells

        cluster_expr = expr[:, cluster_cells]
        other_expr = expr[:, other_cells]

        for i, gene in enumerate(cds.gene_metadata.index):
            # Calculate statistics
            cluster_mean = np.mean(cluster_expr[i, :])
            other_mean = np.mean(other_expr[i, :])

            fc = cluster_mean / (other_mean + 1e-10)

            # Statistical test
            if test_type == "wilcox":
                try:
                    stat, pval = ranksums(cluster_expr[i, :], other_expr[i, :])
                except:
                    pval = 1.0
            elif test_type == "mannwhitney":
                try:
                    stat, pval = mannwhitneyu(
                        cluster_expr[i, :], other_expr[i, :], alternative="two-sided"
                    )
                except:
                    pval = 1.0
            else:
                pval = 1.0

            if cluster_mean > other_mean and pval < 0.05:
                markers.append(
                    {
                        "gene": gene,
                        "cluster": cluster,
                        "cluster_mean": cluster_mean,
                        "other_mean": other_mean,
                        "log2fc": np.log2(fc + 1e-10),
                        "pval": pval,
                    }
                )

    markers_df = pd.DataFrame(markers).sort_values(["cluster", "pval"])

    log.info(f"Found markers for {len(clusters)} clusters")

    return markers_df
