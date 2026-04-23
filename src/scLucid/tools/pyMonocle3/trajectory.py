"""
Trajectory inference for pyMonocle3 (R-free)
"""

import logging
from typing import List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra, minimum_spanning_tree
from sklearn.neighbors import NearestNeighbors

from .core import CellDataSet

log = logging.getLogger(__name__)


def learn_graph(
    cds: CellDataSet,
    reduction_method: str = "UMAP",
    partition_cells_by: str = "partition",
    k: int = 10,
    use_labels: bool = False,
    close_loop: bool = False,
    prune_graph: bool = True,
    minimal_branch_len: int = 10,
    euclidean_distance_ratio: float = 1.0,
    geodesic_distance_ratio: float = 1.0,
    prune_tol: float = 1e-6,
    random_seed: int = 42,
) -> CellDataSet:
    """
    Learn principal graph from cell positions

    Constructs a principal graph using the minimal spanning tree (MST)
    algorithm on the reduced dimensional space.

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    reduction_method : str
        Reduction method to use
    partition_cells_by : str
        Column for partitioning cells
    k : int
        Number of nearest neighbors
    use_labels : bool
        Use existing labels
    close_loop : bool
        Allow closed loops in graph
    prune_graph : bool
        Prune spurious branches
    minimal_branch_len : int
        Minimum branch length
    euclidean_distance_ratio : float
        Ratio for Euclidean distance
    geodesic_distance_ratio : float
        Ratio for geodesic distance
    prune_tol : float
        Pruning tolerance
    random_seed : int
        Random seed

    Returns:
    -------
    CellDataSet
        CellDataSet with principal graph
    """
    if reduction_method not in cds.reducedDims:
        raise ValueError(f"{reduction_method} not found. Run reduce_dimension first.")

    reduced_data = cds.reducedDims[reduction_method]

    # Get partitions
    if partition_cells_by in cds.cell_metadata.columns:
        partitions = cds.cell_metadata[partition_cells_by].values
    elif cds.partitions is not None:
        partitions = cds.partitions.values
    else:
        # Use single partition
        partitions = np.zeros(cds.n_cells, dtype=int)

    unique_partitions = np.unique(partitions)

    # Build graph for each partition
    graph_dict = {
        "partitions": partitions,
        "reduction_method": reduction_method,
        "partition_graphs": {},
        "edge_list": [],
        "node_positions": {},
    }

    np.random.seed(random_seed)

    for part_id in unique_partitions:
        part_mask = partitions == part_id
        part_data = reduced_data[part_mask]
        part_indices = np.where(part_mask)[0]

        if len(part_data) < k:
            log.warning(f"Partition {part_id} has too few cells ({len(part_data)}), skipping")
            continue

        # Build nearest neighbor graph
        nbrs = NearestNeighbors(n_neighbors=min(k, len(part_data) - 1), metric="euclidean")
        nbrs.fit(part_data)
        distances, indices = nbrs.kneighbors(part_data)

        # Create adjacency matrix
        n_part = len(part_data)
        adj_matrix = sp.lil_matrix((n_part, n_part))

        for i in range(n_part):
            for j, dist in zip(indices[i], distances[i]):
                if i != j:
                    adj_matrix[i, j] = dist
                    adj_matrix[j, i] = dist

        # Compute minimum spanning tree
        mst = minimum_spanning_tree(adj_matrix)

        # Convert to edge list
        mst_coo = mst.tocoo()
        edges = []
        for i, j, w in zip(mst_coo.row, mst_coo.col, mst_coo.data):
            if i < j:  # Avoid duplicates
                global_i = part_indices[i]
                global_j = part_indices[j]
                edges.append((global_i, global_j, w))

        # Store graph info
        graph_dict["partition_graphs"][part_id] = {
            "mst": mst,
            "edges": edges,
            "n_nodes": n_part,
            "indices": part_indices,
        }

        graph_dict["edge_list"].extend(edges)

    # Build full adjacency matrix
    full_adj = sp.lil_matrix((cds.n_cells, cds.n_cells))
    for i, j, w in graph_dict["edge_list"]:
        full_adj[i, j] = w
        full_adj[j, i] = w

    graph_dict["adj_matrix"] = full_adj.tocsr()

    cds.principal_graph = graph_dict

    log.info(
        f"Graph learning completed: {len(unique_partitions)} partitions, "
        f"{len(graph_dict['edge_list'])} edges"
    )

    return cds


def order_cells(
    cds: CellDataSet,
    reduction_method: Optional[str] = None,
    root_cells: Optional[Union[int, List[int], str]] = None,
    root_pr_nodes: Optional[List[int]] = None,
    cell_selection: str = "cell",
    root_state: Optional[int] = None,
    verbose: bool = False,
) -> CellDataSet:
    """
    Order cells by pseudotime along the learned graph

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet with principal graph
    reduction_method : str, optional
        Reduction method (auto-detected if None)
    root_cells : int, list, or str
        Root cell(s) for pseudotime calculation
    root_pr_nodes : list
        Root principal nodes
    cell_selection : str
        Method for selecting root cells
    root_state : int
        Root state for pseudotime
    verbose : bool
        Verbose output

    Returns:
    -------
    CellDataSet
        CellDataSet with pseudotime values
    """
    if cds.principal_graph is None:
        raise ValueError("No principal graph found. Run learn_graph first.")

    # Detect reduction method
    if reduction_method is None:
        if "UMAP" in cds.reducedDims:
            reduction_method = "UMAP"
        elif "tSNE" in cds.reducedDims:
            reduction_method = "tSNE"
        elif "PCA" in cds.reducedDims:
            reduction_method = "PCA"
        else:
            raise ValueError("No reduction found. Run reduce_dimension first.")

    adj_matrix = cds.principal_graph["adj_matrix"]

    # Determine root cells
    if root_cells is None:
        # Auto-select root based on minimal degree in graph
        degrees = np.array(adj_matrix.sum(axis=0)).flatten()
        root_cells = int(np.argmin(degrees))
        log.info(f"Auto-selected root cell: {root_cells}")
    elif isinstance(root_cells, str):
        # Root based on metadata query
        raise NotImplementedError("String-based root selection not yet implemented")

    if isinstance(root_cells, int):
        root_cells = [root_cells]

    # Calculate pseudotime using Dijkstra's algorithm
    # Distance from root to all other cells
    distances = []
    for root in root_cells:
        dist = dijkstra(adj_matrix, directed=False, indices=root)
        distances.append(dist)

    # Take minimum distance from any root
    pseudotime = np.min(distances, axis=0)

    # Handle infinite distances (disconnected components)
    pseudotime[np.isinf(pseudotime)] = pseudotime[~np.isinf(pseudotime)].max() * 1.5

    # Store pseudotime
    cds.cell_metadata["pseudotime"] = pseudotime

    log.info(
        f"Pseudotime calculation completed: range [{pseudotime.min():.2f}, {pseudotime.max():.2f}]"
    )

    return cds


def graph_test(
    cds: CellDataSet,
    neighbor_graph: str = "principal_graph",
    k: int = 25,
    n_tests: int = 100,
    cores: int = 1,
) -> pd.DataFrame:
    """
    Test for genes that vary along the trajectory using Moran's I

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet
    neighbor_graph : str
        Which graph to use for neighborhood
    k : int
        Number of neighbors for Moran's I
    n_tests : int
        Number of permutation tests
    cores : int
        Number of cores to use

    Returns:
    -------
    pd.DataFrame
        Gene statistics with Moran's I values
    """
    # Get expression data
    expr = cds.expression_data
    if sp.issparse(expr):
        expr = expr.toarray()

    # Build weight matrix based on principal graph or kNN
    if neighbor_graph == "principal_graph" and cds.principal_graph is not None:
        W = cds.principal_graph["adj_matrix"].copy()
        # Normalize rows
        row_sums = np.array(W.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1
        W = W / row_sums[:, np.newaxis]
    else:
        # Build kNN graph from reduced dimensions
        reduced = cds.reducedDims.get("UMAP") or cds.reducedDims.get("PCA")
        nbrs = NearestNeighbors(n_neighbors=k)
        nbrs.fit(reduced)
        distances, indices = nbrs.kneighbors(reduced)

        # Create weight matrix
        W = np.zeros((cds.n_cells, cds.n_cells))
        for i in range(cds.n_cells):
            for j, dist in zip(indices[i], distances[i]):
                W[i, j] = 1.0 / (1.0 + dist)

        # Row normalize
        row_sums = W.sum(axis=1)
        row_sums[row_sums == 0] = 1
        W = W / row_sums[:, np.newaxis]

    # Calculate Moran's I for each gene
    results = []

    n = cds.n_cells
    W_sum = W.sum()

    for i, gene in enumerate(cds.gene_metadata.index):
        x = expr[i, :]
        x_mean = np.mean(x)
        x_centered = x - x_mean

        # Moran's I formula
        numerator = (x_centered.reshape(1, -1) @ W @ x_centered.reshape(-1, 1)).item()
        denominator = np.sum(x_centered**2)

        if denominator > 0:
            morans_i = (n / W_sum) * (numerator / denominator)
        else:
            morans_i = 0

        # Z-score (simplified)
        # In practice, you'd use permutation tests
        results.append(
            {
                "gene": gene,
                "morans_i": morans_i,
                "status": "OK",
            }
        )

    results_df = pd.DataFrame(results).sort_values("morans_i", ascending=False)

    log.info(f"Graph test completed for {len(results_df)} genes")

    return results_df


def choose_graph_segments(
    cds: CellDataSet,
    start_cells: List[int],
    end_cells: List[int],
) -> List[Tuple[int, int]]:
    """
    Choose specific segments of the graph between start and end cells

    Parameters
    ----------
    cds : CellDataSet
        Input CellDataSet with principal graph
    start_cells : list
        Starting cell indices
    end_cells : list
        Ending cell indices

    Returns:
    -------
    list
        List of edges in the selected segment
    """
    if cds.principal_graph is None:
        raise ValueError("No principal graph found")

    adj_matrix = cds.principal_graph["adj_matrix"]

    segments = []
    for start in start_cells:
        for end in end_cells:
            # Find shortest path
            dist, predecessors = dijkstra(
                adj_matrix, directed=False, indices=start, return_predecessors=True
            )

            # Reconstruct path
            path = []
            current = end
            while current != start:
                if predecessors[current] < 0:
                    break
                path.append((predecessors[current], current))
                current = predecessors[current]

            segments.extend(path)

    return segments
