"""
Core data structures for pyMonocle3 (R-free)
"""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any, Union, List
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData


@dataclass
class CellDataSet:
    """
    Python equivalent of Monocle3's CellDataSet

    A CellDataSet stores single-cell expression data along with
    metadata about cells and genes, dimensionality reduction results,
    cluster assignments, and trajectory graph information.

    Parameters
    ----------
    expression_data : np.ndarray or sp.spmatrix
        Gene expression matrix (genes x cells)
    cell_metadata : pd.DataFrame
        Cell-level metadata (cell x attributes)
    gene_metadata : pd.DataFrame
        Gene-level metadata (gene x attributes)
    reducedDims : dict, optional
        Dictionary of dimensionality reduction results
    clusters : pd.Series, optional
        Cluster assignments for each cell
    partitions : pd.Series, optional
        Partition assignments for trajectory analysis
    principal_graph : dict, optional
        Principal graph for trajectory inference
    """
    expression_data: Union[np.ndarray, sp.spmatrix]
    cell_metadata: pd.DataFrame
    gene_metadata: pd.DataFrame
    reducedDims: Dict[str, np.ndarray] = field(default_factory=dict)
    clusters: Optional[pd.Series] = None
    partitions: Optional[pd.Series] = None
    principal_graph: Optional[Dict] = None
    preprocessing_params: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate dimensions after initialization"""
        n_genes, n_cells = self.expression_data.shape

        if len(self.cell_metadata) != n_cells:
            raise ValueError(
                f"Cell metadata ({len(self.cell_metadata)} rows) doesn't match "
                f"expression data ({n_cells} cells)"
            )
        if len(self.gene_metadata) != n_genes:
            raise ValueError(
                f"Gene metadata ({len(self.gene_metadata)} rows) doesn't match "
                f"expression data ({n_genes} genes)"
            )

        # Normalize purely positional default indices to canonical names.
        if list(self.cell_metadata.index) == list(range(n_cells)):
            self.cell_metadata.index = [f"cell_{i}" for i in range(n_cells)]
        if list(self.gene_metadata.index) == list(range(n_genes)):
            self.gene_metadata.index = [f"gene_{i}" for i in range(n_genes)]

    @property
    def n_cells(self) -> int:
        """Number of cells in the dataset"""
        return self.expression_data.shape[1]

    @property
    def n_genes(self) -> int:
        """Number of genes in the dataset"""
        return self.expression_data.shape[0]

    def __repr__(self) -> str:
        return (
            f"CellDataSet({self.n_cells} cells x {self.n_genes} genes, "
            f"{len(self.reducedDims)} reductions, "
            f"clusters={'yes' if self.clusters is not None else 'no'}, "
            f"graph={'yes' if self.principal_graph is not None else 'no'})"
        )

    def copy(self) -> "CellDataSet":
        """Create a deep copy of the CellDataSet"""
        return CellDataSet(
            expression_data=self.expression_data.copy()
            if sp.issparse(self.expression_data)
            else self.expression_data.copy(),
            cell_metadata=self.cell_metadata.copy(),
            gene_metadata=self.gene_metadata.copy(),
            reducedDims={k: v.copy() for k, v in self.reducedDims.items()},
            clusters=self.clusters.copy() if self.clusters is not None else None,
            partitions=self.partitions.copy() if self.partitions is not None else None,
            principal_graph=self.principal_graph.copy()
            if self.principal_graph is not None
            else None,
            preprocessing_params=self.preprocessing_params.copy(),
        )

    def save(self, filepath: str):
        """Save CellDataSet to file"""
        import pickle

        with open(filepath, "wb") as f:
            pickle.dump(self, f)

    @staticmethod
    def load(filepath: str) -> "CellDataSet":
        """Load CellDataSet from file"""
        import pickle

        with open(filepath, "rb") as f:
            return pickle.load(f)


def new_cell_data_set(
    expression_data: Union[np.ndarray, sp.spmatrix, pd.DataFrame],
    cell_metadata: Optional[pd.DataFrame] = None,
    gene_metadata: Optional[pd.DataFrame] = None,
) -> CellDataSet:
    """
    Create a new CellDataSet from expression data

    Parameters
    ----------
    expression_data : array-like
        Expression matrix. If DataFrame, index=genes, columns=cells.
        If numpy array or sparse matrix, shape is (genes, cells).
    cell_metadata : pd.DataFrame, optional
        Cell metadata. If None, creates default metadata.
    gene_metadata : pd.DataFrame, optional
        Gene metadata. If None, creates default metadata.

    Returns
    -------
    CellDataSet
        New CellDataSet object
    """
    # Convert DataFrame to matrix + metadata
    if isinstance(expression_data, pd.DataFrame):
        gene_names = expression_data.index.tolist()
        cell_names = expression_data.columns.tolist()
        expression_matrix = expression_data.values
    else:
        n_genes, n_cells = expression_data.shape
        gene_names = [f"gene_{i}" for i in range(n_genes)]
        cell_names = [f"cell_{i}" for i in range(n_cells)]
        expression_matrix = expression_data

    # Create default metadata if not provided
    if cell_metadata is None:
        cell_metadata = pd.DataFrame(index=cell_names)
    else:
        cell_metadata = cell_metadata.copy()
        if cell_metadata.index.tolist() != cell_names:
            cell_metadata.index = cell_names

    if gene_metadata is None:
        gene_metadata = pd.DataFrame(index=gene_names)
    else:
        gene_metadata = gene_metadata.copy()
        if gene_metadata.index.tolist() != gene_names:
            gene_metadata.index = gene_names

    return CellDataSet(
        expression_data=expression_matrix
        if not isinstance(expression_matrix, pd.DataFrame)
        else expression_matrix.values,
        cell_metadata=cell_metadata,
        gene_metadata=gene_metadata,
    )


def create_cds_from_scanpy(
    adata: AnnData,
    layer: Optional[str] = None,
) -> CellDataSet:
    """
    Create CellDataSet from Scanpy AnnData

    Parameters
    ----------
    adata : AnnData
        Scanpy AnnData object
    layer : str, optional
        Layer to use. If None, uses .X

    Returns
    -------
    CellDataSet
        New CellDataSet object
    """
    if layer is not None:
        expr = adata.layers[layer].T
    else:
        expr = adata.X.T

    # Convert to dense if sparse
    if sp.issparse(expr):
        expr = expr.toarray()

    return CellDataSet(
        expression_data=expr,
        cell_metadata=adata.obs.copy(),
        gene_metadata=adata.var.copy(),
    )


def export_to_scanpy(cds: CellDataSet) -> AnnData:
    """
    Export CellDataSet to Scanpy AnnData

    Parameters
    ----------
    cds : CellDataSet
        CellDataSet to export

    Returns
    -------
    AnnData
        Scanpy AnnData object
    """
    adata = AnnData(
        X=cds.expression_data.T,
        obs=cds.cell_metadata,
        var=cds.gene_metadata,
    )

    # Add dimensionality reductions
    for key, value in cds.reducedDims.items():
        adata.obsm[f"X_{key}"] = value

    # Add clusters
    if cds.clusters is not None:
        adata.obs["monocle_clusters"] = cds.clusters

    # Add partitions
    if cds.partitions is not None:
        adata.obs["monocle_partitions"] = cds.partitions

    return adata
