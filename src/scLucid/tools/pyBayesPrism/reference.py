"""
Reference data handling for BayesPrism (R-free)

Handles scRNA-seq reference data processing and cell type profile generation.
"""

import numpy as np
import pandas as pd
from scipy import sparse
from typing import Optional, Union, List, Dict, Set
import logging

from .config import ReferenceConfig

log = logging.getLogger(__name__)


class BayesPrismReference:
    """
    BayesPrism reference data class

    Processes scRNA-seq reference data and generates cell type-specific
    expression profiles for deconvolution.

    Parameters
    ----------
    reference : pd.DataFrame or sparse matrix
        Reference expression matrix (genes x cells)
    cell_type_labels : pd.Series
        Cell type labels for each cell
    cell_state_labels : pd.Series, optional
        Cell state labels (defaults to cell_type_labels)
    input_type : str
        Type of input ("count.matrix" or "GEP")
    pseudo_min : float
        Minimum pseudo-count to avoid zeros

    Attributes
    ----------
    phi : np.ndarray
        Cell type-specific expression profiles (genes x cell_types)
    cell_types : List[str]
        Unique cell type names
    cell_states : List[str]
        Unique cell state names

    Examples
    --------
    >>> ref = BayesPrismReference(
    ...     reference=ref_counts,
    ...     cell_type_labels=cell_types,
    ...     input_type="count.matrix"
    ... )
    >>> print(f"Reference contains {len(ref.cell_types)} cell types")
    """

    def __init__(
        self,
        reference: Union[pd.DataFrame, sparse.spmatrix, np.ndarray],
        cell_type_labels: pd.Series,
        cell_state_labels: Optional[pd.Series] = None,
        input_type: str = "count.matrix",
        pseudo_min: float = 1e-8,
    ):
        self.config = ReferenceConfig(
            input_type=input_type,
            pseudo_min=pseudo_min
        )

        # Store labels
        self.cell_type_labels = pd.Series(cell_type_labels)
        self.cell_state_labels = (
            pd.Series(cell_state_labels)
            if cell_state_labels is not None
            else self.cell_type_labels.copy()
        )

        # Validate lengths match
        n_cells = len(self.cell_type_labels)
        if len(self.cell_state_labels) != n_cells:
            raise ValueError(
                f"cell_state_labels length ({len(self.cell_state_labels)}) "
                f"must match cell_type_labels ({n_cells})"
            )

        # Convert to sparse matrix
        self.reference_matrix, self.gene_names, self.cell_names = (
            self._convert_to_sparse(reference, n_cells)
        )

        # Get unique types
        self.cell_types = self.cell_type_labels.unique().tolist()
        self.cell_states = self.cell_state_labels.unique().tolist()

        # Generate reference profiles
        self._generate_reference_profile()

        log.info(
            f"Reference initialized: {self.reference_matrix.shape[0]} genes, "
            f"{n_cells} cells, {len(self.cell_types)} cell types"
        )

    def _convert_to_sparse(
        self,
        reference: Union[pd.DataFrame, sparse.spmatrix, np.ndarray],
        n_cells: int,
    ) -> tuple:
        """Convert reference data to sparse matrix format"""
        if isinstance(reference, pd.DataFrame):
            if reference.shape[1] != n_cells:
                raise ValueError(
                    f"Reference columns ({reference.shape[1]}) must match "
                    f"number of cells ({n_cells})"
                )
            ref_matrix = sparse.csr_matrix(reference.values)
            gene_names = reference.index.tolist()
            cell_names = reference.columns.tolist()
        elif isinstance(reference, np.ndarray):
            if reference.shape[1] != n_cells:
                raise ValueError(
                    f"Reference shape {reference.shape} doesn't match {n_cells} cells"
                )
            ref_matrix = sparse.csr_matrix(reference)
            gene_names = [f"gene_{i}" for i in range(reference.shape[0])]
            cell_names = [f"cell_{i}" for i in range(n_cells)]
        elif sparse.issparse(reference):
            ref_matrix = reference.tocsr()
            gene_names = [f"gene_{i}" for i in range(reference.shape[0])]
            cell_names = [f"cell_{i}" for i in range(n_cells)]
        else:
            raise TypeError(f"Unsupported reference type: {type(reference)}")

        return ref_matrix, gene_names, cell_names

    def _generate_reference_profile(self) -> None:
        """Generate cell type-specific expression profiles (phi matrix)"""
        n_genes = self.reference_matrix.shape[0]
        n_cell_types = len(self.cell_types)

        # Initialize reference profile matrix
        self.phi = np.zeros((n_genes, n_cell_types))

        for i, cell_type in enumerate(self.cell_types):
            # Get cells of this type
            cell_mask = (self.cell_type_labels == cell_type).to_numpy(dtype=bool)
            n_type_cells = int(cell_mask.sum())

            if n_type_cells < self.config.min_cells_per_type:
                log.warning(
                    f"Cell type '{cell_type}' has only {n_type_cells} cells "
                    f"(minimum recommended: {self.config.min_cells_per_type})"
                )

            if sparse.issparse(self.reference_matrix):
                type_expr = self.reference_matrix[:, cell_mask].toarray()
            else:
                type_expr = self.reference_matrix[:, cell_mask]

            # Calculate average expression
            total_counts = type_expr.sum(axis=0, keepdims=True)
            total_counts[total_counts == 0] = 1  # Avoid division by zero

            normalized = type_expr / total_counts
            mean_expr = normalized.mean(axis=1)

            # Add pseudo-count
            mean_expr[mean_expr == 0] = self.config.pseudo_min

            self.phi[:, i] = mean_expr

        # Normalize columns to sum to 1
        col_sums = self.phi.sum(axis=0, keepdims=True)
        col_sums[col_sums == 0] = 1  # Avoid division by zero
        self.phi = self.phi / col_sums

        log.debug(f"Generated phi matrix: {self.phi.shape}")

    def get_cell_state_profile(self) -> Dict[str, np.ndarray]:
        """
        Get cell state-specific expression profiles

        Returns
        -------
        Dict[str, np.ndarray]
            Dictionary mapping state names to expression profiles
        """
        n_genes = self.reference_matrix.shape[0]
        state_profiles = {}

        for state in self.cell_states:
            cell_mask = (self.cell_state_labels == state).to_numpy(dtype=bool)

            if sparse.issparse(self.reference_matrix):
                state_expr = self.reference_matrix[:, cell_mask].toarray()
            else:
                state_expr = self.reference_matrix[:, cell_mask]

            total_counts = state_expr.sum(axis=0, keepdims=True)
            total_counts[total_counts == 0] = 1

            normalized = state_expr / total_counts
            mean_expr = normalized.mean(axis=1)
            mean_expr[mean_expr == 0] = self.config.pseudo_min

            state_profiles[state] = mean_expr

        return state_profiles

    def get_marker_genes(
        self,
        n_markers: int = 100,
        method: str = "fold_change",
        min_fold_change: float = 2.0,
    ) -> Dict[str, List[str]]:
        """
        Identify marker genes for each cell type

        Parameters
        ----------
        n_markers : int
            Number of markers per cell type
        method : str
            Method for marker selection ("fold_change" or "t-test")
        min_fold_change : float
            Minimum fold change for markers

        Returns
        -------
        Dict[str, List[str]]
            Dictionary mapping cell types to marker gene lists
        """
        markers = {}

        for i, cell_type in enumerate(self.cell_types):
            type_expr = self.phi[:, i]
            other_expr = np.delete(self.phi, i, axis=1).mean(axis=1)

            if method == "fold_change":
                fold_change = np.log2(
                    (type_expr + 1e-10) / (other_expr + 1e-10)
                )

                # Filter by minimum fold change
                valid_idx = np.where(fold_change >= np.log2(min_fold_change))[0]

                # Get top indices
                if len(valid_idx) > 0:
                    top_local_idx = np.argsort(fold_change[valid_idx])[-n_markers:]
                    top_idx = valid_idx[top_local_idx]
                else:
                    top_idx = np.argsort(fold_change)[-n_markers:]

            else:
                raise ValueError(f"Unknown method: {method}")

            top_genes = [self.gene_names[idx] for idx in top_idx]
            markers[cell_type] = top_genes[::-1]  # Highest first

        return markers

    def filter_genes(
        self,
        genes_to_keep: List[str],
    ) -> None:
        """
        Filter reference to keep only specified genes

        Parameters
        ----------
        genes_to_keep : List[str]
            List of gene names to retain
        """
        keep_lower = {str(g).lower() for g in genes_to_keep}
        gene_idx = [i for i, g in enumerate(self.gene_names) if str(g).lower() in keep_lower]

        if len(gene_idx) == 0:
            raise ValueError("No genes to keep found in reference")

        self.reference_matrix = self.reference_matrix[gene_idx, :]
        self.phi = self.phi[gene_idx, :]
        self.gene_names = [self.gene_names[i] for i in gene_idx]

        log.info(f"Filtered to {len(gene_idx)} genes")

    def get_cell_type_counts(self) -> pd.Series:
        """Get count of cells per cell type"""
        return self.cell_type_labels.value_counts()

    def summary(self) -> str:
        """Get summary string of reference data"""
        lines = [
            "BayesPrism Reference Summary:",
            f"  Genes: {self.reference_matrix.shape[0]}",
            f"  Cells: {self.reference_matrix.shape[1]}",
            f"  Cell types: {len(self.cell_types)}",
            f"  Cell states: {len(self.cell_states)}",
            "  Cell type counts:",
        ]
        for ct, count in self.get_cell_type_counts().items():
            lines.append(f"    {ct}: {count}")
        return "\n".join(lines)
