"""
Core DWLS class for bulk RNA-seq deconvolution.

This module implements the main DWLS class that orchestrates the
deconvolution workflow, including signature building, gene selection,
and proportion estimation.
"""

import logging
from typing import List, Optional, Union

import numpy as np
import pandas as pd

from .markers import MarkerSelector
from .signature import SignatureBuilder
from .solver import DampenedWLS
from .utils import align_data

log = logging.getLogger(__name__)


class DWLS:
    """
    Dampened Weighted Least Squares (DWLS) for bulk RNA-seq deconvolution.

    DWLS estimates cell type proportions in bulk RNA-seq samples by solving
    a weighted least squares problem with dampening to handle highly
    expressed genes.

    Parameters
    ----------
    signature_matrix : pd.DataFrame, optional
        Gene expression signature matrix (genes x cell types).
    bulk_data : pd.DataFrame, optional
        Bulk RNA-seq data (genes x samples).
    dampen_factor : float, default=1.0
        Factor for dampening highly expressed genes.
    use_nonneg : bool, default=True
        Enforce non-negative proportion constraints.

    Attributes:
    ----------
    signature_matrix_ : pd.DataFrame
        The signature matrix used for deconvolution.
    results_ : pd.DataFrame
        Deconvolution results (samples x cell types).

    Examples:
    --------
    >>> dwls = DWLS()
    >>> signature = dwls.build_signature_matrix(sc_data, cell_labels)
    >>> proportions = dwls.deconvolve(bulk_data)
    >>> print(proportions.head())
    """

    def __init__(
        self,
        signature_matrix: Optional[pd.DataFrame] = None,
        bulk_data: Optional[pd.DataFrame] = None,
        dampen_factor: float = 1.0,
        use_nonneg: bool = True,
    ):
        self.signature_matrix = signature_matrix
        self.bulk_data = bulk_data
        self.dampen_factor = dampen_factor
        self.use_nonneg = use_nonneg
        self.results_: Optional[pd.DataFrame] = None

        # Component objects
        self._signature_builder = SignatureBuilder()
        self._marker_selector = MarkerSelector()
        self._solver = DampenedWLS(
            dampen_factor=dampen_factor,
            use_nonneg=use_nonneg,
        )

    def build_signature_matrix(
        self,
        sc_data: pd.DataFrame,
        cell_type_labels: Union[pd.Series, np.ndarray, List],
        genes_to_use: Optional[List[str]] = None,
        method: str = "mean",
        min_cells: int = 10,
    ) -> pd.DataFrame:
        """
        Build signature matrix from single-cell reference data.

        Parameters
        ----------
        sc_data : pd.DataFrame
            Single-cell expression matrix (genes x cells).
        cell_type_labels : array-like
            Cell type label for each cell.
        genes_to_use : list, optional
            Specific genes to include. If None, uses all genes.
        method : str, default="mean"
            Aggregation method ("mean" or "trimmed_mean").
        min_cells : int, default=10
            Minimum cells required per cell type.

        Returns:
        -------
        pd.DataFrame
            Signature matrix (genes x cell types).

        Examples:
        --------
        >>> signature = dwls.build_signature_matrix(
        ...     sc_data, cell_labels, method="trimmed_mean"
        ... )
        """
        self.signature_matrix = self._signature_builder.build(
            sc_data=sc_data,
            cell_type_labels=cell_type_labels,
            genes_to_use=genes_to_use,
            method=method,
            min_cells=min_cells,
        )
        return self.signature_matrix

    def select_marker_genes(
        self,
        sc_data: pd.DataFrame,
        cell_type_labels: Union[pd.Series, np.ndarray, List],
        n_markers: int = 50,
        method: str = "ratio",
        log_transform: bool = True,
    ) -> List[str]:
        """
        Select marker genes for each cell type.

        Parameters
        ----------
        sc_data : pd.DataFrame
            Single-cell expression data (genes x cells).
        cell_type_labels : array-like
            Cell type labels.
        n_markers : int, default=50
            Number of marker genes per cell type.
        method : str, default="ratio"
            Selection method ("ratio", "difference", or "fold_change").
        log_transform : bool, default=True
            Log-transform data before selection.

        Returns:
        -------
        list
            Selected marker gene names.
        """
        return self._marker_selector.select(
            sc_data=sc_data,
            cell_type_labels=cell_type_labels,
            n_markers=n_markers,
            method=method,
            log_transform=log_transform,
        )

    def deconvolve(
        self,
        bulk_data: Optional[pd.DataFrame] = None,
        verbose: bool = True,
    ) -> pd.DataFrame:
        """
        Deconvolve bulk RNA-seq data into cell type proportions.

        Parameters
        ----------
        bulk_data : pd.DataFrame, optional
            Bulk expression matrix (genes x samples). Uses self.bulk_data
            if not provided.
        verbose : bool, default=True
            Print progress messages.

        Returns:
        -------
        pd.DataFrame
            Cell type proportions (samples x cell types).

        Raises:
        ------
        ValueError
            If signature matrix or bulk data is not provided.

        Examples:
        --------
        >>> proportions = dwls.deconvolve(bulk_data)
        >>> print(proportions.head())
        """
        if bulk_data is None:
            bulk_data = self.bulk_data

        if bulk_data is None:
            raise ValueError("Bulk data must be provided")

        if self.signature_matrix is None:
            raise ValueError("Signature matrix must be built first")

        # Align data to common genes
        sig_aligned, bulk_aligned = align_data(
            self.signature_matrix,
            bulk_data,
        )

        if verbose:
            log.info(f"Deconvolving {bulk_aligned.shape[1]} samples...")

        # Deconvolve each sample
        results_list = []
        n_samples = bulk_aligned.shape[1]

        for i, sample_id in enumerate(bulk_aligned.columns):
            if verbose and (i + 1) % 10 == 0:
                log.info(f"  Processed {i + 1}/{n_samples} samples")

            bulk_sample = bulk_aligned[sample_id].values
            proportions = self._solver.solve(
                sig_aligned.values,
                bulk_sample,
            )

            result_series = pd.Series(
                proportions,
                index=sig_aligned.columns,
                name=sample_id,
            )
            results_list.append(result_series)

        # Combine results
        self.results_ = pd.concat(results_list, axis=1).T

        if verbose:
            log.info("Deconvolution complete!")

        return self.results_

    def solve_single(
        self,
        bulk_sample: Union[pd.Series, np.ndarray],
    ) -> pd.Series:
        """
        Solve for a single bulk sample.

        Parameters
        ----------
        bulk_sample : pd.Series or np.ndarray
            Expression values for one sample.

        Returns:
        -------
        pd.Series
            Cell type proportions.
        """
        if self.signature_matrix is None:
            raise ValueError("Signature matrix must be built first")

        # Align genes
        if isinstance(bulk_sample, pd.Series):
            common_genes = self.signature_matrix.index.intersection(bulk_sample.index)
            sig = self.signature_matrix.loc[common_genes]
            bulk = bulk_sample.loc[common_genes].values
        else:
            sig = self.signature_matrix
            bulk = bulk_sample

        # Solve
        proportions = self._solver.solve(sig.values, bulk)

        return pd.Series(proportions, index=sig.columns)

    def get_results(self) -> Optional[pd.DataFrame]:
        """Get deconvolution results."""
        return self.results_

    def summary(self) -> str:
        """Get summary of DWLS object."""
        lines = ["DWLS Summary:"]

        if self.signature_matrix is not None:
            lines.append(
                f"  Signature: {self.signature_matrix.shape[0]} genes "
                f"x {self.signature_matrix.shape[1]} cell types"
            )
        else:
            lines.append("  Signature: Not built")

        if self.bulk_data is not None:
            lines.append(
                f"  Bulk data: {self.bulk_data.shape[0]} genes "
                f"x {self.bulk_data.shape[1]} samples"
            )

        if self.results_ is not None:
            lines.append(f"  Results: {self.results_.shape}")

        return "\n".join(lines)
