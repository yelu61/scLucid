"""
Signature matrix builder for pyDWLS.

A signature matrix encodes the expected per-gene expression of each cell type.
For DWLS it is the design matrix ``S`` in ``b = S theta``, where ``theta`` is
the cell-type proportion vector to be inferred from the bulk sample ``b``.
"""

import logging
from typing import List, Optional, Union

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

_AggMethod = ("mean", "trimmed_mean")


class SignatureBuilder:
    """
    Build a cell-type signature matrix from single-cell reference data.

    The default ``"mean"`` aggregation is a simple per-cell-type column mean.
    ``"trimmed_mean"`` drops the top and bottom ``trim_percent`` fraction of
    values per gene before averaging, which dampens the influence of
    occasional very-high or zero counts within a cell-type group.

    Parameters
    ----------
    trim_percent : float, default=0.1
        Fraction (between 0 and 0.5) to trim from each tail of the
        per-gene-per-cell-type expression distribution before averaging.
        Only used when ``method="trimmed_mean"``.
    min_cells : int, default=10
        Minimum number of cells required per cell type. Cell types with fewer
        cells in the reference are dropped from the signature matrix.

    Examples:
    --------
    >>> builder = SignatureBuilder(trim_percent=0.1, min_cells=20)
    >>> signature = builder.build(sc_data, cell_type_labels)
    >>> print(signature.shape)  # (n_genes, n_kept_cell_types)
    """

    def __init__(self, trim_percent: float = 0.1, min_cells: int = 10):
        if not 0.0 <= trim_percent < 0.5:
            raise ValueError(
                f"trim_percent must be in [0, 0.5); got {trim_percent}"
            )
        if min_cells < 1:
            raise ValueError(f"min_cells must be >= 1; got {min_cells}")

        self.trim_percent = float(trim_percent)
        self.min_cells = int(min_cells)

    def build(
        self,
        sc_data: pd.DataFrame,
        cell_type_labels: Union[pd.Series, np.ndarray, List],
        genes_to_use: Optional[List[str]] = None,
        method: str = "mean",
        min_cells: Optional[int] = None,
    ) -> pd.DataFrame:
        """
        Aggregate single-cell expression to a cell-type signature matrix.

        Parameters
        ----------
        sc_data : pd.DataFrame
            Expression matrix (genes x cells).
        cell_type_labels : array-like
            Cell type label per cell, length must equal ``sc_data.shape[1]``.
        genes_to_use : list of str, optional
            Restrict the output rows to these genes. If ``None``, all genes
            in ``sc_data`` are kept.
        method : {"mean", "trimmed_mean"}, default="mean"
            Aggregation strategy.
        min_cells : int, optional
            Override the constructor-level ``min_cells`` threshold for this
            build call.

        Returns:
        -------
        pd.DataFrame
            Signature matrix indexed by gene name with one column per kept
            cell type.

        Raises:
        ------
        ValueError
            If ``method`` is unknown or no cell type passes ``min_cells``.

        Examples:
        --------
        >>> sig = builder.build(sc_data, labels, method="trimmed_mean")
        """
        if method not in _AggMethod:
            raise ValueError(
                f"method must be one of {_AggMethod}; got {method!r}"
            )

        effective_min = self.min_cells if min_cells is None else int(min_cells)
        labels = pd.Series(cell_type_labels, index=sc_data.columns)

        gene_index = sc_data.index
        if genes_to_use is not None:
            gene_index = sc_data.index.intersection(pd.Index(genes_to_use))
            if len(gene_index) == 0:
                raise ValueError("None of the requested genes are present in sc_data")
            sc_data = sc_data.loc[gene_index]

        columns: List[str] = []
        signature_values: List[np.ndarray] = []
        dropped: List[str] = []

        for cell_type in labels.unique():
            mask = labels == cell_type
            n_cells = int(mask.sum())
            if n_cells < effective_min:
                dropped.append(f"{cell_type}({n_cells})")
                continue

            subset = sc_data.loc[:, mask.values]
            if method == "mean":
                values = subset.mean(axis=1).to_numpy()
            else:
                values = self._trimmed_mean_per_gene(subset.to_numpy())

            columns.append(str(cell_type))
            signature_values.append(values)

        if not columns:
            log.warning(
                "No cell type met min_cells=%d (available counts: %s); "
                "returning empty signature matrix.",
                effective_min,
                labels.value_counts().to_dict(),
            )
            return pd.DataFrame(index=gene_index)

        if dropped:
            log.info(
                "SignatureBuilder: dropped %d cell types below min_cells=%d: %s",
                len(dropped),
                effective_min,
                ", ".join(dropped),
            )

        signature = pd.DataFrame(
            np.column_stack(signature_values),
            index=gene_index,
            columns=columns,
        )

        log.info(
            "Built signature matrix: %d genes x %d cell types (method=%s)",
            signature.shape[0],
            signature.shape[1],
            method,
        )
        return signature

    def _trimmed_mean_per_gene(self, matrix: np.ndarray) -> np.ndarray:
        """Per-gene trimmed mean across cells."""
        if matrix.shape[1] == 0:
            return np.zeros(matrix.shape[0])

        sorted_matrix = np.sort(matrix, axis=1)
        n_cells = matrix.shape[1]
        k = int(np.floor(n_cells * self.trim_percent))
        if k > 0:
            sorted_matrix = sorted_matrix[:, k:-k] if (n_cells - 2 * k) > 0 else sorted_matrix
        return sorted_matrix.mean(axis=1)
