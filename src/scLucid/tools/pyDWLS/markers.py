"""
Marker-gene selection for pyDWLS.

Choosing a compact set of marker genes per cell type is the most effective
way to improve DWLS deconvolution: it concentrates signal on genes that
actually discriminate between reference cell types and reduces noise from
ubiquitously expressed transcripts.
"""

import logging
from typing import List, Set, Union

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

_SelectionMethod = ("ratio", "difference", "fold_change")


class MarkerSelector:
    """
    Select marker genes per cell type from single-cell reference data.

    Three ranking criteria are available:

    - ``"ratio"`` — ``mean_in_type / (mean_in_others + eps)``. Sensitive to
      genes that are highly specific to a single cell type.
    - ``"difference"`` — ``mean_in_type - mean_in_others``. Favors highly
      expressed genes; less sensitive to low-magnitude specificity.
    - ``"fold_change"`` — ``log2((mean_in_type + 1) / (mean_in_others + 1))``.
      Balances specificity and magnitude; closest to standard scRNA-seq
      marker conventions.

    Examples:
    --------
    >>> selector = MarkerSelector()
    >>> markers = selector.select(sc_data, labels, n_markers=20, method="fold_change")
    >>> print(f"Selected {len(markers)} unique markers")
    """

    def __init__(self):
        pass

    def select(
        self,
        sc_data: pd.DataFrame,
        cell_type_labels: Union[pd.Series, np.ndarray, list],
        n_markers: int = 50,
        method: str = "ratio",
        log_transform: bool = True,
    ) -> List[str]:
        """
        Return up to ``n_markers`` top genes per cell type, deduplicated.

        Parameters
        ----------
        sc_data : pd.DataFrame
            Expression matrix (genes x cells).
        cell_type_labels : array-like
            Cell type label per cell.
        n_markers : int, default=50
            Maximum number of markers to take per cell type before dedup.
        method : {"ratio", "difference", "fold_change"}, default="ratio"
            Ranking criterion.
        log_transform : bool, default=True
            Apply ``log1p`` to ``sc_data`` before ranking. Recommended when
            ``sc_data`` contains raw counts; disable if values are already
            log-normalized.

        Returns:
        -------
        list of str
            Deduplicated marker gene names ordered by first appearance in the
            per-cell-type ranking.

        Raises:
        ------
        ValueError
            If ``method`` is not recognized or fewer than 2 cell types exist.
        """
        if method not in _SelectionMethod:
            raise ValueError(
                f"method must be one of {_SelectionMethod}; got {method!r}"
            )
        if n_markers < 1:
            raise ValueError(f"n_markers must be >= 1; got {n_markers}")

        labels = pd.Series(cell_type_labels, index=sc_data.columns)
        cell_types = labels.unique()
        if len(cell_types) < 2:
            raise ValueError(
                f"Need at least 2 cell types for marker selection; got {len(cell_types)}"
            )

        expression = np.log1p(sc_data) if log_transform else sc_data

        per_type_means = pd.DataFrame(
            {
                str(ct): expression.loc[:, labels == ct].mean(axis=1)
                for ct in cell_types
            }
        )

        ordered_markers: List[str] = []
        seen: Set[str] = set()
        eps = 1e-6

        for cell_type in cell_types:
            ct_name = str(cell_type)
            others_mean = per_type_means.drop(columns=ct_name).mean(axis=1)
            in_type = per_type_means[ct_name]

            if method == "ratio":
                score = in_type / (others_mean + eps)
            elif method == "difference":
                score = in_type - others_mean
            else:  # fold_change
                score = np.log2((in_type + 1.0) / (others_mean + 1.0))

            top_genes = score.sort_values(ascending=False).head(n_markers).index
            for gene in top_genes:
                gene_str = str(gene)
                if gene_str not in seen:
                    ordered_markers.append(gene_str)
                    seen.add(gene_str)

        log.info(
            "MarkerSelector: %d unique markers across %d cell types (method=%s)",
            len(ordered_markers),
            len(cell_types),
            method,
        )
        return ordered_markers
