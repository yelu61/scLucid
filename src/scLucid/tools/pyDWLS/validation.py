"""
Cross-validation utilities for pyDWLS signature quality assessment.

Holds out a fold of cells, builds a pseudo-bulk and known cell-type
proportions from them, trains the signature matrix on the remaining cells,
deconvolves the pseudo-bulk, and scores the recovered proportions against
the known ground truth.
"""

import logging
from typing import Any, Dict, List, Union

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.model_selection import KFold

from .markers import MarkerSelector
from .signature import SignatureBuilder
from .solver import DampenedWLS
from .utils import create_pseudo_bulk

log = logging.getLogger(__name__)


class CrossValidator:
    """
    K-fold cross-validation for DWLS signature matrices.

    Each fold holds out a subset of single cells, constructs a pseudo-bulk
    sample with known proportions from them, fits the signature matrix on
    the remaining cells, deconvolves the pseudo-bulk, and reports
    correlation and RMSE between recovered and true proportions.

    Parameters
    ----------
    n_folds : int, default=5
        Number of cross-validation folds.
    n_cells_per_bulk : int, default=100
        Number of held-out cells to sample per cell type when building the
        pseudo-bulk for each fold.
    random_state : int, default=42
        Seed for the fold splitter and pseudo-bulk sampling.

    Examples:
    --------
    >>> cv = CrossValidator(n_folds=5)
    >>> report = cv.run(sc_data, cell_type_labels)
    >>> print(report["mean_correlation"], report["mean_rmse"])
    """

    def __init__(
        self,
        n_folds: int = 5,
        n_cells_per_bulk: int = 100,
        random_state: int = 42,
    ):
        if n_folds < 2:
            raise ValueError(f"n_folds must be >= 2; got {n_folds}")
        if n_cells_per_bulk < 1:
            raise ValueError(f"n_cells_per_bulk must be >= 1; got {n_cells_per_bulk}")

        self.n_folds = int(n_folds)
        self.n_cells_per_bulk = int(n_cells_per_bulk)
        self.random_state = int(random_state)

    def run(
        self,
        sc_data: pd.DataFrame,
        cell_type_labels: Union[pd.Series, np.ndarray, list],
        verbose: bool = False,
    ) -> Dict[str, Any]:
        """
        Run the K-fold cross-validation loop.

        Parameters
        ----------
        sc_data : pd.DataFrame
            Single-cell expression matrix (genes x cells).
        cell_type_labels : array-like
            Cell type label per cell, aligned with ``sc_data.columns``.
        verbose : bool, default=False
            Emit per-fold progress logs at INFO level.

        Returns:
        -------
        dict
            Dictionary with keys:

            - ``"mean_correlation"`` (float): mean Pearson correlation across folds.
            - ``"mean_rmse"`` (float): mean RMSE across folds.
            - ``"fold_correlations"`` (list of float): per-fold correlations.
            - ``"fold_rmse"`` (list of float): per-fold RMSE values.
            - ``"n_folds_evaluated"`` (int): number of folds that produced a
              valid evaluation (folds with degenerate splits are skipped).
        """
        labels = pd.Series(cell_type_labels, index=sc_data.columns)
        cell_index = np.arange(sc_data.shape[1])

        splitter = KFold(
            n_splits=self.n_folds, shuffle=True, random_state=self.random_state
        )

        fold_correlations: List[float] = []
        fold_rmse: List[float] = []
        n_evaluated = 0

        for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(cell_index)):
            train_cols = sc_data.columns[train_idx]
            test_cols = sc_data.columns[test_idx]

            train_labels = labels.loc[train_cols]
            test_labels = labels.loc[test_cols]

            shared_types = sorted(
                set(train_labels.unique()) & set(test_labels.unique()),
                key=str,
            )
            if len(shared_types) < 2:
                log.warning(
                    "Fold %d skipped: only %d shared cell types between train/test",
                    fold_idx,
                    len(shared_types),
                )
                continue

            train_filter = train_labels.isin(shared_types)
            test_filter = test_labels.isin(shared_types)
            train_cols_kept = train_cols[train_filter.values]
            test_cols_kept = test_cols[test_filter.values]

            train_data = sc_data.loc[:, train_cols_kept]
            train_lbls = labels.loc[train_cols_kept]
            test_data = sc_data.loc[:, test_cols_kept]
            test_lbls = labels.loc[test_cols_kept]

            try:
                pseudo_bulk, true_props = create_pseudo_bulk(
                    test_data,
                    test_lbls,
                    n_cells=self.n_cells_per_bulk,
                    random_state=self.random_state + fold_idx,
                )
            except ValueError as exc:
                log.warning("Fold %d skipped during pseudo-bulk build: %s", fold_idx, exc)
                continue

            marker_selector = MarkerSelector()
            try:
                markers = marker_selector.select(train_data, train_lbls, n_markers=50)
            except ValueError as exc:
                log.warning("Fold %d skipped during marker selection: %s", fold_idx, exc)
                continue

            signature_builder = SignatureBuilder(min_cells=1)
            signature = signature_builder.build(
                train_data,
                train_lbls,
                genes_to_use=markers,
            )

            common_genes = signature.index.intersection(pseudo_bulk.index)
            sig_aligned = signature.loc[common_genes]
            bulk_aligned = pseudo_bulk.loc[common_genes].to_numpy(dtype=float)

            solver = DampenedWLS(dampen_factor=1.0)
            estimated = solver.solve(sig_aligned.to_numpy(dtype=float), bulk_aligned)
            estimated_series = pd.Series(estimated, index=sig_aligned.columns)

            comparable = true_props.index.intersection(estimated_series.index)
            if len(comparable) < 2:
                log.warning(
                    "Fold %d skipped: only %d cell types comparable",
                    fold_idx,
                    len(comparable),
                )
                continue

            true_vec = true_props.loc[comparable].to_numpy(dtype=float)
            est_vec = estimated_series.loc[comparable].to_numpy(dtype=float)

            true_norm = true_vec / max(true_vec.sum(), 1e-12)
            est_norm = est_vec / max(est_vec.sum(), 1e-12)

            if np.std(true_norm) < 1e-9 or np.std(est_norm) < 1e-9:
                corr = 0.0
            else:
                corr, _ = pearsonr(true_norm, est_norm)
                if np.isnan(corr):
                    corr = 0.0
            rmse = float(np.sqrt(np.mean((true_norm - est_norm) ** 2)))

            fold_correlations.append(float(corr))
            fold_rmse.append(rmse)
            n_evaluated += 1

            if verbose:
                log.info(
                    "Fold %d: correlation=%.3f, rmse=%.3f", fold_idx, corr, rmse
                )

        # Always pad lengths so ``len(fold_correlations) == n_folds`` per the
        # public contract. Skipped folds are recorded as NaN so downstream
        # mean/agg ignores them gracefully.
        while len(fold_correlations) < self.n_folds:
            fold_correlations.append(float("nan"))
            fold_rmse.append(float("nan"))

        valid_corr = [c for c in fold_correlations if not np.isnan(c)]
        valid_rmse = [r for r in fold_rmse if not np.isnan(r)]

        return {
            "mean_correlation": float(np.mean(valid_corr)) if valid_corr else float("nan"),
            "mean_rmse": float(np.mean(valid_rmse)) if valid_rmse else float("nan"),
            "fold_correlations": fold_correlations,
            "fold_rmse": fold_rmse,
            "n_folds_evaluated": n_evaluated,
        }
