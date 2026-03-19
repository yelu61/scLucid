"""
Core BayesPrism class for deconvolution (R-free)

Main implementation of Bayesian deconvolution algorithm.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, List, Tuple
from scipy.optimize import nnls
from concurrent.futures import ProcessPoolExecutor
import logging

from .config import PrismConfig, DeconvolutionConfig
from .reference import BayesPrismReference
from .sampling import GibbsSampler
from .utils import cleanup_genes

log = logging.getLogger(__name__)


class BayesPrism:
    """
    Main BayesPrism class for deconvolution

    Implements the complete Bayesian deconvolution workflow:
    1. NNLS initialization for theta
    2. Gibbs sampling for posterior inference
    3. Cell type-specific expression estimation

    Parameters
    ----------
    reference : BayesPrismReference
        Reference data object
    mixture : pd.DataFrame
        Bulk RNA-seq expression (genes x samples)
    config : PrismConfig, optional
        Configuration parameters

    Attributes
    ----------
    theta_initial_ : np.ndarray
        Initial cell type proportions from NNLS
    theta_updated_ : np.ndarray
        Posterior mean proportions from Gibbs sampling
    Z_ : np.ndarray
        Cell type-specific expression (genes x cell_types x samples)
    aligned_genes_ : List[str]
        Genes common to reference and mixture

    Examples
    --------
    >>> ref = BayesPrismReference(ref_data, cell_types)
    >>> bp = BayesPrism(reference=ref, mixture=bulk_data)
    >>> bp.cleanup_genes(remove_ribo=True, remove_mito=True)
    >>> bp.run_deconvolution(n_cores=4)
    >>> fractions = bp.get_fraction()
    """

    def __init__(
        self,
        reference: BayesPrismReference,
        mixture: pd.DataFrame,
        config: Optional[PrismConfig] = None,
    ):
        self.reference = reference
        self.mixture = mixture
        self.config = config if config is not None else PrismConfig()
        self.deconv_config = DeconvolutionConfig()

        # Pydantic configs validate automatically on instantiation

        # Initialize results
        self.theta_initial_: Optional[np.ndarray] = None
        self.theta_updated_: Optional[np.ndarray] = None
        self.Z_: Optional[np.ndarray] = None
        self.aligned_genes_: List[str] = []

        # Align genes
        self._align_genes()

    def _align_genes(self) -> None:
        """Align reference and mixture data by common genes"""
        # Use case-insensitive matching for robust cross-source gene alignment.
        ref_lookup = {str(g).lower(): str(g) for g in self.reference.gene_names}
        mix_lookup = {str(g).lower(): g for g in self.mixture.index}
        common_lower = sorted(set(ref_lookup).intersection(mix_lookup))

        if len(common_lower) == 0:
            raise ValueError("No common genes between reference and mixture!")

        common_ref_genes = [ref_lookup[g] for g in common_lower]
        common_mix_genes = [mix_lookup[g] for g in common_lower]
        log.info(f"Aligned to {len(common_ref_genes)} common genes")

        # Subset reference
        self.reference.filter_genes(common_ref_genes)

        # Subset mixture and normalize index to reference naming.
        self.mixture = self.mixture.loc[common_mix_genes, :].copy()
        self.mixture.index = common_ref_genes
        self.aligned_genes_ = common_ref_genes

    def cleanup_genes(
        self,
        remove_ribo: bool = True,
        remove_mito: bool = True,
        remove_sex: bool = True,
        min_expression: float = 0.0,
    ) -> None:
        """
        Clean up gene list

        Parameters
        ----------
        remove_ribo : bool
            Remove ribosomal genes (RPS*, RPL*)
        remove_mito : bool
            Remove mitochondrial genes (MT-*)
        remove_sex : bool
            Remove sex chromosome genes (XIST, Y*)
        min_expression : float
            Minimum expression threshold
        """
        genes_to_keep = cleanup_genes(
            self.aligned_genes_,
            remove_ribo=remove_ribo,
            remove_mito=remove_mito,
            remove_sex=remove_sex,
        )

        if min_expression > 0:
            expr_means = self.mixture.mean(axis=1)
            genes_to_keep = [g for g in genes_to_keep if expr_means[g] >= min_expression]

        if len(genes_to_keep) == 0:
            raise ValueError("No genes remaining after cleanup!")

        self._update_gene_subset(genes_to_keep)
        log.info(f"Cleaned to {len(genes_to_keep)} genes")

    def _update_gene_subset(self, genes: List[str]) -> None:
        """Update data to gene subset"""
        self.reference.filter_genes(genes)
        self.mixture = self.mixture.loc[genes, :]
        self.aligned_genes_ = genes

    def select_markers(
        self,
        n_markers: int = 500,
        method: str = "fold_change",
    ) -> List[str]:
        """
        Select marker genes for deconvolution

        Parameters
        ----------
        n_markers : int
            Number of markers per cell type
        method : str
            Selection method ("fold_change" or "t-test")

        Returns
        -------
        List[str]
            Selected marker genes
        """
        markers_dict = self.reference.get_marker_genes(
            n_markers=n_markers,
            method=method,
        )

        all_markers = set()
        for genes in markers_dict.values():
            all_markers.update(genes)

        marker_list = sorted(list(all_markers))
        log.info(f"Selected {len(marker_list)} marker genes")

        return marker_list

    def run_deconvolution(
        self,
        n_cores: int = 1,
        verbose: bool = True,
        use_numba: bool = True,
    ) -> None:
        """
        Run deconvolution

        Parameters
        ----------
        n_cores : int
            Number of CPU cores for parallel processing
        verbose : bool
            Whether to display progress
        use_numba : bool
            Use Numba JIT for faster sampling
        """
        n_samples = self.mixture.shape[1]
        n_cell_types = len(self.reference.cell_types)
        n_genes = len(self.aligned_genes_)

        if verbose:
            log.info(f"Starting deconvolution for {n_samples} samples...")

        # Initialize result arrays
        self.theta_initial_ = np.zeros((n_cell_types, n_samples))
        self.theta_updated_ = np.zeros((n_cell_types, n_samples))
        self.Z_ = np.zeros((n_genes, n_cell_types, n_samples))

        # Prepare arguments for parallel processing
        mixture_array = self.mixture.values
        sample_args = [
            (mixture_array[:, i], self.reference.phi, self.config, use_numba)
            for i in range(n_samples)
        ]

        # Run deconvolution
        if n_cores > 1:
            with ProcessPoolExecutor(max_workers=n_cores) as executor:
                results = list(executor.map(_deconvolve_sample_worker, sample_args))
        else:
            results = [_deconvolve_sample_worker(args) for args in sample_args]

        # Collect results
        for i, (theta0, theta_f, Z_i) in enumerate(results):
            self.theta_initial_[:, i] = theta0
            self.theta_updated_[:, i] = theta_f
            self.Z_[:, :, i] = Z_i

        if verbose:
            log.info("Deconvolution complete!")

    def get_fraction(self, updated: bool = True) -> pd.DataFrame:
        """
        Get cell type proportions

        Parameters
        ----------
        updated : bool
            Use updated (Gibbs) or initial (NNLS) estimates

        Returns
        -------
        pd.DataFrame
            Cell type proportions (samples x cell_types)
        """
        if updated:
            if self.theta_updated_ is None:
                raise ValueError("Run deconvolution first")
            theta = self.theta_updated_
        else:
            if self.theta_initial_ is None:
                raise ValueError("Run deconvolution first")
            theta = self.theta_initial_

        return pd.DataFrame(
            theta.T,
            index=self.mixture.columns,
            columns=self.reference.cell_types,
        )

    def get_expression(self, cell_type: Optional[str] = None) -> pd.DataFrame:
        """
        Get cell type-specific expression

        Parameters
        ----------
        cell_type : str, optional
            Specific cell type (None for all)

        Returns
        -------
        pd.DataFrame or dict
            Cell type-specific expression
        """
        if self.Z_ is None:
            raise ValueError("Run deconvolution first")

        if cell_type is None:
            result = {}
            for i, ct in enumerate(self.reference.cell_types):
                result[ct] = pd.DataFrame(
                    self.Z_[:, i, :],
                    index=self.aligned_genes_,
                    columns=self.mixture.columns,
                )
            return result
        else:
            ct_idx = self.reference.cell_types.index(cell_type)
            return pd.DataFrame(
                self.Z_[:, ct_idx, :],
                index=self.aligned_genes_,
                columns=self.mixture.columns,
            )

    def compute_cv(self) -> pd.DataFrame:
        """Compute coefficient of variation for proportions"""
        if self.theta_updated_ is None:
            raise ValueError("Run deconvolution first")

        cv = np.std(self.theta_updated_, axis=1) / (
            np.mean(self.theta_updated_, axis=1) + 1e-10
        )

        return pd.DataFrame({
            'cell_type': self.reference.cell_types,
            'CV': cv,
        })


def _deconvolve_sample_worker(args) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Worker function for parallel deconvolution"""
    mixture_vector, phi, config, use_numba = args

    # Initial NNLS estimate
    theta0 = _initial_estimate_nnls(mixture_vector, phi)

    # Gibbs sampling
    sampler = GibbsSampler(
        n_iter=config.n_iter,
        burnin=config.burnin,
        use_numba=use_numba,
    )

    theta_samples, Z_samples = sampler.sample(
        mixture=mixture_vector,
        theta_init=theta0,
        phi=phi,
        verbose=False,
    )

    theta_f = theta_samples.mean(axis=0)
    Z_mean = Z_samples.mean(axis=0)

    return theta0, theta_f, Z_mean


def _initial_estimate_nnls(mixture: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """NNLS initialization for theta"""
    mixture_norm = mixture / (mixture.sum() + 1e-10)
    theta, _ = nnls(phi, mixture_norm)
    theta = theta / (theta.sum() + 1e-10)
    return theta
