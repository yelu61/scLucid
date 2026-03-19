"""
Embedding learning module for BayesPrism (R-free)

NMF-based gene program learning for tumor expression deconvolution.
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple, List
from sklearn.decomposition import NMF
import logging

from .core import BayesPrism

log = logging.getLogger(__name__)


class BayesPrismEmbedding:
    """
    BayesPrism embedding learning module

    Uses NMF to learn malignant gene programs from tumor-specific expression.
    This is useful for identifying cancer-specific gene expression patterns.

    Parameters
    ----------
    prism : BayesPrism
        BayesPrism object with completed deconvolution
    tumor_key : str
        Cell type identifier for tumor/malignant cells
    n_programs : int
        Number of gene programs to learn

    Attributes
    ----------
    W_ : np.ndarray
        Gene program matrix (genes x programs)
    H_ : np.ndarray
        Program usage matrix (programs x samples)
    W_df_ : pd.DataFrame
        Gene program as DataFrame with gene names
    H_df_ : pd.DataFrame
        Program usage as DataFrame with sample names

    Examples
    --------
    >>> embedding = BayesPrismEmbedding(
    ...     prism=bp,
    ...     tumor_key='Tumor',
    ...     n_programs=5
    ... )
    >>> embedding.run_nmf(max_iter=200)
    >>> programs = embedding.get_gene_programs()
    >>> usage = embedding.get_program_usage()
    """

    def __init__(
        self,
        prism: BayesPrism,
        tumor_key: str,
        n_programs: int = 5,
    ):
        self.prism = prism
        self.tumor_key = tumor_key
        self.n_programs = n_programs

        # Validate tumor key
        if tumor_key not in prism.reference.cell_types:
            raise ValueError(
                f"tumor_key '{tumor_key}' not in cell types: {prism.reference.cell_types}"
            )

        # Extract tumor expression
        self.tumor_expression_ = prism.get_expression(tumor_key)

        # Initialize results
        self.W_: Optional[np.ndarray] = None
        self.H_: Optional[np.ndarray] = None
        self.W_df_: Optional[pd.DataFrame] = None
        self.H_df_: Optional[pd.DataFrame] = None
        self.model_: Optional[NMF] = None
        self.reconstruction_err_: Optional[float] = None

    def run_nmf(
        self,
        max_iter: int = 200,
        tol: float = 1e-4,
        init: str = "random",
        solver: str = "cd",
        verbose: bool = True,
    ) -> None:
        """
        Run NMF to learn gene programs

        Parameters
        ----------
        max_iter : int
            Maximum iterations
        tol : float
            Convergence tolerance
        init : str
            Initialization method
        solver : str
            NMF solver ("cd" for coordinate descent, "mu" for multiplicative update)
        verbose : bool
            Whether to print progress
        """
        if verbose:
            log.info(f"Learning {self.n_programs} gene programs via NMF...")

        # Run NMF
        self.model_ = NMF(
            n_components=self.n_programs,
            init=init,
            max_iter=max_iter,
            tol=tol,
            solver=solver,
            random_state=42,
            verbose=verbose,
        )

        # Fit to tumor expression
        self.W_ = self.model_.fit_transform(self.tumor_expression_.values)
        self.H_ = self.model_.components_
        self.reconstruction_err_ = self.model_.reconstruction_err_

        # Convert to DataFrames
        program_names = [f"Program_{i+1}" for i in range(self.n_programs)]

        self.W_df_ = pd.DataFrame(
            self.W_,
            index=self.tumor_expression_.index,
            columns=program_names,
        )

        self.H_df_ = pd.DataFrame(
            self.H_,
            index=program_names,
            columns=self.tumor_expression_.columns,
        )

        if verbose:
            log.info(
                f"NMF complete! Reconstruction error: {self.reconstruction_err_:.4f}"
            )

    def get_gene_programs(self, top_n: Optional[int] = None) -> pd.DataFrame:
        """
        Get gene program matrix

        Parameters
        ----------
        top_n : int, optional
            If provided, return only top N genes per program

        Returns
        -------
        pd.DataFrame
            Gene program matrix
        """
        if self.W_df_ is None:
            raise ValueError("Run run_nmf() first")

        if top_n is not None:
            # Return top genes for each program
            result = {}
            for col in self.W_df_.columns:
                top_genes = self.W_df_[col].nlargest(top_n)
                result[col] = top_genes
            return pd.DataFrame(result)

        return self.W_df_

    def get_program_usage(self) -> pd.DataFrame:
        """
        Get program usage coefficients

        Returns
        -------
        pd.DataFrame
            Program usage across samples
        """
        if self.H_df_ is None:
            raise ValueError("Run run_nmf() first")

        return self.H_df_

    def get_top_genes(
        self,
        program: str,
        n: int = 50,
    ) -> pd.Series:
        """
        Get top genes for a specific program

        Parameters
        ----------
        program : str
            Program name (e.g., "Program_1")
        n : int
            Number of top genes to return

        Returns
        -------
        pd.Series
            Top genes with weights
        """
        if self.W_df_ is None:
            raise ValueError("Run run_nmf() first")

        if program not in self.W_df_.columns:
            raise ValueError(f"Program '{program}' not found")

        return self.W_df_[program].nlargest(n)

    def score_samples(
        self,
        program: str,
        gene_expression: pd.DataFrame,
    ) -> pd.Series:
        """
        Score samples based on a gene program

        Parameters
        ----------
        program : str
            Program name
        gene_expression : pd.DataFrame
            Gene expression data (genes x samples)

        Returns
        -------
        pd.Series
            Program scores for each sample
        """
        if self.W_df_ is None:
            raise ValueError("Run run_nmf() first")

        program_genes = self.W_df_[program]

        # Get common genes
        common_genes = program_genes.index.intersection(gene_expression.index)

        if len(common_genes) == 0:
            raise ValueError("No common genes between program and expression data")

        # Weighted sum of gene expression
        weights = program_genes[common_genes].values
        expr_subset = gene_expression.loc[common_genes].values

        scores = np.dot(weights, expr_subset)

        return pd.Series(scores, index=gene_expression.columns)

    def get_program_similarity(self) -> pd.DataFrame:
        """
        Compute cosine similarity between gene programs

        Returns
        -------
        pd.DataFrame
            Similarity matrix between programs
        """
        if self.W_df_ is None:
            raise ValueError("Run run_nmf() first")

        from sklearn.metrics.pairwise import cosine_similarity

        similarity = cosine_similarity(self.W_.T)

        return pd.DataFrame(
            similarity,
            index=self.W_df_.columns,
            columns=self.W_df_.columns,
        )

    def get_gene_contributions(self, gene: str) -> pd.Series:
        """
        Get a gene's contribution to each program

        Parameters
        ----------
        gene : str
            Gene name

        Returns
        -------
        pd.Series
            Gene weights across programs
        """
        if self.W_df_ is None:
            raise ValueError("Run run_nmf() first")

        if gene not in self.W_df_.index:
            raise ValueError(f"Gene '{gene}' not found")

        return self.W_df_.loc[gene]

    def summary(self) -> str:
        """Get summary of embedding results"""
        if self.W_df_ is None:
            return "NMF not yet run"

        lines = [
            "BayesPrismEmbedding Summary:",
            f"  Tumor cell type: {self.tumor_key}",
            f"  Number of programs: {self.n_programs}",
            f"  Genes: {self.W_.shape[0]}",
            f"  Samples: {self.H_.shape[1]}",
            f"  Reconstruction error: {self.reconstruction_err_:.4f}",
            "  Top genes per program:",
        ]

        for i in range(self.n_programs):
            program = f"Program_{i+1}"
            top_genes = self.W_df_[program].nlargest(5).index.tolist()
            lines.append(f"    {program}: {', '.join(top_genes)}")

        return "\n".join(lines)


def compare_programs_across_conditions(
    embeddings: List[BayesPrismEmbedding],
    labels: List[str],
) -> pd.DataFrame:
    """
    Compare gene programs across different conditions/datasets

    Parameters
    ----------
    embeddings : List[BayesPrismEmbedding]
        List of embedding objects
    labels : List[str]
        Labels for each embedding

    Returns
    -------
    pd.DataFrame
        Comparison statistics
    """
    results = []

    for emb, label in zip(embeddings, labels):
        if emb.W_df_ is None:
            continue

        for program in emb.W_df_.columns:
            top_genes = set(emb.get_top_genes(program, n=50).index)

            results.append({
                'condition': label,
                'program': program,
                'n_genes': len(top_genes),
                'top_genes': ','.join(list(top_genes)[:10]),
            })

    return pd.DataFrame(results)
