"""
Synthetic data generation for testing scLucid.

Provides fast, reproducible generation of realistic single-cell data
without downloading large datasets.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData


class SyntheticDataGenerator:
    """
    Generator for synthetic single-cell RNA-seq data.

    Creates realistic data with known properties for testing:
    - Cell types with distinct marker expression
    - Batch effects
    - QC metrics variation
    - Known doublets
    """

    def __init__(self, random_state: int = 42):
        self.rng = np.random.RandomState(random_state)

    def generate_adata(
        self,
        n_cells: int = 1000,
        n_genes: int = 2000,
        n_cell_types: int = 4,
        n_batches: int = 2,
        with_qc_metrics: bool = True,
        with_cell_types: bool = True,
        with_batches: bool = True,
        sparsity: float = 0.9,
        marker_fold_change: float = 5.0,
    ) -> AnnData:
        """
        Generate a synthetic AnnData object.

        Args:
            n_cells: Number of cells to generate
            n_genes: Number of genes
            n_cell_types: Number of distinct cell types
            n_batches: Number of batches
            with_qc_metrics: Add realistic QC metrics
            with_cell_types: Add cell type labels
            with_batches: Add batch labels
            sparsity: Proportion of zeros (0-1)
            marker_fold_change: Expression fold change for marker genes

        Returns:
            AnnData with synthetic counts in .X and .layers['counts']
        """
        # Generate base expression matrix
        X = self._generate_expression_matrix(
            n_cells=n_cells,
            n_genes=n_genes,
            sparsity=sparsity,
        )

        # Create AnnData
        obs_names = [f"cell_{i:04d}" for i in range(n_cells)]
        var_names = [f"gene_{i:04d}" for i in range(n_genes)]

        adata = AnnData(
            X=X.copy(),
            obs=pd.DataFrame(index=obs_names),
            var=pd.DataFrame(index=var_names),
        )
        adata.layers["counts"] = X.copy()

        # Add cell type structure
        if with_cell_types:
            self._add_cell_types(
                adata,
                n_cell_types=n_cell_types,
                marker_fold_change=marker_fold_change,
            )

        # Add batch structure
        if with_batches:
            self._add_batches(adata, n_batches=n_batches)

        # Add QC metrics
        if with_qc_metrics:
            self._add_qc_metrics(adata)

        return adata

    def _generate_expression_matrix(
        self,
        n_cells: int,
        n_genes: int,
        sparsity: float,
    ) -> np.ndarray:
        """Generate sparse count matrix using negative binomial distribution."""
        # Negative binomial parameters for realistic counts
        mean_counts = 3.0
        dispersion = 0.5

        # Generate counts
        n = 1.0 / dispersion
        p = n / (n + mean_counts)
        counts = self.rng.negative_binomial(n, p, size=(n_cells, n_genes))

        # Add sparsity
        if sparsity > 0:
            zero_mask = self.rng.random((n_cells, n_genes)) < sparsity
            counts[zero_mask] = 0

        return counts.astype(np.float32)

    def _add_cell_types(
        self,
        adata: AnnData,
        n_cell_types: int,
        marker_fold_change: float,
    ) -> None:
        """Add cell type labels and marker gene expression."""
        n_cells = adata.n_obs
        n_genes = adata.n_vars

        # Assign cells to types
        cell_types = [f"CellType_{i}" for i in range(n_cell_types)]
        assignments = self.rng.choice(cell_types, size=n_cells)
        adata.obs["cell_type"] = pd.Categorical(assignments)

        # Add marker genes for each cell type
        markers_per_type = max(10, n_genes // (n_cell_types * 4))

        for i, ct in enumerate(cell_types):
            # Select marker genes for this cell type
            marker_start = i * markers_per_type
            marker_end = min(marker_start + markers_per_type, n_genes)
            marker_genes = np.arange(marker_start, marker_end)

            # Upregulate in cells of this type
            mask = (adata.obs["cell_type"] == ct).to_numpy()
            adata.X[np.ix_(mask, marker_genes)] *= marker_fold_change
            adata.layers["counts"][np.ix_(mask, marker_genes)] *= marker_fold_change

            # Mark as marker in var
            adata.var[f"marker_{ct}"] = False
            adata.var.loc[adata.var_names[marker_genes], f"marker_{ct}"] = True

    def _add_batches(self, adata: AnnData, n_batches: int) -> None:
        """Add batch labels and batch effects."""
        n_cells = adata.n_obs

        batch_names = [f"batch_{i}" for i in range(n_batches)]
        assignments = self.rng.choice(batch_names, size=n_cells)
        adata.obs["batch"] = pd.Categorical(assignments)
        adata.obs["sampleID"] = adata.obs["batch"]

        # Add batch effects (scaling factors)
        for batch in batch_names:
            mask = adata.obs["batch"] == batch
            scale_factor = self.rng.uniform(0.8, 1.2)
            adata.X[mask] *= scale_factor
            adata.layers["counts"][mask] *= scale_factor

    def _add_qc_metrics(self, adata: AnnData) -> None:
        """Add realistic QC metrics based on expression data."""
        X = adata.X

        # Basic QC metrics
        adata.obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
        adata.obs["total_counts"] = X.sum(axis=1)
        adata.obs["log1p_n_genes_by_counts"] = np.log1p(adata.obs["n_genes_by_counts"])
        adata.obs["log1p_total_counts"] = np.log1p(adata.obs["total_counts"])

        # Mitochondrial genes (simulate ~10% of genes)
        n_mt = max(1, adata.n_vars // 10)
        mt_genes = adata.var_names[:n_mt]
        adata.var["mt"] = False
        adata.var.loc[mt_genes, "mt"] = True

        mt_counts = X[:, :n_mt].sum(axis=1)
        adata.obs["pct_counts_mt"] = 100 * mt_counts / (X.sum(axis=1) + 1e-6)

        # Hemoglobin genes (simulate ~1% of genes)
        n_hb = max(1, adata.n_vars // 100)
        hb_genes = adata.var_names[n_mt : n_mt + n_hb]
        adata.var["hb"] = False
        adata.var.loc[hb_genes, "hb"] = True

        hb_counts = X[:, n_mt : n_mt + n_hb].sum(axis=1)
        adata.obs["pct_counts_hb"] = 100 * hb_counts / (X.sum(axis=1) + 1e-6)

        # Top genes percentage
        for top_n in [20, 50, 100]:
            if top_n < adata.n_vars:
                top_genes_counts = np.partition(X, -top_n, axis=1)[:, -top_n:].sum(axis=1)
                adata.obs[f"pct_counts_in_top_{top_n}_genes"] = (
                    100 * top_genes_counts / (X.sum(axis=1) + 1e-6)
                )

    def generate_with_doublets(
        self,
        n_cells: int = 1000,
        doublet_rate: float = 0.05,
        **kwargs,
    ) -> AnnData:
        """
        Generate data with known doublets.

        Doublets are created by averaging pairs of cells from different cell types.

        Args:
            n_cells: Target number of cells (will be increased for doublets)
            doublet_rate: Proportion of doublets to generate
            **kwargs: Passed to generate_adata()

        Returns:
            AnnData with .obs['is_doublet'] column
        """
        # Generate more cells to account for doublet merging
        n_real_cells = int(n_cells / (1 + doublet_rate))
        adata = self.generate_adata(n_cells=n_real_cells, **kwargs)

        # Create doublets by averaging pairs
        n_doublets = int(n_cells * doublet_rate)
        doublet_indices = []

        cell_types = adata.obs["cell_type"].unique()

        for _ in range(n_doublets):
            # Pick two cells from different cell types
            ct1, ct2 = self.rng.choice(cell_types, size=2, replace=False)
            cells1 = adata.obs_names[adata.obs["cell_type"] == ct1].tolist()
            cells2 = adata.obs_names[adata.obs["cell_type"] == ct2].tolist()

            idx1 = self.rng.choice(cells1)
            idx2 = self.rng.choice(cells2)

            # Create doublet expression
            doublet_expr = (adata[idx1].X + adata[idx2].X) / 2
            doublet_counts = (adata[idx1].layers["counts"] + adata[idx2].layers["counts"]) / 2

            # Add to adata
            new_name = f"doublet_{idx1}_{idx2}"
            doublet_indices.append(new_name)

            # Extend matrices
            adata.obs.loc[new_name] = adata.obs.loc[idx1].copy()
            adata.obs.loc[new_name, "cell_type"] = f"Doublet_{ct1}_{ct2}"
            adata.X = np.vstack([adata.X, doublet_expr])
            adata.layers["counts"] = np.vstack([adata.layers["counts"], doublet_counts])

        # Mark doublets
        adata.obs["is_doublet"] = False
        adata.obs.loc[doublet_indices, "is_doublet"] = True

        return adata

    def generate_batch_effects_data(
        self,
        n_cells: int = 2000,
        n_batches: int = 4,
        batch_effect_strength: float = 0.3,
        **kwargs,
    ) -> AnnData:
        """
        Generate data with strong batch effects for integration testing.

        Args:
            n_cells: Total number of cells
            n_batches: Number of batches
            batch_effect_strength: Strength of batch effect (0-1)
            **kwargs: Passed to generate_adata()

        Returns:
            AnnData with strong batch effects
        """
        adata = self.generate_adata(
            n_cells=n_cells,
            n_batches=n_batches,
            with_batches=True,
            **kwargs,
        )

        # Add strong batch-specific gene expression shifts
        batch_names = adata.obs["batch"].unique()
        n_shift_genes = max(50, adata.n_vars // 20)

        for i, batch in enumerate(batch_names):
            mask = adata.obs["batch"] == batch

            # Select batch-specific genes
            shift_genes = range(i * n_shift_genes, (i + 1) * n_shift_genes)
            shift_genes = [g for g in shift_genes if g < adata.n_vars]

            # Apply batch shift
            shift_factor = 1 + batch_effect_strength * self.rng.uniform(-1, 1)
            adata.X[np.ix_(mask, shift_genes)] *= shift_factor
            adata.layers["counts"][np.ix_(mask, shift_genes)] *= shift_factor

        return adata


def generate_minimal_adata(n_cells: int = 100, n_genes: int = 200) -> AnnData:
    """
    Generate a minimal AnnData for quick tests.

    This is a convenience function for simple test cases.
    """
    gen = SyntheticDataGenerator()
    return gen.generate_adata(
        n_cells=n_cells,
        n_genes=n_genes,
        n_cell_types=2,
        n_batches=1,
        with_qc_metrics=True,
    )


def generate_qc_test_data(n_cells: int = 500) -> AnnData:
    """
    Generate data specifically for QC testing.

    Includes outliers for filtering tests.
    """
    gen = SyntheticDataGenerator()
    adata = gen.generate_adata(n_cells=n_cells, with_qc_metrics=True)

    # Add some low-quality outliers
    n_outliers = n_cells // 20
    outlier_idx = gen.rng.choice(adata.obs_names, size=n_outliers, replace=False)

    # Low gene count outliers
    for idx in outlier_idx[: n_outliers // 2]:
        adata.X[adata.obs_names.get_loc(idx)] *= 0.1
        adata.layers["counts"][adata.obs_names.get_loc(idx)] *= 0.1

    # High mitochondrial outliers
    for idx in outlier_idx[n_outliers // 2 :]:
        mt_genes = adata.var_names[adata.var.get("mt", False)]
        if len(mt_genes) > 0:
            adata[idx, mt_genes].X *= 10
            adata[idx, mt_genes].layers["counts"] *= 10

    # Recalculate QC metrics
    X = adata.X
    adata.obs["n_genes_by_counts"] = (X > 0).sum(axis=1)
    adata.obs["total_counts"] = X.sum(axis=1)
    adata.obs["log1p_n_genes_by_counts"] = np.log1p(adata.obs["n_genes_by_counts"])
    adata.obs["log1p_total_counts"] = np.log1p(adata.obs["total_counts"])

    mt_counts = X[:, adata.var.get("mt", False)].sum(axis=1)
    adata.obs["pct_counts_mt"] = 100 * mt_counts / (X.sum(axis=1) + 1e-6)

    return adata


def generate_integration_test_data(n_cells: int = 1000, n_batches: int = 3) -> AnnData:
    """
    Generate data for batch correction/integration testing.
    """
    gen = SyntheticDataGenerator()
    return gen.generate_batch_effects_data(
        n_cells=n_cells,
        n_batches=n_batches,
        batch_effect_strength=0.5,
    )


# Pytest fixtures
import pytest


@pytest.fixture
def synthetic_generator():
    """Fixture providing a SyntheticDataGenerator instance."""
    return SyntheticDataGenerator()


@pytest.fixture
def minimal_adata():
    """Fixture providing minimal synthetic data."""
    return generate_minimal_adata()


@pytest.fixture
def qc_test_adata():
    """Fixture providing QC test data with outliers."""
    return generate_qc_test_data()


@pytest.fixture
def integration_test_adata():
    """Fixture providing data with batch effects."""
    return generate_integration_test_data()


@pytest.fixture
def doublet_test_adata():
    """Fixture providing data with known doublets."""
    gen = SyntheticDataGenerator()
    return gen.generate_with_doublets(n_cells=500, doublet_rate=0.1)
