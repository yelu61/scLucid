"""
Tumor ecosystem analysis.

This module provides tools for analyzing the tumor ecosystem
as a whole, including cellular composition and spatial organization.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Union
from anndata import AnnData
import logging

log = logging.getLogger(__name__)


class EcosystemAnalyzer:
    """
    Analyze tumor ecosystem composition and properties.

    Parameters
    ----------
    cell_type_key : str
        Column containing cell type annotations

    Attributes
    ----------
    composition_ : pd.DataFrame
        Cellular composition of ecosystem
    """

    def __init__(
        self,
        cell_type_key: str = "cell_type",
    ):
        self.cell_type_key = cell_type_key
        self.composition_: Optional[pd.DataFrame] = None

    def analyze_ecosystem_composition(
        self,
        adata: AnnData,
        groupby: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Analyze cellular composition of tumor ecosystem.

        Parameters
        ----------
        adata : AnnData
            Expression data
        groupby : str, optional
            Column to group by (e.g., sample, patient)

        Returns
        -------
        pd.DataFrame
            Ecosystem composition
        """
        if self.cell_type_key not in adata.obs.columns:
            raise ValueError(f"Cell type column {self.cell_type_key} not found")

        results = []

        if groupby is None:
            # Single analysis
            composition = adata.obs[self.cell_type_key].value_counts(normalize=True)

            result = {"group": "all"}
            for cell_type, freq in composition.items():
                result[f"{cell_type}_freq"] = freq
            result["total_cells"] = len(adata)

            results.append(result)
        else:
            # Grouped analysis
            for group in adata.obs[groupby].unique():
                mask = adata.obs[groupby] == group
                group_adata = adata[mask]

                composition = group_adata.obs[self.cell_type_key].value_counts(normalize=True)

                result = {"group": group}
                for cell_type, freq in composition.items():
                    result[f"{cell_type}_freq"] = freq
                result["total_cells"] = len(group_adata)

                results.append(result)

        self.composition_ = pd.DataFrame(results)
        return self.composition_

    def calculate_ecosystem_diversity(
        self,
        adata: AnnData,
        groupby: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Calculate diversity metrics for tumor ecosystem.

        Parameters
        ----------
        adata : AnnData
            Expression data
        groupby : str, optional
            Column to group by

        Returns
        -------
        pd.DataFrame
            Diversity metrics
        """
        results = []

        if groupby is None:
            groups = [("all", adata)]
        else:
            groups = [
                (g, adata[adata.obs[groupby] == g])
                for g in adata.obs[groupby].unique()
            ]

        for group_name, group_adata in groups:
            composition = group_adata.obs[self.cell_type_key].value_counts(normalize=True)
            proportions = composition.values

            # Shannon diversity
            shannon = -np.sum(proportions * np.log(proportions + 1e-10))

            # Simpson diversity
            simpson = 1 - np.sum(proportions ** 2)

            # Richness
            richness = len(proportions)

            results.append({
                "group": group_name,
                "shannon_diversity": shannon,
                "simpson_diversity": simpson,
                "richness": richness,
                "evenness": shannon / np.log(richness) if richness > 1 else 0,
                "dominant_type": composition.index[0],
                "dominant_freq": proportions[0],
            })

        return pd.DataFrame(results)

    def identify_ecosystem_subtypes(
        self,
        adata: AnnData,
        n_clusters: int = 3,
    ) -> pd.Series:
        """
        Identify ecosystem subtypes based on composition.

        Parameters
        ----------
        adata : AnnData
            Expression data
        n_clusters : int
            Number of ecosystem subtypes

        Returns
        -------
        pd.Series
            Ecosystem subtype assignments
        """
        from sklearn.cluster import KMeans

        # Get composition matrix
        composition = self.analyze_ecosystem_composition(adata)

        # Select frequency columns
        freq_cols = [c for c in composition.columns if c.endswith("_freq")]
        X = composition[freq_cols].fillna(0).values

        # Cluster
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        labels = kmeans.fit_predict(X)

        return pd.Series(labels, index=composition["group"], name="ecosystem_subtype")


def analyze_ecosystem_composition(
    adata: AnnData,
    cell_type_key: str = "cell_type",
    groupby: Optional[str] = None,
) -> pd.DataFrame:
    """
    Analyze tumor ecosystem composition.

    Parameters
    ----------
    adata : AnnData
        Expression data
    cell_type_key : str
        Column with cell types
    groupby : str, optional
        Column to group by

    Returns
    -------
    pd.DataFrame
        Ecosystem composition
    """
    analyzer = EcosystemAnalyzer(cell_type_key=cell_type_key)
    return analyzer.analyze_ecosystem_composition(adata, groupby)


def calculate_tumor_microenvironment_score(
    adata: AnnData,
    cell_type_key: str = "cell_type",
    immune_types: Optional[List[str]] = None,
    stromal_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Calculate tumor microenvironment scores.

    Parameters
    ----------
    adata : AnnData
        Expression data
    cell_type_key : str
        Column with cell types
    immune_types : list, optional
        List of immune cell types
    stromal_types : list, optional
        List of stromal cell types

    Returns
    -------
    pd.DataFrame
        TME scores
    """
    if immune_types is None:
        immune_types = ["T cell", "B cell", "NK cell", "Macrophage", "Monocyte", "DC"]

    if stromal_types is None:
        stromal_types = ["Fibroblast", "Endothelial", "Pericyte"]

    composition = adata.obs[cell_type_key].value_counts(normalize=True)

    # Calculate scores
    immune_score = sum(composition.get(ct, 0) for ct in immune_types)
    stromal_score = sum(composition.get(ct, 0) for ct in stromal_types)

    # Calculate ratios
    tumor_score = 1 - immune_score - stromal_score
    immune_to_tumor = immune_score / (tumor_score + 1e-6)

    return pd.DataFrame([{
        "immune_score": immune_score,
        "stromal_score": stromal_score,
        "tumor_score": tumor_score,
        "immune_to_tumor_ratio": immune_to_tumor,
        "stromal_to_tumor_ratio": stromal_score / (tumor_score + 1e-6),
    }])


def compare_ecosystems(
    adata1: AnnData,
    adata2: AnnData,
    cell_type_key: str = "cell_type",
) -> pd.DataFrame:
    """
    Compare two tumor ecosystems.

    Parameters
    ----------
    adata1 : AnnData
        First ecosystem
    adata2 : AnnData
        Second ecosystem
    cell_type_key : str
        Column with cell types

    Returns
    -------
    pd.DataFrame
        Comparison results
    """
    comp1 = adata1.obs[cell_type_key].value_counts(normalize=True)
    comp2 = adata2.obs[cell_type_key].value_counts(normalize=True)

    # Get all cell types
    all_types = set(comp1.index) | set(comp2.index)

    results = []
    for cell_type in all_types:
        results.append({
            "cell_type": cell_type,
            "freq_ecosystem1": comp1.get(cell_type, 0),
            "freq_ecosystem2": comp2.get(cell_type, 0),
            "difference": comp2.get(cell_type, 0) - comp1.get(cell_type, 0),
        })

    return pd.DataFrame(results).sort_values("difference", ascending=False)
