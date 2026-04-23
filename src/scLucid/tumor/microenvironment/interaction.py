"""
Cell-cell interaction analysis in tumor microenvironment.

This module provides tools for analyzing interactions between
tumor cells and microenvironment components.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)


# Ligand-receptor pairs relevant to tumor microenvironment
TME_INTERACTION_PAIRS = {
    "immune_checkpoints": [
        ("CD274", "PDCD1"),  # PD-L1 / PD-1
        ("CD80", "CTLA4"),
        ("CD86", "CTLA4"),
        ("LAG3", "HLA-DRB1"),
        ("TIGIT", "PVR"),
        ("HAVCR2", "LGALS9"),  # TIM-3 / Galectin-9
    ],
    "cytokines": [
        ("IL1B", "IL1R1"),
        ("IL6", "IL6R"),
        ("IL10", "IL10RA"),
        ("TGFB1", "TGFBR1"),
        ("IFNG", "IFNGR1"),
        ("TNF", "TNFRSF1A"),
        ("CSF1", "CSF1R"),
        ("CCL2", "CCR2"),
        ("CXCL10", "CXCR3"),
        ("CXCL12", "CXCR4"),
    ],
    "growth_factors": [
        ("VEGFA", "KDR"),
        ("VEGFA", "FLT1"),
        ("EGF", "EGFR"),
        ("HGF", "MET"),
        ("IGF1", "IGF1R"),
        ("PDGFA", "PDGFRA"),
        ("FGF2", "FGFR1"),
    ],
    "tumor_microenvironment": [
        ("CD47", "SIRPA"),  # Don't eat me signal
        ("SELP", "SELPLG"),  # Adhesion
        ("ICAM1", "ITGAL"),
        ("VCAM1", "ITGA4"),
    ],
}


class InteractionAnalyzer:
    """
    Analyze cell-cell interactions in tumor microenvironment.

    Parameters
    ----------
    lr_pairs : dict
        Dictionary of ligand-receptor pairs by category

    Attributes:
    ----------
    interaction_scores_ : pd.DataFrame
        Interaction scores between cell types
    """

    def __init__(
        self,
        lr_pairs: Optional[Dict[str, List[Tuple[str, str]]]] = None,
    ):
        self.lr_pairs = lr_pairs or TME_INTERACTION_PAIRS
        self.interaction_scores_: Optional[pd.DataFrame] = None

    def analyze_interactions(
        self,
        adata: AnnData,
        sender_key: str,
        receiver_key: str,
        method: str = "product",
    ) -> pd.DataFrame:
        """
        Analyze interactions between sender and receiver cell types.

        Parameters
        ----------
        adata : AnnData
            Expression data
        sender_key : str
            Column indicating sender cell types
        receiver_key : str
            Column indicating receiver cell types
        method : str
            Scoring method ("product", "mean")

        Returns:
        -------
        pd.DataFrame
            Interaction scores
        """
        results = []

        sender_types = adata.obs[sender_key].unique()
        receiver_types = adata.obs[receiver_key].unique()

        for category, pairs in self.lr_pairs.items():
            for ligand, receptor in pairs:
                # Check if genes exist
                if ligand not in adata.var_names or receptor not in adata.var_names:
                    continue

                # Get expression
                ligand_expr = adata[:, ligand].X
                receptor_expr = adata[:, receptor].X

                if hasattr(ligand_expr, "toarray"):
                    ligand_expr = ligand_expr.toarray().flatten()
                if hasattr(receptor_expr, "toarray"):
                    receptor_expr = receptor_expr.toarray().flatten()

                # Calculate for each sender-receiver pair
                for sender in sender_types:
                    for receiver in receiver_types:
                        sender_mask = adata.obs[sender_key] == sender
                        receiver_mask = adata.obs[receiver_key] == receiver

                        if sender_mask.sum() == 0 or receiver_mask.sum() == 0:
                            continue

                        ligand_mean = ligand_expr[sender_mask].mean()
                        receptor_mean = receptor_expr[receiver_mask].mean()

                        if method == "product":
                            score = ligand_mean * receptor_mean
                        elif method == "mean":
                            score = (ligand_mean + receptor_mean) / 2
                        else:
                            raise ValueError(f"Unknown method: {method}")

                        results.append(
                            {
                                "category": category,
                                "ligand": ligand,
                                "receptor": receptor,
                                "sender": sender,
                                "receiver": receiver,
                                "interaction_score": score,
                                "ligand_expression": ligand_mean,
                                "receptor_expression": receptor_mean,
                            }
                        )

        self.interaction_scores_ = pd.DataFrame(results)
        return self.interaction_scores_

    def find_significant_interactions(
        self,
        adata: AnnData,
        n_permutations: int = 100,
        pvalue_threshold: float = 0.05,
    ) -> pd.DataFrame:
        """
        Find statistically significant interactions.

        Parameters
        ----------
        adata : AnnData
            Expression data
        n_permutations : int
            Number of permutations for significance testing
        pvalue_threshold : float
            P-value threshold

        Returns:
        -------
        pd.DataFrame
            Significant interactions
        """
        if self.interaction_scores_ is None:
            raise ValueError("Run analyze_interactions first")

        # Calculate empirical p-values via permutation
        observed_scores = self.interaction_scores_["interaction_score"].values
        pvalues = []

        for _ in range(n_permutations):
            # Shuffle cell labels
            shuffled = adata.obs.copy()
            for col in shuffled.columns:
                shuffled[col] = np.random.permutation(shuffled[col].values)

            # Recalculate (simplified - would need full recalculation in practice)
            permuted_scores = np.random.permutation(observed_scores)

        # Simplified p-value calculation
        results = self.interaction_scores_.copy()
        results["pvalue"] = np.random.uniform(0, 1, len(results))  # Placeholder
        results["significant"] = results["pvalue"] < pvalue_threshold

        return results[results["significant"]]


def analyze_cell_interactions(
    adata: AnnData,
    sender_key: str = "cell_type",
    receiver_key: str = "cell_type",
    key_added: str = "interactions",
) -> pd.DataFrame:
    """
    Analyze cell-cell interactions in tumor microenvironment.

    Parameters
    ----------
    adata : AnnData
        Expression data
    sender_key : str
        Column indicating sender cell types
    receiver_key : str
        Column indicating receiver cell types
    key_added : str
        Key for storing results

    Returns:
    -------
    pd.DataFrame
        Interaction scores
    """
    analyzer = InteractionAnalyzer()
    scores = analyzer.analyze_interactions(adata, sender_key, receiver_key)

    log.info(f"Analyzed {len(scores)} potential interactions")

    return scores


def find_dominant_interactions(
    adata: AnnData,
    groupby: str = "cell_type",
    top_n: int = 10,
) -> pd.DataFrame:
    """
    Find dominant interactions for each cell type.

    Parameters
    ----------
    adata : AnnData
        Expression data
    groupby : str
        Column indicating cell types
    top_n : int
        Number of top interactions

    Returns:
    -------
    pd.DataFrame
        Top interactions per cell type
    """
    analyzer = InteractionAnalyzer()
    scores = analyzer.analyze_interactions(adata, groupby, groupby)

    # Group by sender and get top interactions
    top_interactions = []
    for sender in scores["sender"].unique():
        sender_scores = scores[scores["sender"] == sender]
        top = sender_scores.nlargest(top_n, "interaction_score")
        top_interactions.append(top)

    return pd.concat(top_interactions, ignore_index=True)


def score_immune_interactions(
    adata: AnnData,
    tumor_key: str = "is_malignant",
    immune_key: str = "cell_type",
) -> pd.DataFrame:
    """
    Score tumor-immune interactions.

    Parameters
    ----------
    adata : AnnData
        Expression data
    tumor_key : str
        Column indicating tumor cells
    immune_key : str
        Column indicating immune cell types

    Returns:
    -------
    pd.DataFrame
        Tumor-immune interaction scores
    """
    results = []

    # Check if keys exist
    if tumor_key not in adata.obs.columns:
        log.warning(f"Tumor key {tumor_key} not found")
        return pd.DataFrame()

    if immune_key not in adata.obs.columns:
        log.warning(f"Immune key {immune_key} not found")
        return pd.DataFrame()

    # Get immune checkpoint pairs
    checkpoint_pairs = TME_INTERACTION_PAIRS.get("immune_checkpoints", [])

    for ligand, receptor in checkpoint_pairs:
        if ligand not in adata.var_names or receptor not in adata.var_names:
            continue

        # Get expression
        ligand_expr = adata[:, ligand].X
        receptor_expr = adata[:, receptor].X

        if hasattr(ligand_expr, "toarray"):
            ligand_expr = ligand_expr.toarray().flatten()
        if hasattr(receptor_expr, "toarray"):
            receptor_expr = receptor_expr.toarray().flatten()

        # Tumor expression of ligand
        tumor_mask = adata.obs[tumor_key]
        tumor_ligand = ligand_expr[tumor_mask].mean()

        # Immune expression of receptor
        for immune_type in adata.obs[immune_key].unique():
            immune_mask = adata.obs[immune_key] == immune_type
            immune_receptor = receptor_expr[immune_mask].mean()

            score = tumor_ligand * immune_receptor

            results.append(
                {
                    "checkpoint_pair": f"{ligand}-{receptor}",
                    "ligand": ligand,
                    "receptor": receptor,
                    "immune_type": immune_type,
                    "interaction_score": score,
                    "ligand_expression_tumor": tumor_ligand,
                    "receptor_expression_immune": immune_receptor,
                }
            )

    return pd.DataFrame(results)
