"""
Drug resistance mechanism identification and scoring.

This module provides tools for identifying and quantifying
drug resistance mechanisms in tumor cells.
"""

import logging
from typing import Dict, List, Optional

import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)

# Known drug resistance signatures
RESISTANCE_SIGNATURES = {
    "chemotherapy": {
        "ABC_transporters": ["ABCB1", "ABCC1", "ABCC2", "ABCC3", "ABCG2"],
        "DNA_repair": ["BRCA1", "BRCA2", "RAD51", "ERCC1"],
        "anti_apoptotic": ["BCL2", "BCL2L1", "MCL1", "BCL2A1"],
        "drug_metabolism": ["CYP1A1", "CYP1A2", "CYP3A4", "CYP3A5", "GSTP1"],
    },
    "targeted_therapy": {
        "RTK_signaling": ["EGFR", "ERBB2", "ERBB3", "MET", "IGF1R"],
        "downstream_kinases": ["MAPK1", "MAPK3", "AKT1", "MTOR", "S6K1"],
        "alternative_pathways": ["PIK3CA", "PTEN", "RAS", "RAF1"],
    },
    "immunotherapy": {
        "IFN_gamma": ["IFNG", "STAT1", "IRF1", "JAK1", "JAK2"],
        "antigen_presentation": ["B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2"],
        "inhibitory_receptors": ["PDCD1", "CTLA4", "LAG3", "TIGIT", "TIM3"],
        "exhaustion_markers": ["TOX", "NR4A1", "NFATC1", "EOMES"],
    },
}

# Drug-specific resistance genes
DRUG_SPECIFIC_RESISTANCE = {
    "cisplatin": {
        "mechanism": "DNA repair upregulation",
        "genes": ["ERCC1", "BRCA1", "BRCA2", "RAD51", "MLH1", "MSH2"],
        "pathways": ["NER", "HR", "MMR"],
    },
    "paclitaxel": {
        "mechanism": "Microtubule dynamics alteration",
        "genes": ["TUBB3", "MAP4", "STMN1", "TP53"],
        "pathways": ["Microtubule stabilization", "Apoptosis"],
    },
    "5FU": {
        "mechanism": "TS overexpression",
        "genes": ["TYMS", "DPD", "TP", "OPRT"],
        "pathways": ["Pyrimidine synthesis", "Drug metabolism"],
    },
    "gemcitabine": {
        "mechanism": "RRM1/2 overexpression",
        "genes": ["RRM1", "RRM2", "dCK", "CDA"],
        "pathways": ["Nucleotide synthesis", "Drug transport"],
    },
    "doxorubicin": {
        "mechanism": "ABC transporter upregulation",
        "genes": ["ABCB1", "ABCC1", "TOP2A"],
        "pathways": ["Drug efflux", "DNA topology"],
    },
    "oxaliplatin": {
        "mechanism": "DNA adduct repair",
        "genes": ["ERCC1", "XRCC1", "GSTP1", "GSTM1"],
        "pathways": ["NER", "MMR", "Glutathione metabolism"],
    },
    "tamoxifen": {
        "mechanism": "ER signaling alteration",
        "genes": ["ESR1", "ESR2", "PGR", "TFF1", "GREB1"],
        "pathways": ["Estrogen signaling", "Growth factor signaling"],
    },
    "trastuzumab": {
        "mechanism": "HER2 bypass signaling",
        "genes": ["ERBB2", "ERBB3", "PIK3CA", "PTEN"],
        "pathways": ["PI3K/AKT", "Alternative RTK"],
    },
    "imatinib": {
        "mechanism": "BCR-ABL mutation/kinesis",
        "genes": ["BCR", "ABL1", "SRC", "LYN"],
        "pathways": ["BCR-ABL signaling", "Src family kinases"],
    },
    "gefitinib": {
        "mechanism": "EGFR T790M mutation / MET amplification",
        "genes": ["EGFR", "MET", "HER2", "PIK3CA"],
        "pathways": ["EGFR signaling", "MET signaling"],
    },
    "vemurafenib": {
        "mechanism": "BRAF bypass / NRAS mutation",
        "genes": ["BRAF", "NRAS", "MAPK1", "MAPK3", "MEK1"],
        "pathways": ["MAPK signaling", "RTK signaling"],
    },
    "pembrolizumab": {
        "mechanism": "IFN-gamma pathway defects",
        "genes": ["IFNG", "JAK1", "JAK2", "B2M", "HLA-A"],
        "pathways": ["IFN-gamma signaling", "Antigen presentation"],
    },
}


class ResistanceAnalyzer:
    """
    Analyze drug resistance mechanisms in tumor cells.

    Parameters
    ----------
    resistance_signatures : dict
        Dictionary of resistance gene signatures
    method : str
        Scoring method ("mean", "sum", "weighted")

    Attributes:
    ----------
    resistance_scores_ : pd.DataFrame
        Resistance scores per cell and mechanism
    """

    def __init__(
        self,
        resistance_signatures: Optional[Dict] = None,
        method: str = "mean",
    ):
        self.resistance_signatures = resistance_signatures or RESISTANCE_SIGNATURES
        self.method = method
        self.resistance_scores_: Optional[pd.DataFrame] = None
        self.drug_scores_: Optional[pd.DataFrame] = None

    def fit(self, adata: AnnData, groupby: Optional[str] = None) -> "ResistanceAnalyzer":
        """
        Calculate resistance scores.

        Parameters
        ----------
        adata : AnnData
            Expression data
        groupby : str, optional
            Column to group cells by (e.g., cell type, clone)

        Returns:
        -------
        ResistanceAnalyzer
            Fitted analyzer
        """
        scores = {}

        # Calculate scores for each resistance mechanism
        for drug_class, mechanisms in self.resistance_signatures.items():
            for mechanism, genes in mechanisms.items():
                score_name = f"{drug_class}_{mechanism}"

                # Get available genes
                available = [g for g in genes if g in adata.var_names]

                if len(available) == 0:
                    log.warning(f"No genes found for {score_name}")
                    continue

                # Calculate score
                expr = adata[:, available].X.mean(axis=1)
                if hasattr(expr, "toarray"):
                    expr = expr.toarray().flatten()

                scores[score_name] = expr

        self.resistance_scores_ = pd.DataFrame(scores, index=adata.obs_names)

        # Calculate per-group means if requested
        if groupby is not None and groupby in adata.obs.columns:
            self.group_scores_ = self.resistance_scores_.groupby(adata.obs[groupby]).mean()

        return self

    def score_drug_resistance(
        self,
        adata: AnnData,
        drug: str,
        custom_genes: Optional[List[str]] = None,
    ) -> pd.Series:
        """
        Calculate resistance score for a specific drug.

        Parameters
        ----------
        adata : AnnData
            Expression data
        drug : str
            Drug name (must be in DRUG_SPECIFIC_RESISTANCE)
        custom_genes : list, optional
            Custom resistance genes to use instead

        Returns:
        -------
        pd.Series
            Resistance scores per cell
        """
        if custom_genes is not None:
            genes = custom_genes
        elif drug in DRUG_SPECIFIC_RESISTANCE:
            genes = DRUG_SPECIFIC_RESISTANCE[drug]["genes"]
        else:
            raise ValueError(f"Unknown drug: {drug}. Use custom_genes or known drug.")

        # Get available genes
        available = [g for g in genes if g in adata.var_names]

        if len(available) == 0:
            log.warning(f"No resistance genes found for {drug}")
            return pd.Series(0, index=adata.obs_names)

        log.info(f"Calculating {drug} resistance score using {len(available)} genes")

        # Calculate score
        expr = adata[:, available].X.mean(axis=1)
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()

        return pd.Series(expr, index=adata.obs_names, name=f"{drug}_resistance")

    def identify_resistance_clones(
        self,
        adata: AnnData,
        drug: str,
        clone_key: str = "clone_id",
        top_n: int = 5,
    ) -> pd.DataFrame:
        """
        Identify clones with highest resistance to a drug.

        Parameters
        ----------
        adata : AnnData
            Expression data
        drug : str
            Drug name
        clone_key : str
            Column containing clone IDs
        top_n : int
            Number of top resistant clones to return

        Returns:
        -------
        pd.DataFrame
            Top resistant clones with scores
        """
        scores = self.score_drug_resistance(adata, drug)

        # Calculate mean score per clone
        clone_scores = scores.groupby(adata.obs[clone_key]).agg(["mean", "std", "count"])
        clone_scores.columns = ["mean_score", "std_score", "n_cells"]

        # Sort by mean score
        clone_scores = clone_scores.sort_values("mean_score", ascending=False)

        return clone_scores.head(top_n)

    def get_resistance_mechanism(self, drug: str) -> Dict:
        """
        Get known resistance mechanism for a drug.

        Parameters
        ----------
        drug : str
            Drug name

        Returns:
        -------
        dict
            Resistance mechanism information
        """
        if drug not in DRUG_SPECIFIC_RESISTANCE:
            return {"mechanism": "Unknown", "genes": [], "pathways": []}

        return DRUG_SPECIFIC_RESISTANCE[drug].copy()


def identify_resistance_mechanisms(
    adata: AnnData,
    drug: str,
    groupby: Optional[str] = None,
    key_added: str = "resistance",
) -> pd.DataFrame:
    """
    Identify drug resistance mechanisms in tumor cells.

    Parameters
    ----------
    adata : AnnData
        Expression data
    drug : str
        Drug name
    groupby : str, optional
        Column to group cells
    key_added : str
        Key prefix for storing results

    Returns:
    -------
    pd.DataFrame
        Resistance mechanism scores
    """
    analyzer = ResistanceAnalyzer().fit(adata, groupby=groupby)

    # Get resistance scores
    scores = analyzer.resistance_scores_

    # Get specific drug resistance
    drug_resistance = analyzer.score_drug_resistance(adata, drug)
    scores[f"{drug}_resistance"] = drug_resistance

    # Store in adata
    for col in scores.columns:
        adata.obs[f"{key_added}_{col}"] = scores[col]

    log.info(f"Identified resistance mechanisms for {drug}")

    return scores


def score_drug_resistance(
    adata: AnnData,
    drug: str,
    custom_genes: Optional[List[str]] = None,
    key_added: str = "resistance",
) -> pd.Series:
    """
    Calculate drug resistance scores for all cells.

    Parameters
    ----------
    adata : AnnData
        Expression data
    drug : str
        Drug name
    custom_genes : list, optional
        Custom resistance genes
    key_added : str
        Key for storing results

    Returns:
    -------
    pd.Series
        Resistance scores
    """
    analyzer = ResistanceAnalyzer()
    scores = analyzer.score_drug_resistance(adata, drug, custom_genes)

    adata.obs[f"{key_added}_{drug}"] = scores

    return scores


def compare_resistance_between_groups(
    adata: AnnData,
    groupby: str,
    drug: str,
    method: str = "t-test",
) -> pd.DataFrame:
    """
    Compare drug resistance between cell groups.

    Parameters
    ----------
    adata : AnnData
        Expression data
    groupby : str
        Column to group cells
    drug : str
        Drug name
    method : str
        Statistical test method

    Returns:
    -------
    pd.DataFrame
        Comparison results
    """
    from scipy import stats

    # Calculate resistance scores
    analyzer = ResistanceAnalyzer()
    scores = analyzer.score_drug_resistance(adata, drug)

    # Get unique groups
    groups = adata.obs[groupby].unique()

    results = []

    # Pairwise comparison
    for i, g1 in enumerate(groups):
        for g2 in groups[i + 1 :]:
            mask1 = adata.obs[groupby] == g1
            mask2 = adata.obs[groupby] == g2

            scores1 = scores[mask1]
            scores2 = scores[mask2]

            if method == "t-test":
                stat, pval = stats.ttest_ind(scores1, scores2)
            elif method == "mann-whitney":
                stat, pval = stats.mannwhitneyu(scores1, scores2)
            else:
                raise ValueError(f"Unknown method: {method}")

            results.append(
                {
                    "group1": g1,
                    "group2": g2,
                    "mean1": scores1.mean(),
                    "mean2": scores2.mean(),
                    "diff": scores1.mean() - scores2.mean(),
                    "statistic": stat,
                    "pvalue": pval,
                }
            )

    return pd.DataFrame(results)
