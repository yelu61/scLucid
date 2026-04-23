"""
Therapeutic target discovery and prioritization.

This module provides tools for discovering and prioritizing
druggable targets in tumor cells.
"""

import logging
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)

# Druggable gene categories
DRUGGABLE_CATEGORIES = {
    "kinase": {
        "description": "Protein kinases",
        "genes": [
            "EGFR",
            "ERBB2",
            "ERBB3",
            "ERBB4",
            "MET",
            "ROS1",
            "ALK",
            "RET",
            "BRAF",
            "RAF1",
            "MAPK1",
            "MAPK3",
            "MAP2K1",
            "MAP2K2",
            "PIK3CA",
            "PIK3CB",
            "PIK3CD",
            "PIK3CG",
            "MTOR",
            "AKT1",
            "AKT2",
            "AKT3",
            "SRC",
            "LCK",
            "FYN",
            "YES1",
            "ABL1",
            "ABL2",
            "KIT",
            "PDGFRA",
            "PDGFRB",
            "FLT1",
            "FLT3",
            "FLT4",
            "KDR",
            "CSF1R",
            "CSF3R",
            "JAK1",
            "JAK2",
            "JAK3",
            "TYK2",
            "STAT1",
            "STAT3",
            "STAT5A",
            "STAT5B",
        ],
    },
    "receptor": {
        "description": "Cell surface receptors",
        "genes": [
            "EGFR",
            "ERBB2",
            "ERBB3",
            "ERBB4",
            "MET",
            "IGF1R",
            "INSR",
            "FGFR1",
            "FGFR2",
            "FGFR3",
            "FGFR4",
            "PDGFRA",
            "PDGFRB",
            "VEGFA",
            "VEGFB",
            "VEGFC",
            "VEGFD",
            "KDR",
            "FLT1",
            "FLT4",
            "TNFRSF1A",
            "TNFRSF1B",
            "FAS",
            "TNFRSF10A",
            "TNFRSF10B",
            "CD40",
            "CD27",
            "TNFRSF4",
            "TNFRSF9",
            "TNFRSF18",
            "CCR1",
            "CCR2",
            "CCR4",
            "CCR5",
            "CXCR4",
            "CXCR2",
            "CXCR3",
        ],
    },
    "immune_checkpoint": {
        "description": "Immune checkpoint molecules",
        "genes": [
            "CD274",
            "PDCD1",
            "CTLA4",
            "CD28",
            "CD80",
            "CD86",
            "LAG3",
            "HAVCR2",
            "TIGIT",
            "SIGLEC15",
            "VTCN1",
            "ICOS",
            "ICOSLG",
            "TNFSF4",
            "TNFRSF4",
            "TNFSF9",
            "TNFRSF9",
            "TNFSF18",
            "TNFRSF18",
            "CD40",
            "CD40LG",
            "CD47",
            "SIRPA",
            "SIRPB1",
            "NECTIN2",
            "CD96",
        ],
    },
    "epigenetic": {
        "description": "Epigenetic regulators",
        "genes": [
            "DNMT1",
            "DNMT3A",
            "DNMT3B",
            "DNMT3L",
            "HDAC1",
            "HDAC2",
            "HDAC3",
            "HDAC6",
            "HDAC8",
            "SIRT1",
            "SIRT2",
            "SIRT3",
            "EZH2",
            "EED",
            "SUZ12",
            "KDM1A",
            "KDM4A",
            "KDM4B",
            "KDM5A",
            "KDM5B",
            "KDM5C",
            "KDM6A",
            "KDM6B",
            "BRD2",
            "BRD3",
            "BRD4",
            "BRDT",
            "IDH1",
            "IDH2",
            "TET1",
            "TET2",
            "TET3",
        ],
    },
    "dna_repair": {
        "description": "DNA repair enzymes",
        "genes": [
            "PARP1",
            "PARP2",
            "PARP3",
            "BRCA1",
            "BRCA2",
            "RAD51",
            "RAD51C",
            "RAD51D",
            "RAD54L",
            "ATM",
            "ATR",
            "CHEK1",
            "CHEK2",
            "TP53",
            "MDM2",
            "MDM4",
            "ERCC1",
            "ERCC2",
            "ERCC3",
            "ERCC4",
            "ERCC5",
            "XPF",
            "XPG",
            "MLH1",
            "MSH2",
            "MSH6",
            "PMS2",
            "MLH3",
        ],
    },
    "metabolism": {
        "description": "Metabolic enzymes",
        "genes": [
            "IDH1",
            "IDH2",
            "FH",
            "SDHA",
            "SDHB",
            "SDHC",
            "SDHD",
            "LDHA",
            "LDHB",
            "PKM",
            "HK1",
            "HK2",
            "HK3",
            "PFKFB3",
            "PFKL",
            "GLS",
            "GLS2",
            "GLUL",
            "ASNS",
            "ACLY",
            "ACC1",
            "FASN",
            "SCD",
            "ACACA",
            "NAMPT",
            "NNMT",
            "IDO1",
            "IDO2",
            "TDO2",
        ],
    },
}

# Known drug-target interactions
DRUG_TARGETS = {
    # Targeted therapies
    "trastuzumab": {"targets": ["ERBB2"], "category": "mAb", "indications": ["breast", "gastric"]},
    "pertuzumab": {"targets": ["ERBB2"], "category": "mAb", "indications": ["breast"]},
    "ado-trastuzumab": {"targets": ["ERBB2"], "category": "ADC", "indications": ["breast"]},
    "lapatinib": {"targets": ["EGFR", "ERBB2"], "category": "TKI", "indications": ["breast"]},
    "neratinib": {
        "targets": ["EGFR", "ERBB2", "ERBB4"],
        "category": "TKI",
        "indications": ["breast"],
    },
    "gefitinib": {"targets": ["EGFR"], "category": "TKI", "indications": ["lung"]},
    "erlotinib": {"targets": ["EGFR"], "category": "TKI", "indications": ["lung", "pancreatic"]},
    "afatinib": {"targets": ["EGFR", "ERBB2", "ERBB4"], "category": "TKI", "indications": ["lung"]},
    "osimertinib": {"targets": ["EGFR"], "category": "TKI", "indications": ["lung"]},
    "crizotinib": {"targets": ["ALK", "ROS1", "MET"], "category": "TKI", "indications": ["lung"]},
    "alectinib": {"targets": ["ALK"], "category": "TKI", "indications": ["lung"]},
    "brigatinib": {"targets": ["ALK", "ROS1"], "category": "TKI", "indications": ["lung"]},
    "lorlatinib": {"targets": ["ALK", "ROS1"], "category": "TKI", "indications": ["lung"]},
    "vemurafenib": {"targets": ["BRAF"], "category": "TKI", "indications": ["melanoma"]},
    "dabrafenib": {"targets": ["BRAF"], "category": "TKI", "indications": ["melanoma", "lung"]},
    "trametinib": {
        "targets": ["MAP2K1", "MAP2K2"],
        "category": "TKI",
        "indications": ["melanoma", "lung"],
    },
    "imatinib": {
        "targets": ["ABL1", "KIT", "PDGFRA"],
        "category": "TKI",
        "indications": ["CML", "GIST"],
    },
    "dasatinib": {"targets": ["ABL1", "SRC"], "category": "TKI", "indications": ["CML", "ALL"]},
    "nilotinib": {"targets": ["ABL1"], "category": "TKI", "indications": ["CML"]},
    "ponatinib": {
        "targets": ["ABL1", "FGFR", "FLT3"],
        "category": "TKI",
        "indications": ["CML", "ALL"],
    },
    "sunitinib": {
        "targets": ["KIT", "PDGFRA", "PDGFRB", "VEGFR", "FLT3"],
        "category": "TKI",
        "indications": ["GIST", "RCC"],
    },
    "sorafenib": {
        "targets": ["RAF", "VEGFR", "PDGFR", "KIT", "FLT3"],
        "category": "TKI",
        "indications": ["HCC", "RCC", "thyroid"],
    },
    "lenvatinib": {
        "targets": ["VEGFR", "FGFR", "PDGFR", "KIT", "RET"],
        "category": "TKI",
        "indications": ["thyroid", "HCC", "RCC"],
    },
    "regorafenib": {
        "targets": ["KIT", "PDGFR", "VEGFR", "RAF"],
        "category": "TKI",
        "indications": ["CRC", "GIST", "HCC"],
    },
    "cabozantinib": {
        "targets": ["MET", "VEGFR", "RET", "KIT", "FLT3"],
        "category": "TKI",
        "indications": ["thyroid", "RCC", "HCC"],
    },
    "vandetanib": {
        "targets": ["RET", "VEGFR", "EGFR"],
        "category": "TKI",
        "indications": ["thyroid"],
    },
    "selpercatinib": {"targets": ["RET"], "category": "TKI", "indications": ["thyroid", "lung"]},
    "pralsetinib": {"targets": ["RET"], "category": "TKI", "indications": ["thyroid", "lung"]},
    "larotrectinib": {
        "targets": ["NTRK1", "NTRK2", "NTRK3"],
        "category": "TKI",
        "indications": ["solid_tumors"],
    },
    "entrectinib": {
        "targets": ["NTRK1", "NTRK2", "NTRK3", "ROS1", "ALK"],
        "category": "TKI",
        "indications": ["solid_tumors", "lung"],
    },
    "alpelisib": {"targets": ["PIK3CA"], "category": "TKI", "indications": ["breast"]},
    "copanlisib": {"targets": ["PIK3CA"], "category": "TKI", "indications": ["lymphoma"]},
    "idelalisib": {"targets": ["PIK3CD"], "category": "TKI", "indications": ["lymphoma", "CLL"]},
    "everolimus": {
        "targets": ["MTOR"],
        "category": "TKI",
        "indications": ["RCC", "breast", "pNET"],
    },
    "temsirolimus": {"targets": ["MTOR"], "category": "TKI", "indications": ["RCC"]},
    "ruxolitinib": {
        "targets": ["JAK1", "JAK2"],
        "category": "TKI",
        "indications": ["myelofibrosis"],
    },
    "fedratinib": {"targets": ["JAK2"], "category": "TKI", "indications": ["myelofibrosis"]},
    "pacritinib": {"targets": ["JAK2"], "category": "TKI", "indications": ["myelofibrosis"]},
    # Immunotherapies
    "pembrolizumab": {"targets": ["PDCD1"], "category": "mAb", "indications": ["multiple"]},
    "nivolumab": {"targets": ["PDCD1"], "category": "mAb", "indications": ["multiple"]},
    "atezolizumab": {"targets": ["CD274"], "category": "mAb", "indications": ["multiple"]},
    "avelumab": {"targets": ["CD274"], "category": "mAb", "indications": ["MCC", "RCC", "UC"]},
    "durvalumab": {"targets": ["CD274"], "category": "mAb", "indications": ["lung", "SCLC"]},
    "cemiplimab": {"targets": ["PDCD1"], "category": "mAb", "indications": ["CSCC", "BCC"]},
    "dostarlimab": {"targets": ["PDCD1"], "category": "mAb", "indications": ["dMMR"]},
    "ipilimumab": {"targets": ["CTLA4"], "category": "mAb", "indications": ["melanoma"]},
    "tremelimumab": {"targets": ["CTLA4"], "category": "mAb", "indications": ["HCC"]},
    # PARP inhibitors
    "olaparib": {
        "targets": ["PARP1", "PARP2"],
        "category": "PARPi",
        "indications": ["ovarian", "breast", "pancreatic", "prostate"],
    },
    "rucaparib": {
        "targets": ["PARP1", "PARP2"],
        "category": "PARPi",
        "indications": ["ovarian", "prostate"],
    },
    "niraparib": {"targets": ["PARP1", "PARP2"], "category": "PARPi", "indications": ["ovarian"]},
    "talazoparib": {
        "targets": ["PARP1", "PARP2"],
        "category": "PARPi",
        "indications": ["breast", "prostate"],
    },
}


class TargetDiscovery:
    """
    Discover and prioritize therapeutic targets.

    Parameters
    ----------
    druggable_categories : dict
        Dictionary of druggable gene categories

    Attributes:
    ----------
    targets_ : pd.DataFrame
        Prioritized targets
    """

    def __init__(
        self,
        druggable_categories: Optional[Dict] = None,
    ):
        self.druggable_categories = druggable_categories or DRUGGABLE_CATEGORIES
        self.targets_: Optional[pd.DataFrame] = None

    def discover_targets(
        self,
        adata: AnnData,
        groupby: Optional[str] = None,
        malignant_key: str = "is_malignant",
        min_expression: float = 0.1,
    ) -> pd.DataFrame:
        """
        Discover therapeutic targets in tumor cells.

        Parameters
        ----------
        adata : AnnData
            Expression data
        groupby : str, optional
            Column to group cells
        malignant_key : str
            Column indicating malignant cells
        min_expression : float
            Minimum expression threshold

        Returns:
        -------
        pd.DataFrame
            Discovered targets with scores
        """
        targets = []

        # Get all druggable genes
        all_druggable = []
        for category, info in self.druggable_categories.items():
            all_druggable.extend(info["genes"])
        all_druggable = list(set(all_druggable))

        # Filter to available genes
        available = [g for g in all_druggable if g in adata.var_names]

        log.info(f"Evaluating {len(available)} druggable genes")

        # Calculate expression metrics
        for gene in available:
            gene_data = {
                "gene": gene,
                "category": self._get_gene_category(gene),
            }

            if malignant_key in adata.obs.columns:
                # Compare malignant vs non-malignant
                malignant_mask = adata.obs[malignant_key]

                mal_expr = adata[malignant_mask, gene].X
                if hasattr(mal_expr, "toarray"):
                    mal_expr = mal_expr.toarray().flatten()

                normal_expr = (
                    adata[~malignant_mask, gene].X if (~malignant_mask).sum() > 0 else np.array([0])
                )
                if hasattr(normal_expr, "toarray"):
                    normal_expr = normal_expr.toarray().flatten()

                gene_data["malignant_mean"] = np.mean(mal_expr)
                gene_data["normal_mean"] = np.mean(normal_expr)
                gene_data["fold_change"] = gene_data["malignant_mean"] / (
                    gene_data["normal_mean"] + 1e-6
                )
                gene_data["pct_expressed_malignant"] = np.mean(mal_expr > min_expression)
            else:
                # Overall expression
                expr = adata[:, gene].X
                if hasattr(expr, "toarray"):
                    expr = expr.toarray().flatten()

                gene_data["mean_expression"] = np.mean(expr)
                gene_data["pct_expressed"] = np.mean(expr > min_expression)

            targets.append(gene_data)

        self.targets_ = pd.DataFrame(targets)

        # Calculate priority score
        if malignant_key in adata.obs.columns and "fold_change" in self.targets_.columns:
            self.targets_["priority_score"] = (
                np.log1p(self.targets_["fold_change"]) * self.targets_["pct_expressed_malignant"]
            )
            self.targets_ = self.targets_.sort_values("priority_score", ascending=False)

        return self.targets_

    def _get_gene_category(self, gene: str) -> List[str]:
        """Get category for a gene."""
        categories = []
        for category, info in self.druggable_categories.items():
            if gene in info["genes"]:
                categories.append(category)
        return ",".join(categories) if categories else "unknown"

    def prioritize_druggable_genes(
        self,
        adata: AnnData,
        n_top: int = 50,
        category_filter: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        Prioritize druggable genes by expression and specificity.

        Parameters
        ----------
        adata : AnnData
            Expression data
        n_top : int
            Number of top genes to return
        category_filter : str, optional
            Filter by category

        Returns:
        -------
        pd.DataFrame
            Prioritized genes
        """
        if self.targets_ is None:
            self.discover_targets(adata)

        targets = self.targets_.copy()

        if category_filter:
            targets = targets[targets["category"].str.contains(category_filter, na=False)]

        return targets.head(n_top)

    def find_drug_combinations(
        self,
        adata: AnnData,
        resistance_genes: List[str],
        n_drugs: int = 3,
    ) -> pd.DataFrame:
        """
        Find drug combinations to overcome resistance.

        Parameters
        ----------
        adata : AnnData
            Expression data
        resistance_genes : list
            Genes conferring resistance
        n_drugs : int
            Number of drugs in combination

        Returns:
        -------
        pd.DataFrame
            Suggested drug combinations
        """
        combinations = []

        # Find drugs targeting resistance mechanisms
        for drug_name, drug_info in DRUG_TARGETS.items():
            overlap = set(drug_info["targets"]) & set(resistance_genes)
            if len(overlap) > 0:
                combinations.append(
                    {
                        "drug": drug_name,
                        "targets": ",".join(drug_info["targets"]),
                        "category": drug_info["category"],
                        "overlap_with_resistance": len(overlap),
                        "resistance_targets": ",".join(overlap),
                    }
                )

        return pd.DataFrame(combinations).sort_values("overlap_with_resistance", ascending=False)


def discover_therapeutic_targets(
    adata: AnnData,
    malignant_key: str = "is_malignant",
    n_top: int = 50,
    key_added: str = "targets",
) -> pd.DataFrame:
    """
    Discover therapeutic targets in tumor cells.

    Parameters
    ----------
    adata : AnnData
        Expression data
    malignant_key : str
        Column indicating malignant cells
    n_top : int
        Number of top targets to return
    key_added : str
        Key for storing results

    Returns:
    -------
    pd.DataFrame
        Top therapeutic targets
    """
    discoverer = TargetDiscovery()
    targets = discoverer.discover_targets(adata, malignant_key=malignant_key)

    top_targets = targets.head(n_top)

    log.info(f"Discovered {len(top_targets)} therapeutic targets")

    return top_targets


def prioritize_druggable_genes(
    adata: AnnData,
    category: Optional[str] = None,
    n_top: int = 50,
) -> pd.DataFrame:
    """
    Prioritize druggable genes by expression.

    Parameters
    ----------
    adata : AnnData
        Expression data
    category : str, optional
        Filter by category ("kinase", "receptor", "immune_checkpoint", etc.)
    n_top : int
        Number of top genes

    Returns:
    -------
    pd.DataFrame
        Prioritized genes
    """
    discoverer = TargetDiscovery()
    prioritized = discoverer.prioritize_druggable_genes(
        adata, n_top=n_top, category_filter=category
    )

    return prioritized


def suggest_targeted_therapies(
    adata: AnnData,
    groupby: Optional[str] = None,
    indication: Optional[str] = None,
) -> pd.DataFrame:
    """
    Suggest targeted therapies based on molecular profile.

    Parameters
    ----------
    adata : AnnData
        Expression data
    groupby : str, optional
        Column to group cells
    indication : str, optional
        Cancer type indication

    Returns:
    -------
    pd.DataFrame
        Suggested therapies
    """
    suggestions = []

    for drug_name, drug_info in DRUG_TARGETS.items():
        # Filter by indication if specified
        if indication and indication not in drug_info["indications"]:
            continue

        # Check target expression
        targets = drug_info["targets"]
        available_targets = [t for t in targets if t in adata.var_names]

        if len(available_targets) == 0:
            continue

        # Calculate mean expression of targets
        expr = adata[:, available_targets].X.mean(axis=1)
        if hasattr(expr, "toarray"):
            expr = expr.toarray().flatten()

        suggestions.append(
            {
                "drug": drug_name,
                "targets": ",".join(targets),
                "category": drug_info["category"],
                "indications": ",".join(drug_info["indications"]),
                "mean_target_expression": np.mean(expr),
                "pct_cells_with_expression": np.mean(expr > 0.1),
            }
        )

    suggestions_df = pd.DataFrame(suggestions)

    if len(suggestions_df) > 0:
        suggestions_df = suggestions_df.sort_values("mean_target_expression", ascending=False)

    return suggestions_df
