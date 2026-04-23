"""
Database interfaces for tumor analysis.

This module provides interfaces to cancer databases like
COSMIC, TCGA, and other cancer genomics resources.
"""

import logging
from typing import Optional

import pandas as pd

log = logging.getLogger(__name__)


# Cancer Gene Census - curated list of cancer genes
CANCER_GENE_CENSUS = {
    "ABL1": {"role": "oncogene", "cancer_types": ["CML", "ALL"]},
    "AKT1": {"role": "oncogene", "cancer_types": ["breast", "colorectal"]},
    "ALK": {"role": "oncogene", "cancer_types": ["lung", "neuroblastoma"]},
    "APC": {"role": "TSG", "cancer_types": ["colorectal"]},
    "ATM": {"role": "TSG", "cancer_types": ["leukemia", "lymphoma"]},
    "BCL2": {"role": "oncogene", "cancer_types": ["lymphoma"]},
    "BRAF": {"role": "oncogene", "cancer_types": ["melanoma", "colorectal", "lung", "thyroid"]},
    "BRCA1": {"role": "TSG", "cancer_types": ["breast", "ovarian"]},
    "BRCA2": {"role": "TSG", "cancer_types": ["breast", "ovarian", "pancreatic"]},
    "CDH1": {"role": "TSG", "cancer_types": ["gastric", "breast"]},
    "CDKN2A": {"role": "TSG", "cancer_types": ["melanoma", "pancreatic", "lung"]},
    "CTNNB1": {"role": "oncogene", "cancer_types": ["colorectal", "liver"]},
    "EGFR": {"role": "oncogene", "cancer_types": ["lung", "glioblastoma"]},
    "ERBB2": {"role": "oncogene", "cancer_types": ["breast", "gastric"]},
    "FBXW7": {"role": "TSG", "cancer_types": ["colorectal", "endometrial"]},
    "FGFR3": {"role": "oncogene", "cancer_types": ["bladder", "myeloma"]},
    "FLT3": {"role": "oncogene", "cancer_types": ["AML"]},
    "GNAQ": {"role": "oncogene", "cancer_types": ["melanoma"]},
    "GNAS": {"role": "oncogene", "cancer_types": ["pituitary", "pancreatic"]},
    "HRAS": {"role": "oncogene", "cancer_types": ["bladder", "thyroid"]},
    "IDH1": {"role": "oncogene", "cancer_types": ["glioma", "AML"]},
    "IDH2": {"role": "oncogene", "cancer_types": ["glioma", "AML"]},
    "JAK2": {"role": "oncogene", "cancer_types": ["myeloproliferative"]},
    "KDR": {"role": "oncogene", "cancer_types": ["angiogenesis"]},
    "KIT": {"role": "oncogene", "cancer_types": ["GIST", "melanoma"]},
    "KRAS": {"role": "oncogene", "cancer_types": ["pancreatic", "colorectal", "lung"]},
    "MAP2K1": {"role": "oncogene", "cancer_types": ["melanoma", "lung"]},
    "MPL": {"role": "oncogene", "cancer_types": ["myeloproliferative"]},
    "MYC": {"role": "oncogene", "cancer_types": ["lymphoma", "Burkitt"]},
    "NF1": {"role": "TSG", "cancer_types": ["neurofibromatosis"]},
    "NF2": {"role": "TSG", "cancer_types": ["schwannoma"]},
    "NOTCH1": {"role": "oncogene", "cancer_types": ["T-ALL"]},
    "NPM1": {"role": "oncogene", "cancer_types": ["AML"]},
    "NRAS": {"role": "oncogene", "cancer_types": ["melanoma", "AML"]},
    "PDGFRA": {"role": "oncogene", "cancer_types": ["GIST"]},
    "PIK3CA": {"role": "oncogene", "cancer_types": ["breast", "colorectal", "endometrial"]},
    "PTEN": {"role": "TSG", "cancer_types": ["glioblastoma", "endometrial", "prostate"]},
    "RB1": {"role": "TSG", "cancer_types": ["retinoblastoma", "osteosarcoma"]},
    "RET": {"role": "oncogene", "cancer_types": ["thyroid"]},
    "SMAD4": {"role": "TSG", "cancer_types": ["pancreatic", "colorectal"]},
    "SMARCB1": {"role": "TSG", "cancer_types": ["rhabdoid"]},
    "SRC": {"role": "oncogene", "cancer_types": ["colon"]},
    "STK11": {"role": "TSG", "cancer_types": ["Peutz-Jeghers"]},
    "TERT": {"role": "oncogene", "cancer_types": ["melanoma", "glioblastoma"]},
    "TP53": {"role": "TSG", "cancer_types": ["multiple"]},
    "VHL": {"role": "TSG", "cancer_types": ["renal"]},
}


# Drug target information
DRUG_TARGETS_DB = {
    "imatinib": {
        "targets": ["ABL1", "KIT", "PDGFRA"],
        "cancer_types": ["CML", "ALL", "GIST"],
        "mechanism": "kinase_inhibitor",
    },
    "trastuzumab": {
        "targets": ["ERBB2"],
        "cancer_types": ["breast", "gastric"],
        "mechanism": "mAb",
    },
    "bevacizumab": {
        "targets": ["VEGFA"],
        "cancer_types": ["colorectal", "lung", "renal", "glioblastoma"],
        "mechanism": "mAb",
    },
    "cetuximab": {
        "targets": ["EGFR"],
        "cancer_types": ["colorectal", "HNSCC"],
        "mechanism": "mAb",
    },
    "gefitinib": {
        "targets": ["EGFR"],
        "cancer_types": ["lung"],
        "mechanism": "TKI",
    },
    "vemurafenib": {
        "targets": ["BRAF"],
        "cancer_types": ["melanoma"],
        "mechanism": "kinase_inhibitor",
    },
    "olaparib": {
        "targets": ["PARP1", "PARP2"],
        "cancer_types": ["ovarian", "breast", "pancreatic", "prostate"],
        "mechanism": "PARP_inhibitor",
    },
    "pembrolizumab": {
        "targets": ["PDCD1"],
        "cancer_types": ["multiple"],
        "mechanism": "checkpoint_inhibitor",
    },
    "nivolumab": {
        "targets": ["PDCD1"],
        "cancer_types": ["multiple"],
        "mechanism": "checkpoint_inhibitor",
    },
}


def query_cancer_gene_census(
    gene: Optional[str] = None,
    cancer_type: Optional[str] = None,
    role: Optional[str] = None,
) -> pd.DataFrame:
    """
    Query the Cancer Gene Census database.

    Parameters
    ----------
    gene : str, optional
        Gene symbol to query
    cancer_type : str, optional
        Filter by cancer type
    role : str, optional
        Filter by gene role ("oncogene" or "TSG")

    Returns:
    -------
    pd.DataFrame
        Query results
    """
    results = []

    for gene_symbol, info in CANCER_GENE_CENSUS.items():
        # Filter by gene
        if gene and gene.upper() != gene_symbol.upper():
            continue

        # Filter by cancer type
        if cancer_type and cancer_type.lower() not in [ct.lower() for ct in info["cancer_types"]]:
            continue

        # Filter by role
        if role and info["role"] != role:
            continue

        results.append(
            {
                "gene": gene_symbol,
                "role": info["role"],
                "cancer_types": ", ".join(info["cancer_types"]),
            }
        )

    return pd.DataFrame(results)


def get_drug_targets(
    drug: Optional[str] = None,
    target: Optional[str] = None,
    cancer_type: Optional[str] = None,
) -> pd.DataFrame:
    """
    Query drug-target information.

    Parameters
    ----------
    drug : str, optional
        Drug name to query
    target : str, optional
        Target gene to query
    cancer_type : str, optional
        Filter by cancer type

    Returns:
    -------
    pd.DataFrame
        Drug-target information
    """
    results = []

    for drug_name, info in DRUG_TARGETS_DB.items():
        # Filter by drug
        if drug and drug.lower() != drug_name.lower():
            continue

        # Filter by target
        if target and target.upper() not in [t.upper() for t in info["targets"]]:
            continue

        # Filter by cancer type
        if cancer_type and cancer_type.lower() not in [ct.lower() for ct in info["cancer_types"]]:
            continue

        results.append(
            {
                "drug": drug_name,
                "targets": ", ".join(info["targets"]),
                "cancer_types": ", ".join(info["cancer_types"]),
                "mechanism": info["mechanism"],
            }
        )

    return pd.DataFrame(results)


def query_tcga_data(
    gene: str,
    cancer_type: str,
    data_type: str = "expression",
) -> pd.DataFrame:
    """
    Query TCGA data (placeholder - requires actual TCGA API access).

    Parameters
    ----------
    gene : str
        Gene symbol
    cancer_type : str
        TCGA cancer type abbreviation (e.g., "BRCA", "LUAD")
    data_type : str
        Type of data ("expression", "mutation", "cnv")

    Returns:
    -------
    pd.DataFrame
        TCGA data (placeholder)
    """
    log.warning("TCGA query requires external data source - returning placeholder")

    # Placeholder return
    return pd.DataFrame(
        {
            "gene": [gene],
            "cancer_type": [cancer_type],
            "data_type": [data_type],
            "note": ["External data source required"],
        }
    )


def is_cancer_gene(
    gene: str,
    cancer_type: Optional[str] = None,
) -> bool:
    """
    Check if a gene is a known cancer gene.

    Parameters
    ----------
    gene : str
        Gene symbol
    cancer_type : str, optional
        Specific cancer type

    Returns:
    -------
    bool
        Whether gene is a known cancer gene
    """
    if gene.upper() in CANCER_GENE_CENSUS:
        if cancer_type:
            info = CANCER_GENE_CENSUS[gene.upper()]
            return cancer_type.lower() in [ct.lower() for ct in info["cancer_types"]]
        return True
    return False


def get_gene_role(
    gene: str,
) -> Optional[str]:
    """
    Get the role of a cancer gene (oncogene or TSG).

    Parameters
    ----------
    gene : str
        Gene symbol

    Returns:
    -------
    str or None
        Gene role or None if not in census
    """
    gene_upper = gene.upper()
    if gene_upper in CANCER_GENE_CENSUS:
        return CANCER_GENE_CENSUS[gene_upper]["role"]
    return None
