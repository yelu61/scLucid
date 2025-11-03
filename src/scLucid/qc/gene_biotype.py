"""
Gene biotype annotation and filtering for single-cell RNA-seq data.

This module provides utilities to annotate genes with their biotypes
(protein-coding, lncRNA, etc.) and enables biotype-aware QC and analysis.
"""

import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Set, Union

import numpy as np
import pandas as pd
from anndata import AnnData

log = logging.getLogger(__name__)

__all__ = [
    "annotate_gene_biotypes",
    "filter_genes_by_biotype",
    "get_biotype_statistics",
    "recommend_biotype_strategy",
]


# --- Predefined Biotype Categories ---
BIOTYPE_CATEGORIES = {
    "protein_coding": {
        "name": "Protein Coding",
        "recommended_for_analysis": True,
        "description": "Genes that code for proteins",
    },
    "lncRNA": {
        "name": "Long Non-coding RNA",
        "recommended_for_analysis": False,
        "description": "Long non-coding RNAs (>200 nucleotides)",
        "aliases": ["lincRNA", "processed_transcript", "antisense"],
    },
    "pseudogene": {
        "name": "Pseudogene",
        "recommended_for_analysis": False,
        "description": "Non-functional gene copies",
        "aliases": [
            "processed_pseudogene",
            "unprocessed_pseudogene",
            "transcribed_pseudogene",
        ],
    },
    "miRNA": {
        "name": "microRNA",
        "recommended_for_analysis": False,
        "description": "MicroRNA genes",
    },
    "snoRNA": {
        "name": "Small Nucleolar RNA",
        "recommended_for_analysis": False,
        "description": "Small nucleolar RNAs",
    },
    "snRNA": {
        "name": "Small Nuclear RNA",
        "recommended_for_analysis": False,
        "description": "Small nuclear RNAs",
    },
    "rRNA": {
        "name": "Ribosomal RNA",
        "recommended_for_analysis": False,
        "description": "Ribosomal RNA genes",
    },
    "Mt_tRNA": {
        "name": "Mitochondrial tRNA",
        "recommended_for_analysis": False,
        "description": "Mitochondrial transfer RNAs",
    },
    "Mt_rRNA": {
        "name": "Mitochondrial rRNA",
        "recommended_for_analysis": False,
        "description": "Mitochondrial ribosomal RNAs",
    },
    "IG_gene": {
        "name": "Immunoglobulin Gene",
        "recommended_for_analysis": True,  # Important for immune cells
        "description": "Immunoglobulin genes (IG_C, IG_V, IG_J, IG_D)",
        "aliases": ["IG_C_gene", "IG_V_gene", "IG_J_gene", "IG_D_gene"],
    },
    "TR_gene": {
        "name": "T-cell Receptor Gene",
        "recommended_for_analysis": True,  # Important for T cells
        "description": "T-cell receptor genes (TR_C, TR_V, TR_J, TR_D)",
        "aliases": ["TR_C_gene", "TR_V_gene", "TR_J_gene", "TR_D_gene"],
    },
}


# --- Helper Functions ---
def _load_ensembl_biotypes(
    species: Literal["human", "mouse", "rat"] = "human",
    ensembl_version: Optional[int] = None,
) -> pd.DataFrame:
    """
    Load gene biotype annotations from Ensembl biomart.

    This function attempts to load from local cache first, then downloads
    from Ensembl if needed.

    Args:
        species: Species name
        ensembl_version: Ensembl version (default: latest)

    Returns:
        DataFrame with columns: gene_id, gene_name, biotype, chromosome
    """
    from io import StringIO

    try:
        import requests
    except ImportError:
        log.error(
            "Package 'requests' is required for downloading Ensembl annotations. "
            "Install with: pip install requests"
        )
        raise

    # Map species to Ensembl dataset names
    dataset_map = {
        "human": "hsapiens_gene_ensembl",
        "mouse": "mmusculus_gene_ensembl",
        "rat": "rnorvegicus_gene_ensembl",
    }

    if species not in dataset_map:
        raise ValueError(
            f"Species '{species}' not supported. Choose from: {list(dataset_map.keys())}"
        )

    dataset = dataset_map[species]

    # Check local cache first
    cache_dir = Path.home() / ".sclucid" / "gene_annotations"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{species}_biotypes_ensembl.csv"

    if cache_file.exists():
        log.info(f"Loading gene biotypes from cache: {cache_file}")
        return pd.read_csv(cache_file)

    # Download from Ensembl BioMart
    log.info(f"Downloading gene biotype annotations for {species} from Ensembl...")

    # Ensembl BioMart REST API
    server = "http://www.ensembl.org"
    if ensembl_version:
        server = f"http://{ensembl_version}.ensembl.org"

    # BioMart XML query
    query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "{dataset}" interface = "default" >
        <Attribute name = "ensembl_gene_id" />
        <Attribute name = "external_gene_name" />
        <Attribute name = "gene_biotype" />
        <Attribute name = "chromosome_name" />
        <Attribute name = "start_position" />
        <Attribute name = "end_position" />
    </Dataset>
</Query>"""

    try:
        response = requests.get(
            f"{server}/biomart/martservice",
            params={"query": query},
            timeout=300,  # 5 minutes timeout
        )
        response.raise_for_status()

        # Parse response
        df = pd.read_csv(StringIO(response.text), sep="\t")
        df.columns = [
            "gene_id",
            "gene_name",
            "biotype",
            "chromosome",
            "start",
            "end",
        ]

        # Remove genes without names
        df = df[df["gene_name"].notna() & (df["gene_name"] != "")]

        # Cache the result
        df.to_csv(cache_file, index=False)
        log.info(f"Downloaded and cached {len(df)} gene annotations")

        return df

    except requests.exceptions.RequestException as e:
        log.error(f"Failed to download Ensembl annotations: {e}")
        raise RuntimeError(
            "Could not download gene annotations. Please check your internet connection "
            "or provide a custom biotype annotation file."
        )


def _match_genes_to_biotypes(
    gene_names: pd.Index, biotype_df: pd.DataFrame, fuzzy_match: bool = True
) -> pd.Series:
    """
    Match gene names to their biotypes.

    Args:
        gene_names: Gene names from AnnData.var_names
        biotype_df: DataFrame from _load_ensembl_biotypes
        fuzzy_match: Whether to try fuzzy matching for unmatched genes

    Returns:
        Series mapping gene names to biotypes
    """
    # Direct match on gene_name
    gene_to_biotype = biotype_df.set_index("gene_name")["biotype"].to_dict()

    biotypes = gene_names.map(gene_to_biotype)

    # Try matching on gene_id for genes with Ensembl IDs
    if fuzzy_match:
        unmatched = biotypes.isna()
        if unmatched.sum() > 0:
            log.info(
                f"Trying Ensembl ID matching for {unmatched.sum()} unmatched genes..."
            )

            gene_id_to_biotype = biotype_df.set_index("gene_id")["biotype"].to_dict()
            biotypes[unmatched] = gene_names[unmatched].map(gene_id_to_biotype)

    matched_count = (~biotypes.isna()).sum()
    match_rate = matched_count / len(gene_names) * 100

    log.info(f"Matched {matched_count}/{len(gene_names)} genes ({match_rate:.1f}%)")

    if match_rate < 50:
        log.warning(
            f"⚠️  Low match rate ({match_rate:.1f}%). "
            "Check if gene symbols match the species annotation. "
            "Common issues: mouse genes should be capitalized (Gapdh), "
            "human genes are all caps (GAPDH)."
        )

    return biotypes


def _categorize_biotype(biotype: str) -> str:
    """
    Map detailed Ensembl biotypes to broader categories.

    Args:
        biotype: Original Ensembl biotype

    Returns:
        Standardized category name
    """
    if pd.isna(biotype):
        return "unknown"

    biotype_lower = biotype.lower()

    # Check each category
    for category, info in BIOTYPE_CATEGORIES.items():
        # Exact match
        if biotype_lower == category.lower():
            return category

        # Check aliases
        if "aliases" in info:
            for alias in info["aliases"]:
                if biotype_lower == alias.lower() or alias.lower() in biotype_lower:
                    return category

    # Special handling for common patterns
    if "pseudogene" in biotype_lower:
        return "pseudogene"
    elif "lincrna" in biotype_lower or "lncrna" in biotype_lower:
        return "lncRNA"
    elif biotype_lower.startswith("ig_"):
        return "IG_gene"
    elif biotype_lower.startswith("tr_"):
        return "TR_gene"
    elif "protein_coding" in biotype_lower:
        return "protein_coding"

    return "other"


# --- Main Functions ---
def annotate_gene_biotypes(
    adata: AnnData,
    species: Literal["human", "mouse", "rat"] = "human",
    biotype_df: Optional[pd.DataFrame] = None,
    method: Literal["ensembl", "custom"] = "ensembl",
    fuzzy_match: bool = True,
    overwrite: bool = False,
) -> AnnData:
    """
    Annotate genes with biotype information.

    This function adds several columns to adata.var:
    - 'biotype': Detailed Ensembl biotype
    - 'biotype_category': Broad category (protein_coding, lncRNA, etc.)
    - 'recommended_for_analysis': Boolean flag based on biotype

    Args:
        adata: AnnData object
        species: Species for annotation
        biotype_df: Custom biotype DataFrame (gene_name, biotype columns required)
        method: Annotation method ('ensembl' or 'custom')
        fuzzy_match: Try to match unmatched genes with Ensembl IDs
        overwrite: Overwrite existing biotype annotations

    Returns:
        AnnData with annotated .var

    Examples:
        >>> # Basic usage with Ensembl
        >>> adata = annotate_gene_biotypes(adata, species='human')
        >>>
        >>> # Using custom annotations
        >>> custom_df = pd.DataFrame({
        ...     'gene_name': ['GAPDH', 'MALAT1', 'MT-CO1'],
        ...     'biotype': ['protein_coding', 'lncRNA', 'protein_coding']
        ... })
        >>> adata = annotate_gene_biotypes(adata, biotype_df=custom_df, method='custom')
    """
    if "biotype" in adata.var.columns and not overwrite:
        log.info(
            "Gene biotypes already annotated. Use overwrite=True to re-annotate."
        )
        return adata

    log.info(f"Annotating gene biotypes for {adata.n_vars} genes...")

    # Load biotype reference
    if method == "ensembl":
        if biotype_df is None:
            biotype_df = _load_ensembl_biotypes(species)
        else:
            log.info("Using provided biotype_df with ensembl method")
    elif method == "custom":
        if biotype_df is None:
            raise ValueError("biotype_df must be provided when method='custom'")

        required_cols = ["gene_name", "biotype"]
        missing_cols = [col for col in required_cols if col not in biotype_df.columns]
        if missing_cols:
            raise ValueError(
                f"Custom biotype_df missing required columns: {missing_cols}"
            )
    else:
        raise ValueError(f"Unknown method: {method}")

    # Match genes to biotypes
    biotypes = _match_genes_to_biotypes(adata.var_names, biotype_df, fuzzy_match)

    # Add to adata.var
    adata.var["biotype"] = biotypes
    adata.var["biotype_category"] = biotypes.apply(_categorize_biotype)

    # Add recommendation flag
    category_recommendations = {
        cat: info["recommended_for_analysis"]
        for cat, info in BIOTYPE_CATEGORIES.items()
    }
    adata.var["recommended_for_analysis"] = adata.var["biotype_category"].map(
        category_recommendations
    )
    adata.var["recommended_for_analysis"].fillna(False, inplace=True)

    # Log statistics
    biotype_counts = adata.var["biotype_category"].value_counts()
    log.info("\nGene biotype distribution:")
    for biotype, count in biotype_counts.items():
        pct = count / adata.n_vars * 100
        recommended = "✓" if category_recommendations.get(biotype, False) else "✗"
        log.info(f"  {recommended} {biotype}: {count} ({pct:.1f}%)")

    # Store metadata
    adata.uns.setdefault("sclucid", {}).setdefault("qc", {})["gene_biotypes"] = {
        "method": method,
        "species": species,
        "n_genes_annotated": (~adata.var["biotype"].isna()).sum(),
        "n_protein_coding": (adata.var["biotype_category"] == "protein_coding").sum(),
        "annotation_rate": (~adata.var["biotype"].isna()).sum() / adata.n_vars,
    }

    return adata


def filter_genes_by_biotype(
    adata: AnnData,
    keep_biotypes: Optional[List[str]] = None,
    use_recommended: bool = True,
    copy: bool = False,
) -> Optional[AnnData]:
    """
    Filter genes based on biotype annotations.

    Args:
        adata: AnnData object with biotype annotations
        keep_biotypes: List of biotype categories to keep (overrides use_recommended)
        use_recommended: If True, keep only recommended biotypes
        copy: Return a copy instead of filtering in-place

    Returns:
        Filtered AnnData if copy=True, else None

    Examples:
        >>> # Keep only protein-coding genes
        >>> adata_pc = filter_genes_by_biotype(adata, keep_biotypes=['protein_coding'], copy=True)
        >>>
        >>> # Keep recommended biotypes (protein-coding, IG, TR genes)
        >>> adata_filtered = filter_genes_by_biotype(adata, use_recommended=True, copy=True)
        >>>
        >>> # Custom selection
        >>> adata_custom = filter_genes_by_biotype(
        ...     adata,
        ...     keep_biotypes=['protein_coding', 'lncRNA'],
        ...     copy=True
        ... )
    """
    if "biotype_category" not in adata.var.columns:
        raise ValueError(
            "Gene biotypes not found. Run annotate_gene_biotypes() first."
        )

    if copy:
        adata = adata.copy()

    initial_genes = adata.n_vars

    if keep_biotypes is not None:
        # Use custom biotype list
        keep_mask = adata.var["biotype_category"].isin(keep_biotypes)
        strategy = f"custom: {', '.join(keep_biotypes)}"
    elif use_recommended:
        # Use recommended biotypes
        keep_mask = adata.var["recommended_for_analysis"] == True
        strategy = "recommended biotypes"
    else:
        log.warning("No filtering criteria specified. Returning original data.")
        return adata if copy else None

    # Apply filter
    adata._inplace_subset_var(keep_mask)

    removed = initial_genes - adata.n_vars
    log.info(
        f"Filtered by biotype ({strategy}): "
        f"kept {adata.n_vars}/{initial_genes} genes, "
        f"removed {removed} ({removed/initial_genes:.1%})"
    )

    # Log what was kept
    kept_biotypes = adata.var["biotype_category"].value_counts()
    log.info("Remaining biotypes:")
    for biotype, count in kept_biotypes.items():
        log.info(f"  - {biotype}: {count}")

    if copy:
        return adata


def get_biotype_statistics(
    adata: AnnData, sample_key: Optional[str] = None
) -> pd.DataFrame:
    """
    Generate comprehensive biotype statistics.

    Args:
        adata: AnnData object with biotype annotations
        sample_key: If provided, calculate per-sample statistics

    Returns:
        DataFrame with biotype statistics
    """
    if "biotype_category" not in adata.var.columns:
        raise ValueError(
            "Gene biotypes not found. Run annotate_gene_biotypes() first."
        )

    stats = []

    # Global statistics
    for biotype in adata.var["biotype_category"].unique():
        if pd.isna(biotype):
            continue

        biotype_mask = adata.var["biotype_category"] == biotype
        biotype_genes = adata.var_names[biotype_mask]

        # Calculate expression statistics
        expr_data = adata[:, biotype_genes].X
        if hasattr(expr_data, "toarray"):
            expr_data = expr_data.toarray()

        stats.append(
            {
                "sample": "global",
                "biotype": biotype,
                "n_genes": len(biotype_genes),
                "pct_genes": len(biotype_genes) / adata.n_vars * 100,
                "mean_expression": np.mean(expr_data),
                "median_expression": np.median(expr_data),
                "detection_rate": (expr_data > 0).sum() / expr_data.size * 100,
            }
        )

    # Per-sample statistics
    if sample_key and sample_key in adata.obs.columns:
        for sample in adata.obs[sample_key].unique():
            sample_mask = adata.obs[sample_key] == sample

            for biotype in adata.var["biotype_category"].unique():
                if pd.isna(biotype):
                    continue

                biotype_mask = adata.var["biotype_category"] == biotype
                biotype_genes = adata.var_names[biotype_mask]

                expr_data = adata[sample_mask, :][:, biotype_genes].X
                if hasattr(expr_data, "toarray"):
                    expr_data = expr_data.toarray()

                stats.append(
                    {
                        "sample": sample,
                        "biotype": biotype,
                        "n_genes": len(biotype_genes),
                        "pct_genes": len(biotype_genes) / adata.n_vars * 100,
                        "mean_expression": np.mean(expr_data),
                        "median_expression": np.median(expr_data),
                        "detection_rate": (expr_data > 0).sum() / expr_data.size * 100,
                    }
                )

    return pd.DataFrame(stats)


def recommend_biotype_strategy(
    adata: AnnData,
    analysis_goal: Literal[
        "cell_typing", "differential_expression", "trajectory", "general"
    ] = "general",
) -> Dict[str, Union[List[str], str]]:
    """
    Recommend biotype filtering strategy based on analysis goals.

    Args:
        adata: AnnData object with biotype annotations
        analysis_goal: Type of downstream analysis

    Returns:
        Dictionary with recommendations
    """
    if "biotype_category" not in adata.var.columns:
        raise ValueError(
            "Gene biotypes not found. Run annotate_gene_biotypes() first."
        )

    recommendations = {
        "cell_typing": {
            "keep_biotypes": ["protein_coding", "IG_gene", "TR_gene"],
            "rationale": (
                "Cell type markers are primarily protein-coding genes. "
                "IG and TR genes are critical for immune cell identification."
            ),
            "expected_gene_retention": 0.85,
        },
        "differential_expression": {
            "keep_biotypes": ["protein_coding", "lncRNA"],
            "rationale": (
                "Include lncRNAs as they can be differentially expressed "
                "and have regulatory functions."
            ),
            "expected_gene_retention": 0.90,
        },
        "trajectory": {
            "keep_biotypes": ["protein_coding"],
            "rationale": (
                "Trajectory inference benefits from stable, high-signal genes. "
                "Protein-coding genes provide the most reliable temporal dynamics."
            ),
            "expected_gene_retention": 0.80,
        },
        "general": {
            "keep_biotypes": ["protein_coding", "IG_gene", "TR_gene"],
            "rationale": (
                "Standard approach for most analyses. "
                "Balances signal quality with biological completeness."
            ),
            "expected_gene_retention": 0.85,
        },
    }

    rec = recommendations[analysis_goal]

    # Calculate actual retention rate
    keep_mask = adata.var["biotype_category"].isin(rec["keep_biotypes"])
    actual_retention = keep_mask.sum() / len(keep_mask)

    rec["actual_gene_retention"] = actual_retention
    rec["n_genes_kept"] = keep_mask.sum()
    rec["n_genes_removed"] = (~keep_mask).sum()

    log.info(f"\n{'='*60}")
    log.info(f"BIOTYPE FILTERING RECOMMENDATION: {analysis_goal.upper()}")
    log.info(f"{'='*60}")
    log.info(f"\nRationale: {rec['rationale']}")
    log.info(f"\nRecommended biotypes to keep: {', '.join(rec['keep_biotypes'])}")
    log.info(
        f"Genes retained: {rec['n_genes_kept']}/{adata.n_vars} ({actual_retention:.1%})"
    )
    log.info(
        f"Genes removed: {rec['n_genes_removed']} ({1-actual_retention:.1%})\n"
    )

    # Show what will be removed
    removed_biotypes = adata.var.loc[~keep_mask, "biotype_category"].value_counts()
    if len(removed_biotypes) > 0:
        log.info("Biotypes to be removed:")
        for biotype, count in removed_biotypes.items():
            log.info(f"  - {biotype}: {count} genes")

    return rec