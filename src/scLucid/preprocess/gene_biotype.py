"""
Gene biotype annotation and filtering for single-cell RNA-seq preprocessing.

This module provides utilities to annotate genes with their biotypes
(protein-coding, lncRNA, etc.) and enables biotype-aware feature filtering
before downstream analysis.
"""

import logging
import tempfile
import time
from functools import wraps
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
from anndata import AnnData

from ..utils.resource_loader import get_resource_path, resource_exists

log = logging.getLogger(__name__)

__all__ = [
    "annotate_gene_biotypes",
    "apply_gene_biotype_strategy",
    "filter_genes_by_biotype",
    "get_biotype_statistics",
    "recommend_biotype_strategy",
    "load_gene_biotypes",
    "list_gene_biotype_resources",
    "get_gene_biotype_cache_dir",
]

GENE_BIOTYPE_RESOURCE_SUBDIR = "gene_biotypes"
GENE_BIOTYPE_RESOURCE_MANIFEST = f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/manifest.json"
SUPPORTED_GENE_BIOTYPE_SPECIES = ("human", "mouse", "rat")
GENE_BIOTYPE_METHODS = ("reference", "ensembl", "custom")


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
def retry(retries=3, delay=5):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            for i in range(retries):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    if i == retries - 1:
                        raise e
                    log.warning(f"Attempt {i+1} failed: {e}. Retrying in {delay}s...")
                    time.sleep(delay)

        return wrapper

    return decorator


def get_gene_biotype_cache_dir() -> Path:
    """Return the user-local cache directory for gene biotype references."""
    cache_dir = Path.home() / ".sclucid" / "gene_annotations"
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        return cache_dir
    except PermissionError:
        fallback = Path(tempfile.gettempdir()) / "sclucid_gene_annotations"
        fallback.mkdir(parents=True, exist_ok=True)
        log.warning(
            f"Gene biotype cache directory '{cache_dir}' is not writable; using temporary fallback '{fallback}'."
        )
        return fallback


def _candidate_biotype_resource_names(
    species: str,
    ensembl_version: Optional[int] = None,
) -> List[str]:
    names: List[str] = []
    if ensembl_version is not None:
        names.extend(
            [
                f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_reference_v{ensembl_version}.csv.gz",
                f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_reference_v{ensembl_version}.csv",
                f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_ensembl_v{ensembl_version}.csv.gz",
                f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_ensembl_v{ensembl_version}.csv",
            ]
        )
    names.extend(
        [
            f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_reference_latest.csv.gz",
            f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_reference_latest.csv",
            f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_ensembl_latest.csv.gz",
            f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_ensembl_latest.csv",
            f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_biotypes_ensembl.csv.gz",
            f"{GENE_BIOTYPE_RESOURCE_SUBDIR}/{species}_biotypes_ensembl.csv",
        ]
    )
    return names


def _candidate_biotype_cache_paths(
    species: str,
    ensembl_version: Optional[int] = None,
    cache_dir: Optional[Path] = None,
) -> List[Path]:
    cache_root = cache_dir or get_gene_biotype_cache_dir()
    paths: List[Path] = []
    if ensembl_version is not None:
        paths.extend(
            [
                cache_root / f"{species}_reference_v{ensembl_version}.csv.gz",
                cache_root / f"{species}_reference_v{ensembl_version}.csv",
                cache_root / f"{species}_ensembl_v{ensembl_version}.csv.gz",
                cache_root / f"{species}_ensembl_v{ensembl_version}.csv",
            ]
        )
    paths.extend(
        [
            cache_root / f"{species}_reference_latest.csv.gz",
            cache_root / f"{species}_reference_latest.csv",
            cache_root / f"{species}_ensembl_latest.csv.gz",
            cache_root / f"{species}_ensembl_latest.csv",
            cache_root / f"{species}_biotypes_ensembl.csv.gz",
            cache_root / f"{species}_biotypes_ensembl.csv",
        ]
    )
    return paths


def _read_biotype_table(path: Path) -> pd.DataFrame:
    """Load and validate a gene biotype reference table."""
    compression = "gzip" if path.suffix == ".gz" else "infer"
    df = pd.read_csv(path, compression=compression)
    required = {"gene_name", "biotype"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Gene biotype resource '{path}' is missing required columns: {sorted(missing)}"
        )
    return df


def _read_gene_biotype_manifest() -> Dict[str, dict]:
    """Return bundled gene biotype manifest data when available."""
    if not resource_exists(GENE_BIOTYPE_RESOURCE_MANIFEST):
        return {}
    manifest_path = get_resource_path(GENE_BIOTYPE_RESOURCE_MANIFEST)
    return pd.read_json(manifest_path).to_dict() if manifest_path.suffix == ".json" else {}


def list_gene_biotype_resources(
    species: Optional[str] = None,
    cache_dir: Optional[Union[str, Path]] = None,
) -> Dict[str, Dict[str, List[str]]]:
    """
    List bundled and cached gene biotype resources discoverable by scLucid.

    Returns a dictionary keyed by species with `bundled` and `cached` file lists.
    """
    requested_species = [species] if species else list(SUPPORTED_GENE_BIOTYPE_SPECIES)
    cache_root = Path(cache_dir) if cache_dir is not None else get_gene_biotype_cache_dir()
    resources_by_species: Dict[str, Dict[str, List[str]]] = {}

    for sp in requested_species:
        bundled = [name for name in _candidate_biotype_resource_names(sp) if resource_exists(name)]
        cached = [
            str(path)
            for path in _candidate_biotype_cache_paths(sp, cache_dir=cache_root)
            if path.exists()
        ]
        resources_by_species[sp] = {"bundled": bundled, "cached": cached}

    return resources_by_species


def _load_reference_biotypes(
    species: Literal["human", "mouse", "rat"] = "human",
    ensembl_version: Optional[int] = None,
    allow_download: bool = True,
    cache_dir: Optional[Union[str, Path]] = None,
    prefer_bundled: bool = True,
    return_source: bool = False,
) -> pd.DataFrame:
    """
    Load gene biotype annotations from the official reference source.

    Resolution order:
    1. bundled package resource
    2. user-local cache
    3. download from an official remote source if allowed

    Args:
        species: Species name
        ensembl_version: Historical version tag for cache/resource selection

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

    # Map species to Ensembl dataset names for the legacy online fallback path.
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

    cache_root = Path(cache_dir) if cache_dir is not None else get_gene_biotype_cache_dir()

    # Prefer packaged resources for fully offline usage when the release bundles them.
    if prefer_bundled:
        for resource_name in _candidate_biotype_resource_names(species, ensembl_version):
            if resource_exists(resource_name):
                resource_path = get_resource_path(resource_name)
                log.info(f"Loading gene biotypes from packaged resource: {resource_name}")
                df = _read_biotype_table(resource_path)
                return (df, f"package:{resource_name}") if return_source else df

    # Fall back to user-local cache populated by a previous download.
    for cache_path in _candidate_biotype_cache_paths(species, ensembl_version, cache_root):
        if cache_path.exists():
            log.info(f"Loading gene biotypes from cache: {cache_path}")
            df = _read_biotype_table(cache_path)
            return (df, f"cache:{cache_path}") if return_source else df

    # Download from Ensembl BioMart as a fallback path when no bundled/cache resource exists.
    if not allow_download:
        raise FileNotFoundError(
            f"No packaged or cached gene biotype reference found for species='{species}'. "
            "Set allow_download=True to fetch once into the local cache."
        )

    log.info(f"Downloading gene biotype annotations for {species} from Ensembl BioMart...")

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

        # Cache the result in the user-local cache; do not write into package resources.
        cache_name = (
            f"{species}_reference_v{ensembl_version}.csv.gz"
            if ensembl_version is not None
            else f"{species}_reference_latest.csv.gz"
        )
        cache_file = cache_root / cache_name
        df.to_csv(cache_file, index=False, compression="gzip")
        log.info(f"Downloaded and cached {len(df)} gene annotations to: {cache_file}")

        return (df, f"download:{cache_file}") if return_source else df

    except requests.exceptions.RequestException as e:
        log.error(f"Failed to download Ensembl annotations: {e}")
        raise RuntimeError(
            "Could not download gene annotations. Please check your internet connection "
            "or provide a custom biotype annotation file."
        )


# Backward-compatible internal alias.
def _load_ensembl_biotypes(*args, **kwargs):
    return _load_reference_biotypes(*args, **kwargs)


def load_gene_biotypes(
    species: Literal["human", "mouse", "rat"] = "human",
    ensembl_version: Optional[int] = None,
    *,
    allow_download: bool = True,
    cache_dir: Optional[Union[str, Path]] = None,
    prefer_bundled: bool = True,
    return_metadata: bool = False,
) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Dict[str, str]]]:
    """
    Load gene biotype annotations with a stable resolution order:
    packaged resource -> user cache -> one-time download.

    Notes:
    - This function does not write into installed package resources.
    - If a release bundles the reference table under `scLucid/resources/gene_biotypes/`,
      loading is fully offline.
    - If not bundled, the first successful download is cached under
      `~/.sclucid/gene_annotations/` and reused offline on later runs.
    """
    df, source = _load_reference_biotypes(
        species=species,
        ensembl_version=ensembl_version,
        allow_download=allow_download,
        cache_dir=cache_dir,
        prefer_bundled=prefer_bundled,
        return_source=True,
    )
    if return_metadata:
        return df, {"source": source, "species": species}
    return df


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

    biotypes = pd.Series(gene_names.map(gene_to_biotype), index=gene_names, dtype="object")

    # Try matching on gene_id for genes with Ensembl IDs when available.
    if fuzzy_match and "gene_id" in biotype_df.columns:
        unmatched = biotypes.isna()
        if unmatched.sum() > 0:
            log.info(f"Trying Ensembl ID matching for {unmatched.sum()} unmatched genes...")

            gene_id_to_biotype = biotype_df.set_index("gene_id")["biotype"].to_dict()
            biotypes.loc[unmatched] = gene_names[unmatched].map(gene_id_to_biotype)

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
    method: Literal["reference", "ensembl", "custom"] = "reference",
    fuzzy_match: bool = True,
    overwrite: bool = False,
    allow_download: bool = True,
    cache_dir: Optional[Union[str, Path]] = None,
    prefer_bundled: bool = True,
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
        method: Annotation method ('reference', legacy 'ensembl', or 'custom')
        fuzzy_match: Try to match unmatched genes with Ensembl IDs
        overwrite: Overwrite existing biotype annotations

    Returns:
        AnnData with annotated .var

    Examples:
        >>> # Basic usage with packaged/cache-backed reference annotations
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
        log.info("Gene biotypes already annotated. Use overwrite=True to re-annotate.")
        return adata

    log.info(f"Annotating gene biotypes for {adata.n_vars} genes...")

    # Load biotype reference
    source_info = {"source": "custom"}
    if method == "ensembl":
        log.info("method='ensembl' is retained as a legacy alias for method='reference'.")
        method = "reference"

    if method == "reference":
        if biotype_df is None:
            biotype_df, source_info = load_gene_biotypes(
                species=species,
                allow_download=allow_download,
                cache_dir=cache_dir,
                prefer_bundled=prefer_bundled,
                return_metadata=True,
            )
        else:
            log.info("Using provided biotype_df with reference method")
            source_info = {"source": "caller_provided_reference_df"}
    elif method == "custom":
        if biotype_df is None:
            raise ValueError("biotype_df must be provided when method='custom'")

        required_cols = ["gene_name", "biotype"]
        missing_cols = [col for col in required_cols if col not in biotype_df.columns]
        if missing_cols:
            raise ValueError(f"Custom biotype_df missing required columns: {missing_cols}")
    else:
        raise ValueError(f"Unknown method: {method}")

    # Match genes to biotypes
    biotypes = _match_genes_to_biotypes(adata.var_names, biotype_df, fuzzy_match)

    # Add to adata.var
    adata.var["biotype"] = biotypes
    adata.var["biotype_category"] = biotypes.apply(_categorize_biotype)

    # Add recommendation flag
    category_recommendations = {
        cat: info["recommended_for_analysis"] for cat, info in BIOTYPE_CATEGORIES.items()
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
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["gene_biotypes"] = {
        "method": method,
        "species": species,
        "reference_source": source_info["source"],
        "n_genes_annotated": (~adata.var["biotype"].isna()).sum(),
        "n_protein_coding": (adata.var["biotype_category"] == "protein_coding").sum(),
        "annotation_rate": (~adata.var["biotype"].isna()).sum() / adata.n_vars,
    }

    return adata


def apply_gene_biotype_strategy(
    adata: AnnData,
    *,
    species: Literal["human", "mouse", "rat"] = "human",
    method: Literal["reference", "ensembl", "custom"] = "reference",
    custom_biotype_df: Optional[pd.DataFrame] = None,
    custom_biotype_path: Optional[Union[str, Path]] = None,
    keep_biotypes: Optional[List[str]] = None,
    use_recommended: bool = False,
    do_filter: bool = True,
    fuzzy_match: bool = True,
    overwrite: bool = True,
    allow_download: bool = True,
    cache_dir: Optional[Union[str, Path]] = None,
    prefer_bundled: bool = True,
    copy: bool = False,
) -> AnnData:
    """
    Convenience helper for notebook/workflow usage:
    annotate gene biotypes, then optionally filter genes in one call.

    Typical use cases:
    - `method="reference"` for package/cache/download-backed official reference loading
    - `method="custom"` with either a DataFrame or a local TSV/CSV path
    """
    if copy:
        adata = adata.copy()

    if method == "ensembl":
        method = "reference"

    if method == "custom":
        if custom_biotype_df is None and custom_biotype_path is not None:
            custom_path = Path(custom_biotype_path)
            sep = "\t" if custom_path.suffix.lower() in {".tsv", ".txt"} else ","
            custom_biotype_df = pd.read_csv(custom_path, sep=sep)
        if custom_biotype_df is None:
            raise ValueError(
                "Provide custom_biotype_df or custom_biotype_path when method='custom'."
            )

        rename_map = {
            "external_gene_name": "gene_name",
            "gene_biotype": "biotype",
            "ensembl_gene_id": "gene_id",
            "chromosome_name": "chromosome",
            "start_position": "start",
            "end_position": "end",
        }
        custom_biotype_df = custom_biotype_df.rename(
            columns={k: v for k, v in rename_map.items() if k in custom_biotype_df.columns}
        )
        annotated = annotate_gene_biotypes(
            adata,
            species=species,
            biotype_df=custom_biotype_df,
            method="custom",
            fuzzy_match=fuzzy_match,
            overwrite=overwrite,
        )
    else:
        annotated = annotate_gene_biotypes(
            adata,
            species=species,
            method="reference",
            fuzzy_match=fuzzy_match,
            overwrite=overwrite,
            allow_download=allow_download,
            cache_dir=cache_dir,
            prefer_bundled=prefer_bundled,
        )

    if not do_filter:
        return annotated

    if keep_biotypes is None and not use_recommended:
        return annotated

    filtered = filter_genes_by_biotype(
        annotated,
        keep_biotypes=keep_biotypes,
        use_recommended=use_recommended,
        copy=True,
    )
    return filtered if filtered is not None else annotated


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
        raise ValueError("Gene biotypes not found. Run annotate_gene_biotypes() first.")

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


def get_biotype_statistics(adata: AnnData, sample_key: Optional[str] = None) -> pd.DataFrame:
    """
    Generate comprehensive biotype statistics.

    Args:
        adata: AnnData object with biotype annotations
        sample_key: If provided, calculate per-sample statistics

    Returns:
        DataFrame with biotype statistics
    """
    if "biotype_category" not in adata.var.columns:
        raise ValueError("Gene biotypes not found. Run annotate_gene_biotypes() first.")

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
        raise ValueError("Gene biotypes not found. Run annotate_gene_biotypes() first.")

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
    log.info(f"Genes retained: {rec['n_genes_kept']}/{adata.n_vars} ({actual_retention:.1%})")
    log.info(f"Genes removed: {rec['n_genes_removed']} ({1-actual_retention:.1%})\n")

    # Show what will be removed
    removed_biotypes = adata.var.loc[~keep_mask, "biotype_category"].value_counts()
    if len(removed_biotypes) > 0:
        log.info("Biotypes to be removed:")
        for biotype, count in removed_biotypes.items():
            log.info(f"  - {biotype}: {count} genes")

    return rec
