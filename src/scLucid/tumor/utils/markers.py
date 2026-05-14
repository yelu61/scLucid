"""
Tumor marker gene definitions and utilities.

This module provides curated lists of marker genes for various
cell types and tumor characteristics.

Note: Marker data is now loaded from resources/ directory via Manager class.
This module provides backward-compatible wrapper functions.
"""

from typing import Dict, List, Optional


def get_tumor_markers(
    cancer_type: str = "pan_cancer",
    custom_markers: Optional[List[str]] = None,
) -> List[str]:
    """
    Get tumor marker genes for a specific cancer type.

    Parameters
    ----------
    cancer_type : str
        Type of cancer ("pan_cancer", "Lung Cancer", "Breast Cancer", etc.)
    custom_markers : list, optional
        Additional custom markers to include

    Returns:
    -------
    list
        Tumor marker gene names
    """
    from ...utils.manager import _get_cancer_markers

    cancer_markers = _get_cancer_markers(species="human")

    # Try exact match first, then case-insensitive
    if cancer_type in cancer_markers:
        markers = cancer_markers[cancer_type]["markers"].copy()
    else:
        # Try case-insensitive lookup
        for key, value in cancer_markers.items():
            if key.lower() == cancer_type.lower():
                markers = value["markers"].copy()
                break
        else:
            # Fall back to pan-cancer
            markers = cancer_markers.get("Pan-Cancer", {}).get("markers", []).copy()

    if custom_markers:
        markers = list(set(markers + custom_markers))

    return markers


def get_immune_markers(
    cell_type: Optional[str] = None,
    include_subtypes: bool = True,
) -> Dict[str, List[str]]:
    """
    Get immune cell marker genes.

    Parameters
    ----------
    cell_type : str, optional
        Specific immune cell type. If None, returns all.
    include_subtypes : bool
        Include subtype markers for T cells

    Returns:
    -------
    dict or list
        Immune marker genes by cell type

    Note:
    ----
    This function now delegates to the global Manager class.
    Consider using `get_marker_manager()` directly for new code.
    """
    from ...utils.manager import get_marker_manager

    mgr = get_marker_manager(species="human")

    # Get all immune cell types
    immune_cells = {}

    # Query from manager
    for name, cell in mgr.CELLS.items():
        # Check if it's an immune cell by looking at hierarchy
        if mgr[name].level == "minor":
            continue  # Skip subtypes initially

        markers = mgr.query("markers", name)[name]
        immune_cells[name] = markers

        if include_subtypes and cell.minor:
            for subtype in cell.minor:
                immune_cells[subtype.name] = subtype.markers

    if cell_type is not None:
        return immune_cells.get(cell_type, [])

    return immune_cells


def get_stromal_markers(
    cell_type: Optional[str] = None,
) -> Dict[str, List[str]]:
    """
    Get stromal cell marker genes.

    Parameters
    ----------
    cell_type : str, optional
        Specific stromal cell type. If None, returns all.

    Returns:
    -------
    dict or list
        Stromal marker genes

    Note:
    ----
    This function now delegates to the global Manager class.
    Consider using `get_marker_manager()` directly for new code.
    """
    from ...utils.manager import get_marker_manager

    mgr = get_marker_manager(species="human")

    # Get stromal cells (Fibroblasts, Endothelial, etc.)
    stromal_cells = {}

    for name, cell in mgr.CELLS.items():
        # Check if it might be stromal by name patterns
        if any(
            keyword in name.lower()
            for keyword in ["fibroblast", "endothelial", "pericyte", "smooth muscle"]
        ):
            markers = mgr.query("markers", name)[name]
            stromal_cells[name] = markers

            # Add subtypes
            for subtype in cell.minor:
                stromal_cells[subtype.name] = subtype.markers

    if cell_type is not None:
        return stromal_cells.get(cell_type, [])

    return stromal_cells


def get_proliferation_markers(
    phase: Optional[str] = None,
) -> Dict[str, List[str]]:
    """
    Get cell cycle/proliferation marker genes.

    Parameters
    ----------
    phase : str, optional
        Cell cycle phase ("proliferation_core", "g1_s_genes", "s_genes", "g2m_genes", "m_genes")

    Returns:
    -------
    dict or list
        Proliferation marker genes

    Note:
    ----
    This function now loads resource gene sets through the unified Manager helpers.
    """
    from ...utils.manager import load_gene_sets

    try:
        genesets = load_gene_sets(species="human", name="cell_cycle")
    except FileNotFoundError:
        # Fallback to minimal set
        return {
            "proliferation_core": ["MKI67", "PCNA", "TOP2A", "AURKA", "AURKB"],
            "g1_s_genes": ["CCNE1", "CCNE2", "E2F1"],
            "s_genes": ["MCM2", "MCM4", "PCNA"],
            "g2m_genes": ["MKI67", "TOP2A", "CCNB1"],
            "m_genes": ["PRC1", "KIF11", "AURKB"],
        }

    if phase is not None:
        # Map old phase names to new keys
        phase_mapping = {
            "core": "proliferation_core",
            "G1_S": "g1_s_genes",
            "S": "s_genes",
            "G2_M": "g2m_genes",
            "M": "m_genes",
        }
        lookup_key = phase_mapping.get(phase, phase)
        return genesets.get(lookup_key, genesets.get("proliferation_core", []))

    return genesets


def get_emt_markers(
    state: str = "all",
) -> Dict[str, List[str]]:
    """
    Get EMT (Epithelial-Mesenchymal Transition) marker genes.

    Parameters
    ----------
    state : str
        EMT state ("epithelial", "mesenchymal", "hybrid", "all")

    Returns:
    -------
    dict or list
        EMT marker genes

    Note:
    ----
    This function now loads resource gene sets through the unified Manager helpers.
    """
    from ...utils.manager import load_gene_sets

    try:
        genesets = load_gene_sets(species="human", name="cancer_signatures")
    except FileNotFoundError:
        # Fallback
        return {
            "epithelial": ["CDH1", "EPCAM", "CLDN1", "CLDN3", "CLDN4"],
            "mesenchymal": ["VIM", "CDH2", "FN1", "ZEB1", "SNAI1"],
            "hybrid": ["CLDN1", "VIM", "CDH1", "ZEB1"],
        }

    emt_states = {
        "epithelial": genesets.get("EMT_epithelial", []),
        "mesenchymal": genesets.get("EMT_mesenchymal", []),
        "hybrid": genesets.get("EMT_hybrid", []),
    }

    if state == "all":
        return emt_states

    return {state: emt_states.get(state, [])}


def get_all_markers() -> Dict[str, Dict[str, List[str]]]:
    """
    Get all marker gene definitions.

    Returns:
    -------
    dict
        All marker gene sets organized by category
    """
    return {
        "tumor": {
            k: v
            for k, v in [
                ("pan_cancer", get_tumor_markers("pan_cancer")),
                ("lung", get_tumor_markers("lung")),
                ("breast", get_tumor_markers("breast")),
                ("colorectal", get_tumor_markers("colorectal")),
                ("liver", get_tumor_markers("liver")),
                ("gastric", get_tumor_markers("gastric")),
            ]
        },
        "immune": get_immune_markers(),
        "stromal": get_stromal_markers(),
        "proliferation": get_proliferation_markers(),
        "emt": get_emt_markers("all"),
    }


def search_markers(
    gene: str,
    return_categories: bool = True,
) -> Dict[str, List[str]]:
    """
    Search for which marker categories a gene belongs to.

    Parameters
    ----------
    gene : str
        Gene symbol to search
    return_categories : bool
        Return categories containing the gene

    Returns:
    -------
    dict
        Categories and specific sets containing the gene
    """
    results = {}
    all_markers = get_all_markers()

    gene_upper = gene.upper()

    for category, marker_sets in all_markers.items():
        if isinstance(marker_sets, dict):
            matches = {
                k: v
                for k, v in marker_sets.items()
                if isinstance(v, list) and gene_upper in [g.upper() for g in v]
            }
            if matches:
                results[category] = matches
        elif isinstance(marker_sets, list):
            if gene_upper in [g.upper() for g in marker_sets]:
                results[category] = marker_sets

    return results


# Deprecated: These constants are kept for backward compatibility
# but will be removed in a future version. Use the functions instead.

TUMOR_MARKERS = {
    "pan_cancer": get_tumor_markers("pan_cancer"),
    "lung": get_tumor_markers("lung"),
    "breast": get_tumor_markers("breast"),
    "colorectal": get_tumor_markers("colorectal"),
    "liver": get_tumor_markers("liver"),
    "gastric": get_tumor_markers("gastric"),
}

IMMUNE_MARKERS = {}  # Now loaded dynamically via get_immune_markers()
STROMAL_MARKERS = {}  # Now loaded dynamically via get_stromal_markers()
PROLIFERATION_MARKERS = {}  # Now loaded dynamically via get_proliferation_markers()
EMT_MARKERS = {}  # Now loaded dynamically via get_emt_markers()
HYPOXIA_MARKERS = []  # Now available via functional signatures
DNA_REPAIR_MARKERS = {}  # Now available via functional signatures
