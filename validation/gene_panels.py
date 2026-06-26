"""Shared marker/program panels for validation benchmarks."""

from __future__ import annotations

MARKER_PANELS: dict[str, tuple[str, ...]] = {
    "immune_t": ("CD3D", "CD3E", "TRAC", "CD4", "CD8A", "NKG7"),
    "myeloid": ("LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7", "LST1"),
    "b_plasma": ("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"),
    "epithelial": ("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1"),
    "stromal": ("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2"),
    "endothelial": ("PECAM1", "VWF", "KDR", "ENG"),
    "proliferation": ("MKI67", "TOP2A", "STMN1", "UBE2C"),
    "hypoxia_stress": ("VEGFA", "CA9", "DDIT3", "HSPA1A", "JUN"),
}

TUMOR_PROGRAM_PANELS: dict[str, tuple[str, ...]] = {
    "epithelial_malignant_like": ("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "TACSTD2"),
    "cell_cycle": ("MKI67", "TOP2A", "STMN1", "UBE2C", "PCNA", "TYMS"),
    "hypoxia_stress": ("VEGFA", "CA9", "DDIT3", "HSPA1A", "JUN", "FOS"),
    "emt_stromal": ("VIM", "FN1", "COL1A1", "COL1A2", "SPARC", "ACTA2"),
    "oxphos": (
        "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6",
        "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8",
        "NDUFS1", "NDUFS2", "SDHA", "UQCRC1", "COX5B", "ATP5F1A",
    ),
    "glycolysis": (
        "HK1", "HK2", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1", "PGAM1",
        "ENO1", "PKM", "LDHA", "LDHB", "SLC2A1", "SLC2A3", "PGM1",
    ),
    "mt_biogenesis": (
        "TFAM", "TFB1M", "TFB2M", "POLRMT", "PPARGC1A", "PPARGC1B",
        "NRF1", "NRF2", "TFB2M", "MTERF1", "MRPL12", "MRPS12",
    ),
}


def present_genes(var_names, genes: tuple[str, ...]) -> list[str]:
    """Return panel genes present in an AnnData var index."""
    names = set(map(str, var_names))
    return [gene for gene in genes if gene in names]
