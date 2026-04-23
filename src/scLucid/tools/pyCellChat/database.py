"""
CellChatDB database implementation (R-free)

Provides ligand-receptor interaction databases for human and mouse.
"""

import logging
from typing import List, Optional

import pandas as pd

log = logging.getLogger(__name__)


# Built-in interaction data (subset of common interactions)
# In production, this would be loaded from files
BUILTIN_INTERACTIONS = {
    "human": [
        # Notch signaling
        {
            "interaction_name": "DLL1_NOTCH1",
            "pathway_name": "NOTCH",
            "ligand": "DLL1",
            "receptor": "NOTCH1",
            "interaction_type": "Cell-Cell Contact",
        },
        {
            "interaction_name": "DLL4_NOTCH1",
            "pathway_name": "NOTCH",
            "ligand": "DLL4",
            "receptor": "NOTCH1",
            "interaction_type": "Cell-Cell Contact",
        },
        {
            "interaction_name": "JAG1_NOTCH1",
            "pathway_name": "NOTCH",
            "ligand": "JAG1",
            "receptor": "NOTCH1",
            "interaction_type": "Cell-Cell Contact",
        },
        {
            "interaction_name": "JAG2_NOTCH1",
            "pathway_name": "NOTCH",
            "ligand": "JAG2",
            "receptor": "NOTCH1",
            "interaction_type": "Cell-Cell Contact",
        },
        # Wnt signaling
        {
            "interaction_name": "WNT3A_FZD4",
            "pathway_name": "WNT",
            "ligand": "WNT3A",
            "receptor": "FZD4",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "WNT5A_FZD5",
            "pathway_name": "WNT",
            "ligand": "WNT5A",
            "receptor": "FZD5",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "WNT7A_FZD7",
            "pathway_name": "WNT",
            "ligand": "WNT7A",
            "receptor": "FZD7",
            "interaction_type": "Secreted",
        },
        # TGFB signaling
        {
            "interaction_name": "TGFB1_TGFBR1",
            "pathway_name": "TGFB",
            "ligand": "TGFB1",
            "receptor": "TGFBR1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "TGFB2_TGFBR1",
            "pathway_name": "TGFB",
            "ligand": "TGFB2",
            "receptor": "TGFBR1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "TGFB3_TGFBR1",
            "pathway_name": "TGFB",
            "ligand": "TGFB3",
            "receptor": "TGFBR1",
            "interaction_type": "Secreted",
        },
        # BMP signaling
        {
            "interaction_name": "BMP2_BMPR1A",
            "pathway_name": "BMP",
            "ligand": "BMP2",
            "receptor": "BMPR1A",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "BMP4_BMPR1A",
            "pathway_name": "BMP",
            "ligand": "BMP4",
            "receptor": "BMPR1A",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "BMP7_BMPR1A",
            "pathway_name": "BMP",
            "ligand": "BMP7",
            "receptor": "BMPR1A",
            "interaction_type": "Secreted",
        },
        # EGF signaling
        {
            "interaction_name": "EGF_EGFR",
            "pathway_name": "EGF",
            "ligand": "EGF",
            "receptor": "EGFR",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "AREG_EGFR",
            "pathway_name": "EGF",
            "ligand": "AREG",
            "receptor": "EGFR",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "EREG_EGFR",
            "pathway_name": "EGF",
            "ligand": "EREG",
            "receptor": "EGFR",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "HBEGF_EGFR",
            "pathway_name": "EGF",
            "ligand": "HBEGF",
            "receptor": "EGFR",
            "interaction_type": "Secreted",
        },
        # FGF signaling
        {
            "interaction_name": "FGF1_FGFR1",
            "pathway_name": "FGF",
            "ligand": "FGF1",
            "receptor": "FGFR1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "FGF2_FGFR1",
            "pathway_name": "FGF",
            "ligand": "FGF2",
            "receptor": "FGFR1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "FGF7_FGFR2",
            "pathway_name": "FGF",
            "ligand": "FGF7",
            "receptor": "FGFR2",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "FGF10_FGFR2",
            "pathway_name": "FGF",
            "ligand": "FGF10",
            "receptor": "FGFR2",
            "interaction_type": "Secreted",
        },
        # VEGF signaling
        {
            "interaction_name": "VEGFA_VEGFR1",
            "pathway_name": "VEGF",
            "ligand": "VEGFA",
            "receptor": "FLT1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "VEGFA_VEGFR2",
            "pathway_name": "VEGF",
            "ligand": "VEGFA",
            "receptor": "KDR",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "VEGFB_VEGFR1",
            "pathway_name": "VEGF",
            "ligand": "VEGFB",
            "receptor": "FLT1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "PGF_VEGFR1",
            "pathway_name": "VEGF",
            "ligand": "PGF",
            "receptor": "FLT1",
            "interaction_type": "Secreted",
        },
        # PDGF signaling
        {
            "interaction_name": "PDGFA_PDGFRA",
            "pathway_name": "PDGF",
            "ligand": "PDGFA",
            "receptor": "PDGFRA",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "PDGFB_PDGFRB",
            "pathway_name": "PDGF",
            "ligand": "PDGFB",
            "receptor": "PDGFRB",
            "interaction_type": "Secreted",
        },
        # IGF signaling
        {
            "interaction_name": "IGF1_IGF1R",
            "pathway_name": "IGF",
            "ligand": "IGF1",
            "receptor": "IGF1R",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IGF2_IGF1R",
            "pathway_name": "IGF",
            "ligand": "IGF2",
            "receptor": "IGF1R",
            "interaction_type": "Secreted",
        },
        # HGF signaling
        {
            "interaction_name": "HGF_MET",
            "pathway_name": "HGF",
            "ligand": "HGF",
            "receptor": "MET",
            "interaction_type": "Secreted",
        },
        # IFN signaling
        {
            "interaction_name": "IFNG_IFNGR1",
            "pathway_name": "IFN-II",
            "ligand": "IFNG",
            "receptor": "IFNGR1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IFNA1_IFNAR1",
            "pathway_name": "IFN-I",
            "ligand": "IFNA1",
            "receptor": "IFNAR1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IFNB1_IFNAR1",
            "pathway_name": "IFN-I",
            "ligand": "IFNB1",
            "receptor": "IFNAR1",
            "interaction_type": "Secreted",
        },
        # IL signaling (subset)
        {
            "interaction_name": "IL1A_IL1R1",
            "pathway_name": "IL1",
            "ligand": "IL1A",
            "receptor": "IL1R1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IL1B_IL1R1",
            "pathway_name": "IL1",
            "ligand": "IL1B",
            "receptor": "IL1R1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IL2_IL2RA",
            "pathway_name": "IL2",
            "ligand": "IL2",
            "receptor": "IL2RA",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IL4_IL4R",
            "pathway_name": "IL4",
            "ligand": "IL4",
            "receptor": "IL4R",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IL6_IL6R",
            "pathway_name": "IL6",
            "ligand": "IL6",
            "receptor": "IL6R",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IL10_IL10RA",
            "pathway_name": "IL10",
            "ligand": "IL10",
            "receptor": "IL10RA",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IL12A_IL12RB1",
            "pathway_name": "IL12",
            "ligand": "IL12A",
            "receptor": "IL12RB1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "IL17A_IL17RA",
            "pathway_name": "IL17",
            "ligand": "IL17A",
            "receptor": "IL17RA",
            "interaction_type": "Secreted",
        },
        # TNF signaling
        {
            "interaction_name": "TNF_TNFRSF1A",
            "pathway_name": "TNF",
            "ligand": "TNF",
            "receptor": "TNFRSF1A",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "TNFSF10_TNFRSF10A",
            "pathway_name": "TRAIL",
            "ligand": "TNFSF10",
            "receptor": "TNFRSF10A",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "FASLG_FAS",
            "pathway_name": "FAS",
            "ligand": "FASLG",
            "receptor": "FAS",
            "interaction_type": "Secreted",
        },
        # CSF signaling
        {
            "interaction_name": "CSF1_CSF1R",
            "pathway_name": "CSF",
            "ligand": "CSF1",
            "receptor": "CSF1R",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "CSF2_CSF2RA",
            "pathway_name": "CSF",
            "ligand": "CSF2",
            "receptor": "CSF2RA",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "CSF3_CSF3R",
            "pathway_name": "CSF",
            "ligand": "CSF3",
            "receptor": "CSF3R",
            "interaction_type": "Secreted",
        },
        # Chemokine signaling
        {
            "interaction_name": "CXCL12_CXCR4",
            "pathway_name": "CXCL",
            "ligand": "CXCL12",
            "receptor": "CXCR4",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "CXCL10_CXCR3",
            "pathway_name": "CXCL",
            "ligand": "CXCL10",
            "receptor": "CXCR3",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "CCL2_CCR2",
            "pathway_name": "CCL",
            "ligand": "CCL2",
            "receptor": "CCR2",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "CCL5_CCR5",
            "pathway_name": "CCL",
            "ligand": "CCL5",
            "receptor": "CCR5",
            "interaction_type": "Secreted",
        },
        # SPP1 signaling
        {
            "interaction_name": "SPP1_CD44",
            "pathway_name": "SPP1",
            "ligand": "SPP1",
            "receptor": "CD44",
            "interaction_type": "ECM-Receptor",
        },
        {
            "interaction_name": "SPP1_ITGAV",
            "pathway_name": "SPP1",
            "ligand": "SPP1",
            "receptor": "ITGAV",
            "interaction_type": "ECM-Receptor",
        },
        # ANGPTL signaling
        {
            "interaction_name": "ANGPT1_TEK",
            "pathway_name": "ANGPTL",
            "ligand": "ANGPT1",
            "receptor": "TEK",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "ANGPT2_TEK",
            "pathway_name": "ANGPTL",
            "ligand": "ANGPT2",
            "receptor": "TEK",
            "interaction_type": "Secreted",
        },
        # APP signaling
        {
            "interaction_name": "APP_CD74",
            "pathway_name": "APP",
            "ligand": "APP",
            "receptor": "CD74",
            "interaction_type": "Cell-Cell Contact",
        },
        # PROS signaling
        {
            "interaction_name": "PROS1_AXL",
            "pathway_name": "GAS",
            "ligand": "PROS1",
            "receptor": "AXL",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "GAS6_AXL",
            "pathway_name": "GAS",
            "ligand": "GAS6",
            "receptor": "AXL",
            "interaction_type": "Secreted",
        },
        # NT signaling
        {
            "interaction_name": "NGF_NTRK1",
            "pathway_name": "NGF",
            "ligand": "NGF",
            "receptor": "NTRK1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "BDNF_NTRK2",
            "pathway_name": "NGF",
            "ligand": "BDNF",
            "receptor": "NTRK2",
            "interaction_type": "Secreted",
        },
        # SEMA signaling
        {
            "interaction_name": "SEMA3A_PLXNA1",
            "pathway_name": "SEMA3",
            "ligand": "SEMA3A",
            "receptor": "PLXNA1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "SEMA4D_PLXNB1",
            "pathway_name": "SEMA4",
            "ligand": "SEMA4D",
            "receptor": "PLXNB1",
            "interaction_type": "Cell-Cell Contact",
        },
        # VISFATIN signaling
        {
            "interaction_name": "NAMPT_ITGAV",
            "pathway_name": "VISFATIN",
            "ligand": "NAMPT",
            "receptor": "ITGAV",
            "interaction_type": "Secreted",
        },
        # COMPLEMENT signaling
        {
            "interaction_name": "C3_C3AR1",
            "pathway_name": "COMPLEMENT",
            "ligand": "C3",
            "receptor": "C3AR1",
            "interaction_type": "Secreted",
        },
    ],
    "mouse": [
        # Mouse orthologs (gene names are capitalized differently in mouse)
        {
            "interaction_name": "Dll1_Notch1",
            "pathway_name": "NOTCH",
            "ligand": "Dll1",
            "receptor": "Notch1",
            "interaction_type": "Cell-Cell Contact",
        },
        {
            "interaction_name": "Dll4_Notch1",
            "pathway_name": "NOTCH",
            "ligand": "Dll4",
            "receptor": "Notch1",
            "interaction_type": "Cell-Cell Contact",
        },
        {
            "interaction_name": "Jag1_Notch1",
            "pathway_name": "NOTCH",
            "ligand": "Jag1",
            "receptor": "Notch1",
            "interaction_type": "Cell-Cell Contact",
        },
        {
            "interaction_name": "Wnt3a_Fzd4",
            "pathway_name": "WNT",
            "ligand": "Wnt3a",
            "receptor": "Fzd4",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Wnt5a_Fzd5",
            "pathway_name": "WNT",
            "ligand": "Wnt5a",
            "receptor": "Fzd5",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Tgfb1_Tgfbr1",
            "pathway_name": "TGFB",
            "ligand": "Tgfb1",
            "receptor": "Tgfbr1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Tgfb2_Tgfbr1",
            "pathway_name": "TGFB",
            "ligand": "Tgfb2",
            "receptor": "Tgfbr1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Bmp2_Bmpr1a",
            "pathway_name": "BMP",
            "ligand": "Bmp2",
            "receptor": "Bmpr1a",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Bmp4_Bmpr1a",
            "pathway_name": "BMP",
            "ligand": "Bmp4",
            "receptor": "Bmpr1a",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Egf_Egfr",
            "pathway_name": "EGF",
            "ligand": "Egf",
            "receptor": "Egfr",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Areg_Egfr",
            "pathway_name": "EGF",
            "ligand": "Areg",
            "receptor": "Egfr",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Hbegf_Egfr",
            "pathway_name": "EGF",
            "ligand": "Hbegf",
            "receptor": "Egfr",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Fgf1_Fgfr1",
            "pathway_name": "FGF",
            "ligand": "Fgf1",
            "receptor": "Fgfr1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Fgf2_Fgfr1",
            "pathway_name": "FGF",
            "ligand": "Fgf2",
            "receptor": "Fgfr1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Vegfa_Kdr",
            "pathway_name": "VEGF",
            "ligand": "Vegfa",
            "receptor": "Kdr",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Vegfa_Flt1",
            "pathway_name": "VEGF",
            "ligand": "Vegfa",
            "receptor": "Flt1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Pgf_Flt1",
            "pathway_name": "VEGF",
            "ligand": "Pgf",
            "receptor": "Flt1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Igf1_Igf1r",
            "pathway_name": "IGF",
            "ligand": "Igf1",
            "receptor": "Igf1r",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Igf2_Igf1r",
            "pathway_name": "IGF",
            "ligand": "Igf2",
            "receptor": "Igf1r",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Hgf_Met",
            "pathway_name": "HGF",
            "ligand": "Hgf",
            "receptor": "Met",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Ifng_Ifngr1",
            "pathway_name": "IFN-II",
            "ligand": "Ifng",
            "receptor": "Ifngr1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Il1a_Il1r1",
            "pathway_name": "IL1",
            "ligand": "Il1a",
            "receptor": "Il1r1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Il1b_Il1r1",
            "pathway_name": "IL1",
            "ligand": "Il1b",
            "receptor": "Il1r1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Il2_Il2ra",
            "pathway_name": "IL2",
            "ligand": "Il2",
            "receptor": "Il2ra",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Il4_Il4r",
            "pathway_name": "IL4",
            "ligand": "Il4",
            "receptor": "Il4r",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Il6_Il6r",
            "pathway_name": "IL6",
            "ligand": "Il6",
            "receptor": "Il6r",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Il10_Il10ra",
            "pathway_name": "IL10",
            "ligand": "Il10",
            "receptor": "Il10ra",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Tnf_Tnfrsf1a",
            "pathway_name": "TNF",
            "ligand": "Tnf",
            "receptor": "Tnfrsf1a",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Tnfsf10_Tnfrsf10a",
            "pathway_name": "TRAIL",
            "ligand": "Tnfsf10",
            "receptor": "Tnfrsf10a",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Fasl_Fas",
            "pathway_name": "FAS",
            "ligand": "Fasl",
            "receptor": "Fas",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Csf1_Csf1r",
            "pathway_name": "CSF",
            "ligand": "Csf1",
            "receptor": "Csf1r",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Csf2_Csf2ra",
            "pathway_name": "CSF",
            "ligand": "Csf2",
            "receptor": "Csf2ra",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Cxcl12_Cxcr4",
            "pathway_name": "CXCL",
            "ligand": "Cxcl12",
            "receptor": "Cxcr4",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Ccl2_Ccr2",
            "pathway_name": "CCL",
            "ligand": "Ccl2",
            "receptor": "Ccr2",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Spp1_Cd44",
            "pathway_name": "SPP1",
            "ligand": "Spp1",
            "receptor": "Cd44",
            "interaction_type": "ECM-Receptor",
        },
        {
            "interaction_name": "Angpt1_Tek",
            "pathway_name": "ANGPTL",
            "ligand": "Angpt1",
            "receptor": "Tek",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Pros1_Axl",
            "pathway_name": "GAS",
            "ligand": "Pros1",
            "receptor": "Axl",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Gas6_Axl",
            "pathway_name": "GAS",
            "ligand": "Gas6",
            "receptor": "Axl",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Ngf_Ntrk1",
            "pathway_name": "NGF",
            "ligand": "Ngf",
            "receptor": "Ntrk1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Bdnf_Ntrk2",
            "pathway_name": "NGF",
            "ligand": "Bdnf",
            "receptor": "Ntrk2",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "Sema3a_Plxna1",
            "pathway_name": "SEMA3",
            "ligand": "Sema3a",
            "receptor": "Plxna1",
            "interaction_type": "Secreted",
        },
        {
            "interaction_name": "C3_C3ar1",
            "pathway_name": "COMPLEMENT",
            "ligand": "C3",
            "receptor": "C3ar1",
            "interaction_type": "Secreted",
        },
    ],
}


class CellChatDB:
    """
    CellChatDB: Database of ligand-receptor interactions (Pure Python)

    Parameters
    ----------
    species : str
        Species ("human" or "mouse")
    version : str
        Database version
    """

    def __init__(self, species: str = "human", version: str = "v2"):
        self.species = species.lower()
        self.version = version

        # Load database
        self.interaction = self._load_interaction_db()
        self.complex = self._load_complex_db()
        self.cofactor = self._load_cofactor_db()
        self.geneInfo = self._load_gene_info()

        log.info(f"Loaded CellChatDB {version} for {species}: {len(self.interaction)} interactions")

    def _load_interaction_db(self) -> pd.DataFrame:
        """Load ligand-receptor interaction database"""
        if self.species in BUILTIN_INTERACTIONS:
            interactions = BUILTIN_INTERACTIONS[self.species]
        else:
            log.warning(f"Species {self.species} not found, using human as default")
            interactions = BUILTIN_INTERACTIONS["human"]

        df = pd.DataFrame(interactions)

        # Add evidence and annotation columns if missing
        if "evidence" not in df.columns:
            df["evidence"] = "KEGG;PubMed"
        if "annotation" not in df.columns:
            df["annotation"] = "activation"

        return df

    def _load_complex_db(self) -> pd.DataFrame:
        """Load complex composition database"""
        # Simplified - no multi-subunit complexes in built-in data
        complex_db = pd.DataFrame(
            {
                "complex_name": [],
                "subunit": [],
                "subunit_type": [],
            }
        )
        return complex_db

    def _load_cofactor_db(self) -> pd.DataFrame:
        """Load cofactor database"""
        cofactor_db = pd.DataFrame(
            {
                "cofactor": [],
                "interaction_name": [],
            }
        )
        return cofactor_db

    def _load_gene_info(self) -> pd.DataFrame:
        """Load gene information"""
        # Create from interaction data
        all_genes = self.get_all_genes()
        gene_info = pd.DataFrame(
            {
                "Symbol": all_genes,
                "GeneName": all_genes,
                "GeneID": range(len(all_genes)),
            }
        )
        return gene_info

    def get_all_genes(self) -> List[str]:
        """Get all genes in database"""
        genes = set()

        # Add ligands
        genes.update(self.interaction["ligand"].dropna().unique())

        # Add receptors
        genes.update(self.interaction["receptor"].dropna().unique())

        # Add complex subunits
        if not self.complex.empty:
            genes.update(self.complex["subunit"].dropna().unique())

        # Add cofactors
        if not self.cofactor.empty:
            genes.update(self.cofactor["cofactor"].dropna().unique())

        return list(genes)

    def get_pathways(self) -> List[str]:
        """Get all pathways in database"""
        return self.interaction["pathway_name"].dropna().unique().tolist()

    def get_interaction_types(self) -> List[str]:
        """Get all interaction types"""
        return self.interaction["interaction_type"].dropna().unique().tolist()

    def subset_db(
        self, interaction_types: Optional[List[str]] = None, pathways: Optional[List[str]] = None
    ) -> "CellChatDB":
        """
        Subset database by interaction types or pathways

        Parameters
        ----------
        interaction_types : Optional[List[str]]
            List of interaction types to keep ("Secreted", "ECM-Receptor", "Cell-Cell Contact")
        pathways : Optional[List[str]]
            List of pathways to keep
        """
        subset_db = CellChatDB(species=self.species, version=self.version)
        subset_db.interaction = self.interaction.copy()

        if interaction_types is not None:
            subset_db.interaction = subset_db.interaction[
                subset_db.interaction["interaction_type"].isin(interaction_types)
            ]

        if pathways is not None:
            subset_db.interaction = subset_db.interaction[
                subset_db.interaction["pathway_name"].isin(pathways)
            ]

        log.info(f"Subset database: {len(subset_db.interaction)} interactions remaining")
        return subset_db

    def update_db(self, custom_interactions: pd.DataFrame):
        """
        Update database with custom interactions

        Parameters
        ----------
        custom_interactions : pd.DataFrame
            Custom interactions to add
        """
        required_cols = ["interaction_name", "pathway_name", "ligand", "receptor"]

        if not all(col in custom_interactions.columns for col in required_cols):
            raise ValueError(f"Custom interactions must contain: {required_cols}")

        self.interaction = pd.concat([self.interaction, custom_interactions], ignore_index=True)

        log.info(f"Added {len(custom_interactions)} custom interactions")

    def search_interaction(
        self,
        ligand: Optional[str] = None,
        receptor: Optional[str] = None,
        pathway: Optional[str] = None,
    ) -> pd.DataFrame:
        """Search interactions by ligand, receptor, or pathway"""
        results = self.interaction.copy()

        if ligand is not None:
            results = results[results["ligand"].str.contains(ligand, case=False, na=False)]

        if receptor is not None:
            results = results[results["receptor"].str.contains(receptor, case=False, na=False)]

        if pathway is not None:
            results = results[results["pathway_name"].str.contains(pathway, case=False, na=False)]

        return results

    def save(self, filepath: str):
        """Save database to file"""
        self.interaction.to_csv(filepath, index=False)
        log.info(f"Saved database to {filepath}")

    @classmethod
    def load(cls, filepath: str, species: str = "human") -> "CellChatDB":
        """Load database from file"""
        db = cls(species=species)
        db.interaction = pd.read_csv(filepath)
        log.info(f"Loaded database from {filepath}: {len(db.interaction)} interactions")
        return db


def get_default_database(species: str = "human") -> CellChatDB:
    """
    Get default CellChatDB database for the specified species

    Parameters
    ----------
    species : str
        Species ("human" or "mouse")

    Returns:
    -------
    CellChatDB
        Default database instance
    """
    return CellChatDB(species=species)


def merge_databases(db1: CellChatDB, db2: CellChatDB) -> CellChatDB:
    """
    Merge two databases

    Parameters
    ----------
    db1 : CellChatDB
        First database
    db2 : CellChatDB
        Second database

    Returns:
    -------
    CellChatDB
        Merged database
    """
    merged = CellChatDB(species=db1.species)
    merged.interaction = pd.concat(
        [db1.interaction, db2.interaction], ignore_index=True
    ).drop_duplicates()

    log.info(f"Merged databases: {len(merged.interaction)} unique interactions")
    return merged
