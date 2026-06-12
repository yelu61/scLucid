# Marker Curation Batch 10: Fibroblast/CAF Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | Single-cell screens identify ADAM12 as a fibroblast checkpoint impeding anti-tumor immunity | 2026 | Cancer Cell | 10.1016/j.ccell.2025.12.018 | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 2 | Single-cell resolution spatial analysis of antigen-presenting CAF niches | 2025 | Cancer Cell | 10.1016/j.ccell.2025.09.001 | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 3 | Conserved spatial subtypes and cellular neighborhoods of CAFs | 2025 | Cancer Cell | 10.1016/j.ccell.2025.03.004 | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 4 | Fibroblast atlas: shared and specific cell types across tissues | 2025 | Science Advances | 10.1126/sciadv.ado0173 | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 5 | Cross-tissue human fibroblast atlas reveals myofibroblast subtypes | 2024 | Cancer Cell | 10.1016/j.ccell.2024.08.020 | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 6 | Classifying cancer-associated fibroblasts-the good, the bad, and the target | 2024 | Cancer Cell | 10.1016/j.ccell.2024.08.011 | review | marker_registry (naming) |
| 7 | Cancer associated fibroblasts in cancer development and therapy | 2025 | J Hematol Oncol | 10.1186/s13045-025-01688-0 | review | marker_registry (validation) |
| 8 | Pan-cancer spatially resolved single-cell analysis reveals CAF-TME crosstalk | 2023 | Mol Cancer | 10.1186/s12943-023-01876-x | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 9 | The advent of immune stimulating CAFs in cancer | 2023 | Nat Rev Cancer | 10.1038/s41568-023-00549-7 | review | marker_registry (naming + framework) |
| 10 | Pan-cancer single-cell analysis reveals heterogeneity and plasticity of CAFs | 2022 | Nat Commun | 10.1038/s41467-022-34395-2 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 11 | Cancer-associated fibroblasts in the single-cell era | 2022 | Nat Cancer | 10.1038/s43018-022-00411-z | review | marker_registry (validation + naming) |
| 12 | Clinical and therapeutic relevance of cancer-associated fibroblasts | 2021 | Nat Rev Clin Oncol | 10.1038/s41571-021-00546-5 | review | marker_registry (framework) |
| 13 | Fibroblasts: origins, definitions, and functions in health and disease | 2021 | Cell | 10.1016/j.cell.2021.06.024 | review | marker_registry (framework) |
| 14 | Molecular features of CAF subtypes | 2021 | Clin Cancer Res | 10.1158/1078-0432.CCR-20-4226 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 15 | Cross-tissue single-cell landscape of human monocytes and macrophages | 2021 | Immunity | 10.1016/j.immuni.2021.07.007 | single_cell_atlas | marker_registry + marker_tumor (covered in Batch 07) |
| 16 | Integration of pan-cancer single-cell and spatial transcriptomics reveals stromal cell features and therapeutic targets in tumor microenvironment | 2024 | Cancer Research | 10.1158/0008-5472.CAN-23-1418 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |

---

## Overview: Fibroblast/CAF Biology Context

### Existing Registry Entries (Current Status)

| Entry | Current Markers | Evidence Tier | Notes |
|-------|----------------|---------------|-------|
| Fibroblasts | FN1, COL1A1, COL3A1, COL1A2, DCN, LUM, PDGFRA, FAP, THY1, VIM | curated_review | Broad identity |
| iCAF | CFD, APOD, IGF1, CXCL12, MGP, FAP | seed | Inflammatory CAF |
| mCAF | POSTN, FN1, VCAN, CTHRC1, COL1A2 | seed | Matrix CAF |
| pCAF | SCG2, IGFBP2, HIST1H4C, PLAU | seed | Proliferative CAF |
| meCAF | NDRG1, BNIP3, TMEM158, ADM | seed | Metabolic CAF |
| EndMTCAF | PLVAP, RAMP2, RGCC | seed | Endothelial-derived |
| pnCAF | S100B, GPM6B, PLP1 | seed | Perineural CAF |
| apCAF | HLA-DRA, LYZ, TYROBP | seed | Antigen-presenting CAF |

### Batch 01 Additions (Cross-tissue atlas)

| Entry | Markers | Source |
|-------|---------|--------|
| PI16+ Universal Fibroblast | PI16, FAP, THY1 | Cross-tissue atlas |
| S01-S16 fibroblast subtypes | PI16, C7, LAMC1, CD9, DPEP1, etc. | Cross-tissue atlas |

---

## Paper-by-Paper Curation

### Article 4: Science Advances 2025 — Fibroblast Atlas

**Dataset**: 249,156 fibroblasts, 73 studies, 10 tissues
**Key finding**: 18 fibroblast subtypes; **TSPAN8+ chromatin remodeling fibroblasts** (novel)

#### New Subtype Entry (`marker_registry`)

**TSPAN8+ Chromatin Remodeling Fibroblast**

```toml
[[subtype]]
name = "TSPAN8+ chromatin remodeling fibroblast"
color = "#8B008B"
markers = ["TSPAN8", "KDM6B", "KDM5B", "EZH2", "KMT2D", "BRD4", "SMARCA4", "CHD7"]
negative_markers = ["PI16", "FAP", "THY1", "ACTA2", "TAGLN"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "fibroblast"
scope = "all"
applies_to = ["all"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Fibroblast atlas: shared and specific cell types across tissues", year = "2025", doi = "10.1126/sciadv.ado0173" }
notes = """
TSPAN8+ chromatin remodeling fibroblasts. Novel population identified in pan-tissue atlas.
Express histone modification and chromatin remodeling genes (KDM6B, KDM5B, EZH2, KMT2D, BRD4, SMARCA4, CHD7).
Higher scores in cell differentiation and resident fibroblast states.
Interact with endothelial cells and T cells via VEGFA-F2R ligand-receptor pair.
Associated with poor prognosis across tissues.
Negative universal fibroblast (PI16, FAP, THY1) and myofibroblast (ACTA2, TAGLN) markers.
"""
```

---

### Article 5: Cancer Cell 2024 — Cross-tissue Fibroblast Atlas

**Dataset**: 517 samples, 11 tissue types
**Key finding**: 4 transcriptionally distinct myofibroblast subpopulations

| Subpopulation | Markers | Function | Prognosis |
|--------------|---------|----------|-----------|
| **LRRC15+ myofibroblast** | LRRC15, ACTA2, TAGLN, MYL9 | Terminally differentiated, pro-tumor | Poor |
| **MMP1+ myofibroblast** | MMP1, MMP3, MMP10, MMP13 | ECM remodeling, pro-tumor | Poor |
| **PI16+ fibroblast** | PI16, FAP, THY1, CD34 | Universal/anti-tumor | Favorable |
| Other myofibroblasts | ACTA2, TAGLN, COL1A1 | General myofibroblast | Context-dependent |

#### New/Updated Subtype Entries (`marker_registry`)

**LRRC15+ Myofibroblast**

```toml
[[subtype]]
name = "LRRC15+ myofibroblast"
color = "#8B0000"
markers = ["LRRC15", "ACTA2", "TAGLN", "MYL9", "TPM2", "DES", "CNN1"]
negative_markers = ["PI16", "FAP", "THY1", "CD34", "CFD", "APOD"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "fibroblast"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Cross-tissue human fibroblast atlas reveals myofibroblast subtypes with distinct roles in immune modulation", year = "2024", doi = "10.1016/j.ccell.2024.08.020" }
notes = """
LRRC15+ myofibroblasts. Terminally differentiated CAF subpopulation.
Express smooth muscle contractile proteins (ACTA2, TAGLN, MYL9, TPM2, DES, CNN1).
Contribute to immune-excluded and immune-suppressive tumor microenvironments.
Strong pro-tumor potential. Associated with poor prognosis.
Negative PI16/FAP/THY1 (universal fibroblast) and CD34 (resting fibroblast).
Negative CFD/APOD (iCAF markers).
"""
```

**MMP1+ Myofibroblast**

```toml
[[subtype]]
name = "MMP1+ myofibroblast"
color = "#FF4500"
markers = ["MMP1", "MMP3", "MMP10", "MMP13", "MMP14", "COL1A1", "COL1A2", "FN1"]
negative_markers = ["PI16", "FAP", "THY1", "CFD", "APOD", "HLA-DRA"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "fibroblast"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Cross-tissue human fibroblast atlas reveals myofibroblast subtypes with distinct roles in immune modulation", year = "2024", doi = "10.1016/j.ccell.2024.08.020" }
notes = "MMP1+ myofibroblasts. ECM remodeling CAF subpopulation. High matrix metalloproteinase expression (MMP1, MMP3, MMP10, MMP13, MMP14). Collagen synthesis (COL1A1, COL1A2) and fibronectin (FN1). Contribute to immune-excluded and immune-suppressive TMEs. Negative PI16/FAP/THY1 and iCAF/apCAF markers."
```

**PI16+ Fibroblast (Anti-tumor)**

```toml
[[subtype]]
name = "PI16+ fibroblast"
color = "#228B22"
markers = ["PI16", "FAP", "THY1", "CD34", "DPP4", "NT5E", "ENG"]
negative_markers = ["ACTA2", "TAGLN", "LRRC15", "MMP1", "MMP3", "POSTN"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "fibroblast"
scope = "all"
applies_to = ["all"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = [
    { title = "Cross-tissue human fibroblast atlas reveals myofibroblast subtypes with distinct roles in immune modulation", year = "2024", doi = "10.1016/j.ccell.2024.08.020" },
    { title = "Cross-tissue single-cell landscape of human monocytes and macrophages in health and disease", year = "2021", doi = "10.1016/j.immuni.2021.07.007" }
]
notes = "PI16+ fibroblasts. Universal/resting fibroblast with anti-tumor potential. Present in adjacent non-cancerous regions. Express universal fibroblast markers (PI16, FAP, THY1, CD34, DPP4, NT5E, ENG). Show anti-tumor functions. Negative myofibroblast markers (ACTA2, TAGLN, LRRC15, MMP1, MMP3, POSTN)."
```

---

### Article 2: Cancer Cell 2025 — Antigen-Presenting CAF Niches

**Dataset**: 15 tissue types and solid tumors
**Key finding**: Two distinct apCAF populations: **mesothelial-like** and **fibrocyte-like**

| apCAF Type | Location | Function | Markers |
|-----------|----------|----------|---------|
| Mesothelial-like apCAF | Near cancer cells | SPP1-mediated tumor promotion | MSLN, WT1, KRT5, KRT7 |
| Fibrocyte-like apCAF | Lymphocyte-enriched niches | SPP1-mediated therapy resistance | COL1A1, COL1A2, CD34, THY1 |

#### Updated Subtype Entry (`marker_registry`)

**Update apCAF entry** to reflect two subtypes:

```toml
[[subtype]]
name = "Antigen-presenting CAF (apCAF)"
color = "#4682B4"
markers = ["HLA-DRA", "HLA-DQA1", "CD74", "LYZ", "TYROBP", "SPP1"]
negative_markers = ["ACTA2", "TAGLN", "PI16", "FAP"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "fibroblast"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = [
    { title = "Single-cell resolution spatial analysis of antigen-presenting cancer-associated fibroblast niches", year = "2025", doi = "10.1016/j.ccell.2025.09.001" },
    { title = "The advent of immune stimulating CAFs in cancer", year = "2023", doi = "10.1038/s41568-023-00549-7" }
]
notes = """
Antigen-presenting CAFs (apCAFs). Two distinct subtypes identified:
1. Mesothelial-like apCAFs: MSLN+, WT1+, KRT5+, KRT7+. Located near cancer cells.
2. Fibrocyte-like apCAFs: COL1A1+, COL1A2+, CD34+, THY1+. Associated with lymphocyte-enriched niches.
Both express MHC-II (HLA-DRA, HLA-DQA1, CD74) and SPP1 (osteopontin).
SPP1 facilitates primary tumor formation, peritoneal metastasis, and therapy resistance.
Can present antigens to T cells (immunostimulatory function per Nat Rev Cancer 2023).
Negative myofibroblast markers (ACTA2, TAGLN, PI16, FAP).
"""
```

---

### Article 1: Cancer Cell 2026 — ADAM12 as Fibroblast Checkpoint

**Key finding**: ADAM12 mediates balance between IFN-I response and TGF-β-driven myofibroblast activation
**Mechanism**: ADAM12 ablation → IFN-I response → progenitor-like state → T cell revitalization → tumor rejection

#### New Subtype Entry (`marker_registry`)

**ADAM12+ Fibroblast Checkpoint**

```toml
[[subtype]]
name = "ADAM12+ fibroblast checkpoint"
color = "#800080"
markers = ["ADAM12", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "ACTA2", "TAGLN", "COL1A1"]
negative_markers = ["IFNAR1", "IFNAR2", "STAT1", "IRF7", "MX1", "OAS1"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "fibroblast"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Single-cell screens identify ADAM12 as a fibroblast checkpoint impeding anti-tumor immunity", year = "2026", doi = "10.1016/j.ccell.2025.12.018" }
notes = """
ADAM12+ fibroblasts represent a molecular checkpoint in the TGF-β-driven myofibroblast activation axis.
ADAM12 mediates the balance between pro-tumoral myofibroblast activation and anti-tumoral IFN-I response.
TGF-β signaling markers (TGFBR1, TGFBR2, SMAD2, SMAD3, ACTA2, TAGLN, COL1A1).
ADAM12 ablation elicits IFN-I-responsive programs, reconfigures myofibroblasts into progenitor-like states, and revitalizes T cell responses.
Negative IFN-I response markers (IFNAR1, IFNAR2, STAT1, IRF7, MX1, OAS1) — these are induced upon ADAM12 loss.
Therapeutic target: ADAM12 inhibition may convert pro-tumor CAFs to anti-tumor phenotype.
"""
```

#### Functional Programs (geneset)

```json
"CAF_ADAM12_checkpoin": {
  "genes": ["ADAM12", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "TGFB1", "TGFB2", "TGFB3", "ACTA2", "TAGLN", "COL1A1", "COL1A2", "FN1", "POSTN"],
  "description": "ADAM12+ fibroblast checkpoint program. TGF-β-driven myofibroblast activation. Pro-tumoral."
},
"CAF_IFN_I_response": {
  "genes": ["IFNAR1", "IFNAR2", "STAT1", "STAT2", "IRF7", "IRF9", "MX1", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT3", "ISG15"],
  "description = "CAF IFN-I response program. Anti-tumoral. Induced by ADAM12 ablation."
},
"CAF_progenitor_state": {
  "genes": ["PI16", "FAP", "THY1", "CD34", "NT5E", "ENG", "DPP4", "PDGFRA"],
  "description": "CAF progenitor-like state program. Anti-tumoral reversion state upon ADAM12 loss."
}
```

---

### Article 3: Cancer Cell 2025 — Conserved Spatial CAF Subtypes

**Dataset**: 14M+ cells, 10 cancer types, 7 spatial platforms
**Key finding**: **4 conserved spatial CAF subtypes**

| Spatial CAF Subtype | Key Features | Spatial Pattern |
|-------------------|-------------|-----------------|
| **ECM-remodeling CAF** | COL1A1, COL3A1, POSTN, FN1 | Tumor-stroma interface |
| **Inflammatory CAF** | CXCL12, CXCL14, CCL2, IL6 | Immune cell neighborhoods |
| **Angiogenic CAF** | VEGFA, ANGPT2, PGF, KDR | Near blood vessels |
| **Antigen-presenting CAF** | HLA-DRA, CD74, SPP1 | Lymphocyte aggregates |

#### Tumor Context (`marker_tumor`)

**ECM-Remodeling Spatial CAF**

```toml
[[cancer_state]]
name = "ECM-remodeling spatial CAF"
color = "#CD853F"
markers = ["COL1A1", "COL3A1", "POSTN", "FN1", "CTHRC1", "MMP11", "LOX", "LOXL2"]
negative_markers = ["CXCL12", "CXCL14", "VEGFA", "HLA-DRA", "CD74"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["all"]
scope = "tumor_context"
applies_to = ["all"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Conserved spatial subtypes and cellular neighborhoods of cancer-associated fibroblasts revealed by single-cell spatial multi-omics", year = "2025", doi = "10.1016/j.ccell.2025.03.004" }
notes = "ECM-remodeling spatial CAF. Located at tumor-stroma interface. Collagen synthesis (COL1A1, COL3A1) and cross-linking (LOX, LOXL2). Matrix organization (POSTN, FN1, CTHRC1, MMP11). Negative inflammatory (CXCL12, CXCL14), angiogenic (VEGFA), and antigen-presenting (HLA-DRA, CD74) markers."
```

**Inflammatory Spatial CAF**

```toml
[[cancer_state]]
name = "Inflammatory spatial CAF"
color = "#FF6347"
markers = ["CXCL12", "CXCL14", "CCL2", "CCL19", "IL6", "IL1B", "PTGS2"]
negative_markers = ["COL1A1", "POSTN", "VEGFA", "HLA-DRA", "ACTA2"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["all"]
scope = "tumor_context"
applies_to = ["all"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Conserved spatial subtypes and cellular neighborhoods of cancer-associated fibroblasts revealed by single-cell spatial multi-omics", year = "2025", doi = "10.1016/j.ccell.2025.03.004" }
notes = "Inflammatory spatial CAF. Located in immune cell neighborhoods. Chemokine production (CXCL12, CXCL14, CCL2, CCL19). Pro-inflammatory cytokines (IL6, IL1B, PTGS2). Shapes immunosuppressive microenvironment. Negative ECM (COL1A1, POSTN), angiogenic (VEGFA), and antigen-presenting (HLA-DRA) markers."
```

---

### Article 8: Mol Cancer 2023 — Pan-cancer Spatial CAF-TME Crosstalk

**Key finding**: mCAFs (matrix) in tumor angiogenesis; iCAFs (inflammatory) in immunosuppression
**iCAF score correlates with immunotherapy response in melanoma**

#### Validation/Update Existing Entries

| Existing Entry | Update | Source |
|---------------|--------|--------|
| mCAF (POSTN, FN1, VCAN, CTHRC1, COL1A2) | Add note: "Spatially located at tumor-stroma interface. Promotes tumor angiogenesis." | Mol Cancer 2023 |
| iCAF (CFD, APOD, IGF1, CXCL12, MGP, FAP) | Add note: "Located in immune cell neighborhoods. Shapes immunosuppressive TME. iCAF score correlates with ICB response in melanoma." | Mol Cancer 2023 |

---

### Article 14: Clin Cancer Res 2021 — Molecular Features of CAF Subtypes

**Key finding**: 6 pan-cancer CAF subtypes; transcriptional drivers: MEF2C, TWIST1, NR1H3, RELB, FOXM1

| Pan-CAF Subtype | Markers | Function |
|----------------|---------|----------|
| CAF-S1 | CXCL12, CXCL14 | Immunosuppressive |
| CAF-S2 | IL6, LIF | Stem cell-promoting |
| CAF-S3 | MEF2C, TWIST1 | EMT-promoting |
| CAF-S4 | NR1H3, RELB | Metabolic |
| CAF-S5 | FOXM1, MKI67 | Proliferative |
| CAF-S6 | ACTA2, TAGLN | Myofibroblast |

#### Functional Programs (geneset)

```json
"CAF_immunosuppressive": {
  "genes": ["CXCL12", "CXCL14", "CCL2", "CCL5", "IL6", "TGFB1", "TGFB2", "TGFB3", "VEGFA", "IDO1"],
  "description": "Immunosuppressive CAF program. CXCL12/CXCL14-driven. Associated with ICB resistance."
},
"CAF_EMT_promoting": {
  "genes": ["MEF2C", "TWIST1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TGFB1", "TGFB2", "WNT5A", "FAP"],
  "description": "EMT-promoting CAF program. MEF2C/TWIST1-driven. Enhances cancer invasiveness."
},
"CAF_metabolic": {
  "genes": ["NR1H3", "RELB", "PPARG", "FABP4", "ADIPOQ", "LEP", "CPT1A", "ACLY", "FASN"],
  "description = "Metabolic CAF program. NR1H3/RELB-driven. Lipid metabolism and energy support."
}
```

---

### Article 16: Cancer Research 2024 — Pan-Cancer Stromal Cell Features (Du et al.)

**Dataset**: 214,972 nonimmune stromal cells from 258 patients across 16 cancer types; spatial transcriptomics from 16 patients across 7 cancer types (including 6 anti–PD-1 treated)
**Key finding**: 39 stromal subsets with distinct functional modules, spatial locations, and clinical features; tumor-associated PGF+ tip cells enriched in immune-depleted TME

#### Stromal Major Types and Subsets

The analysis identified **6 common major stromal types** across cancer types:

| Major Type | Subsets | Key Markers | Function |
|-----------|---------|-------------|----------|
| Fibroblasts | 14 subsets | COL1A1, COL1A2, DCN, LUM, PDGFRA, FAP | ECM production, immune modulation |
| Endothelial cells | 10 subsets | PECAM1, VWF, CDH5, KDR, FLT1 | Angiogenesis, nutrient delivery |
| Mural cells | 6 subsets | RGS5, PDGFRB, ACTA2, TAGLN | Vessel stability, contractility |
| Mesothelial cells | 3 subsets | MSLN, WT1, KRT5, KRT7 | Peritoneal lining, barrier function |
| Glial cells | 3 subsets | CDH19, SOX10, S100B | Neural support (schwannoma) |
| Other | 3 subsets | — | Mixed/transition states |

#### Tumor-Associated PGF+ Tip Endothelial Cells

**Key finding**: PGF+ tip cells are enriched in immune-depleted TME and associated with anti–PD-1 resistance.

```toml
[[cancer_state]]
name = "PGF+ tip endothelial cell"
color = "#FF6347"
markers = ["PGF", "KDR", "ESM1", "COL4A1", "ANGPT2", "DLL4"]
negative_markers = ["ACKR1", "SELP", "NRG1", "ICAM1"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["all"]
scope = "tumor_context"
applies_to = ["all"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Integration of pan-cancer single-cell and spatial transcriptomics reveals stromal cell features and therapeutic targets in tumor microenvironment", year = "2024", doi = "10.1158/0008-5472.CAN-23-1418" }
notes = """
Tumor-associated PGF+ tip endothelial cells. Enriched in immune-depleted TME.
Express tip cell markers (PGF, KDR/VEGFR2, ESM1, COL4A1, ANGPT2, DLL4).
Associated with anti–PD-1 resistance in spatial transcriptomics analysis.
Negative venous markers (ACKR1, SELP, NRG1) and adhesion molecules (ICAM1).
Therapeutic target: PGF inhibition may overcome anti–PD-1 resistance.
"""
```

#### Spatial Stromal Features (from Anti–PD-1 Treated Patients)

| Spatial Feature | Stromal Cell Type | Clinical Association |
|----------------|-------------------|---------------------|
| Immune-depleted TME | PGF+ tip ECs | Anti–PD-1 resistance |
| Immune-enriched TME | iCAF-like fibroblasts | Anti–PD-1 response |
| Tumor-stroma interface | mCAF-like fibroblasts | Invasion/metastasis |
| Perivascular niche | Pericytes, SMCs | Vessel stability |

#### Functional Programs (geneset)

```json
"Stromal_immune_depleted_TME": {
  "genes": ["PGF", "KDR", "ESM1", "ANGPT2", "DLL4", "COL4A1", "VEGFA", "NRP1", "NRP2"],
  "description": "Stromal program associated with immune-depleted TME and anti–PD-1 resistance. PGF+ tip endothelial cell signature."
},
"Stromal_immune_enriched_TME": {
  "genes": ["CXCL12", "CXCL14", "CCL2", "CCL19", "CCL21", "IL6", "IGF1", "CFD", "APOD"],
  "description": "Stromal program associated with immune-enriched TME and anti–PD-1 response. iCAF-like and immune-recruiting signature."
},
"Stromal_angiogenic_tip": {
  "genes": ["PGF", "KDR", "ESM1", "FLT1", "FLT4", "ANGPT1", "ANGPT2", "TEK", "NRP1", "NRP2"],
  "description": "Angiogenic tip endothelial cell program. PGF-driven pro-angiogenic signature."
}
```

#### Validation/Update Existing Entries

| Existing Entry | Update | Source |
|---------------|--------|--------|
| Endothelial cells (PECAM1, VWF, CDH5) | Add note: "Tip cells (PGF+, KDR+, ESM1+) enriched in immune-depleted TME" | Cancer Res 2024 |
| iCAF (CFD, APOD, IGF1, CXCL12, MGP, FAP) | Add note: "iCAF-like fibroblasts enriched in immune-enriched TME; associated with anti–PD-1 response" | Cancer Res 2024 |
| mCAF (POSTN, FN1, VCAN, CTHRC1, COL1A2) | Add note: "mCAF-like fibroblasts at tumor-stroma interface; associated with invasion" | Cancer Res 2024 |

#### Key Biological Insights

1. **PGF+ tip cells are immunosuppressive**: Tumor-associated PGF+ tip endothelial cells create an immune-depleted microenvironment and are associated with anti–PD-1 resistance. scLucid should flag these cells as markers of immunotherapy resistance.

2. **Spatial context determines stromal function**: The same stromal cell type can have different functions depending on spatial location. scLucid should integrate spatial transcriptomics data when available.

3. **Stromal-targeted therapy opportunities**: PGF inhibition, anti-angiogenic therapy, and CAF reprogramming are promising strategies. scLucid tumor interpretation should suggest stromal-targeted therapies based on stromal composition.

---

## Updated Fibroblast/CAF Registry Hierarchy

### Proposed Structure

```
Fibroblasts (compartment) [EXPANDED]
├── Universal/Resting Fibroblast
│   └── PI16+ fibroblast [NEW/UPDATED]
├── Myofibroblast
│   ├── LRRC15+ myofibroblast [NEW]
│   ├── MMP1+ myofibroblast [NEW]
│   └── General myofibroblast [EXISTING]
├── CAF Subtypes
│   ├── mCAF (matrix) [EXISTING]
│   ├── iCAF (inflammatory) [EXISTING]
│   ├── meCAF (metabolic) [EXISTING]
│   ├── pCAF (proliferative) [EXISTING]
│   ├── apCAF (antigen-presenting) [UPDATED]
│   │   ├── Mesothelial-like apCAF [NEW]
│   │   └── Fibrocyte-like apCAF [NEW]
│   ├── EndMTCAF [EXISTING]
│   └── pnCAF [EXISTING]
├── Novel Subtypes
│   ├── TSPAN8+ chromatin remodeling fibroblast [NEW]
│   └── ADAM12+ fibroblast checkpoint [NEW]
├── Spatial CAF States
│   ├── ECM-remodeling spatial CAF [cancer_state]
│   ├── Inflammatory spatial CAF [cancer_state]
│   ├── Angiogenic spatial CAF [cancer_state]
│   └── Antigen-presenting spatial CAF [cancer_state]
└── Tumor contexts
    └── Various cancer_state entries
```

---

## Curation Notes and Conflicts

### 1. CAF Classification Systems

**Issue**: Multiple classification systems:
- scLucid existing: mCAF, iCAF, meCAF, pCAF, apCAF, EndMTCAF, pnCAF
- Cancer Cell 2024: LRRC15+, MMP1+, PI16+ myofibroblasts
- Cancer Cell 2025: 4 spatial subtypes (ECM, inflammatory, angiogenic, antigen-presenting)
- Clin Cancer Res 2021: 6 pan-CAF subtypes (CAF-S1 to S6)
- Mol Cancer 2023: mCAF, iCAF

**Resolution**:
- **Adopt Cancer Cell 2024 framework** (cross-tissue validated) as primary:
  - LRRC15+ myofibroblast (pro-tumor)
  - MMP1+ myofibroblast (pro-tumor)
  - PI16+ fibroblast (anti-tumor)
- Map existing entries:
  - mCAF ≈ MMP1+ myofibroblast (ECM remodeling)
  - iCAF ≈ Inflammatory spatial CAF (chemokine-producing)
  - meCAF ≈ Metabolic CAF (retain existing)
  - apCAF ≈ Antigen-presenting spatial CAF (updated with two subtypes)
- Keep Clin Cancer Res 2021 subtypes as **geneset programs** (not identity markers)

### 2. apCAF Heterogeneity

**Issue**: Two distinct apCAF populations (Cancer Cell 2025) vs single apCAF entry (existing).

**Resolution**:
- Update apCAF entry to include both subtypes
- Use additional markers to distinguish:
  - Mesothelial-like: MSLN, WT1, KRT5, KRT7
  - Fibrocyte-like: COL1A1, COL1A2, CD34, THY1
- Shared: HLA-DRA, HLA-DQA1, CD74, SPP1

### 3. ADAM12 as Therapeutic Target

**Issue**: ADAM12 is a molecular checkpoint, not a stable identity marker.

**Resolution**:
- Add as subtype entry with therapeutic context
- Note: "ADAM12 is a druggable target. Inhibition may convert pro-tumor CAFs to anti-tumor phenotype."
- Include both TGF-β pathway markers (pro-tumor) and IFN-I markers (anti-tumor, negative)

### 4. TSPAN8+ Fibroblasts

**Issue**: TSPAN8+ chromatin remodeling fibroblasts are newly discovered.

**Resolution**:
- Add as new subtype with chromatin remodeling markers
- Note poor prognosis association
- Document interaction with endothelial cells (VEGFA) and T cells (F2R)

### 5. Spatial Context

**Issue**: Spatial CAF subtypes have distinct locations (Cancer Cell 2025).

**Resolution**:
- Add spatial location in notes:
  - ECM-remodeling: tumor-stroma interface
  - Inflammatory: immune cell neighborhoods
  - Angiogenic: near blood vessels
  - Antigen-presenting: lymphocyte aggregates
- Use for spatial transcriptomics interpretation

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New subtype entries | 5 (LRRC15+ myofibroblast, MMP1+ myofibroblast, TSPAN8+ fibroblast, ADAM12+ checkpoint, updated apCAF) | marker_registry |
| New tumor context entries | 5 (ECM-remodeling, inflammatory, angiogenic, antigen-presenting spatial CAFs, PGF+ tip EC) | marker_tumor |
| Existing entry updates | 6 (apCAF, mCAF, iCAF, endothelial, pericyte, SMC notes) | marker_registry |
| New geneset programs | 9 (ADAM12 checkpoint, IFN-I response, CAF progenitor, immunosuppressive, EMT-promoting, metabolic, immune-depleted TME, immune-enriched TME, angiogenic tip) | genesets_cancer_signatures.json |
| Review validations | 6 | marker_registry |
| Conflict resolutions | 5 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **CAFs are not all pro-tumor**: PI16+ fibroblasts show anti-tumor functions (Cancer Cell 2024). scLucid should not assume all CAFs are pro-tumoral.

2. **ADAM12 is a therapeutic target**: ADAM12 inhibition converts pro-tumor CAFs to anti-tumor phenotype via IFN-I response (Cancer Cell 2026). scLucid tumor interpretation should flag ADAM12+ CAFs as druggable targets.

3. **apCAFs are immunostimulatory**: Both mesothelial-like and fibrocyte-like apCAFs can present antigens to T cells (Cancer Cell 2025, Nat Rev Cancer 2023). This challenges the immunosuppressive CAF dogma.

4. **Spatial context matters**: CAF subtypes have distinct spatial locations (Cancer Cell 2025). ECM-remodeling CAFs are at tumor-stroma interface; inflammatory CAFs are near immune cells.

5. **TSPAN8+ fibroblasts are prognostic**: TSPAN8+ chromatin remodeling fibroblasts are associated with poor prognosis (Science Advances 2025). May be a novel therapeutic target.
