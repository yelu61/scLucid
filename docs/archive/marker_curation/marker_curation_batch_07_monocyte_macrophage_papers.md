> **⚠️ ARCHIVED / SUPERSEDED**
>
> This batch curation note is kept for provenance only. The live marker
> resource status is tracked in `docs/marker_curation_literature_index.jsonl`,
> `docs/marker_resource_quality_gaps.jsonl`, `docs/marker_curation_candidates.jsonl`,
> and `docs/MARKER_RESOURCE_CURATION.md`. New curation should follow the current
> contract rather than adding new files to this archive.

---

# Marker Curation Batch 07: Monocyte/Macrophage Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | Using a pan-cancer atlas to investigate tumour associated macrophages as regulators of immunotherapy response | 2024 | Nature Communications | 10.1038/s41467-024-49885-8 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 2 | Timing and location dictate monocyte fate and their transition to tumor-associated macrophages | 2024 | Science Immunology | 10.1126/sciimmunol.adk3981 | single_cell_atlas | marker_registry + marker_tumor |
| 3 | Decoding the spatiotemporal heterogeneity of tumor-associated macrophages | 2024 | Molecular Cancer | 10.1186/s12943-024-02064-1 | review | marker_registry (validation) |
| 4 | Coordinated chemokine expression defines macrophage subsets across tissues | 2024 | Nature Immunology | 10.1038/s41590-024-01826-9 | single_cell_atlas (mouse) | reference_only |
| 5 | Deciphering the performance of macrophages in tumour microenvironment | 2024 | Journal of Hematology and Oncology | 10.1186/s13045-024-01559-0 | review | marker_registry (validation) |
| 6 | An immune cell atlas reveals the dynamics of human macrophage specification during prenatal development | 2023 | Cell | 10.1016/j.cell.2023.08.019 | single_cell_atlas | marker_registry + geneset |
| 7 | Macrophages at the interface of the co-evolving cancer ecosystem | 2023 | Cell | 10.1016/j.cell.2023.02.020 | review | marker_registry (validation) |
| 8 | Tumor macrophage functional heterogeneity can inform the development of novel cancer therapies | 2023 | Trends in Immunology | 10.1016/j.it.2023.10.007 | review | marker_registry (validation) |
| 9 | Macrophage diversity in cancer revisited in the era of single-cell omics | 2022 | Trends in Immunology | 10.1016/j.it.2022.04.008 | review | marker_registry (naming) |
| 10 | Macrophages in health and disease | 2022 | Cell | 10.1016/j.cell.2022.10.007 | review | marker_registry (framework) |
| 11 | Cross-tissue single-cell landscape of human monocytes and macrophages in health and disease | 2021 | Immunity | 10.1016/j.immuni.2021.07.007 | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 12 | Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors | 2017 | Science | 10.1126/science.aah4573 | single_cell_atlas | marker_registry (foundational) |

**Note**: Article 12 was curated in Batch 06 (DC cells). This batch focuses on monocyte/macrophage-specific findings.

---

## Overview: Monocyte/Macrophage Context

### Existing Registry Entries (Current Status)

| Entry | Current Markers | Evidence Tier | Context |
|-------|----------------|---------------|---------|
| Monocytes | FCN1, S100A8, S100A9, CD14, VCAN, LYZ, S100A12, CD36 | curated_review | Blood/tissue |
| CD14+ Mono | FCN1, HLA-DQB1, S100A9, CSF3R, S100A8 | atlas_supported | Classical |
| CD16+ Mono | FCGR3A, LST1, LILRB2, HK3 | atlas_supported | Non-classical |
| CD14+CD16+ Mono | NFKBIA, NFKB1, NLRP3, HLA-DQA1 | seed | Intermediate |
| Macrophages | CD68, CD163, APOE, C1QA, C1QB, MSR1, MRC1, CD14, CSF1R | curated_review | Broad |
| Reg TAM | ARG1, MRC1, CD274, CX3CR1 | atlas_supported | Tumor |
| Inflam TAM | IL1B, CCL3, CXCL1, CXCL2 | atlas_supported | Tumor |
| IFN TAM | IDO1, ISG15, CXCL10, STAT1 | atlas_supported | Tumor |
| Angio TAM | VEGFA, SPP1, FGF2, MMP9 | atlas_supported | Tumor |
| Prolif TAM | MKI67, CDK1, PCNA | atlas_supported | Tumor |
| LA TAM | APOC1, APOE, FABP5 | atlas_supported | Tumor |
| TRM TAM | LYVE1, FOLR2, CX3CR1, MERTK | atlas_supported | Tumor |
| FOLR2+ Resident Mac | FOLR2, LYVE1, MERTK | atlas_supported | Cross-tissue |
| CD5L+ Macrophage | CD5L, MARCO, MERTK | atlas_supported | Cross-tissue |
| TREM2+ Macrophage | TREM2, CD9, LPL, GPNMB, APOE, CST7 | atlas_supported | Tumor |

### Batch 05 Additions

| Entry | Markers | Context |
|-------|---------|---------|
| Tumor-infiltrating classical mono | FCN1, S100A8, S100A9, CD14 | Tumor |
| Tumor-infiltrating non-classical mono | FCGR3A, LST1, LILRB2 | Tumor |
| ICB-responsive myeloid | CD80, CD86, HLA-DRA, CXCL9/10/11 | ICB response |
| ICB-nonresponsive myeloid | CD274, PDCD1LG2, ARG1, MRC1 | ICB resistance |

---

## Paper-by-Paper Curation

### Article 11: Immunity 2021 — MNP-VERSE (Cross-tissue Monocyte/Macrophage Atlas)

**Dataset**: 178,651 MNPs from 13 tissues, 41 datasets
**Key finding**: IL4I1+ CD274+ IDO1+ macrophages (IL4I1_Macs) — tumor periphery, immunosuppressive, T cell-dependent

#### New Subtype Entry (`marker_registry`)

**IL4I1+ Macrophage**

```toml
[[subtype]]
name = "IL4I1+ macrophage"
color = "#8B0000"
markers = ["IL4I1", "CD274", "IDO1", "CD40", "CD83", "CCR7", "FSCN1"]
negative_markers = ["FCN1", "S100A8", "S100A9", "LYZ", "VCAN"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tumor_context"
applies_to = ["Macrophages", "Tumor"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue single-cell landscape of human monocytes and macrophages in health and disease", year = "2021", doi = "10.1016/j.immuni.2021.07.007" }
notes = """
IL4I1+ macrophages (IL4I1_Macs) identified in MNP-VERSE compendium.
Accumulate in tumor periphery in a T cell-dependent manner via IFN-γ and CD40/CD40L signaling.
Immunosuppressive through tryptophan degradation (IDO1+).
Promote Treg entry into tumors.
Co-express CD274 (PD-L1) and IDO1 — dual immunosuppressive mechanism.
Negative FCN1/S100A8/9/LYZ/VCAN distinguishes from inflammatory monocytes.
"""
```

#### Tumor Context (`marker_tumor`)

**Tumor Periphery IL4I1+ Macrophage**

```toml
[[cancer_state]]
name = "Tumor periphery IL4I1+ macrophage"
color = "#8B0000"
markers = ["IL4I1", "CD274", "IDO1", "CD40", "CD83", "CCR7", "FSCN1", "EBI3"]
negative_markers = ["FCN1", "S100A8", "S100A9", "LYZ", "CD14", "VEGFA"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["all"]
scope = "tumor_context"
applies_to = ["all"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue single-cell landscape of human monocytes and macrophages in health and disease", year = "2021", doi = "10.1016/j.immuni.2021.07.007" }
notes = "Tumor periphery IL4I1+ macrophages. IFN-γ and CD40/CD40L-induced maturation from IFN-primed monocytes. T cell-dependent accumulation. Immunosuppressive via tryptophan degradation (IDO1) and PD-L1 (CD274). Promotes Treg infiltration. Spatially restricted to tumor periphery."
```

#### Functional Programs (geneset)

```json
"IL4I1_macrophage_program": {
  "genes": ["IL4I1", "CD274", "IDO1", "CD40", "CD83", "CCR7", "FSCN1", "EBI3", "HLA-DRA", "HLA-DQA1", "CD86", "CD80"],
  "description": "IL4I1+ macrophage program. Tumor periphery, immunosuppressive, T cell-dependent accumulation."
},
"Tryptophan_degradation": {
  "genes": ["IDO1", "IDO2", "TDO2", "KYNU", "KMO", "IL4I1", "AFMID"],
  "description": "Tryptophan degradation program. Immunosuppressive via kynurenine pathway."
}
```

---

### Article 1: Nat Commun 2024 — Pan-Cancer TAM Atlas and ICB Response

**Dataset**: Pan-cancer scRNA-seq atlas of TAMs
**Key findings**:
- TAM composition varies between primary and metastatic tumors
- Macrophage-T cell functional cross-talk
- Two TAM subsets associated with T cell activation
- **Collagen-related gene-upregulating TAM subset associated with ICB response**

#### New Subtype Entry (`marker_registry`)

**Collagen-upregulating TAM**

```toml
[[subtype]]
name = "Collagen-upregulating TAM"
color = "#CD853F"
markers = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "LOX", "LOXL2", "P4HA1", "P4HA2"]
negative_markers = ["VEGFA", "SPP1", "MMP9", "ARG1", "MRC1"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tumor_context"
applies_to = ["Macrophages", "Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Using a pan-cancer atlas to investigate tumour associated macrophages as regulators of immunotherapy response", year = "2024", doi = "10.1038/s41467-024-49885-8" }
notes = """
Collagen-upregulating TAM subset identified in pan-cancer atlas.
Associated with immune checkpoint inhibitor (ICB) response.
Expresses collagen synthesis genes (COL1A1/2, COL3A1, COL5A1) and collagen cross-linking enzymes (LOX, LOXL2).
May promote T cell infiltration via ECM remodeling.
Negative VEGFA/SPP1/MMP9 distinguishes from angiogenic TAM.
Negative ARG1/MRC1 distinguishes from regulatory TAM.
"""
```

#### Tumor Context (`marker_tumor`)

**ICB-Response-Associated Collagen TAM**

```toml
[[cancer_state]]
name = "ICB-response collagen TAM"
color = "#228B22"
markers = ["COL1A1", "COL1A2", "COL3A1", "LOX", "LOXL2", "CD86", "HLA-DRA", "CCL19"]
negative_markers = ["ARG1", "MRC1", "CD274", "VEGFA", "SPP1"]

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
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Using a pan-cancer atlas to investigate tumour associated macrophages as regulators of immunotherapy response", year = "2024", doi = "10.1038/s41467-024-49885-8" }
notes = "Collagen-upregulating TAM associated with ICB response. ECM remodeling may facilitate T cell infiltration. Positive co-stimulatory markers (CD86, HLA-DRA) and T cell recruitment (CCL19). Negative immunosuppressive markers (ARG1, MRC1, CD274)."
```

#### Functional Programs (geneset)

```json
"TAM_collagen_program": {
  "genes": ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL6A1", "LOX", "LOXL2", "P4HA1", "P4HA2", "PLOD1", "PLOD2", "BMP1"],
  "description": "TAM collagen synthesis and ECM remodeling program. Associated with ICB response."
},
"TAM_T_cell_activation": {
  "genes": ["CD86", "HLA-DRA", "HLA-DQA1", "CD80", "CD40", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "IL12A", "IL12B"],
  "description": "TAM subset associated with T cell activation. Co-stimulatory and T cell-recruiting."
}
```

---

### Article 2: Science Immunology 2024 — Monocyte-to-TAM Transition in PDAC

**Dataset**: PDAC mouse model + human validation
**Key finding**: Monocyte → intermediate TAM → two terminally differentiated TAM lineages
**TF**: Maf-dependent differentiation

#### Tumor Context (`marker_tumor`)

**Intermediate TAM**

```toml
[[cancer_state]]
name = "Intermediate TAM"
color = "#DAA520"
markers = ["LYZ", "CD14", "FCN1", "CD68", "CSF1R", "CD163", "APOE"]
negative_markers = ["C1QA", "C1QB", "MRC1", "ARG1", "VEGFA", "SPP1"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["Pancreatic Cancer"]
scope = "cancer_type_specific"
applies_to = ["Pancreas"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Timing and location dictate monocyte fate and their transition to tumor-associated macrophages", year = "2024", doi = "10.1126/sciimmunol.adk3981" }
notes = "Transient intermediate TAM population in PDAC. Expresses monocyte markers (LYZ, CD14, FCN1) and early macrophage markers (CD68, CSF1R, CD163, APOE). Negative terminal differentiation markers (C1QA, C1QB, MRC1, ARG1, VEGFA, SPP1). Gives rise to two terminally differentiated TAM lineages."
```

**Terminal TAM Lineage 1 (Maf-dependent)**

```toml
[[cancer_state]]
name = "Terminal TAM lineage 1"
color = "#8B4513"
markers = ["MAF", "C1QA", "C1QB", "MRC1", "CD163", "LYVE1", "FOLR2"]
negative_markers = ["LYZ", "CD14", "FCN1", "VEGFA", "SPP1"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["Pancreatic Cancer"]
scope = "cancer_type_specific"
applies_to = ["Pancreas"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Timing and location dictate monocyte fate and their transition to tumor-associated macrophages", year = "2024", doi = "10.1126/sciimmunol.adk3981" }
notes = "Terminally differentiated TAM lineage 1 in PDAC. MAF-dependent. Expresses complement components (C1QA, C1QB) and M2-like markers (MRC1, CD163, LYVE1, FOLR2). Tissue-resident-like features. Negative monocyte markers (LYZ, CD14, FCN1). Negative angiogenic markers (VEGFA, SPP1)."
```

**Terminal TAM Lineage 2 (Maf-dependent)**

```toml
[[cancer_state]]
name = "Terminal TAM lineage 2"
color = "#A0522D"
markers = ["MAF", "VEGFA", "SPP1", "MMP9", "FN1", "TIMP1"]
negative_markers = ["LYZ", "CD14", "FCN1", "LYVE1", "FOLR2", "MRC1"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["Pancreatic Cancer"]
scope = "cancer_type_specific"
applies_to = ["Pancreas"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Timing and location dictate monocyte fate and their transition to tumor-associated macrophages", year = "2024", doi = "10.1126/sciimmunol.adk3981" }
notes = "Terminally differentiated TAM lineage 2 in PDAC. MAF-dependent. Angiogenic and pro-tumor (VEGFA+, SPP1+, MMP9+, FN1+, TIMP1+). Negative tissue-resident markers (LYVE1, FOLR2, MRC1). Negative monocyte markers."
```

#### Functional Programs (geneset)

```json
"TAM_differentiation_trajectory": {
  "genes": ["LYZ", "CD14", "FCN1", "CD68", "CSF1R", "CD163", "APOE", "MAF", "C1QA", "C1QB", "MRC1", "LYVE1", "FOLR2", "VEGFA", "SPP1", "MMP9"],
  "description = "TAM differentiation trajectory in PDAC. Monocyte → intermediate → terminal (resident-like or angiogenic)."
},
"TAM_Maf_program": {
  "genes": ["MAF", "NR4A1", "NR4A2", "NR4A3", "KLF2", "KLF4", "IRF4"],
  "description": "Maf-dependent TAM differentiation program. Terminal differentiation of monocyte-derived TAMs."
}
```

---

### Article 6: Cell 2023 — Prenatal Macrophage Atlas

**Dataset**: PCW 4-26, 19 tissues
**Key findings**:
- Microglia-like population in fetal epidermis, testicle, heart
- Proangiogenic macrophages (perivascular, yolk-sac-derived)

#### New Subtype Entry (`marker_registry`)

**Microglia-like Fetal Macrophage**

```toml
[[subtype]]
name = "Microglia-like fetal macrophage"
color = "#4682B4"
markers = ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "HEXB", "CST3", "CSF1R"]
negative_markers = ["CD14", "LYZ", "FCN1", "S100A8", "CD163", "MRC1"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tissue_specific"
applies_to = ["Epidermis", "Testicle", "Heart", "CNS"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "An immune cell atlas reveals the dynamics of human macrophage specification during prenatal development", year = "2023", doi = "10.1016/j.cell.2023.08.019" }
notes = "Microglia-like fetal macrophages identified in prenatal atlas. Present in fetal epidermis, testicle, and heart (not just CNS). Express microglial markers (CX3CR1, P2RY12, TMEM119, AIF1). Interact with neural crest cells, modulating melanocyte differentiation. Negative monocyte/macrophage markers (CD14, LYZ, FCN1, S100A8, CD163, MRC1). Yolk-sac-derived."
```

**Proangiogenic Macrophage**

```toml
[[subtype]]
name = "Proangiogenic macrophage"
color = "#FF6347"
markers = ["VEGFA", "ANGPT1", "ANGPT2", "TEK", "KDR", "FLT1", "PGF", "PECAM1"]
negative_markers = ["CD14", "LYZ", "FCN1", "CD68", "CD163"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tissue_specific"
applies_to = ["all"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "An immune cell atlas reveals the dynamics of human macrophage specification during prenatal development", year = "2023", doi = "10.1016/j.cell.2023.08.019" }
notes = "Proangiogenic fetal macrophages. Perivascular across fetal organs. Yolk-sac-derived. Express angiogenic factors (VEGFA, ANGPT1/2, TEK, KDR, FLT1, PGF). Share endothelial marker PECAM1 (due to perivascular location). Negative monocyte/macrophage markers (CD14, LYZ, FCN1, CD68, CD163). Developmental context; may be rare in adult tissues."
```

---

### Articles 3, 5, 7, 8, 9, 10: Reviews

These reviews provide important context but do not introduce new marker genes. Key concepts:

| Review | Key Concept | Registry Impact |
|--------|-------------|-----------------|
| Mol Cancer 2024 | TAM spatiotemporal heterogeneity | Validates spatial context importance |
| J Hematol Oncol 2024 | M1/M2 dichotomy limitations | Notes: avoid M1/M2 as identity markers |
| Cell 2023 | TAM in cancer evolution | Validates TAM functional diversity |
| Trends Immunol 2023 | TAM functional heterogeneity | Validates single-cell approaches |
| Trends Immunol 2022 | 7 TAM subtypes from pan-cancer | Naming reference (see below) |
| Cell 2022 | Ontogeny framework | Validates embryonic vs BM-derived |

#### Trends Immunology 2022 — Seven TAM Subtypes

From the pan-cancer analysis cited in Trends Immunol 2022:

| TAM Subtype | Markers | Proposed Registry Entry |
|-------------|---------|------------------------|
| INHBA+ TAM | INHBA, ACTA2, TAGLN | New subtype |
| C1QC+ TAM | C1QC, C1QA, C1QB | Existing (broad macrophage) |
| ISG15+ TAM | ISG15, MX1, OAS1 | New state (IFN-response) |
| NLRP3+ TAM | NLRP3, IL1B, PYCARD | Existing (Inflam TAM) |
| LYVE1+ TAM | LYVE1, FOLR2, MERTK | Existing (TRM TAM / FOLR2+ Mac) |
| SPP1+ TAM | SPP1, FN1, VEGFA | Existing (Angio TAM) |

**New entries from this analysis**:

**INHBA+ TAM**

```toml
[[subtype]]
name = "INHBA+ macrophage"
color = "#8B4513"
markers = ["INHBA", "ACTA2", "TAGLN", "MMP11", "FN1", "COL1A1"]
negative_markers = ["LYVE1", "FOLR2", "MRC1", "ARG1"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tumor_context"
applies_to = ["Macrophages", "Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = [
    { title = "Macrophage diversity in cancer revisited in the era of single-cell omics", year = "2022", doi = "10.1016/j.it.2022.04.008" },
    { title = "A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells", year = "2021", doi = "10.1016/j.cell.2021.01.010" }
]
notes = "INHBA+ TAM identified in pan-cancer analysis. Expresses TGF-β family member INHBA (Activin A) and myofibroblast-like markers (ACTA2, TAGLN, FN1, COL1A1). May promote tumor fibrosis and EMT. Negative tissue-resident markers (LYVE1, FOLR2, MRC1, ARG1)."
```

**ISG15+ TAM** (already exists as "IFN TAM", validate)

The existing "IFN TAM" entry (IDO1, ISG15, CXCL10, STAT1) partially overlaps. Add INHBA+ TAM as new, confirm ISG15+ as IFN TAM state.

---

### Article 4: Nature Immunology 2024 — Mouse Macrophage Subsets

**Dataset**: Mouse lung macrophages
**Key finding**: Tissue-resident IMs distinct from recruited macrophages (recMacs)

**Curation note**: This is a **mouse study**. Human orthologs should be used with caution.

| Mouse Population | Mouse Markers | Human Orthologs | Applicability |
|-----------------|---------------|-----------------|---------------|
| Tissue-resident IMs | Lyve1, Folr2, C1qa | LYVE1, FOLR2, C1QA | ✅ Validated |
| Recruited Macs (recMacs) | Fn1, Ccr2, Ly6c2 | FN1, CCR2, — | ⚠️ Partial |
| DC1 | Xcr1, Clec9a | XCR1, CLEC9A | ✅ Validated |
| DC2 (CD301b+) | Cd209a, Mgl2 | CD209A, CLEC10A | ⚠️ Partial |
| Inflammatory DC2 | — | — | ⚠️ Species-specific |

**Human-applicable markers**: LYVE1, FOLR2, C1QA, FN1, CCR2, XCR1, CLEC9A

Add note to existing entries: "Mouse study (Nature Immunology 2024) validates tissue-resident vs recruited distinction."

---

## Summary: Key Conflicts and Resolutions

### 1. M1/M2 Dichotomy

**Issue**: Multiple reviews still reference M1/M2 classification, but single-cell studies show this is oversimplified.

**Resolution**:
- **Do NOT** add M1 or M2 as registry entries
- Use specific subtype/state names instead (Inflam TAM, Reg TAM, Angio TAM, etc.)
- Add note: "M1/M2 dichotomy is oversimplified. Use specific scLucid subtypes and states."

### 2. Monocyte-to-Macrophage Continuum

**Issue**: Monocytes and macrophages exist on a differentiation continuum (Science Immunology 2024).

**Resolution**:
- Add intermediate states ("Tumor-infiltrating classical mono", "Intermediate TAM")
- Use negative markers to distinguish stages:
  - Monocyte: LYZ+, CD14+, FCN1+, S100A8/9+
  - Intermediate: LYZ+, CD14+, CD68+, CSF1R+
  - Terminal: LYZ-, C1QA+, C1QB+, MRC1+/VEGFA+

### 3. Embryonic vs Bone Marrow Origin

**Issue**: Macrophages have dual origins (Cell 2022, Cell 2023).

**Resolution**:
- Add notes to relevant entries about ontogeny
- Embryonic-derived: microglia, some tissue-resident macrophages
- BM-derived: monocyte-derived macrophages, most TAMs
- Do not add separate entries for origin (not RNA-detectable in most datasets)

### 4. IL4I1_Macs vs mregDC

**Issue**: IL4I1_Macs (Immunity 2021) and mregDC (Batch 06) both express CD274 and IDO1.

**Resolution**:
- IL4I1_Macs: macrophage lineage, tumor periphery, tryptophan degradation
- mregDC: DC lineage, mature state, antigen presentation
- Distinguish by lineage markers:
  - IL4I1_Macs: CD68+, CD163+, CSF1R+, IL4I1+
  - mregDC: CD74+, FLT3+, LAMP3+, FSCN1+

### 5. Collagen TAMs and CAFs

**Issue**: Collagen-upregulating TAMs (Nat Commun 2024) express COL1A1, COL1A2 — same as CAFs.

**Resolution**:
- Use immune markers to distinguish: CD68, CSF1R, CD163
- Use negative stromal markers: COL1A1 alone is not sufficient
- Add note: "Collagen TAMs express COL1A1/2 but are CD68+ CSF1R+. Do not confuse with CAFs (CD68-, PTPRC-)."

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New subtype entries | 4 (IL4I1+ Mac, Collagen TAM, INHBA+ TAM, Microglia-like fetal Mac) | marker_registry |
| New tumor context entries | 5 (Tumor periphery IL4I1+, ICB collagen TAM, Intermediate TAM, Terminal TAM L1/L2) | marker_tumor |
| New geneset programs | 6 | genesets_cancer_signatures.json |
| Validation notes | 6 (review validations) | marker_registry |
| Conflict resolutions | 5 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **Monocyte-to-TAM trajectory**: Monocytes do not directly become TAMs. They pass through an intermediate state (Science Immunology 2024). scLucid annotation should flag intermediate states.

2. **IL4I1_Macs are spatially restricted**: They accumulate in tumor periphery, not center (Immunity 2021). Spatial context matters for interpretation.

3. **Collagen TAMs and ICB response**: Collagen-upregulating TAMs are associated with ICB response (Nat Commun 2024). ECM remodeling may facilitate T cell infiltration.

4. **Fetal macrophages are distinct**: Microglia-like and proangiogenic fetal macrophages (Cell 2023) are yolk-sac-derived and developmentally distinct. Rare in adult tissues.

5. **M1/M2 is dead**: Single-cell studies reveal far more complexity than M1/M2. scLucid should use specific subtype/state names.
