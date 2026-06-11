# Marker Curation Batch 08: Neutrophil Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | From complexity to consensus: a roadmap for neutrophil classification | 2025 | Immunity | 10.1016/j.immuni.2025.07.011 | review | marker_registry (naming) |
| 2 | Neutrophil profiling illuminates anti-tumor antigen-presenting potency | 2024 | Cell | 10.1016/j.cell.2024.02.005 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 3 | Deterministic reprogramming of neutrophils within tumors | 2024 | Science | 10.1126/science.adf6493 | single_cell_atlas | marker_registry + marker_tumor |
| 4 | Neutrophil profiling illuminates anti-tumor antigen-presenting potency | 2024 | Cell | 10.1016/j.cell.2024.02.005 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 5 | Liver tumour immune microenvironment subtypes and neutrophil heterogeneity | 2022 | Nature | 10.1038/s41586-022-05400-x | single_cell_atlas | marker_tumor + geneset |
| 6 | High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer | 2022 | Cancer Cell | 10.1016/j.ccell.2022.10.008 | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 7 | Single-cell RNA-seq analysis reveals BHLHE40-driven pro-tumour neutrophils with hyperactivated glycolysis in pancreatic tumour microenvironment | 2022 | Gut | 10.1136/gutjnl-2021-326070 | single_cell_atlas | marker_registry + marker_tumor + geneset |
| 8 | Neutrophils in cancer: heterogeneous and multifaceted | 2022 | Nature Reviews Immunology | 10.1038/s41577-021-00571-6 | review | marker_registry (validation) |
| 9 | Neutrophil diversity and plasticity in tumour progression and therapy | 2020 | Nature Reviews Cancer | 10.1038/s41568-020-0281-y | review | marker_registry (framework) |

---

## Overview: Neutrophil Biology Context

### Existing Registry Entry (Current Status — Minimal)

| Entry | Current Markers | Evidence Tier | Notes |
|-------|----------------|---------------|-------|
| Neutrophils | FCGR3B, CXCR2, CSF3R | curated_review | Very minimal; needs major expansion |

**Critical gap**: The existing neutrophil entry is severely underdeveloped. These 8 papers provide extensive evidence for expansion.

---

## Paper-by-Paper Curation

### Article 2: Cell 2024 — Pan-Cancer Neutrophil Atlas (Antigen-Presenting Neutrophils)

**Dataset**: 17 cancer types, 225 samples, 143 patients
**Key finding**: **10 distinct neutrophil states** including inflammation, angiogenesis, and antigen presentation
**Novel discovery**: Antigen-presenting neutrophils (APNs) — leucine metabolism → H3K27ac → antigen presentation program

#### New Subtype/State Entries (`marker_registry`)

**Antigen-Presenting Neutrophil (APN)**

```toml
[[subtype]]
name = "Antigen-presenting neutrophil"
color = "#228B22"
markers = ["HLA-DRA", "HLA-DQA1", "CD74", "HLA-DPB1", "HLA-DMA", "CD80", "CD86", "CIITA"]
negative_markers = ["FCGR3B", "CXCR2", "S100A8", "S100A9", "VEGFA", "MMP9"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tumor_context"
applies_to = ["Neutrophils", "Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Neutrophil profiling illuminates anti-tumor antigen-presenting potency", year = "2024", doi = "10.1016/j.cell.2024.02.005" }
notes = """
Antigen-presenting neutrophils (APNs) identified in pan-cancer atlas (17 cancer types).
Express MHC-II molecules (HLA-DRA, HLA-DQA1, CD74, HLA-DPB1, HLA-DMA) and co-stimulatory molecules (CD80, CD86).
CIITA+ confirms MHC-II transcriptional regulation.
Can invoke both (neo)antigen-specific and antigen-independent T cell responses.
Associated with favorable survival in most cancers.
Induced by leucine metabolism → histone H3K27ac modification.
Negative FCGR3B/CXCR2 (mature neutrophil markers) may be reduced; negative VEGFA/MMP9 distinguishes from angiogenic neutrophils.
"""
```

#### 10 Neutrophil States (from Cell 2024)

| State | Markers | Function | Prognosis |
|-------|---------|----------|-----------|
| Inflammation | S100A8, S100A9, CXCL8, IL1B | Inflammatory | Context-dependent |
| Angiogenesis | VEGFA, MMP9, CXCL1, CXCL2 | Pro-tumor | Poor |
| Antigen presentation | HLA-DRA, CD74, CD80, CD86 | Anti-tumor | **Favorable** |
| Interferon response | ISG15, MX1, OAS1, IFIT1 | Anti-viral | Context-dependent |
| Leukotriene | ALOX5, LTC4S, GPX1 | Lipid mediator | Unknown |
| Phagocytosis | FCGR3B, CYBB, NCF1, NCF2 | Anti-pathogen | Context-dependent |
| ROS production | CYBB, NCF1, NCF2, MPO | Oxidative burst | Context-dependent |
| Degranulation | MPO, ELANE, PRTN3, AZU1 | Cytotoxic | Context-dependent |
| Maturation | LTF, LCN2, CD177, CEACAM8 | Terminal maturation | Unknown |
| Stress | HSPA1A, HSPB1, DNAJB1 | Cellular stress | Poor |

```toml
# State entries for the 10 neutrophil states:

[[state.minor]]
name = "Neutrophil angiogenic"
color = "#FF6347"
markers = ["VEGFA", "MMP9", "CXCL1", "CXCL2", "CXCL8", "FGF2"]
negative_markers = ["HLA-DRA", "CD74", "CD80", "CD86"]

[state.minor.metadata]
kind = "state"
category = "angiogenesis"
scope = "lineage_restricted"
applies_to = ["Neutrophils"]
alias_of = "Angiogenic_neutrophil"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true

[[state.minor]]
name = "Neutrophil antigen-presenting"
color = "#228B22"
markers = ["HLA-DRA", "HLA-DQA1", "CD74", "CD80", "CD86", "CIITA"]
negative_markers = ["VEGFA", "MMP9", "CXCL8", "S100A8"]

[state.minor.metadata]
kind = "state"
category = "immune_activation"
scope = "lineage_restricted"
applies_to = ["Neutrophils"]
alias_of = "APN"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true

[[state.minor]]
name = "Neutrophil inflammatory"
color = "#B22222"
markers = ["S100A8", "S100A9", "CXCL8", "IL1B", "TNF", "PTGS2"]
negative_markers = ["HLA-DRA", "VEGFA", "HSPA1A"]

[state.minor.metadata]
kind = "state"
category = "inflammation"
scope = "lineage_restricted"
applies_to = ["Neutrophils"]
alias_of = "Inflammatory_neutrophil"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true

[[state.minor]]
name = "Neutrophil interferon-response"
color = "#4169E1"
markers = ["ISG15", "MX1", "OAS1", "IFIT1", "IFIT3", "STAT1"]
negative_markers = ["VEGFA", "S100A8"]

[state.minor.metadata]
kind = "state"
category = "immune_activation"
scope = "lineage_restricted"
applies_to = ["Neutrophils"]
alias_of = "IFN_neutrophil"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

#### Tumor Context (`marker_tumor`)

**Tumor-Associated Antigen-Presenting Neutrophil**

```toml
[[cancer_state]]
name = "Tumor-associated antigen-presenting neutrophil"
color = "#228B22"
markers = ["HLA-DRA", "HLA-DQA1", "CD74", "CD80", "CD86", "CIITA", "CD40"]
negative_markers = ["VEGFA", "MMP9", "ARG1", "CD274"]

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
source = { title = "Neutrophil profiling illuminates anti-tumor antigen-presenting potency", year = "2024", doi = "10.1016/j.cell.2024.02.005" }
notes = "Antigen-presenting neutrophils in tumors. Associated with favorable survival in most cancers. Induced by leucine metabolism and histone H3K27ac modification. Can enhance anti-PD-1 therapy. Negative pro-tumor markers (VEGFA, MMP9, ARG1, CD274)."
```

#### Functional Programs (geneset)

```json
"Neutrophil_antigen_presentation": {
  "genes": ["HLA-DRA", "HLA-DQA1", "HLA-DPB1", "HLA-DMA", "CD74", "CD80", "CD86", "CIITA", "CD40", "HLA-DMB", "HLA-DOA", "HLA-DOB"],
  "description": "Neutrophil antigen presentation program (APN). MHC-II and co-stimulatory molecule expression. Associated with favorable prognosis and enhanced anti-PD-1 therapy."
},
"Neutrophil_angiogenic": {
  "genes": ["VEGFA", "MMP9", "CXCL1", "CXCL2", "CXCL8", "FGF2", "MMP2", "ANGPT1", "ANGPT2"],
  "description": "Neutrophil angiogenic program. Pro-tumor, associated with poor prognosis."
},
"Neutrophil_inflammatory": {
  "genes": ["S100A8", "S100A9", "CXCL8", "IL1B", "TNF", "PTGS2", "CXCL1", "CXCL2", "CCL3"],
  "description = "Neutrophil inflammatory program. Pro-inflammatory cytokine and chemokine production."
}
```

---

### Article 7: Gut 2022 — BHLHE40-Driven Pro-Tumour Neutrophils in PDAC

**Dataset**: 5 PDAC patients
**Key finding**: BHLHE40+ TAN-1 (terminally differentiated, hyperactivated glycolysis, poor prognosis)

#### New Subtype Entry (`marker_registry`)

**BHLHE40+ Pro-Tumor Neutrophil**

```toml
[[subtype]]
name = "BHLHE40+ pro-tumor neutrophil"
color = "#8B0000"
markers = ["BHLHE40", "LDHA", "PKM", "ENO1", "SLC2A1", "HK2", "PFKFB3", "MCT4"]
negative_markers = ["HLA-DRA", "CD74", "CD80", "ISG15", "MX1"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tumor_context"
applies_to = ["Neutrophils", "Pancreatic Cancer"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Single-cell RNA-seq analysis reveals BHLHE40-driven pro-tumour neutrophils with hyperactivated glycolysis in pancreatic tumour microenvironment", year = "2022", doi = "10.1136/gutjnl-2021-326070" }
notes = """
BHLHE40+ pro-tumor neutrophils (TAN-1) in PDAC. Hyperactivated glycolysis (LDHA+, PKM+, ENO1+, SLC2A1+, HK2+, PFKFB3+, MCT4+).
BHLHE40 is downstream of hypoxia and ER stress. Direct transcriptional regulator of TAN-1 marker genes (validated by ChIP).
Immunosuppressive and pro-tumor functions. Unfavorable prognostic value.
Negative antigen-presenting markers (HLA-DRA, CD74, CD80) and IFN-response markers (ISG15, MX1).
"""
```

#### PDAC TAN Subtypes

| Subtype | Markers | Function | Prognosis |
|---------|---------|----------|-----------|
| TAN-1 | BHLHE40, LDHA, PKM, ENO1 | Pro-tumor, glycolytic | **Poor** |
| TAN-2 | S100A8, S100A9, CXCL8 | Inflammatory | Context-dependent |
| TAN-3 | VCAN, FCGR3B, CXCR2 | Transitional | Unknown |
| TAN-4 | ISG15, MX1, OAS1 | IFN-response | Context-dependent |

#### Tumor Context (`marker_tumor`)

**BHLHE40+ TAN (PDAC)**

```toml
[[cancer_state]]
name = "BHLHE40+ TAN"
color = "#8B0000"
markers = ["BHLHE40", "LDHA", "PKM", "ENO1", "SLC2A1", "HK2", "PFKFB3"]
negative_markers = ["HLA-DRA", "CD74", "ISG15", "MX1", "CD80"]

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
source = { title = "Single-cell RNA-seq analysis reveals BHLHE40-driven pro-tumour neutrophils with hyperactivated glycolysis in pancreatic tumour microenvironment", year = "2022", doi = "10.1136/gutjnl-2021-326070" }
notes = "BHLHE40+ tumor-associated neutrophils in PDAC. Hyperactivated glycolysis. Hypoxia and ER stress-driven. Pro-tumor and immunosuppressive. Poor prognosis. Negative antigen-presenting and IFN-response markers."
```

#### Functional Programs (geneset)

```json
"Neutrophil_glycolytic": {
  "genes": ["LDHA", "PKM", "ENO1", "SLC2A1", "HK2", "PFKFB3", "MCT4", "PGK1", "GAPDH", "BHLHE40"],
  "description": "Neutrophil hyperactivated glycolysis program (BHLHE40-driven). Pro-tumor, poor prognosis in PDAC."
},
"Neutrophil_hypoxia_response": {
  "genes": ["BHLHE40", "HIF1A", "EPAS1", "SLC2A1", "VEGFA", "LDHA", "PGK1", "ENO1", "CA9", "BNIP3"],
  "description": "Neutrophil hypoxia response program. BHLHE40/HIF-driven. Associated with pro-tumor phenotype."
}
```

---

### Article 6: Cancer Cell 2022 — Tissue-Resident Neutrophils in NSCLC

**Dataset**: NSCLC, 1,283,972 cells, 556 samples, 318 patients
**Key finding**: Tissue-resident neutrophils (TRNs) with distinct functional properties
**Clinical relevance**: TRN-derived gene signature associated with anti-PD-L1 treatment failure

#### New Subtype Entry (`marker_registry`)

**Tissue-Resident Neutrophil (TRN)**

```toml
[[subtype]]
name = "Tissue-resident neutrophil"
color = "#4682B4"
markers = ["CD69", "ITGA1", "CD44", "CXCR4", "SELPLG", "CD93", "FCGR3B"]
negative_markers = ["CXCR2", "CSF3R", "S100A8", "S100A9", "CXCL8"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tissue_specific"
applies_to = ["Lung", "Tumor"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer", year = "2022", doi = "10.1016/j.ccell.2022.10.008" }
notes = """
Tissue-resident neutrophils (TRNs) in NSCLC. Distinct from circulating neutrophils.
Express tissue-residency markers (CD69, ITGA1/CD49a, CD44, CXCR4, SELPLG, CD93).
Acquire new functional properties in tissue microenvironment (plasticity).
TRN-derived gene signature associated with anti-PD-L1 treatment failure.
Negative CXCR2/CSF3R (circulating neutrophil markers) and inflammatory markers (S100A8/9, CXCL8).
"""
```

#### Tumor Context (`marker_tumor`)

**TRN Signature (ICB Resistance)**

```toml
[[cancer_state]]
name = "TRN ICB-resistant signature"
color = "#FF4500"
markers = ["CD69", "ITGA1", "CD44", "CXCR4", "SELPLG", "CD93", "VEGFA"]
negative_markers = ["CXCR2", "HLA-DRA", "CD80", "CD86"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["Lung Cancer"]
scope = "cancer_type_specific"
applies_to = ["Lung"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer", year = "2022", doi = "10.1016/j.ccell.2022.10.008" }
notes = "Tissue-resident neutrophil signature associated with anti-PD-L1 treatment failure in NSCLC. CD69+ ITGA1+ CD44+ tissue-residency markers. May promote immune exclusion or exhaustion. Negative circulating markers (CXCR2) and antigen-presenting markers (HLA-DRA, CD80, CD86)."
```

---

### Article 5: Nature 2022 — Liver Tumour TIME Subtypes and Neutrophil Heterogeneity

**Dataset**: PLC (primary liver cancer)
**Key finding**: 5 TIME subtypes; neutrophil heterogeneity within liver TME
**Notable**: From Zhang Zemin's group (BIOPIC)

#### Tumor Context (`marker_tumor`)

**Liver Tumor Neutrophil (TIME-ISM)**

The paper identifies neutrophils within the TIME-ISM (Immune Suppressive Myeloid) subtype:

```toml
[[cancer_state]]
name = "Liver tumor neutrophil (TIME-ISM)"
color = "#8B0000"
markers = ["S100A8", "S100A9", "CXCL8", "MMP9", "VEGFA", "ARG1", "CD274"]
negative_markers = ["HLA-DRA", "CD74", "CD80", "CXCL9", "CXCL10"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["Liver Cancer"]
scope = "cancer_type_specific"
applies_to = ["Liver"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Liver tumour immune microenvironment subtypes and neutrophil heterogeneity", year = "2022", doi = "10.1038/s41586-022-05400-x" }
notes = "Neutrophils within TIME-ISM (immune suppressive myeloid) subtype of liver cancer. Express S100A8/9, CXCL8, MMP9, VEGFA, ARG1, CD274. Immunosuppressive and pro-tumor. Negative antigen-presenting markers (HLA-DRA, CD74, CD80) and T cell-recruiting chemokines (CXCL9, CXCL10)."
```

---

### Articles 8, 9: Reviews (Nat Rev Immunol 2022, Nat Rev Cancer 2020)

**Key frameworks for neutrophil classification**:

| Framework Element | Description | Registry Impact |
|-------------------|-------------|-----------------|
| Maturation stages | Myeloblast → promyelocyte → myelocyte → metamyelocyte → band → mature | Add developmental states |
| Tissue localization | Circulating vs tissue-resident vs tumor-infiltrating | Scope classification |
| Functional adaptation | Inflammatory, angiogenic, antigen-presenting, suppressive | State entries |

#### Neutrophil Maturation States

```toml
[[state.minor]]
name = "Neutrophil immature (band)"
color = "#DAA520"
markers = ["CD33", "CD11b", "CD15", "CD66b", "CD101", "CD10", "CD16low"]
negative_markers = ["CD10", "CD16", "CXCR2"]

[state.minor.metadata]
kind = "state"
category = "maturation"
scope = "lineage_restricted"
applies_to = ["Neutrophils"]
alias_of = "Immature_neutrophil"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
notes = "Immature band neutrophil. CD33+ CD11b+ CD15+ CD66b+. Low CD10 and CD16. CD101mid. May be expanded in cancer (left shift)."

[[state.minor]]
name = "Neutrophil mature"
color = "#228B22"
markers = ["CD10", "CD16", "CD66b", "CD15", "CXCR2", "FCGR3B", "CSF3R"]
negative_markers = ["CD33", "CD34", "CD117"]

[state.minor.metadata]
kind = "state"
category = "maturation"
scope = "lineage_restricted"
applies_to = ["Neutrophils"]
alias_of = "Mature_neutrophil"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
notes = "Mature segmented neutrophil. CD10+ CD16+ CD66b+ CXCR2+ FCGR3B+ CSF3R+. Negative progenitor markers (CD33, CD34, CD117)."
```

---

## Updated Neutrophil Registry Hierarchy

### Proposed New Structure

```
Neutrophils (compartment) [EXPANDED]
├── Mature neutrophil [EXISTING - expanded]
│   ├── Circulating neutrophil [state]
│   ├── Tissue-resident neutrophil (TRN) [NEW - subtype]
│   └── Tumor-associated neutrophil (TAN) [state]
├── Immature neutrophil (band) [NEW - state]
├── Antigen-presenting neutrophil (APN) [NEW - subtype]
├── BHLHE40+ pro-tumor neutrophil [NEW - subtype]
├── Neutrophil states:
│   ├── Inflammatory [state]
│   ├── Angiogenic [state]
│   ├── Antigen-presenting [state]
│   ├── Interferon-response [state]
│   ├── Leukotriene [state]
│   ├── Phagocytic [state]
│   ├── ROS production [state]
│   ├── Degranulation [state]
│   ├── Maturation [state]
│   └── Stress [state]
└── Tumor contexts:
    ├── TRN ICB-resistant [cancer_state]
    ├── Liver tumor neutrophil (TIME-ISM) [cancer_state]
    └── BHLHE40+ TAN [cancer_state]
```

---

## Curation Notes and Conflicts

### 1. Neutrophil Short Half-Life

**Issue**: Neutrophils are short-lived (6-48 hours), challenging single-cell analysis.

**Resolution**:
- Add note: "Neutrophil transcriptomes may reflect acute activation states rather than stable identities"
- Distinguish between mature circulating (CXCR2+, FCGR3B+) and tissue-resident (CD69+, ITGA1+) populations
- TAN states are dynamic and context-dependent

### 2. Antigen Presentation by Neutrophils

**Issue**: Neutrophils expressing HLA-DRA, CD80, CD86 is counterintuitive.

**Resolution**:
- Validate with multiple markers (HLA-DRA + CD74 + CIITA + CD80/CD86)
- Single HLA-DRA is insufficient (can be non-specific)
- Add note: "APNs are rare and require co-expression of multiple MHC-II and co-stimulatory markers"

### 3. BHLHE40+ Neutrophils vs TAMs

**Issue**: BHLHE40+ neutrophils (Gut 2022) express glycolytic genes that overlap with metabolic programs in other cells.

**Resolution**:
- Use neutrophil identity markers (FCGR3B, CXCR2) to confirm lineage
- BHLHE40 is transcription factor, not lineage marker
- Add negative markers (HLA-DRA, ISG15) to distinguish from other neutrophil states

### 4. TRN vs TAM

**Issue**: Tissue-resident neutrophils (TRNs) and tumor-associated neutrophils (TANs) may overlap.

**Resolution**:
- TRN: tissue-residency markers (CD69, ITGA1), found in both normal and tumor tissue
- TAN: tumor-specific functional states, may include TRNs or circulating neutrophils
- TRN is **location/identity**; TAN is **context/state**

### 5. Immunity 2025 Naming Framework

**Issue**: New standardized naming framework proposed.

**Resolution**:
- Adopt framework: integrate maturation + tissue localization + functional adaptation
- Example: "Mature tumor-associated angiogenic neutrophil" instead of "TAN"
- Gradually transition naming; keep existing aliases for backward compatibility

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New subtype entries | 3 (APN, TRN, BHLHE40+ pro-tumor) | marker_registry |
| New state entries | 4 (angiogenic, antigen-presenting, inflammatory, IFN-response) | marker_registry |
| New maturation states | 2 (immature band, mature) | marker_registry |
| Existing entry expansion | 1 (Neutrophils: FCGR3B, CXCR2 → expanded) | marker_registry |
| New tumor context entries | 3 (TRN ICB-resistant, liver TIME-ISM, BHLHE40+ TAN) | marker_tumor |
| New geneset programs | 6 | genesets_cancer_signatures.json |
| Naming framework updates | 1 (Immunity 2025) | marker_registry |
| Conflict resolutions | 5 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **Neutrophils are not homogeneous**: The existing registry entry (FCGR3B, CXCR2, CSF3R) is grossly insufficient. Neutrophils exhibit 10+ distinct states (Cell 2024).

2. **APNs are anti-tumor**: Antigen-presenting neutrophils (HLA-DRA+, CD80+) are associated with favorable survival and enhanced anti-PD-1 therapy (Cell 2024).

3. **BHLHE40 drives pro-tumor phenotype**: BHLHE40+ neutrophils in PDAC have hyperactivated glycolysis and poor prognosis (Gut 2022).

4. **TRNs and ICB resistance**: Tissue-resident neutrophil signature is associated with anti-PD-L1 failure in NSCLC (Cancer Cell 2022).

5. **Neutrophil plasticity**: Neutrophils acquire new functional properties in tissue microenvironments (Cancer Cell 2022, Cell 2024). State annotation is critical.

6. **Standardized naming**: Immunity 2025 proposes integrated naming (maturation + location + function). scLucid should adopt this framework.
