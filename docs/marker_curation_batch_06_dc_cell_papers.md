# Marker Curation Batch 06: Dendritic Cell Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | Pan-Cancer Analyses Refine the Single-Cell Portrait of Tumor-Infiltrating Dendritic Cells | 2025 | Cancer Research | 10.1158/0008-5472.CAN-24-3595 | pan_cancer_atlas | marker_registry + marker_tumor |
| 2 | Dendritic cells as orchestrators of anticancer immunity and immunotherapy | 2024 | Nature Reviews Clinical Oncology | 10.1038/s41571-024-00859-1 | review | marker_registry (validation + naming) |
| 3 | Human dendritic cells in cancer | 2022 | Science Immunology | 10.1126/sciimmunol.abm9409 | review | marker_registry (validation + naming) |
| 4 | Tumor-infiltrating dendritic cell states are conserved across solid human cancers | 2021 | Journal of Experimental Medicine | 10.1084/jem.20200264 | meta_analysis | marker_registry + marker_tumor |
| 5 | Dendritic cells in cancer immunology and immunotherapy | 2020 | Nature Reviews Immunology | 10.1038/s41577-019-0210-z | review | marker_registry (validation) |
| 6 | Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors | 2017 | Science | 10.1126/science.aah4573 | single_cell_atlas | marker_registry (foundational) |

---

## Overview: DC Biology Context

### Existing Registry Entries (Current Status)

| Entry | Current Markers | Evidence Tier | Notes |
|-------|----------------|---------------|-------|
| Dendritic cells (compartment) | CD74, FLT3, ITGAX, CD1C | curated_review | Broad identity |
| cDC1 | CLEC9A, XCR1, CADM1 | seed | Cross-presentation |
| cDC2 | CD1C, FCER1A, CLEC10A | seed | CD4+ T helper priming |
| cDC3 | LAMP3, CCR7, FSCN1, CST7 | seed | Migratory mature DC |
| pDC | LILRA4, GZMB, IL3RA, TCF4 | seed | Type I IFN production |

### The Villani Taxonomy (Science 2017) — Foundational Reference

The Science 2017 paper (Villani et al.) established the modern human blood DC taxonomy. It identified **6 DC types** in human blood:

| Villani DC | Markers | Relationship to scLucid Registry | Status |
|-----------|---------|--------------------------------|--------|
| DC1 | CLEC9A, XCR1, CADM1 | = cDC1 | **Consistent** |
| DC2 | CD1C, FCER1A, CLEC10A | = cDC2 | **Consistent** |
| DC3 | AXL, SIGLEC6, CD1C | ≈ cDC2 variant | **New** |
| DC4 | CD1C, S100A8, S100A9, CD14 | CD1C+ monocyte-related | **New** |
| DC5 | LILRA4, GZMB, IL3RA | pDC subtype | ≈ pDC |
| DC6 | LILRA4, TCF4, IL3RA | pDC | **Consistent** |

**Key insight**: DC3 (AXL+SIGLEC6+) is a distinct subset sharing properties with both pDCs and cDCs. DC4 represents a CD1C+ subset related to monocytes.

### The mregDC State (Science Immunology 2022)

"Mature DCs enriched in immunoregulatory molecules" (mregDC) is a **cross-lineage state** (not a subtype) observed across human tumors:

| Feature | mregDC |
|---------|--------|
| Markers | LAMP3, CCR7, CD40, CD274, PDCD1LG2, IDO1, FSCN1 |
| Lineage origin | cDC1, cDC2, or DC3 |
| State | Maturation + immunoregulatory |
| Function | T cell regulation, tolerance induction |
| Location | Tumor, tumor-draining lymph nodes |

**Critical distinction**: mregDC is a **maturation state**, not a separate lineage. The existing `cDC3` entry (LAMP3, CCR7, FSCN1, CST7) likely captures the mregDC population. The question is whether `cDC3` should be reclassified as a **state** rather than a **subtype**.

---

## Paper-by-Paper Curation

### Article 6: Science 2017 — Foundational DC Taxonomy (Villani et al.)

**Role**: Establishes the modern human blood DC taxonomy. All subsequent DC papers build on this foundation.

**Dataset**: ~2,400 blood cells
**Key finding**: 6 DC types, 4 monocyte subtypes, circulating cDC progenitors

#### New Subtype Entries (`marker_registry`)

**1. AXL+SIGLEC6+ DC (DC3)**

```toml
[[compartment.minor.minor.minor]]
name = "AXL+SIGLEC6+ DC"
color = "#DAA520"
markers = ["AXL", "SIGLEC6", "CD1C", "CD2", "IL22RA2", "CD22"]
negative_markers = ["CLEC9A", "XCR1", "CADM1", "FCER1A", "CLEC10A", "LILRA4", "TCF4"]

[compartment.minor.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "all"
applies_to = ["Dendritic cells"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors", year = "2017", doi = "10.1126/science.aah4573" }
notes = """
DC3 from Villani taxonomy. AXL+SIGLEC6+ subset sharing properties with both pDCs and cDCs.
Potently activates T cells (unlike typical pDCs). Also known as AS DC (antigen-stimulating DC).
Negative CLEC9A/XCR1 (cDC1), FCER1A/CLEC10A (cDC2), LILRA4/TCF4 (pDC) confirms distinct identity.
Also identified as rare subset in tumors (Cancer Research 2025).
"""
```

**2. CD1C+ CD14+ DC (DC4 / Monocyte-related DC)**

```toml
[[compartment.minor.minor.minor]]
name = "CD1C+ monocyte-related DC"
color = "#CD853F"
markers = ["CD1C", "CD14", "S100A8", "S100A9", "FCN1", "ITGAX"]
negative_markers = ["CLEC9A", "XCR1", "FCER1A", "CLEC10A", "LILRA4", "AXL", "SIGLEC6"]

[compartment.minor.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "all"
applies_to = ["Dendritic cells", "Monocytes"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors", year = "2017", doi = "10.1126/science.aah4573" }
notes = "DC4 from Villani taxonomy. CD1C+ subset with monocyte-related features (CD14+, S100A8/9+, FCN1+). May represent inflammatory DCs or DCs transitioning from monocytes. Negative cDC1/cDC2/AXL DC/pDC markers."
```

**3. cDC Progenitor**

```toml
[[lineage.minor]]
name = "cDC Progenitor"
color = "#D2691E"
markers = ["FLT3", "CSF2RA", "KLF4", "GATA2", "CD34", "CD123", "IL3RA"]
negative_markers = ["CD14", "CD16", "LYZ", "VCAN"]

[lineage.minor.metadata]
kind = "cell_type"
granularity = "lineage"
compartment = "immune"
lineage = "myeloid"
scope = "all"
applies_to = ["Dendritic cells", "Monocytes"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors", year = "2017", doi = "10.1126/science.aah4573" }
notes = "Circulating common DC progenitor (CDP). FLT3+ CD34+ CD123+ (IL3RA+). Differentiates into cDC1, cDC2, and pDC. Negative monocyte markers (CD14, CD16, LYZ, VCAN). Rare in peripheral blood."
```

#### Evidence Tier Promotions

| Entry | Current Tier | Proposed Tier | Reason |
|-------|-------------|--------------|--------|
| cDC1 (CLEC9A, XCR1, CADM1) | seed | **atlas_supported** | Science 2017 foundational |
| cDC2 (CD1C, FCER1A, CLEC10A) | seed | **atlas_supported** | Science 2017 foundational |
| pDC (LILRA4, GZMB, IL3RA, TCF4) | seed | **atlas_supported** | Science 2017 foundational |
| cDC3 (LAMP3, CCR7, FSCN1, CST7) | seed | **needs_review** | Re-evaluate as state vs subtype |

#### cDC3 Re-evaluation

The existing `cDC3` entry uses LAMP3, CCR7, FSCN1, CST7. However, Science Immunology 2022 and Cancer Research 2025 indicate that LAMP3+ DCs are a **maturation state** (mregDC) rather than a separate lineage.

**Recommended action**:
1. Keep `cDC3` entry but reclassify as **state** (not subtype)
2. OR keep as subtype but add notes: "Represents mature migratory DC state, not distinct lineage"
3. Add `mregDC` as explicit state entry with broader marker set

**Resolution**: Keep `cDC3` as subtype (it is commonly used in the field) but add `mregDC` as a **state** with more comprehensive markers.

---

### Article 3: Science Immunology 2022 — Human DCs in Cancer

**Role**: Review of human DC biology with focus on cancer. Introduces mregDC concept.

**Key concepts**:
- cDC1, cDC2, DC3, pDC lineages
- mregDC: cross-lineage mature state
- DC functions in antitumor immunity and tolerance

#### New State Entry (`marker_registry`)

**mregDC (Mature DC Enriched in Immunoregulatory Molecules)**

```toml
[[state.minor]]
name = "mregDC"
color = "#8A2BE2"
markers = ["LAMP3", "CCR7", "CD40", "CD274", "PDCD1LG2", "IDO1", "FSCN1", "CST7", "EBI3", "CD83"]
negative_markers = ["CLEC9A", "XCR1", "FCER1A", "CLEC10A", "LILRA4"]

[state.minor.metadata]
kind = "state"
category = "immune_activation"
scope = "lineage_restricted"
applies_to = ["Dendritic cells", "cDC1", "cDC2", "DC3"]
alias_of = "mregDC"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
source = [
    { title = "Human dendritic cells in cancer", year = "2022", doi = "10.1126/sciimmunol.abm9409" },
    { title = "Tumor-infiltrating dendritic cell states are conserved across solid human cancers", year = "2021", doi = "10.1084/jem.20200264" }
]
notes = """
mregDC: Mature DCs enriched in immunoregulatory molecules. Cross-lineage state (can arise from cDC1, cDC2, or DC3).
Markers: LAMP3, CCR7, CD40, CD274 (PD-L1), PDCD1LG2 (PD-L2), IDO1, FSCN1, CST7, EBI3, CD83.
Expresses co-stimulatory (CD40), co-inhibitory (PD-L1/PD-L2), and tolerogenic (IDO1) molecules simultaneously.
Found in tumors and tumor-draining lymph nodes.
Negative lineage markers (CLEC9A, XCR1, FCER1A, CLEC10A, LILRA4) confirm loss of lineage identity upon maturation.
"""
```

#### Existing Entry Updates

| Entry | Update | Source |
|-------|--------|--------|
| cDC3 (LAMP3, CCR7, FSCN1, CST7) | Add note: "Overlaps with mregDC state. May represent mature state rather than distinct lineage." | Science Immunology 2022 |
| pDC | Add note: "DC3 (AXL+SIGLEC6+) shares pDC properties but potently activates T cells, unlike typical pDCs." | Science 2017 |

---

### Article 4: JEM 2021 — Conserved DC States Across Cancers

**Role**: Meta-analysis confirming conserved DC states across solid tumors.

**Key finding**: Tumor-infiltrating DC states are **conserved** across:
- Patients
- Cancer types
- Species (human and mouse)

#### Validation of Existing Entries

This meta-analysis validates the stability of existing DC subtype markers across cancer contexts. **Evidence tier promotion recommended**:

| Entry | Current Tier | Proposed Tier | Reason |
|-------|-------------|--------------|--------|
| cDC1 | atlas_supported (from Science 2017) | **consensus** | Conserved across cancers and species |
| cDC2 | atlas_supported | **consensus** | Conserved across cancers and species |
| cDC3 | needs_review | **atlas_supported** | Conserved across cancers (as state) |
| pDC | atlas_supported | **consensus** | Conserved across cancers and species |

---

### Article 1: Cancer Research 2025 — Pan-Cancer DC Atlas

**Dataset**: 2,500+ samples, 33 cancer types
**Key findings**:
- Rare subsets: AXL+SIGLEC6+ DCs, Langerhans cell-like DCs
- Langerhans cell-like → additional origin of tumor-enriched LAMP3+ DCs
- Distinct cellular origins → pleiotropic functional potentials of LAMP3+ DCs
- Machine learning model for DC annotation

#### New Subtype Entries (`marker_registry`)

**1. Langerhans Cell-like DC**

```toml
[[compartment.minor.minor.minor]]
name = "Langerhans cell-like DC"
color = "#DAA520"
markers = ["CD207", "CD1A", "S100B", "IRF4", "S100A9", "FCGR2B"]
negative_markers = ["CLEC9A", "XCR1", "FCER1A", "CLEC10A", "LILRA4", "AXL", "SIGLEC6"]

[compartment.minor.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tissue_specific"
applies_to = ["Skin", "Mucosa", "Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Pan-Cancer Analyses Refine the Single-Cell Portrait of Tumor-Infiltrating Dendritic Cells", year = "2025", doi = "10.1158/0008-5472.CAN-24-3595" }
notes = """
Langerhans cell-like DC identified in pan-cancer atlas. Normally epidermal/mucosal DCs.
Tumor-infiltrating Langerhans-like DCs may represent recruited or tissue-resident population.
Markers: CD207 (Langerin), CD1A, S100B, IRF4, S100A9, FCGR2B.
Can differentiate into LAMP3+ DCs in tumors (additional cellular origin).
Negative cDC1/cDC2/AXL DC/pDC markers.
"""
```

#### Tumor Context (`marker_tumor`)

**Tumor-enriched LAMP3+ DC**

```toml
[[cancer_state]]
name = "Tumor-enriched LAMP3+ DC"
color = "#9932CC"
markers = ["LAMP3", "CCR7", "FSCN1", "CST7", "CD40", "CD274", "IDO1", "EBI3"]
negative_markers = ["CLEC9A", "XCR1", "FCER1A", "CLEC10A", "LILRA4", "CD207", "CD1A"]

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
source = [
    { title = "Pan-Cancer Analyses Refine the Single-Cell Portrait of Tumor-Infiltrating Dendritic Cells", year = "2025", doi = "10.1158/0008-5472.CAN-24-3595" },
    { title = "Human dendritic cells in cancer", year = "2022", doi = "10.1126/sciimmunol.abm9409" }
]
notes = """
Tumor-enriched LAMP3+ DC (mregDC state). Multiple cellular origins: cDC1, cDC2, DC3, and Langerhans cell-like DC.
Expresses immunoregulatory molecules: CD40 (co-stimulatory), CD274/PD-L1 (co-inhibitory), IDO1 (tolerogenic).
Pleiotropic functional potential depending on cellular origin and TME context.
Associated with both antitumor immunity (if CD40-dominated) and immune tolerance (if PD-L1/IDO1-dominated).
Negative lineage markers confirm loss of identity upon maturation.
"""
```

#### Functional Programs (geneset)

```json
"DC_pleiotropic_potential": {
  "genes": ["LAMP3", "CCR7", "CD40", "CD274", "PDCD1LG2", "IDO1", "EBI3", "CD83", "FSCN1", "CST7", "HLA-DRA", "HLA-DQA1"],
  "description": "DC pleiotropic functional potential signature (mregDC/tumor-enriched LAMP3+ DC). Includes co-stimulatory, co-inhibitory, and tolerogenic markers."
},
"DC_cross_presentation": {
  "genes": ["CLEC9A", "XCR1", "CADM1", "WDFY4", "IRF8", "BATF3", "RAB7A", "SEC22B", "HLA-A", "HLA-B", "HLA-C", "B2M", "PSMB9", "PSMB10", "TAP1", "TAP2"],
  "description": "cDC1 cross-presentation program. MHC-I antigen presentation to CD8+ T cells."
},
"DC_Th_priming": {
  "genes": ["CD1C", "FCER1A", "CLEC10A", "HLA-DRA", "HLA-DQA1", "CD80", "CD86", "CD40", "CCL17", "CCL22", "CXCL16", "IL15"],
  "description": "cDC2 Th cell priming program. MHC-II antigen presentation to CD4+ T cells."
},
"DC_type_I_IFN": {
  "genes": ["LILRA4", "GZMB", "IL3RA", "TCF4", "IRF7", "IRF8", "TLR7", "TLR9", "IFNA1", "IFNA2", "IFNB1", "MX1", "OAS1"],
  "description": "pDC type I interferon production program. Antiviral and antitumor immunity."
}
```

---

### Articles 2 & 5: Reviews (Nat Rev Clin Oncol 2024, Nat Rev Immunol 2020)

**Role**: Provide DC biology context and validate therapeutic relevance.

**Key insights for registry**:

| Concept | Implication for scLucid |
|---------|------------------------|
| DC states determine immunotherapy efficacy | scLucid tumor interpretation should score DC functional states |
| DC can promote immunity OR tolerance | State markers (mregDC) are critical for interpretation |
| DC heterogeneity requires single-cell resolution | Manager views should separate DC subtypes and states |
| DC-based vaccines under development | Clinical relevance of accurate DC annotation |

**Registry impact**: No new markers. Validate therapeutic relevance of existing entries.

---

## Updated DC Registry Hierarchy

### Proposed Structure

```
Dendritic cells (compartment)
├── cDC1 (classical DC1) [EXISTING]
│   └── Cross-presentation state [geneset]
├── cDC2 (classical DC2) [EXISTING]
│   └── Th priming state [geneset]
├── AXL+SIGLEC6+ DC (DC3) [NEW - subtype]
│   └── DC3 activation state [state]
├── CD1C+ monocyte-related DC (DC4) [NEW - subtype]
├── pDC (plasmacytoid DC) [EXISTING]
│   └── Type I IFN state [geneset]
├── Langerhans cell-like DC [NEW - subtype]
│   └── Tumor-enriched LAMP3+ DC [cancer_state]
├── mregDC [NEW - state]
│   └── Tumor-enriched LAMP3+ DC [cancer_state]
└── cDC Progenitor [NEW - lineage]
```

---

## Curation Notes and Conflicts

### 1. cDC3 vs mregDC

**Issue**: The existing `cDC3` entry (LAMP3, CCR7, FSCN1, CST7) overlaps with the `mregDC` state (LAMP3, CCR7, FSCN1, CST7, CD40, CD274, IDO1).

**Resolution**:
- **Keep `cDC3` as subtype** (commonly used in field, 2017 taxonomy)
- **Add `mregDC` as state** with more comprehensive markers (including CD40, CD274, IDO1)
- Add cross-reference notes:
  - cDC3: "May represent early mregDC state. Distinct lineage identity in blood but converges to mregDC state in tumors."
  - mregDC: "Cross-lineage maturation state. Can arise from cDC1, cDC2, DC3, or Langerhans cell-like DC."

### 2. DC3 (AXL+SIGLEC6+) vs pDC

**Issue**: DC3 shares some pDC properties but functionally distinct.

**Resolution**:
- DC3 is a **separate subtype** (not a pDC state)
- Markers: AXL, SIGLEC6, CD1C, CD2, IL22RA2
- Function: Potently activates T cells (unlike typical pDCs)
- Negative LILRA4/TCF4 distinguishes from pDC
- Negative CLEC9A/XCR1/FCER1A/CLEC10A distinguishes from cDC1/cDC2

### 3. Langerhans Cell-like DC in Tumors

**Issue**: Langerhans cells (LC) are typically epidermal. Their presence in tumors is unexpected.

**Resolution**:
- Add as subtype with `scope = "tissue_specific"`
- Applies_to: ["Skin", "Mucosa", "Tumor"]
- Note: "May represent tissue-resident LCs recruited to tumors or Langerhans-like differentiation of other DCs"
- Can differentiate into LAMP3+ DCs (Cancer Research 2025)

### 4. cDC Progenitor Rarity

**Issue**: cDC progenitors are rare in peripheral blood and may not be detectable in all datasets.

**Resolution**:
- Add as `lineage` entry (not subtype)
- `use_for_global_annotation = true` but note rarity
- Add note: "Rare circulating population. May not be detectable in all scRNA-seq datasets."

### 5. Evidence Tier Consistency

**Issue**: Multiple papers validate the same DC markers. Need consistent evidence tier assignment.

**Resolution**:
| Entry | Evidence Tier | Sources |
|-------|--------------|---------|
| cDC1 | **consensus** | Science 2017 (foundational) + JEM 2021 (conserved) + Cancer Research 2025 (pan-cancer) |
| cDC2 | **consensus** | Science 2017 + JEM 2021 + Cancer Research 2025 |
| pDC | **consensus** | Science 2017 + JEM 2021 + Cancer Research 2025 |
| cDC3 | **atlas_supported** | Science 2017 + Cancer Research 2025 (but reclassified as state) |
| AXL+SIGLEC6+ DC | **atlas_supported** | Science 2017 + Cancer Research 2025 |
| mregDC | **atlas_supported** | Science Immunology 2022 + Cancer Research 2025 |
| Langerhans cell-like DC | **atlas_supported** | Cancer Research 2025 |

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New subtype entries | 3 (AXL+SIGLEC6+ DC, CD1C+ mono DC, Langerhans cell-like DC) | marker_registry |
| New state entries | 1 (mregDC) | marker_registry |
| New lineage entries | 1 (cDC progenitor) | marker_registry |
| Evidence tier promotions | 4 (cDC1, cDC2, cDC3, pDC → consensus/atlas_supported) | marker_registry |
| New tumor context entries | 1 (Tumor-enriched LAMP3+ DC) | marker_tumor |
| New geneset programs | 4 (DC pleiotropic, cross-presentation, Th priming, type I IFN) | genesets_cancer_signatures.json |
| Existing entry updates | 2 (cDC3 notes, pDC notes) | marker_registry |
| Conflict resolutions | 5 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **DC maturation state matters more than lineage in tumors**: mregDC (LAMP3+ CD274+ IDO1+) is a convergent state from multiple lineages. scLucid tumor interpretation should score mregDC state rather than relying solely on lineage markers.

2. **LAMP3+ DC functional ambiguity**: LAMP3+ DCs express both co-stimulatory (CD40) and co-inhibitory (PD-L1, IDO1) molecules. Their net effect depends on TME context. scLucid should not assume pro- or anti-tumor function based on LAMP3 alone.

3. **DC3 (AXL+SIGLEC6+)**: A rare but functionally important subset that bridges pDC and cDC properties. Potent T cell activation. scLucid should flag this subset when detected.

4. **Conserved DC states**: JEM 2021 meta-analysis confirms DC states are conserved across cancer types and species. This validates scLucid's universal DC markers.

5. **Therapeutic relevance**: DCs are key determinants of immunotherapy efficacy (Nat Rev Clin Oncol 2024). Accurate DC annotation in scLucid directly impacts clinical interpretation.
