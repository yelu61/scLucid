# Marker Curation Batch 05: Myeloid Cell Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells | 2021 | Cell | 10.1016/j.cell.2021.01.010 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 2 | A single-cell pan-cancer analysis to show the variability of tumor-infiltrating myeloid cells in immune checkpoint blockade | 2024 | Nature Communications | 10.1038/s41467-024-50478-8 | pan_cancer_atlas | marker_tumor + geneset |
| 3 | Single-cell resolution characterization of myeloid-derived cell states with implication in cancer outcome | 2024 | Nature Communications | 10.1038/s41467-024-49916-4 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |

---

## Overview: Myeloid Cell Landscape in Cancer

### Existing Registry Entries (Current Status)

| Entry | Current Markers | Evidence Tier | Context |
|-------|----------------|---------------|---------|
| Monocytes | FCN1, S100A8, S100A9, CD14, VCAN, LYZ, S100A12, CD36 | curated_review | Blood/tissue |
| CD14+ Mono | FCN1, HLA-DQB1, S100A9, CSF3R, S100A8 | seed | Classical |
| CD16+ Mono | FCGR3A, LST1, LILRB2, HK3 | seed | Non-classical |
| CD14+CD16+ Mono | NFKBIA, NFKB1, NLRP3, HLA-DQA1 | seed | Intermediate |
| Macrophages | CD68, CD163, APOE, C1QA, C1QB, MSR1, MRC1, CD14, CSF1R | curated_review | Broad |
| Reg TAM | ARG1, MRC1, CD274, CX3CR1 | seed | Tumor |
| Inflam TAM | IL1B, CCL3, CXCL1, CXCL2 | seed | Tumor |
| IFN TAM | IDO1, ISG15, CXCL10, STAT1 | seed | Tumor |
| Angio TAM | VEGFA, SPP1, FGF2, MMP9 | seed | Tumor |
| Prolif TAM | MKI67, CDK1, PCNA | seed | Tumor |
| LA TAM | APOC1, APOE, FABP5 | seed | Tumor |
| TRM TAM | LYVE1, FOLR2, CX3CR1, MERTK | seed | Tumor |
| cDC1 | CLEC9A, XCR1, CADM1 | seed | Cross-tissue |
| cDC2 | CD1C, FCER1A, CLEC10A | seed | Cross-tissue |
| cDC3 | LAMP3, CCR7, FSCN1, CST7 | seed | Cross-tissue |
| pDC | LILRA4, GZMB, IL3RA, TCF4 | seed | Cross-tissue |

### Batch 01 Additions (Cross-tissue atlas)

| Entry | Markers | Source |
|-------|---------|--------|
| FOLR2+ Resident Macrophage | FOLR2, LYVE1, MERTK | Cross-tissue atlas (Nature 2025) |
| CD5L+ Macrophage | CD5L, MARCO, MERTK | Cross-tissue atlas |
| CD14+ Mono (classical) | CD14 | Cross-tissue atlas |
| CD16+ Mono (non-classical) | FCGR3A | Cross-tissue atlas |

---

## Paper-by-Paper Curation

### Article 1: Cell Pan-Cancer Myeloid Atlas (Zhang et al., 2021)

**Dataset**: Pan-cancer single-cell atlas of tumor-infiltrating myeloid cells
**Key methods**: SmartSeq2 + droplet-based scRNA-seq; batch correction; DIG (dissociation-induced gene) filtering
**Key findings**: Comprehensive myeloid heterogeneity across cancer types; tissue dissociation artifacts identified

#### Key Markers from Methods (Directly from Paper)

The Cell 2021 paper explicitly identifies these marker genes in their methods:

| Cell Type | Marker Genes | Context |
|-----------|-------------|---------|
| Macrophage | LYZ, C1QA, C1QB | General identity |
| AT2 cell | SFTPC, SFTPA1, SFTPA2 | Lung contamination |
| Melanocyte | GPNMB, TYR, PMEL, MLPH | Melanoma contamination |

**DIG (Dissociation-Induced Genes)**: Heat shock protein-encoding genes identified as artifacts. These overlap with:
- `artifact.Stress-high` (HSPA1A, HSPB1, etc.)
- `state.TSTR` (T cell stress response)
- Various stress states in B cells, NK cells

**Resolution**: The Cell 2021 DIG list validates scLucid's artifact filtering approach. Confirm that DIG genes are excluded from identity markers.

#### Subtype Refinement (marker_registry)

The Cell 2021 atlas provides pan-cancer validation for existing myeloid subtypes. **Evidence tier promotion recommended**:

| Entry | Current Tier | Proposed Tier | Reason |
|-------|-------------|--------------|--------|
| Monocytes | curated_review | consensus | Validated across cancer types |
| CD14+ Mono | seed | atlas_supported | Cell 2021 + Batch 01 |
| CD16+ Mono | seed | atlas_supported | Cell 2021 + Batch 01 |
| Macrophages | curated_review | consensus | Validated across cancer types |
| Reg TAM | seed | atlas_supported | Cell 2021 pan-cancer |
| Inflam TAM | seed | atlas_supported | Cell 2021 pan-cancer |
| IFN TAM | seed | atlas_supported | Cell 2021 pan-cancer |
| Angio TAM | seed | atlas_supported | Cell 2021 pan-cancer |
| cDC1 | seed | atlas_supported | Cell 2021 + Batch 01 |
| cDC2 | seed | atlas_supported | Cell 2021 + Batch 01 |
| cDC3 | seed | atlas_supported | Cell 2021 + Batch 01 |
| pDC | seed | atlas_supported | Cell 2021 + Batch 01 |

#### Tumor Context (marker_tumor)

**Tumor-Infiltrating Monocyte States**

Cell 2021 identifies distinct monocyte states in tumors vs normal tissues:

```toml
[[cancer_state]]
name = "Tumor-infiltrating classical monocyte"
color = "#1E90FF"
markers = ["FCN1", "S100A8", "S100A9", "CD14", "VCAN", "LYZ", "S100A12"]
negative_markers = ["CD68", "CD163", "MRC1", "CSF1R", "FCGR3A"]

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
source = { title = "A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells", year = "2021", doi = "10.1016/j.cell.2021.01.010" }
notes = "Classical monocytes infiltrating tumors. FCN1+ S100A8/9+ CD14+. Negative for macrophage markers (CD68, CD163, MRC1) indicating non-differentiated state. Precursors to TAM differentiation."
```

```toml
[[cancer_state]]
name = "Tumor-infiltrating non-classical monocyte"
color = "#87CEFA"
markers = ["FCGR3A", "LST1", "LILRB2", "HK3", "FCER1G", "PTPRC"]
negative_markers = ["FCN1", "S100A8", "S100A9", "CD14", "CD68"]

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
source = { title = "A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells", year = "2021", doi = "10.1016/j.cell.2021.01.010" }
notes = "Non-classical monocytes in tumors. FCGR3A+ (CD16+). Patrolling function. Negative FCN1/S100A8/9 confirms non-classical identity."
```

#### Functional Programs (geneset)

```json
"Myeloid_pan_cancer_atlas": {
  "genes": ["LYZ", "CD68", "CSF1R", "CD163", "MRC1", "CD14", "FCGR3A", "ITGAM", "ITGAX", "HLA-DRA", "HLA-DQA1", "CD74", "C1QA", "C1QB", "APOE", "MSR1", "FCN1", "S100A8", "S100A9", "VCAN"],
  "description": "Pan-cancer tumor-infiltrating myeloid cell core signature (Cell 2021)"
}
```

---

### Article 2: Nat Commun - ICB Variability (2024)

**Dataset**: 8 cancer types, 192 tumor samples, 129 patients
**Key finding**: Myeloid cell variability before/after ICB; treatment response-associated myeloid states

#### Tumor Context (marker_tumor)

**ICB-Responsive Myeloid State**

```toml
[[cancer_state]]
name = "ICB-responsive myeloid"
color = "#228B22"
markers = ["CD80", "CD86", "HLA-DRA", "HLA-DQA1", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11"]
negative_markers = ["CD274", "PDCD1LG2", "ARG1", "MRC1", "IL10", "TGFB1"]

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
source = { title = "A single-cell pan-cancer analysis to show the variability of tumor-infiltrating myeloid cells in immune checkpoint blockade", year = "2024", doi = "10.1038/s41467-024-50478-8" }
notes = "Myeloid state associated with ICB response. High co-stimulation (CD80/CD86), antigen presentation (HLA-DR), and T cell recruitment (CCL19/21, CXCL9/10/11). Negative immunosuppressive markers (CD274/PD-L1, PDCD1LG2/PD-L2, ARG1, MRC1)."
```

**ICB-Nonresponsive Myeloid State**

```toml
[[cancer_state]]
name = "ICB-nonresponsive myeloid"
color = "#B22222"
markers = ["CD274", "PDCD1LG2", "ARG1", "MRC1", "IL10", "TGFB1", "VEGFA", "SPP1"]
negative_markers = ["CD80", "CD86", "CXCL9", "CXCL10", "CXCL11", "CCL19", "CCL21"]

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
source = { title = "A single-cell pan-cancer analysis to show the variability of tumor-infiltrating myeloid cells in immune checkpoint blockade", year = "2024", doi = "10.1038/s41467-024-50478-8" }
notes = "Myeloid state associated with ICB resistance. Immunosuppressive phenotype: PD-L1/PD-L2+, ARG1+, MRC1+, IL10+, TGFB1+. Angiogenic (VEGFA+) and pro-tumor (SPP1+). Negative T cell recruitment chemokines."
```

#### Functional Programs (geneset)

```json
"Myeloid_ICB_response": {
  "genes": ["CD80", "CD86", "HLA-DRA", "HLA-DQA1", "HLA-DPA1", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "IDO1", "CD40"],
  "description": "Myeloid program associated with ICB response. Co-stimulatory, antigen-presenting, T cell-recruiting."
},
"Myeloid_ICB_resistance": {
  "genes": ["CD274", "PDCD1LG2", "ARG1", "MRC1", "IL10", "TGFB1", "VEGFA", "SPP1", "MMP9", "CCL22", "CCL17"],
  "description": "Myeloid program associated with ICB resistance. Immunosuppressive, angiogenic, pro-tumor."
}
```

---

### Article 3: Nat Commun - MDC States (2024)

**Dataset**: Integrated scRNA-seq across cancer types; 29 MDC subpopulations
**Key findings**: 
- TREM2+ and FOLR2+ subpopulations distinguished
- TREM2+ PD-1+ and FOLR2+ PDL-2+ as independent prognostic markers
- FOLR2+ macrophages correlate with poor outcomes in ovarian and TNBC

#### Critical Finding: TREM2+ vs FOLR2+ Macrophages

This paper provides **critical validation** for Batch 01 entries:

| Feature | FOLR2+ Macrophage (Batch 01) | TREM2+ Macrophage |
|---------|------------------------------|-------------------|
| Source | Cross-tissue atlas (Nature 2025) | Nat Commun 2024 |
| Markers | FOLR2, LYVE1, MERTK | TREM2, CD9, LPL, GPNMB |
| Context | Resident macrophage | Tumor-associated |
| Prognosis | Not specified | Context-dependent (TREM2 alone unreliable) |
| ICB | Not specified | TREM2+ PD-1+ associated with resistance |

**Key insight**: TREM2 **alone** does not reliably predict prognosis. The combination matters:
- TREM2+ PD-1+ (CD274+): poor prognosis
- TREM2+ without PD-1: variable
- FOLR2+ PDL-2+ (PDCD1LG2+): poor prognosis (ovarian, TNBC)

#### Registry Updates

**Update FOLR2+ Resident Macrophage** (from Batch 01):

Add prognostic context:

```toml
# UPDATE to existing Batch 01 entry:
# Add to notes:
notes = """
Pan-cancer validated (Nat Commun 2024). 
FOLR2+ macrophages correlate with poor clinical outcomes in ovarian and triple-negative breast cancers.
Independent prognostic marker when co-expressed with PDCD1LG2 (PD-L2).
Resident tissue macrophage distinct from inflammatory monocytes.
"""
```

**New Entry: TREM2+ Macrophage**

```toml
[[subtype]]
name = "TREM2+ Macrophage"
color = "#8B4513"
markers = ["TREM2", "CD9", "LPL", "GPNMB", "APOE", "CST7"]
negative_markers = ["FOLR2", "LYVE1", "FCN1", "S100A8", "S100A9"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "tumor_context"
aplies_to = ["Macrophages"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = [
    { title = "Single-cell resolution characterization of myeloid-derived cell states with implication in cancer outcome", year = "2024", doi = "10.1038/s41467-024-49916-4" },
    { title = "A pan-cancer single-cell transcriptional atlas of tumor infiltrating myeloid cells", year = "2021", doi = "10.1016/j.cell.2021.01.010" }
]
notes = """
TREM2+ macrophages are tumor-associated myeloid cells distinct from FOLR2+ resident macrophages.
Markers: TREM2, CD9, LPL, GPNMB, APOE, CST7.
Negative FOLR2/LYVE1 distinguishes from resident macrophages.
Negative FCN1/S100A8 distinguishes from inflammatory monocytes.
CRITICAL: TREM2 alone is NOT a reliable prognostic marker.
Context-dependent prognosis: TREM2+ PD-1+ (CD274+) associated with poor outcomes.
"""
```

#### Tumor Context (marker_tumor)

**TREM2+ PD-1+ Macrophage (Poor Prognosis)**

```toml
[[cancer_state]]
name = "TREM2+ PD-1+ macrophage"
color = "#8B0000"
markers = ["TREM2", "CD9", "LPL", "GPNMB", "CD274", "PDCD1"]
negative_markers = ["FOLR2", "LYVE1", "FCN1", "S100A8", "CD80", "CD86"]

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
source = { title = "Single-cell resolution characterization of myeloid-derived cell states with implication in cancer outcome", year = "2024", doi = "10.1038/s41467-024-49916-4" }
notes = "TREM2+ macrophages co-expressing PD-1 (CD274/PD-L1). Independent prognostic marker for poor outcomes. Immunosuppressive phenotype. Negative FOLR2/LYVE1 confirms non-resident. Negative co-stimulatory markers (CD80/CD86)."
```

**FOLR2+ PDL-2+ Macrophage (Poor Prognosis in Ovarian/TNBC)**

```toml
[[cancer_state]]
name = "FOLR2+ PDL-2+ macrophage"
color = "#800080"
markers = ["FOLR2", "LYVE1", "MERTK", "PDCD1LG2", "CD163", "MRC1"]
negative_markers = ["TREM2", "CD9", "FCN1", "S100A8", "CD80"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["Ovarian Cancer", "Triple-Negative Breast Cancer"]
scope = "cancer_type_specific"
applies_to = ["Ovarian", "Breast"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Single-cell resolution characterization of myeloid-derived cell states with implication in cancer outcome", year = "2024", doi = "10.1038/s41467-024-49916-4" }
notes = "FOLR2+ macrophages co-expressing PD-L2 (PDCD1LG2). Specific poor prognosis in ovarian and triple-negative breast cancers. M2-like phenotype (CD163+, MRC1+). Negative TREM2/CD9 distinguishes from TREM2+ macrophages."
```

#### Functional Programs (geneset)

```json
"TREM2_macrophage_signature": {
  "genes": ["TREM2", "CD9", "LPL", "GPNMB", "APOE", "CST7", "SPP1", "CTSB", "CTSD", "C1QA", "C1QB", "C1QC"],
  "description": "TREM2+ macrophage signature. Tumor-associated, lipid metabolism, phagocytic."
},
"FOLR2_macrophage_signature": {
  "genes": ["FOLR2", "LYVE1", "MERTK", "CD163", "MRC1", "COLEC12", "STAB1", "STAB2", "C1QA", "C1QB"],
  "description": "FOLR2+ resident macrophage signature. Tissue-resident, homeostatic."
},
"Myeloid_immunosuppressive": {
  "genes": ["CD274", "PDCD1LG2", "ARG1", "MRC1", "IL10", "TGFB1", "VEGFA", "IDO1", "CCL22", "CCL17"],
  "description": "Myeloid immunosuppressive program. ICB resistance-associated."
},
"Myeloid_immune_activated": {
  "genes": ["CD80", "CD86", "HLA-DRA", "HLA-DQA1", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CD40", "IDO1"],
  "description": "Myeloid immune-activated program. ICB response-associated."
}
```

---

## Summary: Key Conflicts and Resolutions

### 1. TREM2+ vs FOLR2+ Macrophages

**Issue**: Two distinct tumor macrophage populations with different markers and prognosis.

**Resolution**:
- Create separate entries for TREM2+ and FOLR2+ macrophages
- Use negative markers to distinguish:
  - TREM2+: negative FOLR2, LYVE1
  - FOLR2+: negative TREM2, CD9
- Both are tumor-context (`cancer_state`), not primary identity

### 2. TREM2 as Prognostic Marker

**Issue**: TREM2 alone is NOT a reliable prognostic marker (Nat Commun 2024).

**Resolution**:
- Document in notes: "TREM2 alone unreliable; requires context (PD-1 co-expression)"
- Add TREM2+ PD-1+ as specific cancer_state with poor prognosis
- Do NOT promote TREM2 as standalone prognostic marker

### 3. FOLR2+ Macrophage Prognosis

**Issue**: FOLR2+ was initially described as "resident" (Batch 01) but now shown to correlate with poor outcomes in specific cancers.

**Resolution**:
- FOLR2+ macrophages are **resident identity** but their **state** (PD-L2 co-expression) determines prognosis
- Keep FOLR2+ as subtype (resident identity)
- Add FOLR2+ PDL-2+ as cancer_state (poor prognosis in ovarian/TNBC)
- This reflects the principle: identity vs state separation

### 4. Batch Effects and DIG

**Issue**: Cell 2021 identifies heat shock proteins (HSPA1A, etc.) as dissociation artifacts.

**Resolution**:
- Validates existing `artifact.Stress-high` entry
- Confirms HSPA1A should NOT be used as identity marker
- TSTR (T cell stress response) and NK stress states are biological, not artifact
- Document distinction: artifact = dissociation-induced; state = TME-induced

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| Evidence tier promotions | 11 | marker_registry |
| New subtype entries | 1 (TREM2+ macrophage) | marker_registry |
| New tumor context entries | 6 | marker_tumor |
| Existing entry updates | 1 (FOLR2+ notes) | marker_registry |
| New geneset programs | 6 | genesets_cancer_signatures.json |
| Conflict resolutions | 4 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **Myeloid heterogeneity is context-dependent**: TREM2+ vs FOLR2+ macrophages have different functions and prognostic implications. scLucid should not treat all macrophages as equivalent.

2. **Combination markers for prognosis**: Single markers (TREM2 alone) are unreliable. Combinations (TREM2+ PD-1+, FOLR2+ PDL-2+) are required for prognostic prediction.

3. **ICB response prediction**: Myeloid states (CD80/CD86+ vs CD274/PDCD1LG2+) are stronger predictors than myeloid abundance alone.

4. **Artifact awareness**: Cell 2021's DIG list validates scLucid's artifact filtering. Heat shock genes should be excluded from identity annotation but can be used for state annotation (TSTR, stress response).
