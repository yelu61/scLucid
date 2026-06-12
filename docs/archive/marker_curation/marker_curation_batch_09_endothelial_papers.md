# Marker Curation Batch 09: Endothelial Cell Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | Microenvironmental determinants of endothelial cell heterogeneity | 2025 | Nature Reviews Molecular Cell Biology | 10.1038/s41580-024-00825-w | review | marker_registry (naming + framework) |
| 2 | An organotypic atlas of human vascular cells | 2024 | Nature Medicine | 10.1038/s41591-024-03376-x | single_cell_atlas | marker_registry + marker_tissue + geneset |
| 3 | A systems view of the vascular endothelium in health and disease | 2024 | Cell | 10.1016/j.cell.2024.07.012 | review | marker_registry (framework) |
| 4 | Pan-cancer integrative analyses dissect the remodeling of endothelial cells in human cancers | 2024 | National Science Review | 10.1093/nsr/nwae231 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 5 | Pan-cancer landscape of tumour endothelial cells pinpoints insulin receptor as a novel antiangiogenic target and predicts immunotherapy response | 2023 | Clinical and Translational Medicine | 10.1002/ctm2.1501 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 6 | Understanding tumour endothelial cell heterogeneity and function from single-cell omics | 2023 | Nature Reviews Cancer | 10.1038/s41568-023-00591-5 | review | marker_registry (naming + validation) |
| 7 | Tumour vasculature at single-cell resolution | 2024 | Nature | 10.1038/s41586-024-07698-1 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |

---

## Overview: Endothelial Cell Biology Context

### Existing Registry Entries (Current Status)

| Entry | Current Markers | Evidence Tier | Notes |
|-------|----------------|---------------|-------|
| Endothelial cells | PECAM1, VWF, CDH5, RAMP2, ENG, CLDN5, SELE, ICAM1, ACKR1 | curated_review | Broad identity |
| Arterial EC | SEMA3G, GJA5, BMX | seed | Batch 01 added IGFBP3 |
| Venous EC | ACKR1, SELP, NRG1 | seed | Confirmed |
| Capillary EC | RGCC, CA4, CD36 | seed | Batch 01 added VIPR1 |
| Lymphatic EC | PROX1, LYVE1, CCL21 | seed | Confirmed |
| TipEC | COL4A1, KDR, ESM1 | seed | Needs expansion |
| EndoMT | PDGFRB | seed | Partial |

### Batch 01 Additions (Cross-tissue atlas)

| Entry | Markers | Source |
|-------|---------|--------|
| Arterial EC | IGFBP3 | Cross-tissue atlas |
| Capillary EC | CA4, VIPR1 | Cross-tissue atlas |
| Venous EC | ACKR1 | Cross-tissue atlas |
| Lymphatic EC | LYVE1 | Cross-tissue atlas |

---

## Paper-by-Paper Curation

### Article 2: Nature Medicine 2024 — Organotypic Atlas of Human Vascular Cells

**Dataset**: 19 human organs, ~67,000 cells, 62 donors
**Key finding**: 42 vascular cell states; organotypic and angiotypic signatures

#### New Subtype Entries (`marker_registry`)

**1. Splenic Littoral Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "Splenic littoral endothelial cell"
color = "#800080"
markers = ["CD8A", "CD8B", "LYVE1", "STAB2", "CLEC1B", "FCER1G", "TIMP3"]
negative_markers = ["PECAM1", "VWF", "CDH5", "ENG", "CLDN5"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tissue_specific"
applies_to = ["Spleen"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "An organotypic atlas of human vascular cells", year = "2024", doi = "10.1038/s41591-024-03376-x" }
notes = """
Splenic littoral endothelial cells. Specialized for blood filtration and immune surveillance.
Express CD8 (not T cell CD8 but sialoadhesin CD8), LYVE1, STAB2 (scavenger receptor), CLEC1B (CLEC-2 ligand).
Negative conventional endothelial markers (PECAM1, VWF, CDH5, ENG, CLDN5) — note: this is CD8+ littoral cell identity, not PECAM1+ sinusoidal.
Tissue-specific to spleen.
"""
```

**2. Blood-Brain Barrier Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "Blood-brain barrier endothelial cell"
color = "#4682B4"
markers = ["CLDN5", "OCLN", "SLC2A1", "ABCB1", "ABCG2", "SLCO1C1", "MFSD2A", "TJP1"]
negative_markers = ["PLVAP", "ACKR1", "SELE", "ICAM1", "VCAM1"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tissue_specific"
applies_to = ["Brain", "CNS"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "An organotypic atlas of human vascular cells", year = "2024", doi = "10.1038/s41591-024-03376-x" }
notes = """
Blood-brain barrier (BBB) endothelial cells. Tight junction proteins (CLDN5, OCLN, TJP1).
Transporters: SLC2A1 (GLUT1), ABCB1 (P-gp), ABCG2 (BCRP), SLCO1C1, MFSD2A.
Negative fenestration markers (PLVAP) and inflammatory adhesion molecules (ACKR1, SELE, ICAM1, VCAM1).
Highly specialized for CNS homeostasis.
"""
```

**3. Large Artery Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "Large artery endothelial cell"
color = "#B22222"
markers = ["IGFBP3", "SEMA3G", "GJA5", "BMX", "HEY2", "EFNB2", "CXCL12"]
negative_markers = ["ACKR1", "SELP", "NRG1", "PLVAP", "CA4"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "all"
applies_to = ["Artery"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "An organotypic atlas of human vascular cells", year = "2024", doi = "10.1038/s41591-024-03376-x" }
notes = "Large artery endothelial cells (aorta, elastic arteries). IGFBP3+ SEMA3G+ GJA5+ BMX+ HEY2+ EFNB2+. High CXCL12. Negative venous markers (ACKR1, SELP, NRG1) and capillary markers (PLVAP, CA4). Angiotypic transitional signature from large to small caliber vessels."
```

**4. Small Artery/Arteriole Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "Small artery endothelial cell"
color = "#CD5C5C"
markers = ["IGFBP3", "SEMA3G", "GJA5", "S100A4", "ITGA4", "CD36"]
negative_markers = ["BMX", "HEY2", "ACKR1", "SELP", "CA4"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "all"
applies_to = ["Artery"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "An organotypic atlas of human vascular cells", year = "2024", doi = "10.1038/s41591-024-03376-x" }
notes = "Small artery/arteriole endothelial cells. Transitional between large artery and capillary. IGFBP3+ SEMA3G+ GJA5+ S100A4+ ITGA4+ CD36+. Negative large artery (BMX, HEY2) and venous (ACKR1, SELP) markers."
```

#### Functional Programs (geneset)

```json
"Endothelial_organotypic_signatures": {
  "genes": ["CLDN5", "OCLN", "TJP1", "SLC2A1", "ABCB1", "ABCG2", "MFSD2A", "CD8A", "LYVE1", "STAB2", "CLEC1B", "IGFBP3", "SEMA3G", "GJA5", "BMX", "HEY2", "PROX1", "LYVE1", "CCL21", "ACKR1", "SELP", "RGCC", "CA4", "CD36", "VIPR1"],
  "description = "Organotypic endothelial cell signatures across 19 human tissues. Includes BBB, splenic littoral, arterial, venous, capillary, and lymphatic markers."
},
"Endothelial_FOXF1_lung": {
  "genes": ["FOXF1", "SOX17", "HEY2", "EFNB2", "CXCL12", "IGFBP3", "SEMA3G"],
  "description": "FOXF1-driven lung vascular subpopulation program. Tissue-specific transcriptional regulation."
},
"Endothelial_Notch_Wnt": {
  "genes": ["NOTCH1", "NOTCH4", "DLL4", "JAG1", "JAG2", "WNT2", "WNT7B", "LEF1", "TCF7", "RSPO3", "FZD4"],
  "description": "Endothelial Notch-Wnt signaling program. Angiotypic and organotypic communication."
},
"Endothelial_retinoic_acid": {
  "genes": ["RBP1", "ALDH1A2", "ALDH1A3", "RDH10", "STRA6", "RAR", "RXR"],
  "description": "Endothelial retinoic acid signaling program. Organotypic vascular niche maintenance."
}
```

---

### Article 4: National Science Review 2024 — Pan-Cancer EC Remodeling

**Dataset**: 575 cancer patients, 19 solid tumor types
**Key findings**:
- CXCR4+ tip cells: prominent angiogenic phenotype
- SELE+ veins: proinflammatory phenotype
- Tumors: increased CXCR4+ tip cells, depleted SELE+ veins
- Associated with anti-angiogenic therapy and immunotherapy response

#### New Subtype Entries (`marker_registry`)

**CXCR4+ Tip Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "CXCR4+ tip endothelial cell"
color = "#FF4500"
markers = ["CXCR4", "KDR", "FLT1", "ESM1", "PGF", "COL4A1", "NRP1", "ANGPT2"]
negative_markers = ["ACKR1", "SELP", "PROX1", "LYVE1", "CLDN5"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Pan-cancer integrative analyses dissect the remodeling of endothelial cells in human cancers", year = "2024", doi = "10.1093/nsr/nwae231" }
notes = """
CXCR4+ tip endothelial cells. Prominent angiogenic phenotype across 19 solid tumor types.
Express tip cell markers (CXCR4, KDR, FLT1, ESM1, PGF, COL4A1, NRP1, ANGPT2).
Increased prevalence in tumor vs adjacent non-tumor tissue.
Associated with anti-angiogenic therapy targets.
Negative venous (ACKR1, SELP), lymphatic (PROX1, LYVE1), and BBB (CLDN5) markers.
"""
```

**SELE+ Venous Endothelial Cell (Proinflammatory)**

```toml
[[compartment.minor.minor]]
name = "SELE+ proinflammatory venous endothelial cell"
color = "#4169E1"
markers = ["SELE", "SELP", "VCAM1", "ICAM1", "ACKR1", "CCL2", "CXCL1", "CXCL8"]
negative_markers = ["CXCR4", "KDR", "ESM1", "PGF", "SEMA3G", "GJA5"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Pan-cancer integrative analyses dissect the remodeling of endothelial cells in human cancers", year = "2024", doi = "10.1093/nsr/nwae231" }
notes = """
SELE+ proinflammatory venous endothelial cells. Depleted in tumor vs adjacent non-tumor tissue.
Express adhesion molecules (SELE, SELP, VCAM1, ICAM1, ACKR1) and inflammatory chemokines (CCL2, CXCL1, CXCL8).
Proinflammatory phenotype. May facilitate immune cell infiltration.
Negative tip cell markers (CXCR4, KDR, ESM1, PGF) and arterial markers (SEMA3G, GJA5).
"""
```

#### Tumor Context (`marker_tumor`)

**Angiogenic Tumor Endothelial Cell**

```toml
[[cancer_state]]
name = "Angiogenic tumor endothelial cell"
color = "#FF4500"
markers = ["CXCR4", "KDR", "FLT1", "ESM1", "PGF", "COL4A1", "NRP1", "ANGPT2", "VEGFA"]
negative_markers = ["ACKR1", "SELP", "SELE", "VCAM1", "CLDN5"]

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
source = { title = "Pan-cancer integrative analyses dissect the remodeling of endothelial cells in human cancers", year = "2024", doi = "10.1093/nsr/nwae231" }
notes = "Angiogenic tumor endothelial cells (TEC). CXCR4+ tip phenotype. Target of anti-angiogenic therapies. Negative venous adhesion molecules (ACKR1, SELP, SELE, VCAM1) and BBB tight junctions (CLDN5)."
```

---

### Article 5: CTM 2023 — Pan-Cancer TEC Landscape

**Key findings**:
- INSR+ tip ECs (C1): 67.6% of all tip ECs, elevated in 8/12 cancer types
- PGF+ tip ECs (C2): secondary cluster
- INSR is critical metabolic gene
- FCGR2B+ capillary ECs: diminished in liver and lung tumors
- Lymphatic tip cells increased in tumors

#### New Subtype Entries (`marker_registry`)

**INSR+ Tip Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "INSR+ tip endothelial cell"
color = "#DC143C"
markers = ["INSR", "KDR", "FLT1", "ESM1", "PGF", "VEGFA", "NRP1", "ANGPT2"]
negative_markers = ["PGF", "FCGR2B", "CA4", "RGCC"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Pan-cancer landscape of tumour endothelial cells pinpoints insulin receptor as a novel antiangiogenic target and predicts immunotherapy response", year = "2023", doi = "10.1002/ctm2.1501" }
notes = """
INSR+ tip endothelial cells. Novel tip EC cluster identified in pan-cancer TEC atlas.
Comprising 67.6% of all tip ECs. Elevated in 8/12 cancer types.
INSR is a critical metabolic gene. Expresses VEGFR1/VEGFR2 (KDR/FLT1) — angiogenic.
Potential target for anti-angiogenic therapy.
Negative capillary markers (FCGR2B, CA4, RGCC).
"""
```

**PGF+ Tip Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "PGF+ tip endothelial cell"
color = "#FF6347"
markers = ["PGF", "KDR", "FLT1", "ESM1", "NRP1", "ANGPT2", "VEGFA"]
negative_markers = ["INSR", "FCGR2B", "CA4", "RGCC"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Pan-cancer landscape of tumour endothelial cells pinpoints insulin receptor as a novel antiangiogenic target and predicts immunotherapy response", year = "2023", doi = "10.1002/ctm2.1501" }
notes = "PGF+ tip endothelial cells. Secondary tip EC cluster. PGF (placental growth factor) is a VEGF family member. Negative INSR distinguishes from INSR+ tip ECs."
```

**FCGR2B+ Capillary Endothelial Cell (Normal Tissue)**

```toml
[[compartment.minor.minor]]
name = "FCGR2B+ scavenger capillary endothelial cell"
color = "#20B2AA"
markers = ["FCGR2B", "CA4", "RGCC", "VIPR1", "CD36", "PLVAP"]
negative_markers = ["CXCR4", "INSR", "SEMA3G", "ACKR1", "PROX1"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tissue_specific"
applies_to = ["Liver", "Lung"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Pan-cancer landscape of tumour endothelial cells pinpoints insulin receptor as a novel antiangiogenic target and predicts immunotherapy response", year = "2023", doi = "10.1002/ctm2.1501" }
notes = "FCGR2B+ scavenger capillary endothelial cells. Abundant in normal liver and lung but substantially diminished in tumors. Scavenger function (FCGR2B = Fc gamma receptor IIb). Capillary markers (CA4, RGCC, VIPR1, CD36, PLVAP). Negative tip (CXCR4, INSR), arterial (SEMA3G), venous (ACKR1), and lymphatic (PROX1) markers."
```

#### Tumor Context (`marker_tumor`)

**Lymphatic Tip Endothelial Cell (Tumor)**

```toml
[[cancer_state]]
name = "Lymphatic tip endothelial cell (tumor)"
color = "#9370DB"
markers = ["PROX1", "LYVE1", "CCL21", "VEGFR3", "FLT4", "NRP2", "ANGPT2"]
negative_markers = ["PECAM1hi", "VWF", "CLDN5"]

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
source = { title = "Pan-cancer landscape of tumour endothelial cells pinpoints insulin receptor as a novel antiangiogenic target and predicts immunotherapy response", year = "2023", doi = "10.1002/ctm2.1501" }
notes = "Lymphatic tip endothelial cells increased in tumors. Lymph-angiogenesis marker (VEGFR3/FLT4, NRP2, ANGPT2). Express lymphatic identity (PROX1, LYVE1, CCL21). Negative blood vascular markers (PECAM1hi, VWF, CLDN5). Associated with tumor metastasis."
```

---

### Article 7: Nature 2024 — Tumour Vasculature at Single-Cell Resolution

**Dataset**: ~200,000 cells, 372 donors, 31 cancer types
**Key findings**:
- Tumour angiogenesis initiates from **venous ECs** and extends toward **arterial ECs**
- Angiogenic stages: **SI → SII → SIII**
- **APLN+ TipSI** cells: associated with disease progression, poor prognosis, predict anti-VEGF therapy response
- **TipSIII** cells: increased Notch signalling
- **Stalk cells**: transition from high chemokine to elevated **TEK (Tie2)** expression
- **Lymphatic ECs**: two lineages — lymphangiogenesis and **antigen presentation**
- **BASP1+ matrix-producing pericytes**: ER stress-associated, proangiogenic
- Neovascular ECs shape immunosuppressive microenvironment

#### New Subtype Entries (`marker_registry`)

**APLN+ Tip Endothelial Cell (TipSI)**

```toml
[[compartment.minor.minor]]
name = "APLN+ tip endothelial cell (TipSI)"
color = "#FF4500"
markers = ["APLN", "KDR", "FLT1", "ESM1", "PGF", "NRP1", "ANGPT2", "COL4A1"]
negative_markers = ["TEK", "ANGPT1", "SEMA3G", "ACKR1", "SELP"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Tumour vasculature at single-cell resolution", year = "2024", doi = "10.1038/s41586-024-07698-1" }
notes = """
APLN+ tip endothelial cells at angiogenic stage I (TipSI).
Leading cells at the tips of vascular angiogenic sprouts.
Express APLN (apelin), KDR, FLT1, ESM1, PGF, NRP1, ANGPT2, COL4A1.
Associated with disease progression and poor prognosis.
Hold promise for predicting response to anti-VEGF therapy.
Negative stalk cell markers (TEK/Tie2, ANGPT1) and arterial/venous markers (SEMA3G, ACKR1, SELP).
"""
```

**Tip Endothelial Cell (TipSIII)**

```toml
[[compartment.minor.minor]]
name = "Tip endothelial cell (TipSIII)"
color = "#FF6347"
markers = ["KDR", "FLT1", "ESM1", "NRP1", "ANGPT2", "DLL4", "JAG1", "NOTCH1", "HEY1"]
negative_markers = ["APLN", "TEK", "ANGPT1", "ACKR1"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Tumour vasculature at single-cell resolution", year = "2024", doi = "10.1038/s41586-024-07698-1" }
notes = "Tip endothelial cells at angiogenic stage III (TipSIII). Advanced angiogenic stage with increased Notch signalling (DLL4, JAG1, NOTCH1, HEY1). More mature tip phenotype than TipSI. Negative APLN (TipSI marker) and stalk/venous markers (TEK, ANGPT1, ACKR1)."
```

**Stalk Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "Stalk endothelial cell"
color = "#DAA520"
markers = ["TEK", "ANGPT1", "KDR", "FLT1", "ESM1", "CCL2", "CCL8", "CXCL12"]
negative_markers = ["APLN", "NRP1", "ANGPT2", "SEMA3G", "ACKR1"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Tumour vasculature at single-cell resolution", year = "2024", doi = "10.1038/s41586-024-07698-1" }
notes = "Stalk endothelial cells following tip cells in angiogenic sprouts. Express TEK (Tie2) and ANGPT1 (Ang1). Transition from high chemokine expression (CCL2, CCL8, CXCL12) to elevated TEK expression. Negative tip markers (APLN, NRP1, ANGPT2) and arterial/venous markers (SEMA3G, ACKR1)."
```

**Lymphatic Antigen-Presenting Endothelial Cell**

```toml
[[compartment.minor.minor]]
name = "Lymphatic antigen-presenting endothelial cell"
color = "#9370DB"
markers = ["PROX1", "LYVE1", "CCL21", "HLA-DRA", "HLA-DQA1", "CD74", "CIITA", "FLT4"]
negative_markers = ["PECAM1hi", "VWF", "CLDN5", "KDR", "FLT1"]

[compartment.minor.minor.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "stromal"
lineage = "endothelial"
scope = "tumor_context"
applies_to = ["Tumor"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Tumour vasculature at single-cell resolution", year = "2024", doi = "10.1038/s41586-024-07698-1" }
notes = "Lymphatic endothelial cells with antigen-presenting capacity. Express MHC-II molecules (HLA-DRA, HLA-DQA1, CD74, CIITA) alongside lymphatic identity (PROX1, LYVE1, CCL21, FLT4/VEGFR3). One of two lymphatic EC lineages (the other is lymphangiogenesis). Negative blood vascular markers (PECAM1hi, VWF, CLDN5, KDR, FLT1)."
```

#### Tumor Context (`marker_tumor`)

**Neovascular Endothelial Cell (Immunosuppressive)**

```toml
[[cancer_state]]
name = "Neovascular endothelial cell (immunosuppressive)"
color = "#8B0000"
markers = ["KDR", "FLT1", "ESM1", "CXCL12", "CCL2", "CD274", "PDCD1LG2", "IDO1"]
negative_markers = ["APLN", "TEK", "ANGPT1", "CLDN5", "HLA-DRA"]

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
source = { title = "Tumour vasculature at single-cell resolution", year = "2024", doi = "10.1038/s41586-024-07698-1" }
notes = "Neovascular endothelial cells shaping immunosuppressive microenvironment conducive to angiogenesis. Express immunosuppressive molecules (CXCL12, CCL2, CD274/PD-L1, PDCD1LG2/PD-L2, IDO1). Negative APLN (TipSI), TEK/ANGPT1 (stalk), CLDN5 (BBB), and HLA-DRA (antigen-presenting)."
```

#### Angiogenic Trajectory

```
Venous EC (ACKR1+, SELP+)
    ↓ (initiation)
Angiogenic SI
    ├── APLN+ TipSI (leading cells)
    └── Stalk (chemokine-high)
    ↓
Angiogenic SII
    ├── TipSII (intermediate)
    └── Stalk (transitional)
    ↓
Angiogenic SIII
    ├── TipSIII (Notch-high, DLL4+, JAG1+)
    └── Stalk (TEK+/Tie2+high)
    ↓
Arterial EC (SEMA3G+, GJA5+)
```

#### Functional Programs (geneset)

```json
"Tumour_angiogenesis_trajectory": {
  "genes": ["ACKR1", "SELP", "NRG1", "APLN", "KDR", "FLT1", "ESM1", "PGF", "NRP1", "ANGPT2", "DLL4", "JAG1", "NOTCH1", "HEY1", "TEK", "ANGPT1", "SEMA3G", "GJA5", "BMX", "IGFBP3"],
  "description": "Tumour angiogenesis trajectory from venous to arterial. Stages: venous → TipSI (APLN+) → TipSIII (Notch+) → stalk (TEK+/Tie2+) → arterial."
},
"TipSI_antiVEGF_prediction": {
  "genes": ["APLN", "KDR", "FLT1", "ESM1", "PGF", "NRP1", "ANGPT2", "COL4A1"],
  "description": "APLN+ TipSI signature for predicting anti-VEGF therapy response. Associated with disease progression and poor prognosis."
},
"Stalk_cell_signature": {
  "genes": ["TEK", "ANGPT1", "KDR", "FLT1", "ESM1", "CCL2", "CCL8", "CXCL12", "SEMA3G"],
  "description": "Stalk endothelial cell signature. TEK+/Tie2+ following tip cells. Chemokine-to-Tie2 transition."
},
"Lymphatic_antigen_presentation": {
  "genes": ["PROX1", "LYVE1", "CCL21", "FLT4", "HLA-DRA", "HLA-DQA1", "CD74", "CIITA", "CD80", "CD86"],
  "description": "Lymphatic endothelial cell antigen presentation program. MHC-II and co-stimulatory molecule expression."
},
"Pericyte_BASP1_proangiogenic": {
  "genes": ["BASP1", "PDGFRB", "RGS5", "COL1A1", "COL1A2", "FN1", "POSTN", "MMP2", "HSPA5", "DDIT3", "ATF4"],
  "description": "BASP1+ matrix-producing pericyte program. ER stress-associated (HSPA5, DDIT3, ATF4), proangiogenic."
}
```

---

### Article 6: Nature Reviews Cancer 2023 — TEC Heterogeneity Review

**Key findings**:
- Universal TEC markers: ACKR1, PLVAP, IGFBP3
- Immunomodulation and ECM organization are common enriched signatures
- Lack of uniform nomenclature

#### Validation of Existing Entries

| Marker | Context | Status |
|--------|---------|--------|
| ACKR1 | Universal TEC marker | Confirmed (existing in broad endothelial) |
| PLVAP | Universal TEC marker | **Add to broad endothelial** |
| IGFBP3 | Universal TEC marker | Confirmed (existing in arterial) |

**Recommended update** to broad Endothelial cells entry:

```toml
# UPDATE existing entry:
# Current: PECAM1, VWF, CDH5, RAMP2, ENG, CLDN5, SELE, ICAM1, ACKR1
# Add: PLVAP
markers = ["PECAM1", "VWF", "CDH5", "RAMP2", "ENG", "CLDN5", "SELE", "ICAM1", "ACKR1", "PLVAP"]
```

---

### Articles 1, 3: Reviews (Nat Rev Mol Cell Biol 2025, Cell 2024)

**Key frameworks**:

| Framework Element | Description | Registry Impact |
|-------------------|-------------|-----------------|
| Organotypic specialization | ECs acquire organ-specific functions | Scope classification |
| Capillary zonation | Arterial → capillary → venous zonation | Subtype hierarchy |
| Perivascular niche | Pericytes, smooth muscle cells influence ECs | Interaction notes |
| Angiocrine factors | Tissue-specific factors from ECs | Functional programs |
| Disease-related erasure | Loss of organotypic signatures in disease | State entries |

---

## Updated Endothelial Cell Registry Hierarchy

### Proposed New Structure

```
Endothelial cells (compartment) [EXPANDED]
├── Arterial EC [EXISTING - expanded]
│   ├── Large artery EC [NEW]
│   ├── Small artery/arteriole EC [NEW]
│   └── Artery-capillary transitional [NEW]
├── Capillary EC [EXISTING - expanded]
│   ├── General capillary [EXISTING]
│   ├── FCGR2B+ scavenger capillary [NEW]
│   └── TipEC [EXISTING - expanded]
│       ├── APLN+ TipSI [NEW - Nature 2024]
│       ├── TipSIII (Notch-high) [NEW - Nature 2024]
│       ├── CXCR4+ tip EC [NEW]
│       ├── INSR+ tip EC [NEW]
│       └── PGF+ tip EC [NEW]
├── Venous EC [EXISTING]
│   └── SELE+ proinflammatory venous [NEW]
├── Lymphatic EC [EXISTING - expanded]
│   ├── General lymphatic [EXISTING]
│   ├── Lymphatic tip EC (tumor) [NEW]
│   └── Lymphatic antigen-presenting EC [NEW - Nature 2024]
├── Organotypic ECs [NEW]
│   ├── Blood-brain barrier EC [NEW]
│   └── Splenic littoral EC [NEW]
├── Stalk cells [NEW - Nature 2024]
│   └── Stalk endothelial cell (TEK+/Tie2+high) [NEW]
├── Pericytes [NEW - Nature 2024]
│   └── BASP1+ matrix-producing pericyte [NEW]
└── Tumor contexts [NEW]
    ├── Angiogenic TEC [cancer_state]
    ├── Lymphatic tip TEC [cancer_state]
    └── Neovascular immunosuppressive EC [cancer_state - Nature 2024]
```

---

## Curation Notes and Conflicts

### 1. PLVAP as Universal TEC Marker

**Issue**: PLVAP is added to broad endothelial markers but is also a fenestration marker.

**Resolution**:
- PLVAP marks fenestrated/sinusoidal endothelium
- Present in many TECs (universal TEC marker per Nat Rev Cancer 2023)
- Absent in BBB (non-fenestrated) — use as negative marker for BBB
- Add to broad endothelial entry with note

### 2. Tip EC Heterogeneity

**Issue**: Multiple tip EC subtypes identified:
- Existing: COL4A1, KDR, ESM1
- Batch 01: KDR, ESM1 (confirmed)
- Nat Sci Rev 2024: CXCR4+ tip
- CTM 2023: INSR+ tip (67.6%), PGF+ tip

**Resolution**:
- TipEC is a **functional state**, not a single subtype
- Multiple tip EC clusters exist with distinct marker profiles
- Keep existing TipEC as broad entry
- Add CXCR4+, INSR+, PGF+ as specific subtypes
- Note: "Tip ECs are heterogeneous. CXCR4+, INSR+, and PGF+ represent distinct angiogenic programs."

### 3. BBB vs Blood Vessel EC

**Issue**: BBB ECs express CLDN5 but lack PLVAP, ACKR1, SELE.

**Resolution**:
- BBB EC: CLDN5+ OCLN+ TJP1+ SLC2A1+ ABCB1+ ABCG2+
- General EC: PLVAP+ ACKR1+ SELE+
- Use negative markers to distinguish:
  - BBB negative: PLVAP, ACKR1, SELE, ICAM1
  - General EC negative: tight junctions (incomplete — many ECs express CLDN5 at low levels)

### 4. Splenic Littoral vs Sinusoidal EC

**Issue**: Splenic littoral cells express CD8 (sialoadhesin) but not PECAM1.

**Resolution**:
- Splenic littoral cells are a specialized endothelial population
- CD8+ STAB2+ CLEC1B+ confirms identity
- PECAM1- VWF- distinguishes from conventional sinusoidal EC
- Add note: "Splenic littoral cells are PECAM1-negative. Do not confuse with conventional endothelial cells."

### 5. EndoMT Confusion

**Issue**: EndoMT (endothelial-to-mesenchymal transition) overlaps with CAF markers.

**Resolution**:
- EndoMT: PDGFRB+ (existing entry)
- CAF: COL1A1+ COL1A2+ ACTA2+
- Co-expression required for EndoMT identity:
  - Positive: PECAM1 (or residual), VWF, PDGFRB
  - Negative: full CAF program
- Add note: "EndoMT requires residual endothelial identity (PECAM1/VWF) plus mesenchymal transition (PDGFRB). Pure CAFs are PECAM1-negative."

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New subtype entries | 11 (BBB EC, splenic littoral, large/small artery, CXCR4+ tip, INSR+ tip, PGF+ tip, FCGR2B+ capillary, **APLN+ TipSI**, **TipSIII**, **Stalk EC**, **Lymphatic antigen-presenting EC**) | marker_registry |
| New tumor context entries | 3 (angiogenic TEC, lymphatic tip TEC, **neovascular immunosuppressive EC**) | marker_tumor |
| Existing entry expansion | 1 (broad endothelial: +PLVAP) | marker_registry |
| New geneset programs | 8 (organotypic, FOXF1 lung, Notch-Wnt, retinoic acid, **tumour angiogenesis trajectory**, **TipSI anti-VEGF prediction**, **stalk cell signature**, **lymphatic antigen presentation**, **BASP1+ pericyte**) | genesets_cancer_signatures.json |
| Review validations | 3 (naming frameworks) | marker_registry |
| Conflict resolutions | 5 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **Endothelial cells are organotypic**: ECs in different organs have distinct molecular profiles (Nat Med 2024). scLucid should consider tissue context when interpreting endothelial markers.

2. **Tip EC heterogeneity**: Tip ECs are not homogeneous. CXCR4+, INSR+, PGF+, and **APLN+ TipSI** represent distinct angiogenic programs (Nat Sci Rev 2024, CTM 2023, Nature 2024). **APLN+ TipSI cells predict anti-VEGF therapy response** and are associated with poor prognosis.

3. **Angiogenic trajectory**: Tumour angiogenesis initiates from venous ECs and extends toward arterial ECs through stages SI→SII→SIII (Nature 2024). Tip cells advance from APLN+ (SI) to Notch-high (SIII), while stalk cells transition from chemokine-high to TEK+/Tie2+-high.

4. **TEC vs normal EC**: Tumor endothelial cells show reduced SELE+ veins and increased CXCR4+/APLN+ tip cells (Nat Sci Rev 2024, Nature 2024). This remodeling is associated with therapy response.

5. **BBB is highly specialized**: BBB ECs express tight junctions and transporters but lack fenestration (PLVAP-) and inflammatory adhesion molecules (Nat Med 2024).

6. **Universal TEC markers**: ACKR1, PLVAP, and IGFBP3 are universally expressed across TECs (Nat Rev Cancer 2023). These can be used as pan-cancer endothelial markers.

7. **Lymphatic EC antigen presentation**: Lymphatic endothelial cells have two lineages — one for lymphangiogenesis and one for **antigen presentation** (HLA-DRA+, CD74+, CIITA+) (Nature 2024). This challenges the view of lymphatic ECs as purely structural.

8. **Neovascular ECs are immunosuppressive**: Newly formed tumour endothelial cells express CD274/PD-L1, PDCD1LG2/PD-L2, and IDO1, shaping an immunosuppressive microenvironment (Nature 2024).
