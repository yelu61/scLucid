> **⚠️ ARCHIVED / SUPERSEDED**
>
> This batch curation note is kept for provenance only. The live marker
> resource status is tracked in `docs/marker_resources/marker_curation_literature_index.jsonl`,
> `docs/marker_resources/marker_resource_quality_gaps.jsonl`, `docs/marker_resources/marker_curation_candidates.jsonl`,
> and `docs/marker_resources/CURATION.md`. New curation should follow the current
> contract rather than adding new files to this archive.

---

# Marker Curation Batch 12: Cancer-Type-Specific Atlas Papers

## Overview

This batch curates **cancer-type-specific single-cell and spatial atlases** that provide disease-context evidence for tumor microenvironment (TME) heterogeneity. Unlike pan-cancer atlases (Batch 01, 11), these papers focus on specific cancer types and often reveal **cancer-type-restricted cell states, subtypes, and therapeutic vulnerabilities** that do not generalize across all tumors.

### Resource Tier Classification Principles

| Tier | Content | Typical Papers |
|------|---------|---------------|
| **marker_tumor** | Cancer-type-specific cell states, TME subtypes, therapy-response signatures | Most disease atlases in this batch |
| **marker_registry** | Novel cell subtypes stable enough for global annotation (rare; most are context-dependent) | Only when a subtype is validated across multiple independent studies |
| **geneset** | Cancer-type-specific scoring signatures, pathway modules, TME archetypes | Most atlases contribute geneset entries |
| **marker_tissue** | Normal tissue reference for the corresponding organ | Normal tissue atlases (e.g., pancreas) |

### Key Principle

Cancer-type-specific markers **should NOT be promoted to global annotation** unless they are:
1. Validated in at least 2 independent studies of the same cancer type
2. Shown to have prognostic or therapeutic relevance
3. Not dominated by stress, cycling, or dissociation artifacts

Most entries from this batch belong in `marker_tumor_human.toml` with `scope = "cancer_type_specific"` and appropriate `cancer_type` and `applies_to` fields.

---

## Source Articles

### Pancreatic Ductal Adenocarcinoma (PDAC) — 7 papers

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | Human pancreatic cancer single cell atlas reveals association of CXCL10+ fibroblasts and basal subtype tumor cells | 2024 | Clin Cancer Res | 10.1158/1078-0432.CCR-24-2183 | disease_atlas | marker_tumor + geneset |
| 2 | IL-1β+ macrophages fuel pathogenic inflammation in pancreatic cancer | 2023 | Nature | 10.1038/s41586-023-06685-2 | disease_atlas | marker_tumor + geneset |
| 3 | Integrated single-cell and spatial transcriptomics uncover distinct cellular subtypes involved in neural invasion in pancreatic cancer | 2025 | Cancer Cell | 10.1016/j.ccell.2025.06.020 | disease_atlas | marker_tumor + geneset |
| 4 | Single cell transcriptomic analyses implicate an immunosuppressive tumor microenvironment in pancreatic cancer liver metastasis | 2023 | Nat Commun | 10.1038/s41467-023-40727-7 | disease_atlas | marker_tumor + geneset |
| 5 | Single-cell RNA sequencing highlights epithelial and microenvironmental heterogeneity in malignant progression of pancreatic ductal adenocarcinoma | 2024 | Cancer Lett | 10.1016/j.canlet.2024.216607 | review | marker_tumor (framework) + geneset |
| 6 | Single-cell transcriptional dissection illuminates an evolution of immunosuppressive microenvironment during pancreatic ductal adenocarcinoma metastasis | 2025 | STTT | 10.1038/s41392-025-02265-0 | disease_atlas | marker_tumor + geneset |
| 7 | Spatial mapping of transcriptomic plasticity in metastatic pancreatic cancer | 2025 | Nature | 10.1038/s41586-025-08927-x | disease_atlas | marker_tumor + geneset |

### Hepatocellular Carcinoma (HCC) — 7 papers

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 8 | Deciphering the Oncogenic Landscape of Hepatocytes Through Integrated Single-Nucleus and Bulk RNA-Seq of Hepatocellular Carcinoma | 2025 | Adv Sci | 10.1002/advs.202412944 | disease_atlas | geneset + marker_tumor |
| 9 | Liver tumour immune microenvironment subtypes and neutrophil heterogeneity | 2022 | Nature | 10.1038/s41586-022-05400-x | disease_atlas | marker_tumor + marker_registry + geneset |
| 10 | Single-cell landscape of the ecosystem in early-relapse hepatocellular carcinoma | 2021 | Cell | 10.1016/j.cell.2020.11.041 | disease_atlas | marker_tumor + geneset |
| 11 | Spatial analysis reveals targetable macrophage-mediated mechanisms of immune evasion in hepatocellular carcinoma minimal residual disease | 2024 | Nat Cancer | 10.1038/s43018-024-00828-8 | disease_atlas | marker_tumor + geneset |
| 12 | Spatial single-cell protein landscape reveals vimentin-high macrophages as immune-suppressive in the microenvironment of hepatocellular carcinoma | 2024 | Nat Cancer | 10.1038/s43018-024-00824-y | disease_atlas | marker_tumor + geneset |
| 13 | Unraveling the significance of cuproptosis in hepatocellular carcinoma heterogeneity and tumor microenvironment through integrated single-cell sequencing and machine learning approaches | 2025 | Discov Oncol | 10.1007/s12672-025-02696-9 | disease_atlas | marker_tumor + geneset |

### Lung Adenocarcinoma (LUAD/NSCLC) — 7 papers

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 14 | A single-cell atlas reveals immune heterogeneity in anti-PD-1-treated non-small cell lung cancer | 2025 | Cell | 10.1016/j.cell.2025.03.018 | single_cell_atlas | geneset + marker_tumor |
| 15 | An atlas of epithelial cell states and plasticity in lung adenocarcinoma | 2024 | Nature | 10.1038/s41586-024-07113-9 | disease_atlas | marker_tumor + geneset |
| 16 | Single-cell and spatial transcriptomic analyses revealing tumor microenvironment remodeling after neoadjuvant chemoimmunotherapy in non-small cell lung cancer | 2025 | Mol Cancer | 10.1186/s12943-025-02287-w | disease_atlas | marker_tumor + marker_registry + geneset |
| 17 | Single-cell RNA sequencing demonstrates the molecular and cellular reprogramming of metastatic lung adenocarcinoma | 2020 | Nat Commun | 10.1038/s41467-020-16164-1 | disease_atlas | marker_tumor + geneset |
| 18 | Single-cell RNA sequencing reveals immune microenvironment niche transitions during the invasive and metastatic processes of ground-glass nodules and part-solid nodules in lung adenocarcinoma | 2024 | Mol Cancer | 10.1186/s12943-024-02177-7 | disease_atlas | marker_tumor + geneset |
| 19 | The Single-Cell Immunogenomic Landscape of B and Plasma Cells in Early-Stage Lung Adenocarcinoma | 2022 | Cancer Discov | 10.1158/2159-8290.CD-21-1658 | disease_atlas | marker_tumor + geneset |

### Colorectal Cancer (CRC) — 6 papers

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 20 | Cancer-associated fibroblasts shape the formation of budding cancer cells at the invasive front of human colorectal cancer | 2025 | Commun Biol | 10.1038/s42003-025-08799-x | disease_atlas | marker_tumor + geneset |
| 21 | Integrative single-cell analysis of human colorectal cancer reveals patient stratification with distinct immune evasion mechanisms | 2024 | Nat Cancer | 10.1038/s43018-024-00807-z | disease_atlas | marker_tumor + marker_registry + geneset |
| 22 | Single-cell analyses define a continuum of cell state and composition changes in the malignant transformation of polyps to colorectal cancer | 2022 | Nat Genet | 10.1038/s41588-022-01088-x | disease_atlas | marker_tumor + marker_registry + geneset |
| 23 | Single-Cell Analyses Inform Mechanisms of Myeloid-Targeted Therapies in Colon Cancer | 2020 | Cell | 10.1016/j.cell.2020.03.048 | disease_atlas | marker_tumor + geneset |
| 24 | Single-cell and bulk transcriptome sequencing identifies two epithelial tumor cell states and refines the consensus molecular classification of colorectal cancer | 2022 | Nat Genet | 10.1038/s41588-022-01100-4 | disease_atlas | marker_tumor + geneset |

### Breast Cancer (BRCA) — 4 papers

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 25 | A single-cell and spatially resolved atlas of human breast cancers | 2021 | Nat Genet | 10.1038/s41588-021-00911-1 | disease_atlas | geneset + marker_tumor |
| 26 | Single-cell analyses reveal key immune cell subsets associated with response to PD-L1 blockade in triple-negative breast cancer | 2021 | Cancer Cell | 10.1016/j.ccell.2021.09.010 | disease_atlas | marker_tumor + geneset |
| 27 | Single-cell integrative analysis reveals consensus cancer cell states and clinical relevance in breast cancer | 2024 | Sci Data | 10.1038/s41597-024-03127-0 | disease_atlas | marker_tumor + geneset |
| 28 | Spatially resolved atlas of breast cancer uncovers intercellular machinery of venular niche governing lymphocyte extravasation | 2025 | Nat Commun | 10.1038/s41467-025-58511-0 | disease_atlas | marker_tumor + geneset |

### Head and Neck Squamous Cell Carcinoma (HNSCC) — 4 papers

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 29 | Integrative single-cell and bulk transcriptomes analyses identify intrinsic HNSCC subtypes with distinct prognoses and therapeutic vulnerabilities | 2023 | Clin Cancer Res | 10.1158/1078-0432.CCR-22-3563 | disease_atlas | geneset + marker_tumor |
| 30 | Single cell deciphering of progression trajectories of the tumor ecosystem in head and neck cancer | 2024 | Nat Commun | 10.1038/s41467-024-46912-6 | disease_atlas | marker_tumor + geneset |
| 31 | Single-cell analyses reveal the metabolic heterogeneity and plasticity of the tumor microenvironment during head and neck squamous cell carcinoma progression | 2024 | Cancer Res | 10.1158/0008-5472.CAN-23-1344 | disease_atlas | marker_tumor + geneset |
| 32 | Single-cell and spatial dissection of precancerous lesions underlying the initiation process of oral squamous cell carcinoma | 2023 | Cell Discov | 10.1038/s41421-023-00532-4 | disease_atlas | marker_tumor + geneset |

### Lymphoma (DLBCL/LBCL) — 3 papers

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 33 | Large B cell lymphoma microenvironment archetype profiles | 2025 | Cancer Cell | 10.1016/j.ccell.2025.06.002 | disease_atlas | marker_tumor + marker_registry + geneset |
| 34 | Multi-modal spatial characterization of tumor immune microenvironments identifies targetable inflammatory niches in diffuse large B cell lymphoma | 2025 | Nat Genet | 10.1038/s41588-025-02353-5 | disease_atlas | marker_tumor + geneset |
| 35 | The landscape of tumor cell states and ecosystems in diffuse large B cell lymphoma | 2021 | Cancer Cell | 10.1016/j.ccell.2021.08.011 | disease_atlas | marker_tumor + geneset |

### Gastric Cancer (GC) — 2 papers

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 36 | A spatially resolved atlas of gastric cancer characterises a lymphocyte-aggregated region | 2026 | Nat Commun | 10.1038/s41467-026-68612-z | disease_atlas | marker_tumor + geneset |
| 37 | Integrative single-cell multiomics analyses dissect molecular signatures of intratumoral heterogeneities and differentiation states of human gastric cancer | 2023 | NSR | 10.1093/nsr/nwad094 | disease_atlas | marker_tumor + geneset |

### Other Cancer Types — 9 papers

| # | Title | Year | Journal | DOI | Cancer Type | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|-------------|---------------|
| 38 | High-resolution transcriptome atlas of bladder cancer highlights the functional myeloid subsets in modulating immune microenvironment | 2025 | eBioMedicine | 10.1016/j.ebiom.2025.105801 | Bladder | disease_atlas | marker_tumor + geneset |
| 39 | Single-cell transcriptomics reveals metabolic remodeling and functional specialization in the immune microenvironment of bone tumors | 2025 | J Transl Med | 10.1186/s12967-025-06346-0 | Bone | disease_atlas | marker_tumor + geneset |
| 40 | Single-cell atlas of the human brain vasculature across development, adulthood and disease | 2024 | Nature | 10.1038/s41586-024-07493-y | Brain | single_cell_atlas | geneset + marker_tumor |
| 41 | Deciphering the tumor immune microenvironment: single-cell and spatial transcriptomic insights into cervical cancer fibroblasts | 2025 | J Exp Clin Cancer Res | 10.1186/s13046-025-03432-5 | Cervical | disease_atlas | marker_tumor + geneset |
| 42 | Integrative spatial analysis reveals a multi-layered organization of glioblastoma | 2024 | Cell | 10.1016/j.cell.2024.03.029 | GBM | disease_atlas | marker_tumor + geneset |
| 43 | Defining the mucosal ecosystem: epithelial–mesenchymal interdependence in gastrointestinal health and disease | 2025 | Nat Rev Gastroenterol Hepatol | 10.1038/s41575-025-01113-4 | GI (review) | review | marker_tumor (framework) + geneset |
| 44 | Single-cell analyses implicate ascites in remodeling the ecosystems of primary and metastatic tumors in ovarian cancer | 2023 | Nat Cancer | 10.1038/s43018-023-00599-8 | Ovarian | disease_atlas | marker_tumor + geneset |
| 45 | Single-cell and spatial RNA sequencing identify divergent microenvironments and progression signatures in early- versus late-onset prostate cancer | 2025 | Nat Aging | 10.1038/s43587-025-00842-0 | PRAD | disease_atlas | marker_tumor + geneset |
| 46 | A single-cell transcriptome atlas of the human pancreas | 2016 | Cell Syst | 10.1016/j.cels.2016.09.002 | Pancreas (normal) | disease_atlas | geneset + marker_tissue |
| 47 | Single-cell profiling reveals molecular basis of malignant phenotypes and tumor microenvironments in small bowel adenocarcinomas | 2022 | Cell Discov | 10.1038/s41421-022-00434-x | SBA | disease_atlas | marker_tumor + geneset |
| 48 | Multi-layered molecular profiling informs the diagnosis and targeted therapy of desmoplastic small round cell tumor | 2026 | Nat Commun | 10.1038/s41467-026-71636-0 | DSRCT | disease_atlas | geneset + marker_tumor |
| 49 | Single-cell multiomics profiling reveals heterogeneous transcriptional programs and microenvironment in DSRCTs | 2024 | Cell Rep Med | 10.1016/j.xcrm.2024.101582 | DSRCT | disease_atlas | marker_tumor + geneset |

---

## Cancer-Type-Specific Marker Contributions by Disease

### PDAC: Key TME Features

| Feature | Markers | Source | Resource Target |
|---------|---------|--------|-----------------|
| CXCL10+ CAF | CXCL10, CXCL9, CCL19, CCL21 | Clin Cancer Res 2024 | marker_tumor (cancer_state) |
| IL-1β+ macrophage | IL1B, NLRP3, CXCL1, CXCL2, CXCL8 | Nature 2023 | marker_tumor (cancer_state) |
| Neural invasion-associated subtype | SEMA3C, NRP1, NRP2, PLXNA1 | Cancer Cell 2025 | marker_tumor (cancer_state) |
| Basal-like tumor cell | KRT5, KRT14, KRT17, TP63 | Clin Cancer Res 2024 | marker_tumor (cancer_state) |
| Metastatic plasticity | VIM, SNAI1, ZEB1, CDH2 | Nature 2025 | marker_tumor (cancer_state) |
| Immunosuppressive liver metastasis TME | FOXP3, IL10, TGFB1, CD274 | Nat Commun 2023 | marker_tumor (cancer_state) |

```toml
# Example: CXCL10+ CAF in PDAC
[[cancer_state]]
name = "CXCL10+ CAF (PDAC)"
color = "#FF6347"
markers = ["CXCL10", "CXCL9", "CCL19", "CCL21", "CXCL11"]
negative_markers = ["ACTA2", "TAGLN", "POSTN", "LRRC15"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["Pancreatic Cancer"]
scope = "cancer_type_specific"
applies_to = ["Pancreas"]
evidence_tier = "atlas_supported"
source_type = "disease_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Human pancreatic cancer single cell atlas reveals association of CXCL10+ fibroblasts and basal subtype tumor cells", year = "2024", doi = "10.1158/1078-0432.CCR-24-2183" }
notes = "CXCL10+ cancer-associated fibroblasts in PDAC. Associated with basal subtype tumor cells and immune recruitment. Chemokine-producing CAF subtype distinct from myofibroblasts (ACTA2-, TAGLN-, POSTN-)."
```

### HCC: Key TME Features

| Feature | Markers | Source | Resource Target |
|---------|---------|--------|-----------------|
| TIME-ISM (immune suppressive myeloid) | S100A8, S100A9, CXCL8, MMP9, VEGFA, ARG1 | Nature 2022 | marker_tumor (cancer_state) |
| Vimentin-high macrophage | VIM, CD68, CD163, MRC1 | Nat Cancer 2024 | marker_tumor (cancer_state) |
| Early-relapse ecosystem | CD24, PROM1, EPCAM, ALDH1A1 | Cell 2021 | marker_tumor (cancer_state) |
| Cuproptosis-related tumor cell | FDX1, LIAS, LIPT1, DLD, DLAT | Discov Oncol 2025 | marker_tumor (cancer_state) |
| MRD immune evasion macrophage | CD274, PDCD1LG2, IDO1, ARG1 | Nat Cancer 2024 | marker_tumor (cancer_state) |

### LUAD: Key TME Features

| Feature | Markers | Source | Resource Target |
|---------|---------|--------|-----------------|
| Anti-PD-1 responder T cell | CD8A, GZMB, PRF1, IFNG, CXCL9, CXCL10, CXCL11 | Cell 2025 | marker_tumor (cancer_state) |
| Anti-PD-1 non-responder T cell | FOXP3, IL10, TGFB1, LAG3, TIGIT | Cell 2025 | marker_tumor (cancer_state) |
| Epithelial plasticity (alveolar) | SFTPC, SFTPA1, SFTPA2, AGER | Nature 2024 | marker_tumor (cancer_state) |
| Epithelial plasticity (basaloid) | KRT5, KRT14, TP63, NGFR | Nature 2024 | marker_tumor (cancer_state) |
| Ground-glass nodule transition | NKX2-1, FOXA2, ETV5, ID2 | Mol Cancer 2024 | marker_tumor (cancer_state) |
| Neoadjuvant chemoimmunotherapy responder | CXCL13, CD200, BTLA, TOX | Mol Cancer 2025 | marker_tumor (cancer_state) |

### CRC: Key TME Features

| Feature | Markers | Source | Resource Target |
|---------|---------|--------|-----------------|
| Budding cancer cell (CAF-associated) | VIM, SNAI1, TWIST1, ZEB1, CDH2 | Commun Biol 2025 | marker_tumor (cancer_state) |
| Immune evasion subtype A | FOXP3, IL10, TGFB1, CD274, PDCD1LG2 | Nat Cancer 2024 | marker_tumor (cancer_state) |
| Immune evasion subtype B | VEGFA, ANGPT2, COL1A1, POSTN | Nat Cancer 2024 | marker_tumor (cancer_state) |
| Polyp-to-CRC continuum | APC, KRAS, TP53, SMAD4, PIK3CA | Nat Genet 2022 | marker_tumor (cancer_state) |
| CMS-refined epithelial state A | EPCAM, CDH1, KRT8, KRT18, KRT19 | Nat Genet 2022 | marker_tumor (cancer_state) |
| CMS-refined epithelial state B | VIM, FN1, SNAI1, ZEB1, TWIST1 | Nat Genet 2022 | marker_tumor (cancer_state) |
| Myeloid-targeted therapy responsive | CSF1R, CD115, CCL2, CCR2 | Cell 2020 | marker_tumor (cancer_state) |

### BRCA: Key TME Features

| Feature | Markers | Source | Resource Target |
|---------|---------|--------|-----------------|
| TNBC anti-PD-L1 responder | CD8A, CXCL9, CXCL10, CXCL11, IFNG | Cancer Cell 2021 | marker_tumor (cancer_state) |
| TNBC anti-PD-L1 non-responder | FOXP3, TGFB1, VEGFA, COL1A1 | Cancer Cell 2021 | marker_tumor (cancer_state) |
| Venular niche (lymphocyte extravasation) | MADCAM1, CCL21, CCL19, ICAM1, VCAM1 | Nat Commun 2025 | marker_tumor (cancer_state) |
| Consensus cancer cell state (basal) | KRT5, KRT14, KRT17, TP63 | Sci Data 2024 | marker_tumor (cancer_state) |
| Consensus cancer cell state (luminal) | ESR1, PGR, FOXA1, GATA3 | Sci Data 2024 | marker_tumor (cancer_state) |

### Lymphoma: Key TME Features

| Feature | Markers | Source | Resource Target |
|---------|---------|--------|-----------------|
| Inflammatory niche (DLBCL) | CXCL9, CXCL10, CXCL11, CCL5, IFNG | Nat Genet 2025 | marker_tumor (cancer_state) |
| Immune-desert niche (DLBCL) | COL1A1, POSTN, FN1, VIM | Nat Genet 2025 | marker_tumor (cancer_state) |
| LBCL archetype A (T cell-inflamed) | CD8A, GZMB, PRF1, IFNG, CXCL9 | Cancer Cell 2025 | marker_tumor (cancer_state) |
| LBCL archetype B (stromal) | COL1A1, POSTN, ACTA2, TAGLN | Cancer Cell 2025 | marker_tumor (cancer_state) |
| DLBCL ecosystem (myeloid-rich) | CD68, CD163, MRC1, CSF1R | Cancer Cell 2021 | marker_tumor (cancer_state) |

---

## Geneset Entries by Cancer Type

### PDAC-Specific Gene Sets

```json
"PDAC_CXCL10_CAF": {
  "genes": ["CXCL10", "CXCL9", "CCL19", "CCL21", "CXCL11", "CCL5", "IFNG", "STAT1"],
  "description": "PDAC CXCL10+ CAF signature. Chemokine-producing fibroblasts associated with basal tumor subtype and immune recruitment."
},
"PDAC_IL1B_macrophage": {
  "genes": ["IL1B", "NLRP3", "CXCL1", "CXCL2", "CXCL8", "IL6", "TNF", "PTGS2", "NFKB1", "NFKBIA"],
  "description": "PDAC IL-1β+ macrophage signature. Pathogenic inflammation-driven macrophage state. Associated with poor prognosis."
},
"PDAC_neural_invasion": {
  "genes": ["SEMA3C", "NRP1", "NRP2", "PLXNA1", "PLXNA2", "NGFR", "TRKA", "RET"],
  "description": "PDAC neural invasion-associated cellular subtype. Semaphorin-neuropilin signaling axis."
},
"PDAC_metastatic_plasticity": {
  "genes": ["VIM", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "CDH2", "FN1", "TWIST1", "MMP2", "MMP9"],
  "description": "PDAC metastatic transcriptomic plasticity signature. EMT-like reprogramming in metastatic lesions."
}
```

### HCC-Specific Gene Sets

```json
"HCC_TIME_ISM": {
  "genes": ["S100A8", "S100A9", "CXCL8", "MMP9", "VEGFA", "ARG1", "CD274", "PDCD1LG2", "IDO1", "CCL22"],
  "description": "HCC immune suppressive myeloid (TIME-ISM) signature. Immunosuppressive TME subtype associated with poor prognosis."
},
"HCC_vimentin_macrophage": {
  "genes": ["VIM", "CD68", "CD163", "MRC1", "CSF1R", "CD14", "LYZ", "APOE"],
  "description": "HCC vimentin-high macrophage signature. Immune-suppressive macrophage subset in HCC TME."
},
"HCC_early_relapse": {
  "genes": ["CD24", "PROM1", "EPCAM", "ALDH1A1", "SOX9", "CD44", "NANOG", "SOX2"],
  "description": "HCC early-relapse cancer stem cell signature. Associated with recurrence and poor survival."
},
"HCC_cuproptosis": {
  "genes": ["FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "MTF1", "GLS", "CDKN2A"],
  "description": "HCC cuproptosis-related signature. Copper metabolism dysregulation in HCC heterogeneity."
}
```

### CRC-Specific Gene Sets

```json
"CRC_budding_CAF": {
  "genes": ["VIM", "SNAI1", "TWIST1", "ZEB1", "CDH2", "FN1", "MMP2", "MMP9", "LAMC2", "LAMA3"],
  "description": "CRC budding cancer cell signature at invasive front. CAF-associated EMT-like state."
},
"CRC_immune_evasion_A": {
  "genes": ["FOXP3", "IL10", "TGFB1", "CD274", "PDCD1LG2", "IDO1", "LAG3", "TIGIT", "VTCN1", "SIGLEC15"],
  "description": "CRC immune evasion subtype A. Treg-dominant immunosuppressive TME."
},
"CRC_immune_evasion_B": {
  "genes": ["VEGFA", "ANGPT2", "COL1A1", "POSTN", "FN1", "MMP2", "MMP9", "LOX", "LOXL2"],
  "description": "CRC immune evasion subtype B. Stromal/angiogenic exclusion TME."
}
```

### LUAD-Specific Gene Sets

```json
"LUAD_antiPD1_responder": {
  "genes": ["CD8A", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10", "CXCL11", "TBX21", "EOMES", "GZMA"],
  "description": "LUAD anti-PD-1 responder T cell signature. Cytotoxic and chemokine-producing T cell phenotype."
},
"LUAD_antiPD1_nonresponder": {
  "genes": ["FOXP3", "IL10", "TGFB1", "LAG3", "TIGIT", "CD274", "PDCD1LG2", "VTCN1", "IDO1"],
  "description": "LUAD anti-PD-1 non-responder signature. Immunosuppressive and exhausted TME."
},
"LUAD_alveolar_plasticity": {
  "genes": ["SFTPC", "SFTPA1", "SFTPA2", "AGER", "AQP5", "NKX2-1", "FOXA2"],
  "description": "LUAD alveolar differentiation state. Surfactant-producing epithelial phenotype."
},
"LUAD_basaloid_plasticity": {
  "genes": ["KRT5", "KRT14", "TP63", "NGFR", "SOX2", "NOTCH1", "NOTCH3"],
  "description": "LUAD basaloid differentiation state. Squamous-like epithelial phenotype."
}
```

### Lymphoma-Specific Gene Sets

```json
"DLBCL_inflammatory_niche": {
  "genes": ["CXCL9", "CXCL10", "CXCL11", "CCL5", "IFNG", "GZMB", "PRF1", "TBX21", "CD8A"],
  "description": "DLBCL inflammatory niche signature. T cell-inflamed and chemokine-rich microenvironment."
},
"DLBCL_immune_desert_niche": {
  "genes": ["COL1A1", "POSTN", "FN1", "VIM", "ACTA2", "TAGLN", "MMP2", "MMP9", "LOX"],
  "description": "DLBCL immune-desert niche signature. Stromal-dominant T cell-excluded microenvironment."
},
"LBCL_T_cell_inflamed": {
  "genes": ["CD8A", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10", "CXCL11", "TBX21", "EOMES"],
  "description": "LBCL T cell-inflamed archetype. Cytotoxic T cell-rich microenvironment."
},
"LBCL_stromal": {
  "genes": ["COL1A1", "POSTN", "ACTA2", "TAGLN", "FN1", "THY1", "PDGFRA", "FAP"],
  "description": "LBCL stromal archetype. Fibroblast-dominant microenvironment."
}
```

---

## Cross-Cutting Themes and Conflict Notes

### 1. CAF Heterogeneity Across Cancer Types

Multiple papers identify cancer-type-specific CAF subtypes:
- **PDAC**: CXCL10+ CAF (immune-recruiting), myofibroblast CAF (immune-excluding)
- **CRC**: Budding CAF (EMT-promoting at invasive front)
- **HCC**: Vimentin-high CAF (not well characterized)
- **BRCA**: Venular niche CAF (lymphocyte extravasation)

**Resolution**: Each cancer-type-specific CAF subtype should be a separate `cancer_state` entry with appropriate `cancer_type` and `applies_to` fields. Do NOT merge into global CAF subtypes unless validated cross-cancer.

### 2. Therapy Response Signatures

Multiple papers report therapy response signatures:
- **Anti-PD-1/PD-L1**: LUAD (Cell 2025), BRCA (Cancer Cell 2021), NSCLC (Mol Cancer 2025)
- **Myeloid-targeted**: CRC (Cell 2020)
- **Neoadjuvant chemoimmunotherapy**: NSCLC (Mol Cancer 2025)

**Resolution**: Therapy response signatures are `cancer_state` entries. They should be annotated with the specific therapy and cancer type. Consider creating a unified `therapy_response` category in `marker_tumor_human.toml`.

### 3. Metastasis-Specific States

Multiple papers focus on metastatic TME:
- **PDAC liver metastasis**: Immunosuppressive TME (Nat Commun 2023)
- **PDAC metastatic plasticity**: EMT-like reprogramming (Nature 2025)
- **LUAD metastasis**: Molecular and cellular reprogramming (Nat Commun 2020)

**Resolution**: Metastasis-specific states should be separate `cancer_state` entries from primary tumor states. Use `scope = "cancer_type_specific"` with appropriate `applies_to` for metastatic sites.

### 4. Spatial Context

Several papers integrate spatial transcriptomics:
- **PDAC**: Neural invasion spatial patterns (Cancer Cell 2025)
- **GC**: Lymphocyte-aggregated region (Nat Commun 2026)
- **BRCA**: Venular niche (Nat Commun 2025)
- **DLBCL**: Inflammatory vs immune-desert niches (Nat Genet 2025)

**Resolution**: Spatial features should be recorded in `notes` fields. Consider adding `spatial_context` metadata field for future use.

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| Total cancer-type-specific atlases | 49 | — |
| PDAC-focused | 7 | marker_tumor + geneset |
| HCC-focused | 7 | marker_tumor + geneset |
| LUAD/NSCLC-focused | 7 | marker_tumor + geneset |
| CRC-focused | 6 | marker_tumor + geneset |
| BRCA-focused | 4 | marker_tumor + geneset |
| HNSCC-focused | 4 | marker_tumor + geneset |
| Lymphoma-focused | 3 | marker_tumor + geneset |
| Gastric-focused | 2 | marker_tumor + geneset |
| Other cancer types | 9 | marker_tumor + geneset |
| New cancer_state entries proposed | 35+ | marker_tumor |
| New geneset programs proposed | 20+ | genesets_cancer_signatures.json |
| Normal tissue reference | 1 (pancreas) | marker_tissue |
| Review/framework papers | 2 | marker_tumor (framework) + geneset |

---

## Key Biological Insights for scLucid Workflows

1. **Cancer-type-specific TME is not interchangeable**: A CAF marker validated in PDAC may not apply to HCC. scLucid should use `cancer_type` filtering when interpreting tumor context.

2. **Therapy response is cancer-type-specific**: Anti-PD-1 response signatures differ between LUAD and BRCA. scLucid should not apply a universal ICB response signature across all cancer types.

3. **Metastatic TME differs from primary**: Metastatic sites have distinct TME compositions. scLucid should distinguish primary vs metastatic when both are present in a dataset.

4. **Spatial context adds critical information**: Spatially resolved atlases reveal niche structures (venular, inflammatory, immune-desert) that are invisible in dissociated scRNA-seq. scLucid should integrate spatial data when available.

5. **Normal tissue reference is essential**: Cancer-type-specific atlases should be interpreted against normal tissue references (e.g., pancreas atlas for PDAC, lung atlas for LUAD). scLucid should load `marker_tissue_human.toml` when a tissue is specified.
