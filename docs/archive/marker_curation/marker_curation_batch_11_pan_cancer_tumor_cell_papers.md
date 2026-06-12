# Marker Curation Batch 11: Pan-Cancer Tumor Cell Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | The Curated Cancer Cell Atlas provides a comprehensive characterization of tumors at single-cell resolution | 2025 | Nature Cancer | 10.1038/s43018-025-00957-8 | pan_cancer_atlas | geneset + marker_tumor |
| 2 | A gene set enrichment analysis for cancer hallmarks | 2025 | J Pharm Anal | 10.1016/j.jpha.2024.101065 | computational | geneset |
| 3 | Single-cell spatial transcriptomics unravels cell states and ecosystems associated with clinical response to immunotherapy | 2025 | JITC | 10.1136/jitc-2024-011308 | single_cell_atlas | marker_tumor + geneset |
| 4 | Cancer cell states: lessons from ten years of single-cell RNA-sequencing of human tumors | 2024 | Cancer Cell | 10.1016/j.ccell.2024.08.005 | review | marker_tumor (framework) |
| 5 | Embracing cancer complexity: hallmarks of systemic disease | 2024 | Cell | 10.1016/j.cell.2024.02.009 | review | geneset (framework) |
| 6 | How chemokines organize the tumour microenvironment | 2024 | Nat Rev Cancer | 10.1038/s41568-023-00635-w | review | marker_tumor (validation) |
| 7 | Single cell multi-omics reveal intra-cell-line heterogeneity across human cancer cell lines | 2023 | Nat Commun | 10.1038/s41467-023-43991-9 | single_cell_atlas | geneset |
| 8 | Precise identification of cell states altered in disease using healthy single-cell references | 2023 | Nat Genet | 10.1038/s41588-023-01523-7 | computational | marker_tumor |
| 9 | Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours | 2023 | Nature | 10.1038/s41586-023-06130-4 | pan_cancer_atlas | geneset + marker_tumor |
| 10 | Mechanisms driving the immunoregulatory function of cancer cells | 2023 | Nat Rev Cancer | 10.1038/s41568-022-00544-4 | review | marker_tumor (validation) |
| 11 | Metabolic programming and immune suppression in the tumor microenvironment | 2023 | Cancer Cell | 10.1016/j.ccell.2023.01.009 | review | geneset |
| 12 | Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment | 2022 | Nat Genet | 10.1038/s41588-022-01141-9 | single_cell_atlas | geneset + marker_tumor |
| 13 | Atlas of clinically distinct cell states and ecosystems across human solid tumors | 2021 | Cell | 10.1016/j.cell.2021.09.014 | computational | marker_tumor + geneset |
| 14 | Hallmarks of response, resistance, and toxicity to immune checkpoint blockade | 2021 | Cell | 10.1016/j.cell.2021.09.020 | review | geneset + marker_tumor |
| 15 | Pan-cancer single-cell RNA-seq identifies recurring programs of cellular heterogeneity | 2020 | Nat Genet | 10.1038/s41588-020-00726-6 | pan_cancer_atlas | geneset + marker_tumor |
| 16 | Deciphering Human Tumor Biology by Single-Cell Expression Profiling | 2019 | Annu Rev Cancer Biol | 10.1146/annurev-cancerbio-030518-055609 | review | marker_tumor (framework) |
| 17 | CancerSEA: a cancer single-cell state atlas | 2019 | NAR | 10.1093/nar/gky939 | database | geneset |
| 18 | Single cell atlas reveals multilayered metabolic heterogeneity across tumour types | 2024 | eBioMedicine | 10.1016/j.ebiom.2024.105389 | pan_cancer_atlas | marker_tumor + geneset |

---

## Overview: Pan-Cancer Tumor Cell Context

This batch focuses on **cancer cell states, meta-programs, and cancer hallmarks** that recur across tumor types. Unlike previous batches that curated specific lineage markers, this batch provides:

1. **Meta-programs (MPs)**: Recurring gene expression programs in malignant cells
2. **Cancer hallmarks gene sets**: Functional signatures for tumor characterization
3. **Tumor cell states**: Stress, hypoxia, EMT, proliferation, etc.
4. **ICB response/resistance/toxicity signatures**

### Resource Tier Classification

| Tier | Content | Applicable Papers |
|------|---------|-------------------|
| **geneset** | MPs, hallmark signatures, functional programs | 1, 2, 5, 9, 12, 14, 15, 17 |
| **marker_tumor** | Cancer state signatures | 1, 3, 8, 9, 12, 13, 14, 15 |
| **marker_registry** | Not applicable (no new identity markers) | — |

**Key principle**: Tumor cell states are **functional programs**, not identity markers. They belong in `genesets_cancer_signatures.json` or `marker_tumor_human.toml`, NOT in `marker_registry_human.toml`.

---

## Paper-by-Paper Curation

### Article 15: Nature Genetics 2020 — Recurring Programs in CCLE (Tirosh Lab)

**Dataset**: CCLE cell lines, multiplexed scRNA-seq
**Key finding**: 12 **meta-programs (MPs)** in cancer cell lines

#### Meta-Programs Defined

| MP | Key Genes | Biological Interpretation |
|----|-----------|--------------------------|
| MP1 | cell cycle | Proliferation |
| MP2 | pEMT | Partial EMT |
| MP3 | stress | Stress response |
| MP4 | hypoxia | Hypoxia response |
| MP5 | IFN | Interferon response |
| MP6 | protein secretion | Secretory phenotype |
| MP7 | MYC targets | MYC activation |
| MP8 | oxidative phosphorylation | Metabolic |
| MP9 | unfolded protein response | ER stress |
| MP10 | TNF-alpha signaling | Inflammation |
| MP11 | DNA repair | Genomic maintenance |
| MP12 | glycolysis | Warburg effect |

#### Geneset Entries

```json
"MP_cell_cycle": {
  "genes": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "CDK1", "AURKB", "BIRC5", "UBE2C", "CENPF"],
  "description": "Meta-program 1: Cell cycle. Recurring proliferation program across cancer cell lines."
},
"MP_pEMT": {
  "genes": ["VIM", "S100A4", "S100A6", "CDH2", "FN1", "TGFBI", "LAMC2", "LAMA3", "SPARC", "MMP9"],
  "description": "Meta-program 2: Partial EMT. Recurring epithelial-mesenchymal transition program."
},
"MP_stress": {
  "genes": ["HSPA1A", "HSPB1", "HSPA6", "HSP90AA1", "DNAJB1", "HSPD1", "FOS", "JUN", "ATF3", "DDIT3"],
  "description": "Meta-program 3: Stress response. Recurring cellular stress program."
},
"MP_hypoxia": {
  "genes": ["VEGFA", "CA9", "BNIP3", "PGK1", "ENO1", "LDHA", "SLC2A1", "ALDOA", "P4HA1", "P4HA2"],
  "description": "Meta-program 4: Hypoxia response. Recurring hypoxia-induced program."
},
"MP_IFN": {
  "genes": ["ISG15", "MX1", "OAS1", "OAS2", "OAS3", "IFIT1", "IFIT3", "IRF7", "STAT1", "IFI27"],
  "description": "Meta-program 5: Interferon response. Recurring IFN-stimulated gene program."
},
"MP_secretory": {
  "genes": ["SEC61B", "SEC61G", "SEC11C", "SPCS1", "SPCS2", "SPCS3", "SRP14", "SRP9", "RPN1", "RPN2"],
  "description": "Meta-program 6: Protein secretion. Recurring secretory pathway program."
},
"MP_MYC": {
  "genes": ["MYC", "NCL", "NPM1", "EIF4A1", "EIF4G2", "NOP56", "NOP58", "FBL", "DDX21", "DKC1"],
  "description": "Meta-program 7: MYC targets. Recurring MYC-driven ribosomal biogenesis program."
},
"MP_oxphos": {
  "genes": ["NDUFA4", "COX4I1", "COX5B", "ATP5F1A", "ATP5F1B", "UQCRB", "SDHB", "CYC1", "ATP5MG", "ATP5ME"],
  "description": "Meta-program 8: Oxidative phosphorylation. Recurring mitochondrial respiration program."
},
"MP_UPR": {
  "genes": ["HSPA5", "HSP90B1", "PDIA3", "PDIA4", "PDIA6", "CALR", "CANX", "HYOU1", "HERPUD1", "DNAJB11"],
  "description": "Meta-program 9: Unfolded protein response. Recurring ER stress program."
},
"MP_TNF": {
  "genes": ["NFKBIA", "NFKB1", "RELB", "TNFAIP3", "BCL2A1", "CCL2", "CXCL1", "CXCL2", "ICAM1", "VCAM1"],
  "description": "Meta-program 10: TNF-alpha signaling. Recurring inflammatory NF-kB program."
},
"MP_DNA_repair": {
  "genes": ["BRCA1", "BRCA2", "RAD51", "CHEK1", "CHEK2", "ATM", "ATR", "PARP1", "FEN1", "PCNA"],
  "description": "Meta-program 11: DNA repair. Recurring genomic maintenance program."
},
"MP_glycolysis": {
  "genes": ["HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "PGAM1", "TPI1"],
  "description": "Meta-program 12: Glycolysis. Recurring Warburg effect program."
}
```

---

### Article 9: Nature 2023 — Hallmarks of Transcriptional ITH (Tirosh Lab)

**Dataset**: 77 studies, 1,456 samples, 24 cancer types, 2.6M cells
**Key finding**: **149 MPs across 8 cell types**; MPs explain large fraction of ITH

#### New Geneset Entries

```json
"MP_malignant_stemness": {
  "genes": ["SOX2", "NANOG", "PROM1", "NES", "CD44", "ALDH1A1", "SOX4", "SOX9", "ID1", "ID2"],
  "description": "Malignant cell meta-program: Stemness. Cancer stem cell-like state recurring across tumors."
},
"MP_malignant_differentiation": {
  "genes": ["KRT5", "KRT14", "KRT8", "KRT18", "EPCAM", "CDH1", "CLDN4", "MUC1", "KRT19", "DSP"],
  "description": "Malignant cell meta-program: Differentiation. Epithelial differentiation state recurring across carcinomas."
},
"MP_malignant_invasion": {
  "genes": ["VIM", "FN1", "SPARC", "S100A4", "S100A6", "MMP2", "MMP9", "TIMP1", "LAMC2", "LAMA3"],
  "description": "Malignant cell meta-program: Invasion. Invasive phenotype recurring across tumors."
},
"MP_immune_regulatory": {
  "genes": ["TGFB1", "TGFB2", "TGFB3", "IL10", "VEGFA", "IDO1", "ARG1", "CD274", "PDCD1LG2", "CCL22"],
  "description": "Malignant cell meta-program: Immune regulation. Immunosuppressive phenotype recurring across tumors."
}
```

---

### Article 12: Nature Genetics 2022 — Cancer Cell States Recur (Yanai/Tirosh Lab)

**Dataset**: 19 primary tumors, 9 cancer types
**Key finding**: Malignant cell states interact with specific TME cell types

| Cancer Cell State | TME Interaction | Clinical Relevance |
|------------------|----------------|-------------------|
| pEMT | CAFs at leading edge | Invasion |
| Stress | Macrophages | Immune evasion |
| Hypoxia | Endothelial cells | Angiogenesis |
| IFN response | T cells | Immune activation |
| Stemness | Perivascular niche | Therapy resistance |

#### Tumor Context (`marker_tumor`)

**pEMT Cancer Cell State**

```toml
[[cancer_state]]
name = "Partial EMT cancer cell state"
color = "#FF6347"
markers = ["VIM", "S100A4", "S100A6", "CDH2", "FN1", "TGFBI", "LAMC2", "LAMA3", "SPARC", "MMP9"]
negative_markers = ["CDH1", "EPCAM", "CLDN4", "OCLN", "TJP1"]

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
source = { title = "Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment", year = "2022", doi = "10.1038/s41588-022-01141-9" }
notes = "Partial EMT cancer cell state. Recurring across 9 cancer types. Located at tumor leading edge. Interacts with CAFs to mediate invasion. Negative epithelial markers (CDH1, EPCAM, CLDN4, OCLN, TJP1)."
```

**Stemness Cancer Cell State**

```toml
[[cancer_state]]
name = "Stemness cancer cell state"
color = "#9400D3"
markers = ["SOX2", "NANOG", "PROM1", "NES", "CD44", "ALDH1A1", "SOX4", "SOX9", "ID1", "ID2"]
negative_markers = ["CD24", "MUC1", "KRT8", "KRT18", "KRT19"]

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
source = { title = "Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment", year = "2022", doi = "10.1038/s41588-022-01141-9" }
notes = "Stemness cancer cell state. Cancer stem cell-like. Associated with perivascular niche. Therapy resistance. Recurring across tumor types. Negative differentiation markers (CD24, MUC1, KRT8, KRT18, KRT19)."
```

---

### Article 14: Cell 2021 — Hallmarks of ICB Response, Resistance, Toxicity

**Key finding**: Framework for ICB outcomes

#### Geneset Entries

```json
"ICB_response_signature": {
  "genes": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10", "CXCL11", "TBX21", "EOMES", "PDCD1", "CTLA4", "LAG3", "TIGIT"],
  "description": "Immune checkpoint blockade (ICB) response signature. T cell-inflamed phenotype."
},
"ICB_resistance_signature": {
  "genes": ["FOXP3", "IL10", "TGFB1", "VEGFA", "ARG1", "IDO1", "CD274", "PDCD1LG2", "SIGLEC15", "VTCN1", "WNT5A", "AXL", "TGFB2"],
  "description = "Immune checkpoint blockade (ICB) resistance signature. Immunosuppressive and T cell-excluded phenotype."
},
"ICB_toxicity_signature": {
  "genes": ["CTLA4", "PDCD1", "LAG3", "TIGIT", "CD28", "ICOS", "CD40LG", "TNFRSF9", "TNFRSF4", "TNFRSF18"],
  "description": "Immune checkpoint blockade (ICB) toxicity signature. Co-stimulatory and checkpoint molecule overexpression associated with irAEs."
}
```

---

### Article 5: Cell 2024 — Embracing Cancer Complexity (Hanahan)

**Key finding**: Updated cancer hallmarks incorporating systemic disease

#### Updated Hallmarks Gene Sets

```json
"Hallmarks_angiogenesis": {
  "genes": ["VEGFA", "VEGFB", "VEGFC", "KDR", "FLT1", "FLT4", "ANGPT1", "ANGPT2", "TEK", "PGF", "NRP1", "NRP2", "EFNB2", "DLL4"],
  "description": "Cancer hallmark: Inducing angiogenesis. Updated 2024."
},
"Hallmarks_immune_evasion": {
  "genes": ["CD274", "PDCD1LG2", "CTLA4", "IDO1", "FOXP3", "IL10", "TGFB1", "SIGLEC15", "VTCN1", "HAVCR2", "LAG3", "TIGIT"],
  "description": "Cancer hallmark: Avoiding immune destruction. Updated 2024 with new checkpoint molecules."
},
"Hallmarks_metabolic_reprogramming": {
  "genes": ["HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "SLC2A1", "MYC", "HIF1A", "EPAS1", "PKM2", "PDK1"],
  "description": "Cancer hallmark: Reprogramming cellular metabolism. Warburg effect and beyond."
},
"Hallmarks_inflammation": {
  "genes": ["IL6", "IL1B", "TNF", "CXCL8", "CXCL1", "CXCL2", "CCL2", "CCL5", "PTGS2", "NFKB1", "RELA", "STAT3", "JAK2"],
  "description": "Cancer hallmark: Tumor-promoting inflammation. Updated 2024."
},
"Hallmarks_genomic_instability": {
  "genes": ["TP53", "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2", "RAD51", "PARP1", "MLH1", "MSH2", "MSH6", "PMS2"],
  "description": "Cancer hallmark: Genome instability and mutation. DNA repair deficiency."
}
```

---

### Article 11: Cancer Cell 2023 — Metabolic Programming and Immune Suppression

**Key finding**: Metabolic programs in TME drive immune suppression

#### Geneset Entries

```json
"TME_metabolic_immune_suppression": {
  "genes": ["IDO1", "TDO2", "ARG1", "ARG2", "NT5E", "CD39", "CD73", "LGALS1", "LGALS3", "LGALS9", "CEACAM1", "VISTA", "B7H3", "B7H4"],
  "description": "TME metabolic immune suppression program. Tryptophan degradation, arginine depletion, adenosine generation, galectin-mediated suppression."
},
"TME_nutrient_depletion": {
  "genes": ["SLC7A11", "SLC1A5", "SLC38A2", "SLC7A5", "GLS", "GLUL", "ASCT2", "LAT1", "SNAT2"],
  "description": "TME nutrient depletion program. Glutamine, tryptophan, arginine deprivation."
}
```

---

### Article 13: Cell 2021 — EcoTyper (Cell States and Ecosystems)

**Dataset**: 16 types of human carcinoma
**Key finding**: 69 cell states, 10 multicellular communities (ecotypes)

| Ecotype | Composition | Prognosis |
|---------|-------------|-----------|
| EC1-3 | Myeloid + stromal | Adverse |
| EC4 | Normal tissue-enriched | Favorable |
| EC5-6 | Early cancer development | Intermediate |
| EC7-10 | Other distinct patterns | Variable |

#### Tumor Context (`marker_tumor`)

**Ecotype Adverse (Myeloid-Stromal)**

```toml
[[cancer_state]]
name = "Ecotype adverse (myeloid-stromal)"
color = "#8B0000"
markers = ["CD68", "CD163", "MRC1", "CD274", "COL1A1", "POSTN", "FN1", "ACTA2", "TAGLN"]
negative_markers = ["CD8A", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10"]

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
source = { title = "Atlas of clinically distinct cell states and ecosystems across human solid tumors", year = "2021", doi = "10.1016/j.cell.2021.09.014" }
notes = "EcoTyper ecotype: adverse prognosis. Myeloid and stromal enrichment. Macrophage markers (CD68, CD163, MRC1, CD274) and CAF markers (COL1A1, POSTN, FN1, ACTA2, TAGLN). Negative cytotoxic T cell markers (CD8A, GZMB, PRF1, IFNG, CXCL9, CXCL10)."
```

---

### Article 17: NAR 2019 — CancerSEA Database

**Resource**: Cancer Single-cell State Atlas
**Key feature**: 14 functional states for 41,900 cancer cells from 25 cancer types

| State | Description |
|-------|-------------|
| Angiogenesis | VEGF signaling |
| Apoptosis | Programmed cell death |
| Cell cycle | Proliferation |
| Differentiation | Cell fate determination |
| DNA damage | Genomic instability |
| DNA repair | Damage response |
| EMT | Epithelial-mesenchymal transition |
| Hypoxia | Low oxygen response |
| Inflammation | Immune response |
| Invasion | Metastatic potential |
| Metastasis | Distant spread |
| Proliferation | Cell growth |
| Quiescence | Dormant state |
| Stemness | Cancer stem cells |

#### Geneset Entries

```json
"CancerSEA_angiogenesis": {
  "genes": ["VEGFA", "VEGFB", "VEGFC", "KDR", "FLT1", "FLT4", "ANGPT1", "ANGPT2", "TEK", "PGF", "NRP1", "NRP2", "EFNB2"],
  "description": "CancerSEA angiogenesis state signature."
},
"CancerSEA_apoptosis": {
  "genes": ["BAX", "BAK1", "BCL2", "BCL2L1", "BIRC5", "CASP3", "CASP8", "CASP9", "FAS", "TNFRSF10A", "PARP1", "BAD", "BID"],
  "description": "CancerSEA apoptosis state signature."
},
"CancerSEA_EMT": {
  "genes": ["CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "FN1", "MMP2", "MMP9", "TGFB1"],
  "description": "CancerSEA EMT state signature."
},
"CancerSEA_stemness": {
  "genes": ["SOX2", "NANOG", "PROM1", "NES", "CD44", "ALDH1A1", "POU5F1", "MYC", "KLF4", "LIN28A", "SALL4", "DPPA4"],
  "description": "CancerSEA stemness state signature."
},
"CancerSEA_invasion": {
  "genes": ["MMP1", "MMP2", "MMP3", "MMP9", "MMP13", "PLAU", "PLAUR", "CTSD", "CTSB", "SERPINE1", "SERPINE2"],
  "description": "CancerSEA invasion state signature."
},
"CancerSEA_metastasis": {
  "genes": ["S100A4", "S100A6", "VIM", "CDH2", "TGFB1", "TGFB2", "MMP2", "MMP9", "CXCR4", "CXCL12", "LOX", "LOXL2"],
  "description": "CancerSEA metastasis state signature."
}
```

---

### Article 18: eBioMedicine 2024 — Multilayered Metabolic Heterogeneity (Zhou et al.)

**Dataset**: 296 tumour and normal samples spanning 6 common cancer types (BRCA, CRC, HCC, LUAD, OV, PDAC)
**Key methods**: NMF-based identification of metabolic meta-programs (MMPs); SCENIC for metabolic regulons
**Key finding**: Shared glycolysis upregulation and divergent citric acid cycle regulation across cell types; CRC-specific MMP7 associated with cuproptosis resistance

#### Metabolic Meta-Programs (MMPs)

The paper identified **7 metabolic meta-programs** showing intratumour heterogeneity across cancer types:

| MMP | Key Genes | Biological Interpretation | Cancer Type Association |
|-----|-----------|--------------------------|------------------------|
| MMP1 | Glycolysis core (HK2, PFKP, ALDOA, GAPDH, PGK1, ENO1, PKM, LDHA) | Shared glycolysis upregulation | All cancer types |
| MMP2 | Citric acid cycle (IDH1, IDH2, SUCLA2, SDHB, FH, MDH2) | TCA cycle regulation | Context-dependent |
| MMP3 | Oxidative phosphorylation (NDUFA4, COX4I1, ATP5F1A, UQCRB, SDHB) | Mitochondrial respiration | Context-dependent |
| MMP4 | Fatty acid metabolism (ACLY, FASN, SCD, ELOVL6, ACACA) | Lipid synthesis | HCC, BRCA enriched |
| MMP5 | Amino acid metabolism (GLS, GCLC, PSAT1, ASNS, GPT2) | Glutamine/asparagine metabolism | PDAC enriched |
| MMP6 | Nucleotide metabolism (TYMS, RRM1, RRM2, DTYMK, IMPDH1) | DNA/RNA synthesis | Proliferating cells |
| MMP7 | Copper metabolism (FDX1, LIAS, LIPT1, DLD, DLAT) | Cuproptosis-related | CRC-specific |

#### Geneset Entries

```json
"MMP_glycolysis": {
  "genes": ["HK2", "PFKP", "ALDOA", "GAPDH", "PGK1", "ENO1", "PKM", "LDHA", "SLC2A1", "PGAM1", "TPI1"],
  "description": "Metabolic meta-program 1: Glycolysis. Shared upregulation across cell types and cancer types. Core Warburg effect genes."
},
"MMP_TCA_cycle": {
  "genes": ["IDH1", "IDH2", "SUCLA2", "SDHB", "FH", "MDH2", "CS", "ACO2", "OGDH", "DLST"],
  "description": "Metabolic meta-program 2: Citric acid cycle. Divergent regulation across cell types. Context-dependent activity."
},
"MMP_oxphos": {
  "genes": ["NDUFA4", "COX4I1", "ATP5F1A", "ATP5F1B", "UQCRB", "SDHB", "CYC1", "ATP5MG"],
  "description": "Metabolic meta-program 3: Oxidative phosphorylation. Mitochondrial respiration program."
},
"MMP_fatty_acid": {
  "genes": ["ACLY", "FASN", "SCD", "ELOVL6", "ACACA", "SREBF1", "ACSS2", "GPAM"],
  "description": "Metabolic meta-program 4: Fatty acid metabolism. Lipid synthesis and storage. Enriched in HCC and BRCA."
},
"MMP_amino_acid": {
  "genes": ["GLS", "GCLC", "PSAT1", "ASNS", "GPT2", "GLUD1", "GLUL", "ASS1"],
  "description": "Metabolic meta-program 5: Amino acid metabolism. Glutamine and asparagine metabolism. Enriched in PDAC."
},
"MMP_nucleotide": {
  "genes": ["TYMS", "RRM1", "RRM2", "DTYMK", "IMPDH1", "IMPDH2", "DHFR", "TK1"],
  "description": "Metabolic meta-program 6: Nucleotide metabolism. DNA/RNA synthesis support. Associated with proliferating cells."
},
"MMP_copper_cuproptosis": {
  "genes": ["FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "MTF1"],
  "description": "Metabolic meta-program 7: Copper metabolism and cuproptosis. CRC-specific. Associated with resistance to elesclomol (cuproptosis inducer)."
}
```

#### Tumor Context (`marker_tumor`)

**CRC-Specific Cuproptosis-Resistant Metabolic State**

```toml
[[cancer_state]]
name = "CRC cuproptosis-resistant metabolic state"
color = "#B87333"
markers = ["FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "MTF1"]
negative_markers = ["FDX1_low", "LIAS_low", "LIPT1_low"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
cancer_type = ["Colorectal Cancer"]
scope = "cancer_type_specific"
applies_to = ["Colon", "Rectum"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = true
species = "human"
source = { title = "Single cell atlas reveals multilayered metabolic heterogeneity across tumour types", year = "2024", doi = "10.1016/j.ebiom.2024.105389" }
notes = "CRC-specific metabolic meta-program 7 (MMP7) associated with cuproptosis resistance. High expression of FDX1, LIAS, LIPT1, DLD, DLAT. Cells with high MMP7 scores show resistance to elesclomol (cuproptosis inducer). Copper metabolism dysregulation in CRC."
```

#### Key Biological Insights

1. **Metabolic heterogeneity is multilayered**: Both inter-tumour (across cancer types) and intra-tumour (across cell types within a tumour) metabolic heterogeneity exists. scLucid should score metabolic programs at single-cell resolution.

2. **Glycolysis is universally upregulated**: MMP1 (glycolysis) is shared across all cancer types and cell types, confirming the Warburg effect as a pan-cancer metabolic hallmark.

3. **Cuproptosis resistance is CRC-specific**: MMP7 identifies a CRC-specific metabolic vulnerability. scLucid tumor interpretation should flag high MMP7 scores as potential elesclomol resistance in CRC.

4. **Cell type-specific metabolic signatures**: Malignant cells, immune cells, and stromal cells have distinct metabolic programs. scLucid should score metabolic programs separately by cell type.

---

## Summary: Key Conflicts and Resolutions

### 1. MPs vs Cancer Hallmarks

**Issue**: Meta-programs (MPs, Tirosh lab) and cancer hallmarks (Hanahan) are overlapping but distinct frameworks.

**Resolution**:
- **MPs**: Data-driven, recurring gene expression programs from scRNA-seq (gene modules)
- **Hallmarks**: Conceptual framework for cancer biology (biological processes)
- Both can be represented as **genesets** for scoring
- Cross-reference MPs to hallmarks where applicable (e.g., MP_hypoxia ↔ Hallmarks_hypoxia)

### 2. Cell Line vs Primary Tumor MPs

**Issue**: MPs defined in cell lines (Nature Genetics 2020) may not fully recapitulate primary tumor MPs (Nature 2023).

**Resolution**:
- Cell line MPs: 12 programs (Nature Genetics 2020)
- Primary tumor MPs: 149 programs across 8 cell types (Nature 2023)
- Use primary tumor MPs as primary reference
- Cell line MPs as supplementary validation

### 3. Identity vs State in Tumor Cells

**Issue**: Tumor cell states (pEMT, stemness, stress) overlap with normal cell states.

**Resolution**:
- Tumor cell states are **cancer_state** entries in `marker_tumor_human.toml`
- They are NOT identity markers (do not add to `marker_registry`)
- Use `use_for_malignancy_interpretation = true` and `use_for_global_annotation = false`
- Distinguish from normal counterparts by context (tumor vs normal tissue)

### 4. Geneset Size

**Issue**: Some hallmark/MP signatures may exceed the recommended 5-50 gene limit for genesets.

**Resolution**:
- Hallmark signatures: ~50-200 genes (standard MSigDB size)
- MP signatures: ~20-50 genes (meta-program cores)
- CancerSEA signatures: ~20-50 genes
- All are acceptable for module scoring
- For very large signatures, provide "core" subset of top 20-30 genes

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New geneset programs | 47+ | genesets_cancer_signatures.json |
| New tumor context entries | 4 (pEMT, stemness, ecotype adverse, CRC cuproptosis-resistant) | marker_tumor |
| Review validations | 6 (hallmarks framework, chemokines, immunoregulation) | marker_tumor |
| Conflict resolutions | 4 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **Meta-programs explain ITH**: Tirosh lab's MPs explain a large fraction of intratumor heterogeneity. scLucid should score MPs for tumor cell characterization.

2. **pEMT is distinct from full EMT**: Partial EMT (pEMT) is a recurring cancer cell state associated with invasion, not a terminal differentiation. scLucid should distinguish pEMT from full EMT (which is more like a CAF/endothelial transition).

3. **Stemness state is dynamic**: Cancer stem cell states are plastic, not fixed identities. scLucid should score stemness as a continuous state rather than a binary classification.

4. **Hallmarks are functional, not identity**: Hanahan's hallmarks describe biological capabilities, not cell types. All belong in genesets.

5. **ICB response/resistance/toxicity are distinct**: Response, resistance, and toxicity to ICB have different molecular signatures. scLucid tumor interpretation should score all three separately.
