> **⚠️ ARCHIVED / SUPERSEDED**
>
> This batch curation note is kept for provenance only. The live marker
> resource status is tracked in `docs/marker_resources/marker_curation_literature_index.jsonl`,
> `docs/marker_resources/marker_resource_quality_gaps.jsonl`, `docs/marker_resources/marker_curation_candidates.jsonl`,
> and `docs/marker_resources/CURATION.md`. New curation should follow the current
> contract rather than adding new files to this archive.

---

# Marker Curation Batch 01: Cross-Tissue Single-Cell Atlas Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type |
|---|-------|------|---------|-----|-------------|
| 1 | Cross-tissue multicellular coordination and its rewiring in cancer | 2025 | Nature | 10.1038/s41586-025-09053-4 | single_cell_atlas |
| 2 | Single-cell atlases: shared and tissue-specific cell types across human organs | 2022 | Nature Reviews Genetics | 10.1038/s41576-022-00449-w | review |
| 3 | The tabula sapiens: a multiple-organ, single-cell transcriptomic atlas of humans | 2022 | Science | 10.1126/science.abl4896 | single_cell_atlas |
| 4 | A single–cell type transcriptomics map of human tissues | 2021 | Science Advances | 10.1126/sciadv.abh2169 | single_cell_atlas |
| 5 | Construction of a human cell landscape at single-cell level | 2020 | Nature | 10.1038/s41586-020-2157-4 | single_cell_atlas |
| 6 | Cross-tissue immune cell analysis reveals tissue-specific features in humans | 2022 | Science | 10.1126/science.abl5197 | single_cell_atlas |

## Atlas Scope Summary

- **Article 1** (2025, Nature): Cross-tissue atlas from BIOPIC/PKU. ~45% immune compartment. 76 non-epithelial subsets annotated with marker genes across 30+ tissues.
- **Article 3** (2022, Science): Tabula Sapiens. 483,152 cells from 24 tissues, 475 cell types. Same-donor multi-tissue design.
- **Article 5** (2020, Nature): Human Cell Landscape (HCL). Microwell-seq based atlas covering major organs.
- **Article 4** (2021, Sci Adv): 192 cell type clusters with spatial antibody-based protein profiling.
- **Article 2** (2022, Nat Rev Genet): Review summarizing cross-tissue insights for epithelial, fibroblast, vascular, and immune cells.
- **Article 6** (2022, Science): Cross-tissue immune cell atlas from 16 tissues, 12 donors. ~360,000 cells. Developed CellTypist automated annotation framework. Identified 101 immune cell types/states with tissue-specific features.

## Curation Principles Applied

- **Species**: human
- **Scope**: all (cross-tissue shared markers) or tissue_specific (tissue-restricted subtypes)
- **Evidence tier**: `atlas_supported` (supported by at least one well-described atlas)
- **Review status**: `needs_review` (default for atlas-extracted markers pending manual verification)
- **Source type**: `single_cell_atlas` for primary atlas papers; `review` for the Nat Rev Genet summary
- **Routing flags**:
  - Compartment/lineage entries: `use_for_global_annotation = true`
  - State/program entries: `use_for_global_annotation = false`
  - All normal tissue entries: `use_for_malignancy_interpretation = false`

## Extracted Marker Evidence

### Compartment Level (Cross-tissue validated)

These confirm existing compartment markers and add tissue-contextual nuance:

| Compartment | Core Markers | Cross-tissue Atlas Support | Notes |
|-------------|-------------|---------------------------|-------|
| Immune | PTPRC (CD45) | All 5 atlases | ~45% of cells in cross-tissue atlas (Article 1) |
| Epithelial | EPCAM, KRT8/18/19, CDH1 | All 5 atlases | Tissue-specific variation noted in Articles 2, 3 |
| Stromal | COL1A1, PDGFRA, VIM | All 5 atlases | Includes fibroblast, pericyte, SMC subtypes |
| Endothelial | PECAM1, VWF, CDH5 | All 5 atlases | Shared across tissues with subtle differences (Article 3) |
| Neural | RBFOX3, GFAP, PLP1 | Articles 3, 5 | Tissue-restricted to neural tissues |

### Lineage/Subtype Level (New or Validated Markers)

#### 1. Fibroblast Subtypes (from Article 1, 76-subset annotation)

Article 1 provides the most comprehensive cross-tissue fibroblast subtype annotation:

| Subtype | Marker | Tissue Context | Evidence |
|---------|--------|---------------|----------|
| Universal fibroblast | PI16 | Cross-tissue | S01_Fb_PI16 - pan-tissue universal |
| Fibroblast - C7+ | C7 | Cross-tissue | S02_Fb_C7 |
| Fibroblast - LAMC1+ | LAMC1 | Cross-tissue | S03_Fb_LAMC1 |
| Fibroblast - CD9+ | CD9 | Cross-tissue | S04_Fb_CD9 |
| Fibroblast - DPEP1+ | DPEP1 | Cross-tissue | S05_Fb_DPEP1 |
| Fibroblast - ARID5B+ | ARID5B | Cross-tissue | S06_Fb_ARID5B |
| Fibroblast - TCF21+ | TCF21 | Cross-tissue | S07_Fb_TCF21 (mesothelial-associated) |
| Fibroblast - PTGES+ | PTGES | Cross-tissue | S08_Fb_PTGES |
| Fibroblast - CEBPB+ | CEBPB | Cross-tissue | S09_Fb_CEBPB |
| Fibroblast - MXRA5+ | MXRA5 | Cross-tissue | S10_Fb_MXRA5 |
| Fibroblast - IGFBP2+ | IGFBP2 | Cross-tissue | S11_Fb_IGFBP2 |
| Fibroblast - MMP11+ | MMP11 | Cross-tissue | S12_Fb_MMP11 |
| Pericyte | CD36, RGS5, PDGFRB | Cross-tissue | S13_Pericyte_CD36 |
| Smooth muscle (RERGL+) | RERGL | Cross-tissue | S14_SMC_RERGL |
| Smooth muscle (ACTG2+) | ACTA2, ACTG2 | Cross-tissue | S15_SMC_ACTG2 |
| Smooth muscle (FBLN5+) | FBLN5 | Cross-tissue | S16_SMC_FBLN5 |

**Curation note**: Many of these fibroblast subtypes overlap with existing iCAF/mCAF/pCAF annotations in the current registry. The Article 1 annotations are more granular and tissue-contextual. Recommended approach: add as `subtype` entries with `scope = "tissue_specific"` or merge with existing CAF subtypes after manual review.

#### 2. Endothelial Subtypes (from Article 1)

| Subtype | Marker | Tissue Context | Evidence |
|---------|--------|---------------|----------|
| Arterial EC | IGFBP3 | Cross-tissue | E01_Artery_IGFBP3 |
| Capillary EC | CA4, VIPR1 | Cross-tissue | E02_Capillary_CA4, E03_Capillary_VIPR1 |
| Venous EC | ACKR1 | Cross-tissue | E04_Vein_ACKR1 |
| Lymphatic EC | LYVE1 | Cross-tissue | E05_Lymph_LYVE1 |

**Validation against existing registry**:
- Arterial: existing uses SEMA3G, GJA5, BMX → IGFBP3 is complementary
- Capillary: existing uses RGCC, CA4, CD36 → CA4 confirmed; VIPR1 is new
- Venous: existing uses ACKR1, SELP, NRG1 → ACKR1 confirmed
- Lymphatic: existing uses PROX1, LYVE1, CCL21 → LYVE1 confirmed

#### 3. Myeloid Subtypes (from Article 1)

| Subtype | Marker | Tissue Context | Evidence |
|---------|--------|---------------|----------|
| cDC1 | CLEC9A | Cross-tissue | M01_cDC1_CLEC9A |
| cDC2 | CD1C | Cross-tissue | M02_cDC2_CD1C |
| Migratory DC | LAMP3 | Cross-tissue | M03_cDC_LAMP3 |
| pDC | LILRA4 | Cross-tissue | M04_pDC_LILRA4 |
| Classical Monocyte | CD14 | Cross-tissue | M05_Mo_CD14 |
| Non-classical Monocyte | FCGR3A | Cross-tissue | M06_Mo_FCGR3A |
| Inflammatory Macrophage | FCN1 | Cross-tissue | M07_Mph_FCN1 |
| NLRP3+ Macrophage | NLRP3 | Cross-tissue | M08_Mph_NLRP3 |
| Resident Macrophage | FOLR2 | Cross-tissue | M09_Mph_FOLR2 |
| CD5L+ Macrophage | CD5L | Cross-tissue | M10_Mph_CD5L |
| PPARG+ Macrophage | PPARG | Cross-tissue | M11_Mph_PPARG |
| MT1X+ Macrophage | MT1X | Cross-tissue | M12_Mph_MT1X |
| Immature Neutrophil | MMP8 | Cross-tissue | M13_immNeu_MMP8 |
| Mature Neutrophil | CXCR2 | Cross-tissue | M14_mNeu_CXCR2 |
| Mast cell | CPA3 | Cross-tissue | M15_Mast_CPA3 |

**Validation against existing registry**:
- cDC1: CLEC9A confirmed (existing: CLEC9A, XCR1, CADM1)
- cDC2: CD1C confirmed (existing: CD1C, FCER1A, CLEC10A)
- cDC3/migratory: LAMP3 confirmed (existing: LAMP3, CCR7, FSCN1, CST7)
- pDC: LILRA4 confirmed (existing: LILRA4, GZMB, IL3RA, TCF4)
- Classical Mono: CD14 confirmed (existing: FCN1, S100A8, S100A9, CD14, VCAN)
- Non-classical Mono: FCGR3A confirmed (existing: FCGR3A, LST1, LILRB2, HK3)
- Mast: CPA3 confirmed (existing: TPSAB1, KIT, CPA3, TPSB2)

**New additions suggested**:
- FOLR2+ resident macrophage (distinct from existing TRM TAM: LYVE1, FOLR2, CX3CR1, MERTK)
- CD5L+ macrophage subtype
- PPARG+ macrophage subtype
- MT1X+ macrophage subtype
- MMP8+ immature neutrophil (more specific than existing FCGR3B, CXCR2, CSF3R)

#### 4. NK/ILC Subtypes (from Article 1)

| Subtype | Marker | Tissue Context | Evidence |
|---------|--------|---------------|----------|
| CD16hi NK (CREM+) | CREM | Cross-tissue | I01_CD16hiNK_CREM |
| CD16hi NK (SYNE2+) | SYNE2 | Cross-tissue | I02_CD16hiNK_SYNE2 |
| CD16hi NK (GZMB+) | GZMB | Cross-tissue | I03_CD16hiNK_GZMB |
| CD16hi NK (HSPA1A+) | HSPA1A | Cross-tissue | I04_CD16hiNK_HSPA1A |
| CD16lo NK (SELL+) | SELL | Cross-tissue | I05_CD16loNK_SELL |
| CD16lo NK (NR4A2+) | NR4A2 | Cross-tissue | I06_CD16loNK_NR4A2 |
| CD16lo NK (CXCR6+) | CXCR6 | Cross-tissue | I07_CD16loNK_CXCR6 |
| ILC3 | KIT | Cross-tissue | I08_ILC3_KIT |
| γδ T | ITGA1 | Cross-tissue | I09_gdT_ITGA1 |

**Validation against existing registry**:
- CD56+CD16- NK: existing uses NCAM1, KLRK1, KIT, GZMK → overlaps with CD16lo
- CD56-CD16+ NK: existing uses FCGR3A, NCR3, GZMB, PRF1 → overlaps with CD16hi
- ILC: existing uses ID2, IL7R, RORC, GATA3 → KIT is consistent
- Gamma delta T: existing uses TRDC, TRGC1, TRGC2, NKG7 → ITGA1 is new marker

**New additions suggested**:
- CD16hi NK subtypes: CREM, SYNE2, HSPA1A (stress response)
- CD16lo NK subtypes: SELL (naive-like), NR4A2, CXCR6
- γδ T: ITGA1 as additional marker

#### 5. T Cell Subtypes (from Article 1)

| Subtype | Marker | Tissue Context | Evidence |
|---------|--------|---------------|----------|
| CD4 Tn (SOX4+) | SOX4 | Cross-tissue | CD4T01_Tn_SOX4 |
| CD4 Tn (CCR7+) | CCR7 | Cross-tissue | CD4T02_Tn_CCR7 |
| CD4 Tn (NR4A1+) | NR4A1 | Cross-tissue | CD4T03_Tn_NR4A1 |
| CD4 Tfh | IL6ST | Cross-tissue | CD4T04_Tfh_IL6ST |
| CD4 Tm (LTB+) | LTB | Cross-tissue | CD4T05_Tm_LTB |
| CD4 Tm (ANXA1+) | ANXA1 | Cross-tissue | CD4T06_Tm_ANXA1 |
| CD4 Tem | GZMK | Cross-tissue | CD4T07_Tem_GZMK |
| Treg | FOXP3 | Cross-tissue | CD4T08_Treg_FOXP3 |
| CD8 Tn | CCR7 | Cross-tissue | CD8T01_Tn_CCR7 |
| CD8 Tem | GZMK | Cross-tissue | CD8T02_Tem_GZMK |
| CD8 Trm (ITGA1+) | ITGA1 | Cross-tissue | CD8T03_Trm_ITGA1 |
| CD8 Trm (HSPA1A+) | HSPA1A | Cross-tissue | CD8T04_Trm_HSPA1A |
| CD8 Temra | GZMH | Cross-tissue | CD8T05_Temra_GZMH |
| MAIT | SLC4A10 | Cross-tissue | CD8T06_MAIT_SLC4A10 |

**Validation against existing registry**:
- CD4 Tn: existing uses TCF7, SELL, CCR7, LEF1, MAL → CCR7 confirmed; SOX4 and NR4A1 are new
- CD4 Tm: existing uses IL7R, S100A4, GZMK, CCL4, CAPG → LTB and ANXA1 are complementary
- CD4 Tem: existing uses GZMK, KLRB1, CCL5, CXCR4 → GZMK confirmed
- Treg: existing uses FOXP3, IL2RA, CTLA4, TNFRSF9, RTKN2 → FOXP3 confirmed
- CD8 Tn: existing uses TCF7, SELL, CCR7, LEF1, MAL → CCR7 confirmed
- CD8 Tem: existing uses GZMK, CCL3, CXCR4, KLRG1 → GZMK confirmed
- CD8 Trm: existing uses ZNF683, CD52, HOPX, S100A4 → ITGA1 and HSPA1A are new
- CD8 Temra: existing uses CX3CR1, GZMH, TBX21, KLRD1 → GZMH confirmed
- MAIT: existing uses TRAV1-2, SLC4A10, KLRB1, ZBTB16 → SLC4A10 confirmed
- Tfh: existing uses CXCR5, BCL6, IL21, MAF, PDCD1, TOX, TCF7 → IL6ST is new

#### 6. B Cell Subtypes (from Article 1)

| Subtype | Marker | Tissue Context | Evidence |
|---------|--------|---------------|----------|
| Naive B (IGHM+) | IGHM | Cross-tissue | B01_Bn_IGHM |
| Naive B (TCL1A+) | TCL1A | Cross-tissue | B02_Bn_TCL1A |
| Naive B (NR4A2+) | NR4A2 | Cross-tissue | B03_Bn_NR4A2 |
| Memory B (CD27+) | CD27 | Cross-tissue | B04_Bm_CD27 |
| Memory B (NR4A2+) | NR4A2 | Cross-tissue | B05_Bm_NR4A2 |
| Memory B (ITGB1+) | ITGB1 | Cross-tissue | B06_Bm_ITGB1 |
| Memory B (HSPA1A+) | HSPA1A | Cross-tissue | B07_Bm_HSPA1A |
| Atypical memory B | FCRL5 | Cross-tissue | B08_ABC_FCRL5 |
| GC B | RGS13 | Cross-tissue | B09_GCB_RGS13 |
| Plasmablast | MKI67 | Cross-tissue | B10_Plasmablast_MKI67 |
| Plasma (IGHG2+) | IGHG2 | Cross-tissue | B11_Plasma_IGHG2 |
| Plasma (IGHA2+) | IGHA2 | Cross-tissue | B12_Plasma_IGHA2 |

**Validation against existing registry**:
- Naive B: existing uses TCL1A, FCER2, IGHD, YBX3 → TCL1A confirmed; IGHM and NR4A2 are complementary
- Memory B: existing uses CD27, IGHG1, AIM2, ITGB1 → CD27 and ITGB1 confirmed
- GC B: existing uses LMO2, CXCR4, MKI67, MEF2B → RGS13 is new
- Plasmablast: existing uses TXNDC5, MYDGF → MKI67 (proliferation marker) is state indicator
- Plasma: existing uses JCHAIN, XBP1, SDC1, MZB1 → IGHG2 and IGHA2 are isotype-specific

#### 7. Other Tissue-Specific Cell Types (from Article 1)

| Cell Type | Marker | Tissue | Evidence |
|-----------|--------|--------|----------|
| Cardiomyocyte | TNNT2 | Heart | S17_CMC_TNNT2 |
| Skeletal muscle | MYOZ1 | Muscle | S18_SkMC_MYOZ1 |
| Satellite cell | CXCL14 | Muscle | S19_Satellite_CXCL14 |
| Glia | CDH19 | Neural | S20_Glia_CDH19 |
| Melanocyte | DCT | Skin | S21_Melanocyte_DCT |

**Curation note**: These belong in `marker_tissue_human.toml` under their respective tissue sections, not in the global registry.

---

### Cross-Tissue Immune Cell Atlas (Article 6: Domínguez Conde et al., Science 2022)

**Dataset**: 16 tissues from 12 adult donors, ~360,000 single-cell transcriptomes (~330,000 immune cells)
**Key methods**: CellTypist automated annotation (machine learning with logistic regression); VDJ sequencing for clonal dynamics
**Key finding**: 101 immune cell types/states with tissue-specific expression modules; CellTypist framework for automated immune cell annotation

#### CellTypist Annotation Framework

CellTypist provides a two-level hierarchy for immune cell annotation:

| Hierarchy Level | Cell Types | Accuracy (F1) |
|----------------|-----------|---------------|
| High (low-resolution) | 32 cell types | ~0.9 |
| Low (high-resolution) | 91 cell types + 10 novel = 101 | ~0.9 |

**Key immune cell populations identified**:

| Compartment | Cell Types | Tissue Distribution |
|------------|-----------|---------------------|
| Mononuclear phagocytes | 8 populations | Lung, liver, spleen, lymph nodes, gut |
| Dendritic cells | 3 subsets (DC1, DC2, migDC) | All tissues; migDCs enriched in lymph nodes |
| B cells | 12 states | Spleen, lymph nodes, gut, bone marrow |
| T cells | 18 subtypes | All tissues; tissue-resident in gut, lung |
| NK/ILC | 3 clusters | All tissues; CD16hi vs CD56bright |

#### Tissue-Specific Myeloid Markers (from Article 6)

| Cell Type | Tissue-Specific Markers | Tissue | Notes |
|-----------|------------------------|--------|-------|
| Alveolar macrophages | GPNMB, TREM2 | Lung | Disease-associated macrophage markers |
| Intermediate macrophages | TNIP3 | Lung | A20-binding protein; inflammatory regulation |
| Erythrophagocytic macrophages | CD5L, SLC40A1, SPIC | Spleen, liver, lymph nodes | Iron recycling; red pulp/Kupffer cells |
| Gut macrophages | CD209 (DC-SIGN), IGF1 | Gut | M2-like mature intestinal macrophages |
| Classical monocytes | CD14, FCGR3A, CX3CR1 | Blood, lung, liver | CD14hi FCGR3Alo |
| Non-classical monocytes | FCGR3A, CX3CR1 | Blood, lung, liver | CD14lo FCGR3Ahi |
| DC1 | XCR1, CLEC9A | All tissues | Cross-presentation specialized |
| DC2 | CD1C, CLEC10A | All tissues | Conventional DC2 |
| migDCs | CCR7, LAMP3, AIRE, PDLIM4, EBI3 | Lymph nodes | Migratory; extrathymic AIRE expression |

#### Tissue-Specific Lymphoid Markers (from Article 6)

| Cell Type | Tissue-Specific Markers | Tissue | Notes |
|-----------|------------------------|--------|-------|
| Treg (gut) | CCR9, ITGAE, ITGA1 | Gut, lung | Tissue-resident regulatory T cells |
| Trm_gut_CD8 | CCR9, ITGAE (CD103), ITGA1 (CD49a) | Gut | Tissue-resident memory T cells |
| Trm/em_CD8 | CCR9, CX3CR1, CRTAM | Gut, lung | Effector memory / resident memory |
| Tem/emra_CD8 | CX3CR1, CRTAM | Blood, spleen | Effector memory re-expressing CD45RA |
| MAIT | TRAV1-2, TRAJ33/12/29/36 | Spleen, liver, gut | Tissue-specific TRAJ segment usage |
| gd T cells (ITGAD+) | ITGAD, CD52-CD127- | Gut, lung | Integrin alpha subunit specific |
| NK_CD16+ | FCGR3A | Blood, spleen | Cytotoxic NK subset |
| NK_CD56bright | NCAM1 | All tissues | Less mature/regulatory NK subset |
| ILC3 | KIT, RORC, IL7R | Gut | Innate lymphoid cell type 3 |
| Atypical memory B | ITGAX, TBX21, FCRL2 | Lymph nodes, spleen | CD11c+ T-bet+ B cells |
| Plasma cells (gut) | IGHA2 | Gut (jejunum LP) | Mucosal IgA2 dominance |
| Plasma cells (bone marrow) | IGHG2 | Bone marrow, liver, spleen | Systemic IgG2 |

#### Novel Cross-Tissue Insights

1. **Tissue-resident memory T cells show highest clonal expansion**: TRM cells harbor the most frequent clonal sharing between resident and effector memory populations.

2. **Macrophages show the most prominent tissue specificity**: Erythrophagocytic macrophages in liver/spleen share iron-recycling features with mesenteric lymph nodes.

3. **Plasma cell isotype bias is tissue-dependent**: IgA2 dominates in gut; IgG2 in bone marrow/liver/spleen; IgA1 in mesenteric lymph nodes.

4. **Migratory DCs express extrathymic AIRE**: AIRE+ migDCs in lymph nodes may play a role in peripheral tolerance.

#### Validation against Existing Registry

| Existing Entry | Atlas Support (Article 6) | Notes |
|---------------|--------------------------|-------|
| cDC1 (CLEC9A, XCR1) | Confirmed; XCR1 and CLEC9A co-expression | Cross-presentation specialized |
| cDC2 (CD1C, CLEC10A) | Confirmed; CD1C and CLEC10A co-expression | Conventional DC2 |
| migDC (LAMP3, CCR7) | Confirmed; added AIRE, PDLIM4, EBI3 | Lymph node enriched |
| Classical monocytes (CD14+) | Confirmed; CD14, FCGR3A, CX3CR1 | Blood, lung, liver |
| Non-classical monocytes (FCGR3A+) | Confirmed; FCGR3A, CX3CR1 | Blood, lung, liver |
| CD8 Trm (ZNF683, CD52, HOPX) | Confirmed; added ITGA1, ITGAE, CCR9 | Gut and lung specific |
| MAIT (TRAV1-2, SLC4A10) | Confirmed; tissue-specific TRAJ usage | Spleen, liver, gut |
| NK_CD16+ (FCGR3A) | Confirmed | Cytotoxic subset |
| NK_CD56bright (NCAM1) | Confirmed | Less mature subset |

#### Recommended TOML Additions (from Article 6)

```toml
# --- Tissue-specific myeloid states ---

[[state.minor]]
name = "Alveolar macrophage-like"
color = "#9370DB"
markers = ["GPNMB", "TREM2", "MERTK"]
negative_markers = ["FCN1", "S100A8", "S100A9", "CD14"]

[state.minor.metadata]
kind = "state"
category = "tissue_resident"
scope = "tissue_specific"
applies_to = ["Macrophages", "Lung"]
alias_of = "Alveolar_macrophage"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
notes = "Alveolar macrophage-like state. GPNMB+ TREM2+. Disease-associated macrophage markers. Lung-specific. Negative inflammatory monocyte markers (FCN1, S100A8/9, CD14)."

[[state.minor]]
name = "Erythrophagocytic macrophage"
color = "#8B4513"
markers = ["CD5L", "SLC40A1", "SPIC", "CCL2"]
negative_markers = ["FCN1", "IL1B", "CXCL8"]

[state.minor.metadata]
kind = "state"
category = "tissue_function"
scope = "tissue_specific"
applies_to = ["Macrophages", "Spleen", "Liver", "Lymph nodes"]
alias_of = "Red_pulp_macrophage"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
notes = "Erythrophagocytic macrophage state. CD5L+ SLC40A1+ SPIC+. Iron recycling function. Enriched in spleen (red pulp), liver (Kupffer cells), and lymph nodes. Negative inflammatory markers (FCN1, IL1B, CXCL8)."

[[state.minor]]
name = "Gut-resident macrophage"
color = "#228B22"
markers = ["CD209", "IGF1", "MRC1", "CD163"]
negative_markers = ["FCN1", "S100A8", "S100A9", "CXCL8"]

[state.minor.metadata]
kind = "state"
category = "tissue_resident"
scope = "tissue_specific"
applies_to = ["Macrophages", "Gut"]
alias_of = "Intestinal_macrophage"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
notes = "Gut-resident macrophage state. CD209 (DC-SIGN)+ IGF1+ M2-like mature intestinal macrophage. Negative inflammatory monocyte markers (FCN1, S100A8/9, CXCL8)."
```

## Key Cross-Tissue Insights from Article 2 (Nat Rev Genet Review)

The review highlights several important principles for marker curation:

1. **Shared vs tissue-specific cell types**: Immune cells (T, B, myeloid, NK) are broadly shared across tissues with subtle expression differences. Epithelial and stromal cells show more tissue-specific diversity.

2. **Endothelial heterogeneity**: Endothelial cells are shared across tissues but show clear tissue-specific expression patterns (confirmed in Article 3).

3. **Fibroblast diversity**: Cross-tissue fibroblasts show remarkable diversity with both universal and tissue-restricted subtypes.

4. **Immune cell recirculation**: T cell clones shared between organs; B cell hypermutation rates vary by organ (Article 3).

## Recommended TOML Additions

### For `marker_registry_human.toml`

#### New Subtype Entries (from Article 1)

```toml
# --- Fibroblast cross-tissue subtypes ---

[[subtype]]
name = "PI16+ Universal Fibroblast"
color = "#D2691E"
markers = ["PI16", "FAP", "THY1"]
negative_markers = ["PTPRC", "EPCAM", "PECAM1", "ACTA2"]

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
source = { title = "Cross-tissue multicellular coordination and its rewiring in cancer", year = "2025", doi = "10.1038/s41586-025-09053-4" }
notes = "Pan-tissue universal fibroblast marker from cross-tissue atlas (S01_Fb_PI16). Validates across 30+ tissues."

[[subtype]]
name = "FOLR2+ Resident Macrophage"
color = "#9370DB"
markers = ["FOLR2", "LYVE1", "MERTK"]
negative_markers = ["FCN1", "S100A8", "S100A9", "CD14"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "all"
applies_to = ["Macrophages"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue multicellular coordination and its rewiring in cancer", year = "2025", doi = "10.1038/s41586-025-09053-4" }
notes = "Resident tissue macrophage subtype distinct from inflammatory monocytes (M09_Mph_FOLR2). Low FCN1/S100A8 confirms non-inflammatory state."

[[subtype]]
name = "CD5L+ Macrophage"
color = "#8A2BE2"
markers = ["CD5L", "MARCO", "MERTK"]
negative_markers = ["FCN1", "IL1B"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "myeloid"
scope = "lineage_restricted"
applies_to = ["Macrophages", "Liver", "Lung"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue multicellular coordination and its rewiring in cancer", year = "2025", doi = "10.1038/s41586-025-09053-4" }
notes = "M10_Mph_CD5L from cross-tissue atlas. CD5L is a marker of alternatively activated/resident macrophages."

[[subtype]]
name = "CD16hi NK (GZMB+)"
color = "#228B22"
markers = ["FCGR3A", "GZMB", "PRF1", "NCR3"]
negative_markers = ["CD3D", "CD3E", "NCAM1"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "all"
applies_to = ["NK cells"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue multicellular coordination and its rewiring in cancer", year = "2025", doi = "10.1038/s41586-025-09053-4" }
notes = "Cytotoxic CD16high NK subset (I03_CD16hiNK_GZMB). High GZMB/PRF1 indicates strong cytotoxic potential. Negative NCAM1 reflects CD56-low subset."

[[subtype]]
name = "CD16lo NK (SELL+)"
color = "#90EE90"
markers = ["SELL", "KIT", "GZMK", "NCAM1"]
negative_markers = ["CD3D", "CD3E", "FCGR3A"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "all"
applies_to = ["NK cells"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue multicellular coordination and its rewiring in cancer", year = "2025", doi = "10.1038/s41586-025-09053-4" }
notes = "Less mature/regulatory CD16low NK subset (I05_CD16loNK_SELL). SELL (CD62L) indicates naive/central memory-like phenotype. Positive NCAM1 reflects CD56-bright subset."

[[subtype]]
name = "GC B Cell (RGS13+)"
color = "#4682B4"
markers = ["RGS13", "LMO2", "CXCR4", "MEF2B"]
negative_markers = ["TCL1A", "CD27", "IGHM"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "lineage_restricted"
applies_to = ["B cells"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue multicellular coordination and its rewiring in cancer", year = "2025", doi = "10.1038/s41586-025-09053-4" }
notes = "Germinal center B cell subset (B09_GCB_RGS13). RGS13 is a GC-specific regulator of G-protein signaling. Negative naive markers (TCL1A, IGHM) confirm GC identity."

[[subtype]]
name = "Atypical Memory B (FCRL5+)"
color = "#5F9EA0"
markers = ["FCRL5", "ITGB1", "CD27"]
negative_markers = ["TCL1A", "IGHM", "RGS13"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "lineage_restricted"
applies_to = ["B cells"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue multicellular coordination and its rewiring in cancer", year = "2025", doi = "10.1038/s41586-025-09053-4" }
notes = "Atypical/activated memory B cell subset (B08_ABC_FCRL5). FCRL5 marks autoreactive/atypical memory B cells. Distinct from GC B cells."

[[subtype]]
name = "CD8 Trm (ITGA1+)"
color = "#DAA520"
markers = ["ITGA1", "ZNF683", "CD52", "HOPX"]
negative_markers = ["SELL", "CCR7", "TCF7"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "all"
applies_to = ["CD8+ T"]
evidence_tier = "atlas_supported"
source_type = "single_cell_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Cross-tissue multicellular coordination and its rewiring in cancer", year = "2025", doi = "10.1038/s41586-025-09053-4" }
notes = "Tissue-resident memory CD8 T cell (CD8T03_Trm_ITGA1). ITGA1 (CD49a) is a canonical TRM marker. Negative SELL/CCR7 confirms tissue residency."
```

#### Validation Notes for Existing Entries

The following existing entries in `marker_registry_human.toml` are **strongly supported** by these 5 cross-tissue atlases and could have their `evidence_tier` promoted to `consensus`:

| Existing Entry | Current Evidence Tier | Atlas Support | Recommended Action |
|----------------|----------------------|---------------|-------------------|
| Immune (PTPRC) | `curated_review` | All 5 atlases | Promote to `consensus` |
| Epithelial (EPCAM, KRT8/18/19) | `curated_review` | All 5 atlases | Promote to `consensus` |
| Stromal (COL1A1) | `curated_review` | All 5 atlases | Promote to `consensus` |
| Endothelial (PECAM1, VWF, CDH5) | `curated_review` | All 5 atlases | Promote to `consensus` |
| cDC1 (CLEC9A, XCR1, CADM1) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| cDC2 (CD1C, FCER1A, CLEC10A) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| cDC3 (LAMP3, CCR7, FSCN1, CST7) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| pDC (LILRA4, GZMB, IL3RA, TCF4) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| Classical Mono (FCN1, S100A8, CD14) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| Non-classical Mono (FCGR3A, LST1) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| Treg (FOXP3, IL2RA, CTLA4) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| MAIT (TRAV1-2, SLC4A10) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| Naive B (TCL1A, FCER2, IGHD) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |
| Memory B (CD27, IGHG1) | `seed` | Articles 1, 3 | Promote to `atlas_supported` |

### For `marker_tissue_human.toml`

The following tissue-specific cell types from Article 1 should be added to their respective tissue sections:

```toml
# Add to existing Heart section:
# { name = "Cardiomyocyte", color = "#B22222", markers = ["TNNT2", "MYH6", "MYH7"] }

# Add to existing Muscle section (new):
[["Muscle"]]
name = "Muscle Tissue"
color = "#C0C0C0"
markers = []
metadata = { kind = "tissue_context", granularity = "tissue", tissue = "muscle", compartment = "mixed", scope = "tissue_specific", applies_to = ["Muscle"], evidence_tier = "atlas_supported", source_type = "single_cell_atlas", review_status = "needs_review", use_for_global_annotation = false, use_for_state_annotation = false, use_for_malignancy_interpretation = false }
minor = [
    { name = "Skeletal muscle cell", color = "#CD853F", markers = ["MYOZ1", "MYH1", "MYH2"] },
    { name = "Satellite cell", color = "#DEB887", markers = ["CXCL14", "PAX7", "NCAM1"] }
]

# Add to Skin section:
# { name = "Melanocyte", color = "#8B4513", markers = ["DCT", "MLANA", "TYRP1"] }

# Add to Neural tissue section:
# { name = "Peripheral glia", color = "#708090", markers = ["CDH19", "SOX10", "S100B"] }
```

## Conflict and Specificity Notes

1. **HSPA1A as stress marker**: HSPA1A appears in CD16hi NK (I04), CD8 Trm (CD8T04), and Memory B (B07) subtypes. This is likely a stress/dissociation artifact rather than a true identity marker. Recommend **excluding** HSPA1A from identity markers or flagging as `artifact` context.

2. **NR4A2 in multiple subtypes**: NR4A2 appears in CD4 Tn (CD4T03), CD16lo NK (I06), Naive B (B03), and Memory B (B05). NR4A2 is an immediate-early response gene activated by TCR/BCR signaling. It may indicate activation state rather than identity. Recommend classifying as `state` marker rather than primary identity.

3. **MKI67 in Plasmablast**: MKI67 is a proliferation marker. The existing Plasmablast entry uses TXNDC5, MYDGF. MKI67 marks proliferating plasmablasts but is not identity-specific. Recommend keeping in subtype but noting proliferative state.

4. **IGGH2/IGHA2 in Plasma**: These are isotype-specific constant region genes. They mark IgG2+ and IgA2+ plasma cells respectively. Useful for isotype classification but may not be detectable in all scRNA-seq protocols (3' end bias may miss constant region). Recommend adding as `subtype` with caveat about protocol dependency.

## Summary Statistics

| Category | Count |
|----------|-------|
| Total articles curated | 6 |
| New subtype entries proposed | 8 |
| New tissue-specific state entries proposed | 3 (alveolar macrophage, erythrophagocytic macrophage, gut-resident macrophage) |
| Existing entries recommended for promotion | 20 |
| Tissue-specific additions proposed | 4 |
| Conflict/flagged markers | 4 |
| Cross-tissue validated compartments | 5 |

## Next Steps

1. **Manual review**: Each proposed entry needs manual verification against original source figures/tables
2. **Gene symbol validation**: Confirm all gene symbols are official HGNC symbols
3. **Manager view testing**: Test new entries with `get_marker_manager("human", view="subtype_annotation")`
4. **Promote evidence tiers**: After verification, promote validated entries from `atlas_supported` to `consensus`
