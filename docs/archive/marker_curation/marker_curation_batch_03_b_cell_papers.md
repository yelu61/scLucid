> **⚠️ ARCHIVED / SUPERSEDED**
>
> This batch curation note is kept for provenance only. The live marker
> resource status is tracked in `docs/marker_resources/marker_curation_literature_index.jsonl`,
> `docs/marker_resources/marker_resource_quality_gaps.jsonl`, `docs/marker_resources/marker_curation_candidates.jsonl`,
> and `docs/marker_resources/CURATION.md`. New curation should follow the current
> contract rather than adding new files to this archive.

---

# Marker Curation Batch 03: B Cell Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | A blueprint for tumor-infiltrating B cells across human cancers | 2024 | Science | 10.1126/science.adj4857 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 2 | Pan-cancer single-cell dissection reveals phenotypically distinct B cell subtypes | 2024 | Cell | 10.1016/j.cell.2024.06.038 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 3 | A pan-cancer single-cell RNA-seq atlas of intratumoral B cells | 2024 | Cancer Cell | 10.1016/j.ccell.2024.09.011 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |

---

## Overview: Pan-Cancer B Cell Biology

All 3 papers are 2024 pan-cancer B cell atlases with overlapping but complementary findings. Key shared discovery:

### Core Finding: TAAB / AtM B Cells

| Paper | Term | Key Characteristics |
|-------|------|---------------------|
| Science (Blueprint) | **AtM** (Atypical Memory) | Exhausted, bystander, EF-derived, T-bet+/BATF+, α-KG dependent, immunosuppressive |
| Cell (Dissection) | **TAAB** (Tumor-Associated Atypical B) | High clonal expansion, proliferative, close CD4+ T interactions, ICB predictive |
| Cancer Cell (Atlas) | Atypical B subset | ICB response-associated, B-T crosstalk |

**Assessment**: TAAB (Cell) and AtM (Science) likely describe overlapping populations. AtM is the developmental progenitor in the extrafollicular pathway; TAAB is the tumor-enriched functional state. Both share markers: FCRL4/5, TBX21 (T-bet), CD11c (ITGAX).

### Two Developmental Pathways to ASCs

```
Naive B → Activated B → Memory B
                        ↓
            ┌─────────────────────────┐
            │  Germinal Center (GC)   │  ← Canonical, favorable prognosis
            │  DZ: CXCR4, MKI67       │
            │  LZ: CD86, AICDA        │
            └───────────┬─────────────┘
                        ↓
            GC B → Memory B → Plasma (IGHG)
            (LZ/DZ cycling)

Naive B → Activated B → AtM/TAAB ─────┐
                        (EF pathway)   │  ← Alternative, worse prognosis
                        TBX21+, BATF+  │     ICB resistant
                        FCRL4+, ITGAX+ │
                        ↓              │
            EF-derived Plasma (IGHA/IGHE)
            (Short-lived, low affinity)
```

---

## Resource Tier Classification

### Tier 1a: `marker_registry_human.toml` (Subtype Identity)

**What goes here**: Stable B cell subtype identity markers (3-10 genes), cross-tissue applicable.

- TAAB/AtM B cell identity (new subtype)
- Stress-response memory B (new subtype or state)
- Refined GC/EF B cell identity

### Tier 1b: `marker_registry_human.toml` (State)

**What goes here**: B cell activation/differentiation states.

- EF pathway activation state
- GC reaction state
- B cell exhaustion-like state
- Immunoglobulin class-switch states

### Tier 2: `marker_tumor_human.toml` (Tumor Context)

**What goes here**: Tumor-enriched B cell states and malignancy-associated programs.

- TAAB in tumor microenvironment
- EF-dominant vs GC-dominant TME
- TLS-associated B cell states
- Immunotherapy response-associated B cells

### Tier 3: `genesets_cancer_signatures.json` (Functional Programs)

**What goes here**: Broader gene signatures for module scoring.

- B-T cell crosstalk program
- EF differentiation program
- GC reaction program
- Humoral immunity program
- Immunoglobulin isotype signatures (IgG, IgA, IgE)
- Stress response program in B cells

---

## Detailed Curation by Paper

### Article 1: Science Blueprint (Ma et al., 2024)

**Dataset**: 269 patients, 20 cancer types, 474,718 B cells
**Data types**: scRNA-seq, BCR-seq, scATAC-seq
**Key innovation**: Integrated transcriptome + repertoire + chromatin accessibility

#### Subtype Identity (marker_registry)

**1. Atypical Memory B (AtM) — Refined**

The existing registry has an "Atypical Memory B (FCRL5+)" entry from Batch 01 (Cross-tissue atlas). Science 2024 provides more detailed characterization:

| Feature | Batch 01 (Cross-tissue) | Science 2024 (Pan-cancer) |
|---------|------------------------|---------------------------|
| Name | Atypical Memory B | AtM B cell (Atypical Memory) |
| Markers | FCRL5, ITGB1, CD27 | TBX21, FCRL4, ITGAX, CD27, ITGB1 |
| Context | Cross-tissue normal | Tumor, EF pathway |
| TF | Not specified | TBX21 (T-bet), BATF |
| Function | Not specified | Exhausted, bystander, immunosuppressive |
| Location | Cross-tissue | Immature TLS center → periphery |

**Proposed update** to existing Atypical Memory B entry:

```toml
[[subtype]]
name = "Atypical Memory B (AtM/TAAB)"
color = "#5F9EA0"
markers = ["FCRL4", "TBX21", "ITGAX", "CD27", "ITGB1", "FCRL5"]
negative_markers = ["TCL1A", "IGHM", "RGS13", "CXCR4", "MKI67"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "lineage_restricted"
applies_to = ["B cells", "Memory B"]
evidence_tier = "consensus"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = [
    { title = "A blueprint for tumor-infiltrating B cells across human cancers", year = "2024", doi = "10.1126/science.adj4857" },
    { title = "Pan-cancer single-cell dissection reveals phenotypically distinct B cell subtypes", year = "2024", doi = "10.1016/j.cell.2024.06.038" }
]
notes = """
AtM/TAAB B cells are a distinct atypical memory population. 
Pan-cancer studies (Science 2024, Cell 2024) confirm TBX21 (T-bet), FCRL4, ITGAX (CD11c) as core markers.
Negative: naive (TCL1A, IGHM), GC (RGS13, CXCR4, MKI67) markers.
Tumor-enriched; absent or rare in normal tissues. Present in immature TLS center.
Cross-tissue atlas (Batch 01) confirms FCRL5, ITGB1, CD27 as additional markers.
"""
```

**2. EF-derived Plasma Cell**

```toml
[[subtype]]
name = "EF-derived Plasma Cell"
color = "#8B4513"
markers = ["XBP1", "JCHAIN", "IGHA1", "IGHA2", "IGHE", "MZB1"]
negative_markers = ["IGHG1", "IGHG2", "IGHG3", "IGHG4", "MKI67", "MEF2B"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "tumor_context"
applies_to = ["Plasma", "B cells"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "A blueprint for tumor-infiltrating B cells across human cancers", year = "2024", doi = "10.1126/science.adj4857" }
notes = "Short-lived plasma cells from extrafollicular pathway. IgA/IgE-skewed. Negative for GC markers (IGHG, MKI67, MEF2B). Associated with worse prognosis."
```

**3. GC-derived Plasma Cell**

```toml
[[subtype]]
name = "GC-derived Plasma Cell"
color = "#CD853F"
markers = ["XBP1", "JCHAIN", "IGHG1", "IGHG2", "IGHG3", "MZB1", "SDC1"]
negative_markers = ["IGHA1", "IGHA2", "IGHE", "TBX21", "FCRL4"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "lineage_restricted"
applies_to = ["Plasma", "B cells", "GC B"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "A blueprint for tumor-infiltrating B cells across human cancers", year = "2024", doi = "10.1126/science.adj4857" }
notes = "Long-lived plasma cells from germinal center pathway. IgG-skewed. High-affinity antibodies. Associated with favorable prognosis."
```

#### State Entries (marker_registry)

**1. GC Dark Zone (DZ) State**

```toml
[[state.minor]]
name = "GC Dark Zone"
color = "#6495ED"
markers = ["CXCR4", "MKI67", "AICDA", "BCL6", "STMN1"]
negative_markers = ["CD86", "FCER2", "IL4R"]

[state.minor.metadata]
kind = "state"
category = "immune_activation"
scope = "lineage_restricted"
applies_to = ["B cells", "GC B"]
alias_of = "GC_DZ"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

**2. GC Light Zone (LZ) State**

```toml
[[state.minor]]
name = "GC Light Zone"
color = "#87CEEB"
markers = ["CD86", "FCER2", "IL4R", "CD69", "BCL6"]
negative_markers = ["CXCR4", "MKI67", "AICDA"]

[state.minor.metadata]
kind = "state"
category = "immune_activation"
scope = "lineage_restricted"
applies_to = ["B cells", "GC B"]
alias_of = "GC_LZ"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

**3. B cell exhaustion-like**

```toml
[[state.minor]]
name = "B cell exhaustion-like"
color = "#A9A9A9"
markers = ["FCRL4", "TBX21", "ITGAX", "CD72", "IGHG"]
negative_markers = ["TCL1A", "IGHM", "CD27"]

[state.minor.metadata]
kind = "state"
category = "immune_dysfunction"
scope = "lineage_restricted"
applies_to = ["B cells", "Memory B", "AtM"]
alias_of = "Exhausted_B"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

#### Tumor Context (marker_tumor_human.toml)

**1. EF-dominant Tumor Microenvironment**

```toml
[[malignancy_program]]
name = "EF-dominant B cell response"
color = "#B22222"
markers = ["TBX21", "FCRL4", "ITGAX", "IGHA1", "IGHA2", "IGHE"]
negative_markers = ["BCL6", "AICDA", "CXCR4", "IGHG1", "IGHG2"]

[malignancy_program.metadata]
kind = "functional_program"
granularity = "malignancy_program"
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
source = { title = "A blueprint for tumor-infiltrating B cells across human cancers", year = "2024", doi = "10.1126/science.adj4857" }
notes = "Extrafollicular-dominant B cell response. Associated with worse clinical outcomes and ICB resistance. IgA/E-skewed plasma cells. Negative GC markers (BCL6, AICDA, CXCR4, IGHG)."
```

**2. GC-dominant Tumor Microenvironment**

```toml
[[malignancy_program]]
name = "GC-dominant B cell response"
color = "#228B22"
markers = ["BCL6", "AICDA", "CXCR4", "IGHG1", "IGHG2", "IGHG3"]
negative_markers = ["TBX21", "FCRL4", "ITGAX", "IGHA1", "IGHA2"]

[malignancy_program.metadata]
kind = "functional_program"
granularity = "malignancy_program"
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
source = { title = "A blueprint for tumor-infiltrating B cells across human cancers", year = "2024", doi = "10.1126/science.adj4857" }
notes = "Germinal center-dominant B cell response. Associated with favorable prognosis. IgG-skewed, high-affinity antibodies. Negative EF markers (TBX21, FCRL4, ITGAX, IGHA)."
```

#### Functional Programs (geneset)

```json
"B_cell_EF_differentiation": {
  "genes": ["TBX21", "BATF", "IRF4", "PRDM1", "XBP1", "FCRL4", "ITGAX", "CD27", "ITGB1", "IGHA1", "IGHA2", "IGHE", "JCHAIN", "MZB1"],
  "description": "Extrafollicular B cell differentiation program (AtM → EF plasma). Associated with ICB resistance."
},
"B_cell_GC_differentiation": {
  "genes": ["BCL6", "AICDA", "CXCR4", "CD86", "FCER2", "IL4R", "MEF2B", "LMO2", "MKI67", "IGHG1", "IGHG2", "IGHG3", "JCHAIN", "SDC1"],
  "description": "Germinal center B cell differentiation program (GC B → memory → plasma). Associated with favorable prognosis."
},
"B_cell_TBAF_metabolic": {
  "genes": ["TBX21", "BATF", "MTOR", "RPTOR", "AKT1", "AKT2", "GOT1", "GLUD1", "IDH1", "IDH2", "OGDH"],
  "description": "T-bet/BATF-mTOR metabolic program driving AtM differentiation. Glutamine-αKG dependent."
}
```

---

### Article 2: Cell Dissection (Zhang et al., 2024)

**Dataset**: 649 patients, 19 cancer types
**Key innovations**: TAAB identification, stress-response memory B, IgG-skewness

#### Subtype Identity (marker_registry)

**1. Tumor-Associated Atypical B (TAAB)**

Cell 2024 identifies TAABs as a tumor-enriched atypical B cell population. Overlaps with AtM from Science 2024.

| Feature | Description |
|---------|-------------|
| Markers | FCRL4, TBX21, ITGAX, CD27, ITGB1, FCRL5 |
| Characteristics | High clonal expansion, proliferative, CD4+ T interactions |
| Prognosis | Pan-cancer prognostic potential |
| ICB | Predictive of immunotherapy response |
| Location | Tumor, TLS |

TAAB is essentially the same population as AtM (Science 2024). **Do not create separate entry**. Use unified "Atypical Memory B (AtM/TAAB)" entry (see above).

**2. Stress-Response Memory B**

A new population identified in Cell 2024:

```toml
[[subtype]]
name = "Stress-response Memory B"
color = "#FF6347"
markers = ["HSPA1A", "HSPB1", "DNAJB1", "HSPD1", "HSPE1", "CD27", "ITGB1"]
negative_markers = ["FCRL4", "TBX21", "ITGAX", "TCL1A", "IGHM"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "tumor_context"
applies_to = ["B cells", "Memory B"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Pan-cancer single-cell dissection reveals phenotypically distinct B cell subtypes", year = "2024", doi = "10.1016/j.cell.2024.06.038" }
notes = "Stress-response memory B cells identified in pan-cancer atlas. Express heat shock proteins but distinct from TSTR (T cell stress response). May represent bystander-activated memory B cells. Negative for atypical markers (FCRL4, TBX21, ITGAX)."
```

**Curation note**: Stress-response memory B cells share HSPA1A/HSPB1 with TSTR (T cell stress response) and artifact stress markers. However, they are a distinct B cell subtype, not a dissociation artifact. The CD27+ ITGB1+ memory context distinguishes them.

#### Tumor Context (marker_tumor_human.toml)

**1. TAAB in Tumor**

```toml
[[cancer_state]]
name = "TAAB-enriched tumor"
color = "#8B0000"
markers = ["FCRL4", "TBX21", "ITGAX", "MKI67", "TOP2A"]
negative_markers = ["TCL1A", "IGHM", "BCL6", "AICDA"]

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
source = { title = "Pan-cancer single-cell dissection reveals phenotypically distinct B cell subtypes", year = "2024", doi = "10.1016/j.cell.2024.06.038" }
notes = "TAAB-enriched tumor microenvironment. TAABs are highly clonally expanded, proliferative (MKI67+), and interact with CD4+ T cells. Predictive of immunotherapy response (may indicate active immune response)."
```

**2. IgG-skewed ASC Tumor**

```toml
[[cancer_state]]
name = "IgG-skewed humoral response"
color = "#4682B4"
markers = ["IGHG1", "IGHG2", "IGHG3", "IGHG4", "JCHAIN", "XBP1"]
negative_markers = ["IGHA1", "IGHA2", "IGHE", "IGHD", "IGHM"]

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
source = { title = "Pan-cancer single-cell dissection reveals phenotypically distinct B cell subtypes", year = "2024", doi = "10.1016/j.cell.2024.06.038" }
notes = "IgG-skewed antibody-secreting cell response in tumors. Associated with GC pathway and high-affinity antibodies."
```

#### Functional Programs (geneset)

```json
"B_cell_TAAB_signature": {
  "genes": ["FCRL4", "TBX21", "ITGAX", "CD27", "ITGB1", "MKI67", "TOP2A", "PCNA", "CD72", "IGHG"],
  "description": "TAAB (tumor-associated atypical B) signature. Proliferative, clonally expanded, CD4+ T-interacting."
},
"B_cell_stress_response": {
  "genes": ["HSPA1A", "HSPB1", "HSPA6", "HSP90AA1", "DNAJB1", "HSPD1", "HSPE1", "CD27", "ITGB1"],
  "description": "B cell stress response signature (distinct from T cell TSTR)"
},
"B_IgG_skewness": {
  "genes": ["IGHG1", "IGHG2", "IGHG3", "IGHG4", "AICDA", "BCL6", "CD86"],
  "description": "IgG-skewed humoral immunity program (GC-derived)"
},
"B_IgA_skewness": {
  "genes": ["IGHA1", "IGHA2", "IGHE", "TBX21", "BATF", "FCRL4"],
  "description": "IgA-skewed humoral immunity program (EF-derived)"
}
```

---

### Article 3: Cancer Cell Atlas (Chen et al., 2024)

**Dataset**: Pan-cancer B and plasma cell atlas
**Key innovations**: Subset-specific ICB responses, B-T crosstalk, spatial validation

#### Subtype Identity (marker_registry)

Cancer Cell 2024 confirms B cell subset signatures identified in the other two papers. No entirely new subtypes, but provides:
- Independent validation of TAAB/AtM markers
- B cell subset-specific checkpoint inhibitor responses

#### Tumor Context (marker_tumor_human.toml)

**1. B cell subset-specific ICB response**

Cancer Cell 2024 finds that different B cell subsets have distinct effects on ICB response:

| B Cell Subset | ICB Effect | Mechanism |
|--------------|-----------|-----------|
| GC B / GC-derived plasma | Favorable | High-affinity antibodies, TLS maturity |
| AtM / TAAB | Unfavorable | Immunosuppressive, EF pathway |
| Naive B | Neutral/early | Antigen presentation |
| Memory B | Context-dependent | Recall response |

**2. B-T Crosstalk in Tumor**

Cancer Cell 2024 identifies ligand-receptor pairs between B and T cells:

| Pair | B Cell | T Cell | Function |
|------|--------|--------|----------|
| CD40-CD40LG | B (all) | CD4+ T (activated) | Co-stimulation |
| ICOS-ICOSL | B (GC) | Tfh | GC reaction |
| PDCD1-PDCD1LG2 | B (AtM) | CD8+ T (exhausted) | Inhibition |
| CD80/CD86-CTLA4 | B (activated) | Treg | Inhibition |
| LTA-LTB-TNFRSF1A | B | T | Inflammation |

These are **functional interactions**, not identity markers. Document in notes of relevant entries.

#### Functional Programs (geneset)

```json
"B_T_crosstalk": {
  "genes": ["CD40", "CD40LG", "ICOS", "ICOSL", "PDCD1", "PDCD1LG2", "CD80", "CD86", "CTLA4", "LTA", "LTB", "TNFRSF1A", "CXCL13", "CCL19", "CCL21"],
  "description": "B-T cell crosstalk program (ligand-receptor pairs). Spatially validated in TLS."
},
"TLS_mature": {
  "genes": ["CCL19", "CCL21", "CXCL13", "LTB", "LTA", "CD40", "ICOSL", "BCL6", "AICDA"],
  "description": "Mature tertiary lymphoid structure (TLS) program. Associated with organized GC reactions."
},
"TLS_immature": {
  "genes": ["CXCL13", "CCL19", "CCL21", "TBX21", "FCRL4", "ITGAX", "CD27", "ITGB1"],
  "description": "Immature tertiary lymphoid structure (TLS) program. AtM/TAAB-enriched, EF pathway."
}
```

---

## Unified Registry Entries Summary

### New Subtype Entries (`marker_registry`)

| Entry | Markers | Evidence Tier | Source |
|-------|---------|--------------|--------|
| Atypical Memory B (AtM/TAAB) | FCRL4, TBX21, ITGAX, CD27, ITGB1, FCRL5 | consensus | Science 2024, Cell 2024 |
| EF-derived Plasma | XBP1, JCHAIN, IGHA1, IGHA2, IGHE, MZB1 | atlas_supported | Science 2024 |
| GC-derived Plasma | XBP1, JCHAIN, IGHG1, IGHG2, IGHG3, MZB1, SDC1 | atlas_supported | Science 2024 |
| Stress-response Memory B | HSPA1A, HSPB1, DNAJB1, HSPD1, HSPE1, CD27, ITGB1 | atlas_supported | Cell 2024 |

### New State Entries (`marker_registry`)

| Entry | Markers | Category | Source |
|-------|---------|----------|--------|
| GC Dark Zone | CXCR4, MKI67, AICDA, BCL6, STMN1 | immune_activation | Science 2024 |
| GC Light Zone | CD86, FCER2, IL4R, CD69, BCL6 | immune_activation | Science 2024 |
| B cell exhaustion-like | FCRL4, TBX21, ITGAX, CD72, IGHG | immune_dysfunction | Science 2024 |

### New Tumor Context Entries (`marker_tumor`)

| Entry | Markers | Source |
|-------|---------|--------|
| EF-dominant B cell response | TBX21, FCRL4, ITGAX, IGHA1, IGHA2, IGHE | Science 2024 |
| GC-dominant B cell response | BCL6, AICDA, CXCR4, IGHG1, IGHG2, IGHG3 | Science 2024 |
| TAAB-enriched tumor | FCRL4, TBX21, ITGAX, MKI67, TOP2A | Cell 2024 |
| IgG-skewed humoral response | IGHG1, IGHG2, IGHG3, IGHG4, JCHAIN, XBP1 | Cell 2024 |

### New Geneset Programs

| Program | Genes (approx.) | Source |
|---------|----------------|--------|
| B_cell_EF_differentiation | 15 | Science 2024 |
| B_cell_GC_differentiation | 15 | Science 2024 |
| B_cell_TBAF_metabolic | 11 | Science 2024 |
| B_cell_TAAB_signature | 10 | Cell 2024 |
| B_cell_stress_response | 10 | Cell 2024 |
| B_IgG_skewness | 7 | Cell 2024 |
| B_IgA_skewness | 6 | Cell 2024 |
| B_T_crosstalk | 16 | Cancer Cell 2024 |
| TLS_mature | 10 | Cancer Cell 2024 |
| TLS_immature | 9 | Cancer Cell 2024 |

---

## Curation Notes and Conflicts

### 1. TAAB vs AtM vs Atypical Memory B (FCRL5+)

**Issue**: Three different names for overlapping populations:
- Batch 01 (Cross-tissue atlas): "Atypical Memory B (FCRL5+)"
- Cell 2024: "TAAB" (Tumor-Associated Atypical B)
- Science 2024: "AtM" (Atypical Memory)

**Resolution**:
- Unify under single entry: **"Atypical Memory B (AtM/TAAB)"**
- Include both FCRL4 (primary, from pan-cancer) and FCRL5 (secondary, from cross-tissue) as markers
- Note that FCRL4 is the pan-cancer consensus marker; FCRL5 is additional validation
- Keep both source citations

### 2. HSPA1A in Stress-Response Memory B

**Issue**: HSPA1A appears in:
- `artifact.Stress-high` (QC)
- `state.TSTR` (T cell stress, Batch 02)
- `subtype.Stress-response Memory B` (this batch)

**Resolution**:
- `Stress-high` artifact: dissociation stress, broad across cell types → keep as artifact
- `TSTR` state: tumor microenvironment stress in T cells → keep as state
- `Stress-response Memory B`: B cell-specific stress response, memory context (CD27+ ITGB1+) → add as subtype
- The CD27+ ITGB1+ context distinguishes stress-response memory B from artifact

### 3. EF vs GC Pathway Markers

**Issue**: Some markers are shared between EF and GC pathways (e.g., IRF4, PRDM1 in plasma differentiation).

**Resolution**:
- Use **negative markers** to distinguish:
  - EF: TBX21+, FCRL4+, IGHA/IGHE+, BCL6-, CXCR4-
  - GC: BCL6+, AICDA+, CXCR4+, IGHG+, TBX21-, FCRL4-
- Pathway assignment requires multiple markers, not single genes

### 4. IgG vs IgA Skewness

**Issue**: Immunoglobulin constant region genes (IGHG, IGHA, IGHE, IGHD, IGHM) are used as markers but may not be detectable in all scRNA-seq protocols.

**Resolution**:
- Include in subtype markers but add caveat in notes
- For protocols with 3' bias, constant region genes may be under-detected
- Use additional markers (XBP1, JCHAIN, MZB1, SDC1) as proxy for plasma identity
- Ig class-switch recombination markers (AICDA for GC, TBX21/BATF for EF) are more reliably detected

### 5. B-T Crosstalk as Functional Program

**Issue**: B-T cell interactions (CD40-CD40LG, ICOS-ICOSL, etc.) are important but not identity markers.

**Resolution**:
- Add to `genesets_cancer_signatures.json` as "B_T_crosstalk" program
- Document ligand-receptor pairs in notes of relevant entries
- Do not add to marker_registry (these are interactions, not identity markers)

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New subtype entries | 4 | marker_registry |
| New state entries | 3 | marker_registry |
| New tumor context entries | 4 | marker_tumor |
| New geneset programs | 10 | genesets_cancer_signatures.json |
| Existing entry updates | 1 (Atypical Memory B) | marker_registry |
| Conflict resolutions | 5 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **B cell prognostic value**: EF-dominant tumors have worse outcomes than GC-dominant tumors. scLucid tumor interpretation should consider B cell pathway balance.

2. **ICB response prediction**: TAAB/AtM enrichment may indicate either:
   - Active immune response (Cell 2024: TAAB predicts ICB response)
   - Immunosuppressive microenvironment (Science 2024: AtM dampens T cell responses)
   Context-dependent interpretation needed.

3. **TLS maturity**: Mature TLS (GC-dominant) vs immature TLS (EF-dominant) have different marker profiles. Spatial context matters.

4. **Metabolic regulation**: Glutamine-αKG-T-bet-BATF-mTORC1 axis drives AtM differentiation. Metabolic programs may be useful for functional scoring.
