# Marker Curation Batch 04: NK Cell Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | High-dimensional single-cell analysis of human natural killer cell heterogeneity | 2024 | Nature Immunology | 10.1038/s41590-024-01883-0 | single_cell_atlas | marker_registry + geneset |
| 2 | Pan-cancer profiling of tumor-infiltrating natural killer cells through transcriptional reference mapping | 2024 | Nature Immunology | 10.1038/s41590-024-01884-z | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 3 | A pan-cancer single-cell panorama of human natural killer cells | 2023 | Cell | 10.1016/j.cell.2023.07.034 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |
| 4 | Roles of natural killer cells in immunity to cancer, and applications to immunotherapy | 2023 | Nature Reviews Immunology | 10.1038/s41577-022-00732-1 | review | marker_registry (validation) |
| 5 | Natural killer cells in antitumour adoptive cell immunotherapy | 2022 | Nature Reviews Cancer | 10.1038/s41568-022-00491-0 | review | marker_registry (validation) |

---

## NK Cell Biology Context

### Existing Registry Entries (Current Status)

| Entry | Current Markers | Evidence Tier | Notes |
|-------|----------------|---------------|-------|
| NK cells (compartment) | NCAM1, NKG7, KLRF1, KLRD1, GNLY, FCGR3A, GZMB, PRF1, CD244 | curated_review | Broad identity |
| CD56+CD16- NK | NCAM1, KLRK1, KIT, GZMK | seed | CD56bright ( cytokine-producing) |
| CD56-CD16+ NK | FCGR3A, NCR3, GZMB, PRF1 | seed | CD56dim (cytotoxic) |

### Traditional vs New Classification

| System | Subsets | Markers | Context |
|--------|---------|---------|---------|
| Traditional | CD56bright CD16- | NCAM1, KIT, CCL5, XCL1 | Cytokine-producing, low cytotoxicity |
| Traditional | CD56dim CD16+ | FCGR3A, GZMB, PRF1, NCR3, GNLY | Cytotoxic, high perforin/granzyme |
| Nature Immunology 2024 | NK1, NK2, NK3 | See below | Blood-based, transcriptional + CITE-seq |
| Nature Immunology 2024 | 6 subgroups | See below | Refined NK1/2/3 |
| Pan-cancer atlases | Tumor NK states | See below | Tumor microenvironment |

---

## Paper-by-Paper Curation

### Article 1: High-dimensional NK analysis (Nat Immunol, 2024)

**Methods**: scRNA-seq + CITE-seq of healthy blood, lung, tonsils, intraepithelial lymphocytes, and 22 tumor types
**Key finding**: 3 major subsets (NK1, NK2, NK3) → 6 subgroups; two ontogenetic origins

#### New Subtype Entries (`marker_registry`)

**NK1 (Cytotoxic/Adaptive-like)**

```toml
[[compartment.minor.minor]]
name = "NK1 (cytotoxic)"
color = "#228B22"
markers = ["FCGR3A", "NCR3", "PRF1", "GZMB", "GNLY", "KLRD1", "SPON2", "CCL4", "CCL5"]
negative_markers = ["KIT", "XCL1", "XCL2", "NCAM1hi", "SELL"]

[compartment.minor.minor.metadata]
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
source = { title = "High-dimensional single-cell analysis of human natural killer cell heterogeneity", year = "2024", doi = "10.1038/s41590-024-01883-0" }
notes = "NK1: cytotoxic/adaptive-like subset. FCGR3A+ (CD16+), NCR3+ (NKp30+), high PRF1/GZMB/GNLY. Low NCAM1 (CD56dim). Negative KIT/XCL1/XCL2 distinguishes from NK2. Low SELL indicates effector differentiation. SPON2 marks mature cytotoxic NK cells."
```

**NK2 (Cytokine-producing/Tissue-resident-like)**

```toml
[[compartment.minor.minor]]
name = "NK2 (cytokine)"
color = "#4682B4"
markers = ["NCAM1", "KIT", "XCL1", "XCL2", "CCL5", "KLRC1", "SELL"]
negative_markers = ["FCGR3A", "NCR3", "PRF1", "GZMB", "GNLY", "SPON2"]

[compartment.minor.minor.metadata]
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
source = { title = "High-dimensional single-cell analysis of human natural killer cell heterogeneity", year = "2024", doi = "10.1038/s41590-024-01883-0" }
notes = "NK2: cytokine-producing/tissue-resident-like subset. NCAM1hi (CD56bright), KIT+ (CD117+), XCL1/XCL2+. Produces cytokines but low cytotoxicity. SELL (CD62L) indicates naive/central memory potential. Negative FCGR3A/PRF1/GZMB/GNLY distinguishes from NK1."
```

**NK3 (Immature/Regulatory-like)**

```toml
[[compartment.minor.minor]]
name = "NK3 (immature)"
color = "#DAA520"
markers = ["IL7R", "KIT", "CD27", "NCAM1", "SELL", "TCF7"]
negative_markers = ["FCGR3A", "PRF1", "GZMB", "GNLY", "NCR3", "SPON2"]

[compartment.minor.minor.metadata]
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
source = { title = "High-dimensional single-cell analysis of human natural killer cell heterogeneity", year = "2024", doi = "10.1038/s41590-024-01883-0" }
notes = "NK3: immature/regulatory-like subset. IL7R+ (CD127+), KIT+, CD27+, TCF7+ indicates stem-like/immature potential. High SELL (CD62L) indicates naive state. Low cytotoxicity (negative PRF1/GZMB/GNLY). May represent NK progenitors or ILC1-like cells."
```

#### 6 Subgroups (Refined)

The 6 subgroups represent further subdivision of NK1/2/3. These are **state** or **transitional** rather than stable identities:

| Subgroup | Parent | Markers | Interpretation |
|----------|--------|---------|---------------|
| NK1a | NK1 | FCGR3Ahi, SPON2hi, PRF1hi | Mature cytotoxic |
| NK1b | NK1 | FCGR3Ahi, SPON2lo, GZMBhi | Effector cytotoxic |
| NK1c | NK1 | FCGR3Aint, NCR3hi, KLRC2hi | Adaptive-like (NKG2C+) |
| NK2a | NK2 | NCAM1hi, KIT+, XCL1+ | Cytokine-producing |
| NK2b | NK2 | NCAM1hi, KIT+, CCL5+ | Tissue-resident-like |
| NK3a | NK3 | IL7R+, TCF7+, SELL+ | Immature/naive |

**Curation decision**: The 6 subgroups are too granular for primary annotation. Add as **state** entries with `use_for_global_annotation = false`.

```toml
[[state.minor]]
name = "NK cytotoxic mature"
color = "#006400"
markers = ["FCGR3A", "SPON2", "PRF1", "GZMB", "GNLY"]
negative_markers = ["KIT", "XCL1", "SELL", "IL7R"]

[state.minor.metadata]
kind = "state"
category = "cytotoxicity"
scope = "lineage_restricted"
applies_to = ["NK cells", "NK1"]
alias_of = "NK_mature_cytotoxic"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

```toml
[[state.minor]]
name = "NK adaptive-like"
color = "#32CD32"
markers = ["FCGR3A", "NCR3", "KLRC2", "ZEB2", "FCER1G"]
negative_markers = ["KIT", "SELL", "IL7R"]

[state.minor.metadata]
kind = "state"
category = "immune_activation"
scope = "lineage_restricted"
applies_to = ["NK cells", "NK1"]
alias_of = "NK_adaptive"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
notes = "Adaptive-like NK cells expressing NKG2C (KLRC2), characteristic of CMV/HCMV exposure. ZEB2 and FCER1G mark mature adaptive NK cells."
```

```toml
[[state.minor]]
name = "NK cytokine-producing"
color = "#4169E1"
markers = ["NCAM1", "KIT", "XCL1", "XCL2", "CCL5", "IFNG"]
negative_markers = ["FCGR3A", "PRF1", "GZMB", "SPON2"]

[state.minor.metadata]
kind = "state"
category = "cytokine_production"
scope = "lineage_restricted"
applies_to = ["NK cells", "NK2"]
alias_of = "NK_cytokine"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

```toml
[[state.minor]]
name = "NK immature"
color = "#B8860B"
markers = ["IL7R", "KIT", "CD27", "TCF7", "SELL"]
negative_markers = ["FCGR3A", "PRF1", "GZMB", "GNLY", "SPON2"]

[state.minor.metadata]
kind = "state"
category = "immaturity"
scope = "lineage_restricted"
applies_to = ["NK cells", "NK3"]
alias_of = "NK_immature"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

#### Two Ontogenetic Origins

Nature Immunology 2024 identifies two developmental origins:
1. **Bone marrow-derived**: conventional NK (cNK) → NK1, NK2
2. **Tissue-resident**: thymic or liver-derived → NK3, ILC1-like

This is biological context, not directly a marker distinction. Document in notes.

---

### Article 2: Pan-cancer NK profiling (Nat Immunol, 2024)

**Methods**: Reference map + transfer learning (39 datasets, 7 solid tumors, 427 patients)
**Key finding**: 6 functional states; CD56bright stressed vs CD56dim cytotoxic dichotomy

#### Tumor Context (`marker_tumor`)

**CD56bright Dysfunctional/Stressed NK**

```toml
[[cancer_state]]
name = "NK CD56bright stressed"
color = "#FF6347"
markers = ["NCAM1", "HSPA1A", "HSPB1", "HSPD1", "DNAJB1", "HSPA6", "HSP90AA1"]
negative_markers = ["FCGR3A", "PRF1", "GZMB", "GNLY", "SPON2", "NCR3"]

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
source = { title = "Pan-cancer profiling of tumor-infiltrating natural killer cells through transcriptional reference mapping", year = "2024", doi = "10.1038/s41590-024-01884-z" }
notes = "CD56bright stressed NK cells in tumor microenvironment. Express heat shock proteins (HSPA1A, HSPB1) indicating cellular stress. Susceptible to TME-induced immunosuppression. Negative cytotoxic markers (FCGR3A, PRF1, GZMB) indicate functional impairment. Common across tumor types. Associated with unfavorable prognosis when enriched."
```

**CD56dim Cytotoxic/TME-Resistant NK**

```toml
[[cancer_state]]
name = "NK CD56dim cytotoxic TME-resistant"
color = "#006400"
markers = ["FCGR3A", "PRF1", "GZMB", "GNLY", "SPON2", "NCR3", "KLRD1"]
negative_markers = ["HSPA1A", "HSPB1", "HSPD1", "NCAM1hi", "DNAJB1"]

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
source = { title = "Pan-cancer profiling of tumor-infiltrating natural killer cells through transcriptional reference mapping", year = "2024", doi = "10.1038/s41590-024-01884-z" }
notes = "CD56dim cytotoxic NK cells resistant to TME immunosuppression. High PRF1/GZMB/GNLY/SPON2. Negative heat shock proteins. TME-resistant effector state. Ratio to CD56bright stressed NK predictive of patient outcome in melanoma and osteosarcoma."
```

#### Functional Programs (geneset)

```json
"NK_stress_response": {
  "genes": ["HSPA1A", "HSPB1", "HSPA6", "HSP90AA1", "DNAJB1", "HSPD1", "HSPE1", "HSPA8", "HSPA1B", "HSPH1", "HSP90AB1", "HSPA4", "HSPA5"],
  "description": "NK cell stress response signature. CD56bright stressed state in TME. Associated with immunosuppression susceptibility and unfavorable prognosis."
},
"NK_cytotoxic_effector": {
  "genes": ["PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "KLRD1", "KLRC1", "NCR3", "FCGR3A", "SPON2", "CCL4", "CCL5", "IFNG"],
  "description": "NK cell cytotoxic effector program. CD56dim TME-resistant state. Associated with anti-tumor immunity."
},
"NK_bright_dim_ratio": {
  "genes": ["NCAM1", "FCGR3A", "KIT", "HSPA1A", "HSPB1", "PRF1", "GZMB", "GNLY", "SPON2", "XCL1", "XCL2"],
  "description": "NK cell CD56bright/CD56dim ratio signature. Predictive of patient outcome in melanoma and osteosarcoma."
}
```

---

### Article 3: Pan-cancer NK panorama (Cell, 2023)

**Dataset**: 716 patients, 24 cancer types
**Key finding**: Tumor-associated NK cells; LAMP3+ DC-mediated regulation

#### Tumor Context (`marker_tumor`)

**Tumor-Associated NK (TaNK)**

```toml
[[cancer_state]]
name = "Tumor-associated NK"
color = "#8B0000"
markers = ["HSPA1A", "HSPB1", "HSPD1", "DNAJB1", "VEGFA", "TGFB1", "IL10"]
negative_markers = ["PRF1", "GZMB", "GNLY", "NCR3", "SPON2", "IFNG"]

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
source = { title = "A pan-cancer single-cell panorama of human natural killer cells", year = "2023", doi = "10.1016/j.cell.2023.07.034" }
notes = "Tumor-associated NK cells (TaNK) enriched in tumors across 24 cancer types. Show impaired anti-tumor functions (negative PRF1/GZMB/GNLY/IFNG). Express immunosuppressive cytokines (TGFB1, IL10) and angiogenic factors (VEGFA). Associated with unfavorable prognosis and ICB resistance. Distinct from CD56bright stressed state by additional immunosuppressive markers."
```

**LAMP3+ DC-regulated NK**

Cell 2023 identifies LAMP3+ dendritic cells as regulators of NK cell anti-tumor immunity. This is an **interaction**, not an identity:

```json
"NK_DC_LAMP3_interaction": {
  "genes": ["LAMP3", "CCR7", "FSCN1", "CST7", "KLRC1", "KLRC2", "NCR3", "KLRK1", "KLRD1", "IFNG", "PRF1", "GZMB"],
  "description": "LAMP3+ DC - NK cell interaction program. LAMP3+ DCs regulate NK anti-tumor immunity."
}
```

#### Functional Programs (geneset)

```json
"NK_tumor_associated": {
  "genes": ["HSPA1A", "HSPB1", "HSPD1", "DNAJB1", "VEGFA", "TGFB1", "IL10", "NCAM1", "KIT", "XCL1", "XCL2"],
  "description": "Tumor-associated NK (TaNK) signature. Impaired anti-tumor function, immunosuppressive, poor prognosis."
},
"NK_anti_tumor": {
  "genes": ["PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "IFNG", "CCL4", "CCL5", "NCR3", "KLRK1", "KLRD1", "KLRC2"],
  "description": "NK cell anti-tumor effector program. Cytotoxic, cytokine-producing, anti-tumor immunity."
}
```

---

### Article 4: NK cells in cancer immunity (Nat Rev Immunol, 2023)

**Role**: Review providing NK receptor biology and immunotherapy context.

**Key information for marker registry**:

| Receptor | Gene | Function | Clinical Relevance |
|----------|------|----------|-------------------|
| NKG2D | KLRK1 | Activating, recognizes MIC-A/B, ULBPs | Target for NK cell therapy |
| NKG2A | KLRC1 | Inhibitory, recognizes HLA-E | Monalizumab target |
| KIRs | KIR3DL1, KIR2DL1, etc. | Inhibitory, recognizes HLA | KIR mismatch in allo-HSCT |
| DNAM-1 | CD226 | Activating | Checkpoint in NK exhaustion |
| NCRs | NCR1, NCR3, NCR2 | Activating (NKp46, NKp30, NKp44) | Tumor recognition |
| CD16 | FCGR3A | ADCC | Rituximab, trastuzumab mechanism |
| CD57 | B3GAT1 | Terminal maturation | Marker of mature NK |
| LAG3 | LAG3 | Inhibitory | Dual checkpoint with PD-1 |

**Registry updates**:
- Add KLRK1 (NKG2D) to existing NK cell markers if not already present
- Add KLRC1 (NKG2A) as inhibitory receptor marker
- Add KIR genes to mature NK identity (context-dependent, requires HLA matching)
- Note CD226 (DNAM-1) as activation marker
- Add B3GAT1 (CD57) as terminal maturation marker

**Curation note**: KIR genes are highly polymorphic and HLA-restricted. Use with caution for annotation.

---

### Article 5: NK cells in adoptive immunotherapy (Nat Rev Cancer, 2022)

**Role**: Review on NK cell therapy strategies.

**Key information**:
- CAR-NK development
- Cytokine armouring: IL-15, IL-21 for persistence
- Checkpoint inhibition: anti-NKG2A, anti-KIR, anti-TIGIT
- ADCC enhancement

This is primarily **therapeutic context**, not identity markers. Document in notes of relevant entries.

---

## Updated NK Cell Registry Hierarchy

### Proposed New Structure

```
NK cells (compartment)
├── NK1 (cytotoxic/adaptive-like) [NEW]
│   ├── NK1a: Mature cytotoxic [state]
│   ├── NK1b: Effector cytotoxic [state]
│   └── NK1c: Adaptive-like (NKG2C+) [state]
├── NK2 (cytokine-producing) [NEW]
│   ├── NK2a: Cytokine-producing [state]
│   └── NK2b: Tissue-resident-like [state]
├── NK3 (immature/ILC1-like) [NEW]
│   └── NK3a: Immature [state]
├── CD56bright stressed (tumor) [NEW - tumor]
├── CD56dim cytotoxic TME-resistant (tumor) [NEW - tumor]
└── Tumor-associated NK (TaNK) [NEW - tumor]
```

### Existing Entries: Updates Recommended

| Existing Entry | Current Markers | Proposed Update | Source |
|---------------|----------------|-----------------|--------|
| NK cells | NCAM1, NKG7, KLRF1, KLRD1, GNLY, FCGR3A, GZMB, PRF1, CD244 | Add KLRK1 (NKG2D), KLRC1 (NKG2A), KIR3DL1, CD226 (DNAM-1) | Nat Rev Immunol 2023 |
| CD56+CD16- NK | NCAM1, KLRK1, KIT, GZMK | Rename to "NK2 (cytokine)"; update markers | Nat Immunol 2024 |
| CD56-CD16+ NK | FCGR3A, NCR3, GZMB, PRF1 | Rename to "NK1 (cytotoxic)"; add SPON2, GNLY | Nat Immunol 2024 |

---

## Functional Programs (geneset)

### New Geneset Entries

```json
{
  "_categories": {
    "NK_Cell_Function": [
      "NK_cytotoxic_effector",
      "NK_cytokine_production",
      "NK_adaptive",
      "NK_stress_response",
      "NK_tumor_associated",
      "NK_anti_tumor",
      "NK_bright_dim_ratio"
    ],
    "NK_Receptor_Signaling": [
      "NK_activating_receptors",
      "NK_inhibitory_receptors"
    ]
  },
  "NK_cytotoxic_effector": {
    "genes": ["PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "KLRD1", "KLRC1", "NCR3", "FCGR3A", "SPON2", "CCL4", "CCL5", "IFNG"],
    "description": "NK cell cytotoxic effector program (CD56dim, TME-resistant)"
  },
  "NK_cytokine_production": {
    "genes": ["NCAM1", "KIT", "XCL1", "XCL2", "CCL5", "IFNG", "TNF", "IL10", "CSF1", "CSF2"],
    "description": "NK cell cytokine production program (CD56bright)"
  },
  "NK_adaptive": {
    "genes": ["FCGR3A", "NCR3", "KLRC2", "ZEB2", "FCER1G", "PRF1", "GZMB", "SPON2", "IFNG"],
    "description": "Adaptive-like NK cell program (NKG2C+, CMV/HCMV-exposed)"
  },
  "NK_stress_response": {
    "genes": ["HSPA1A", "HSPB1", "HSPA6", "HSP90AA1", "DNAJB1", "HSPD1", "HSPE1", "HSPA8", "HSPA1B", "HSPH1", "HSP90AB1", "HSPA4", "HSPA5"],
    "description": "NK cell stress response signature. CD56bright stressed state in TME."
  },
  "NK_tumor_associated": {
    "genes": ["HSPA1A", "HSPB1", "HSPD1", "DNAJB1", "VEGFA", "TGFB1", "IL10", "NCAM1", "KIT", "XCL1", "XCL2"],
    "description": "Tumor-associated NK (TaNK) signature. Impaired anti-tumor function."
  },
  "NK_anti_tumor": {
    "genes": ["PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "IFNG", "CCL4", "CCL5", "NCR3", "KLRK1", "KLRD1", "KLRC2"],
    "description": "NK cell anti-tumor effector program"
  },
  "NK_activating_receptors": {
    "genes": ["KLRK1", "NCR3", "NCR1", "NCR2", "CD226", "KLRC2", "FCGR3A"],
    "description": "NK cell activating receptor signature"
  },
  "NK_inhibitory_receptors": {
    "genes": ["KLRC1", "KIR3DL1", "KIR2DL1", "KIR2DL3", "KIR2DL4", "KIR3DL2", "LAG3", "TIGIT", "PDCD1"],
    "description": "NK cell inhibitory receptor signature"
  }
}
```

---

## Curation Notes and Conflicts

### 1. NK1/NK2/NK3 vs CD56dim/CD56bright

**Issue**: Two classification systems:
- Traditional: CD56bright (NCAM1hi, FCGR3A-) vs CD56dim (NCAM1lo, FCGR3A+)
- New (Nat Immunol 2024): NK1 (cytotoxic), NK2 (cytokine), NK3 (immature)

**Resolution**:
- NK1 ≈ CD56dim (cytotoxic): FCGR3A+, PRF1+, GZMB+
- NK2 ≈ CD56bright (cytokine): NCAM1hi, KIT+, XCL1+
- NK3 = novel immature subset: IL7R+, TCF7+, SELL+
- Map traditional to new:
  - Rename "CD56+CD16- NK" → "NK2 (cytokine)"
  - Rename "CD56-CD16+ NK" → "NK1 (cytotoxic)"
  - Add "NK3 (immature)" as new

### 2. HSPA1A in NK Stress Response

**Issue**: HSPA1A appears in:
- `artifact.Stress-high` (QC)
- `state.TSTR` (T cell stress, Batch 02)
- `cancer_state.NK CD56bright stressed` (this batch)
- `cancer_state.Tumor-associated NK` (this batch)

**Resolution**:
- `Stress-high` artifact: dissociation stress → keep as artifact
- `TSTR` state: T cell tumor stress → keep as state
- `NK CD56bright stressed`: NK-specific TME stress → add as cancer_state
- `TaNK`: NK tumor-associated (includes stress + immunosuppression) → add as cancer_state
- All are context-specific and do not conflict

### 3. NKG2C (KLRC2) as Adaptive NK Marker

**Issue**: KLRC2 (NKG2C) marks adaptive-like NK cells (CMV-exposed) but is also expressed on some cNK cells.

**Resolution**:
- Use KLRC2 in combination with ZEB2 and FCER1G for adaptive NK identity
- Add as state (not primary identity) since CMV exposure is not universal
- Document in notes: "KLRC2 requires CMV/HCMV exposure context"

### 4. KIR Polymorphism

**Issue**: KIR genes (KIR3DL1, KIR2DL1, etc.) are highly polymorphic.

**Resolution**:
- Include KIR3DL1 as representative inhibitory receptor in broad NK entry
- Add note: "KIR expression is highly polymorphic and HLA-restricted. Use for population-level analysis, not single-cell identity annotation."
- Do not add individual KIRs as primary identity markers

### 5. ILC1 Overlap with NK3

**Issue**: NK3 (IL7R+, TCF7+, SELL+) overlaps with ILC1 (ID2+, IL7R+, RORC-).

**Resolution**:
- NK3 and ILC1 are closely related but distinct
- ILC1: ID2+, IL7R+, GATA3-/RORC- (tissue-resident, non-cytotoxic)
- NK3: IL7R+, TCF7+, SELL+, NCAM1+ (circulating immature NK)
- Use NCAM1+ and ID2- to distinguish NK3 from ILC1

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New subtype entries | 3 (NK1, NK2, NK3) | marker_registry |
| New state entries | 4 (NK mature, adaptive, cytokine, immature) | marker_registry |
| New tumor context entries | 3 (CD56bright stressed, CD56dim cytotoxic, TaNK) | marker_tumor |
| Existing entry renames | 2 (CD56+CD16- → NK2, CD56-CD16+ → NK1) | marker_registry |
| New geneset programs | 8 | genesets_cancer_signatures.json |
| Conflict resolutions | 5 | cross-reference |

---

## Key Biological Insights for scLucid Workflows

1. **NK cell prognostic value**: CD56bright stressed / CD56dim cytotoxic ratio predicts outcome in melanoma and osteosarcoma. scLucid tumor interpretation should calculate this ratio.

2. **NK cell therapy design**: 
   - Target CD56dim cytotoxic NK for adoptive transfer
   - Avoid CD56bright stressed NK (dysfunctional, immunosuppressed)
   - NKG2A blockade (monalizumab) can enhance NK function

3. **Tumor microenvironment effects**: TME induces NK stress (HSPA1A+) and dysfunction (VEGFA+, TGFB1+, IL10+). scLucid should flag these as state changes, not identity changes.

4. **Ontogenetic diversity**: NK1/NK2 (bone marrow-derived) vs NK3 (tissue-resident/ILC1-like). Different developmental origins may respond differently to therapies.
