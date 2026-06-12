# Marker Curation Batch 02: T Cell Lineage Papers

## Source Articles

| # | Title | Year | Journal | DOI | Source Type | Resource Tier |
|---|-------|------|---------|-----|-------------|---------------|
| 1 | Guidelines for T cell nomenclature | 2026 | Nature Reviews Immunology | 10.1038/s41577-025-01238-2 | review | naming_reference |
| 2 | Asian diversity in human immune cells | 2025 | Cell | 10.1016/j.cell.2025.02.017 | single_cell_atlas | marker_registry + geneset |
| 3 | An automatic annotation tool and reference database for T cell subtypes and states at single-cell resolution | 2025 | Science Bulletin | 10.1016/j.scib.2025.02.043 | single_cell_atlas | marker_registry (primary) |
| 4 | Integrative mapping of human CD8+ T cells in inflammation and cancer | 2025 | Nature Methods | 10.1038/s41592-024-02530-0 | pan_cancer_atlas | marker_registry + geneset |
| 5 | TCellSI: A novel method for T cell state assessment and its applications in immune environment prediction | 2024 | iMeta | 10.1002/imt2.231 | computational_tool | **geneset only** |
| 6 | Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance | 2023 | Nature Medicine | 10.1038/s41591-023-02371-y | pan_cancer_atlas | marker_registry (state) + geneset |
| 7 | Pan-cancer single-cell landscape of tumor-infiltrating T cells | 2021 | Science | 10.1126/science.abe6474 | pan_cancer_atlas | marker_registry + marker_tumor + geneset |

---

## Resource Tier Classification Principles

Before extracting markers, we apply scLucid's tier classification to each paper:

### Tier 1: `marker_registry_human.toml` (Identity Annotation)
**What goes here**: Concise marker sets (3-12 genes) that define cell identity at compartment, lineage, or subtype level. These are used for **global annotation** (`use_for_global_annotation = true`).

- **Compartment/Lineage entries** (3-12 markers): Broad T cell identity markers
- **Subtype entries** (3-10 markers): T cell subtypes with stable identity across contexts
- **State entries** (3-15 markers, `use_for_global_annotation = false`): Transient cell states used for **state annotation** only

**Applicable papers**:
- Article 3 (STCAT): 33 subtypes → identity markers for stable subtypes
- Article 4 (CD8 atlas): CD8+ subtypes → identity markers
- Article 2 (Asian diversity): Population-validated identity markers
- Article 1 (Nomenclature): Standardized naming → informs registry naming

### Tier 2: `marker_tumor_human.toml` (Tumor Context)
**What goes here**: Tumor-specific T cell states, malignancy-associated programs, and cancer-context markers. These require tumor context to interpret.

**Applicable papers**:
- Article 7 (Pan-cancer T cell landscape, Science 2021): Tumor-infiltrating T cell composition
- Article 6 (TSTR, Nat Med 2023): Stress response in tumor microenvironment
- Article 4 (CD8 atlas): Exhaustion subtypes in cancer vs inflammation

### Tier 3: `genesets_cancer_signatures.json` (Functional Scoring)
**What goes here**: Broader gene signatures (5-50+ genes) used for **module scoring** and functional interpretation, not direct cell identity annotation.

- **Functional programs**: Pathway-level signatures for scoring cell states
- **Malignancy programs**: Tumor-context functional modules

**Applicable papers**:
- Article 5 (TCellSI): 8 T cell state scoring programs → **geneset only**
- Article 6 (TSTR): Heat shock stress response program → **geneset**
- Article 7 (Pan-cancer): Exhaustion trajectory programs → **geneset**
- Article 4 (CD8 atlas): Cytotoxicity, memory, exhaustion gene modules → **geneset**

---

## Paper-by-Paper Curation

### Article 1: Guidelines for T cell nomenclature (Nat Rev Immunol, 2026)

**Role**: Naming standardization reference. Does not provide new markers but standardizes existing terminology.

**Key principles for scLucid naming**:
1. **Modular nomenclature**: Instead of fixed subset labels, describe individual biological properties
   - Example: "CD8+ CCR7+ SELL+ TCF7+" rather than "CD8+ Tn"
   - Implication: scLucid's hierarchical naming (CD8+ T → CD8+ Tn) is consistent but should support modular descriptors

2. **Standardized definitions** for existing subsets:
   - Naive T: CCR7+ SELL+ (or CD45RA+)
   - Central memory: CCR7+ SELL- (or CD45RO+)
   - Effector memory: CCR7- SELL- CD45RO+
   - Effector memory RA+: CCR7- SELL- CD45RA+
   - Tissue-resident memory: CD69+ CD103+ (or ITGA1+)
   - Treg: FOXP3+ CD25+ (IL2RA+)
   - Th1: TBX21+ IFNG+
   - Th2: GATA3+ IL4+
   - Th17: RORC+ IL17A+
   - Tfh: BCL6+ CXCR5+

3. **Recommendation**: scLucid's existing T cell hierarchy is largely consistent. Minor alignment needed:
   - Add `CD45RA`/`PTPRC` as alternative naive marker alongside `SELL`
   - Add `CD45RO` reference in notes (not RNA-detectable, protein-level)
   - Confirm modular naming support in annotation output

**Registry impact**: No new markers; naming alignment notes only.

---

### Article 2: Asian diversity in human immune cells (Cell, 2025)

**Role**: Population-specific validation of immune cell markers in Asian cohorts.

**Scope**: Single-cell atlas of Asian donors. Validates cross-population stability of T cell markers.

**Expected contributions** (abstract not available in Zotero; based on title and Cell scope):
- Validation that existing T cell markers are stable across ancestral backgrounds
- Potential population-specific expression differences in T cell activation/exhaustion states
- May reveal population-enriched T cell subtypes or states

**Curation approach**:
- Use as **validation source** for existing entries rather than new marker discovery
- If population-specific subtypes are identified, add with `scope = "population_specific"` and `applies_to = ["Asian"]`
- Add source citation to existing entries that this atlas validates

**Registry entries affected**:
- All existing T cell entries → add validation note: "Supported by Asian diversity atlas (Cell 2025)"
- Any new population-specific findings → `subtype` with restricted scope

---

### Article 3: STCAT — T cell reference and annotation tool (Sci Bull, 2025)

**Role**: Primary source for T cell **subtype identity markers**. Largest human T cell reference (1.35M cells, 33 subtypes, 68 categories).

**Key findings from abstract**:
- 33 subtypes classified across 35 conditions and 16 tissues
- 68 categories stratified by subtype + state
- Th17 enriched in late-stage lung cancer
- MAIT prevalent in milder-stage COVID-19
- Treg cytotoxicity decreased in post-treatment ovarian cancer
- CD4+ Treg enriched in tumor samples
- CD8+ naive-related cells abundant in healthy individuals

**Resource tier mapping**:

#### Tier 1a: Subtype Identity Markers (marker_registry)

The 33 subtypes from STCAT include known subtypes. We extract those that provide **new or validated identity markers**:

**CD4+ T subtypes** (from STCAT reference):

| Subtype | Proposed Markers | Evidence | Registry Action |
|---------|-----------------|----------|-----------------|
| CD4+ Tn (naive) | TCF7, SELL, CCR7, LEF1, MAL | Confirmed across 16 tissues | Validate existing |
| CD4+ Tcm (central memory) | CCR7, IL7R, S100A4 | Confirmed | Validate existing |
| CD4+ Tem (effector memory) | GZMK, KLRB1, CCL5 | Confirmed | Validate existing |
| CD4+ Temra | CX3CR1, GZMH, TBX21 | Confirmed | Validate existing |
| CD4+ Treg | FOXP3, IL2RA, CTLA4, TNFRSF9, RTKN2 | Tumor-enriched | Validate + add tumor context |
| CD4+ Th17 | RORC, IL17A, IL17F, CCR6 | Late-stage lung cancer enriched | Validate + add disease context |
| CD4+ Tfh | CXCR5, BCL6, IL21, MAF | Confirmed | Validate existing |
| CD4+ Th1 | TBX21, IFNG, CXCR3 | Confirmed | Validate existing |
| CD4+ Th2 | GATA3, IL4, IL5, IL13 | Confirmed | Validate existing |
| CD4+ Th9 | IL9, PU.1 (SPI1) | Less common | Add as new subtype? |
| CD4+ Tfr (follicular regulatory) | FOXP3, CXCR5, BCL6 | Confirmed | Add as new subtype |
| CD4+ Treg cytotoxic | FOXP3, GZMB, PRF1 | Decreased post-treatment ovarian cancer | Add as state variant |

**CD8+ T subtypes** (from STCAT reference):

| Subtype | Proposed Markers | Evidence | Registry Action |
|---------|-----------------|----------|-----------------|
| CD8+ Tn | TCF7, SELL, CCR7, LEF1 | Confirmed | Validate existing |
| CD8+ Tcm | CCR7, IL7R, S100A4 | Confirmed | Validate existing |
| CD8+ Tem | GZMK, CCL3, CXCR4, KLRG1 | Confirmed | Validate existing |
| CD8+ Temra | CX3CR1, GZMH, TBX21, KLRD1 | Confirmed | Validate existing |
| CD8+ Trm | ZNF683, CD69, ITGA1 (CD103), CD52 | Confirmed | Validate existing |
| CD8+ Tex (exhausted) | PDCD1, TOX, LAG3, HAVCR2 | Confirmed | Validate existing |
| CD8+ Tex-progenitor | TCF7, SELL, PDCD1, TOX | Progenitor exhaustion | Add as new state |
| CD8+ Tex-terminal | PDCD1, HAVCR2, TOX, ENTPD1, CXCL13 | Terminal exhaustion | Validate existing |
| MAIT | TRAV1-2, SLC4A10, KLRB1, ZBTB16 | COVID-19 mild stage enriched | Validate existing |
| Gamma-delta T | TRDC, TRGC1, TRGC2, NKG7 | Confirmed | Validate existing |
| NK-like T | KLRD1, NKG7, TYROBP, CD160 | Confirmed | Add as new subtype |

**New subtype proposals from STCAT** (not in existing registry):

1. **Tfr (Follicular Regulatory T)** — CD4+ FOXP3+ CXCR5+ BCL6+
   - Distinct from Treg (non-follicular) and Tfh (non-regulatory)
   - Markers: FOXP3, CXCR5, BCL6, IL2RA, CTLA4
   - Negative: IL21 (low), IL17A (low)
   - Registry: `subtype`, compartment=immune, lineage=lymphoid
   - Scope: lineage_restricted (applies_to: ["T cells", "B cell zones"])

2. **Treg cytotoxic** — FOXP3+ GZMB+ PRF1+
   - Cytotoxic Treg subset with suppressive + killing functions
   - Markers: FOXP3, GZMB, PRF1, IL2RA, CTLA4
   - Negative: IL10 (may vary)
   - Registry: `state` (not primary identity; cytotoxicity is a state)
   - Note: Post-treatment ovarian cancer shows decreased cytotoxicity

3. **CD8+ Tex-progenitor** — TCF7+ SELL+ PDCD1+ TOX+
   - Progenitor exhausted T cells with stem-like potential
   - Markers: TCF7, SELL, PDCD1, TOX, IL7R
   - Negative: GZMB, PRF1 (low cytotoxicity)
   - Registry: `state`, applies_to: ["CD8+ T"]
   - Evidence from Articles 3, 4, 7

#### Tier 1b: State Markers (marker_registry, state section)

States identified by STCAT (from 68 categories):

| State | Markers | Evidence | Registry Section |
|-------|---------|----------|------------------|
| Quiescence | TCF7, SELL, CCR7, LEF1 | Naive-like, low activation | state → Immune Functional States |
| Activation | CD69, IL2RA, FOS, JUNB, CD38 | Early activation | state → Immune Functional States (existing) |
| Cytotoxicity | GZMA, GZMB, PRF1, NKG7, GNLY | Effector function | state → Immune Functional States |
| Proliferation | MKI67, TOP2A, PCNA, CCNB1 | Cycling | state → Cell Cycle and Growth States (existing) |
| Exhaustion | PDCD1, LAG3, HAVCR2, TIGIT, TOX | Dysfunctional | state → Immune Functional States (existing) |
| Senescence | CDKN1A, CDKN2A, GLB1, TP53 | Aging/dysfunction | state → Metabolic and Plasticity States |
| Helper function | IL2, IL4, IL17A, IL21, IFNG | Cytokine production | state → Immune Functional States |

**Note**: STCAT's 68 categories combine subtype + state. The pure state markers should be extracted and added to `marker_registry` state section.

#### Tier 3: Functional Programs (geneset)

STCAT provides reference expression profiles for all 33 subtypes. These can be used to derive:
- Subtype-specific gene signatures for scoring
- State transition gene modules

**Not directly extracted** into geneset (STCAT is primarily an identity reference, not a functional scoring resource).

---

### Article 4: Integrative mapping of human CD8+ T cells (Nat Methods, 2025)

**Role**: scAtlasVAE integration of 1.15M CD8+ T cells across 42 diseases. Primary source for **CD8+ subtype refinement** and **exhaustion subtypes**.

**Key findings**:
- Three distinct exhausted T cell subtypes
- TCR clonal expansion connects cell subtypes
- Autoimmune vs cancer inflammation patterns
- Automatic annotation framework

**Resource tier mapping**:

#### Tier 1a: Subtype Identity (marker_registry)

**Three exhaustion subtypes** (from scAtlasVAE):

| Exhaustion Subtype | Proposed Markers | Biological Role | Registry Action |
|-------------------|-----------------|-----------------|-----------------|
| Tex-progenitor (Tex-p) | TCF7, SELL, IL7R, PDCD1, TOX | Stem-like, self-renewal, respond to ICB | Add to state |
| Tex-transitory (Tex-t) | GZMK, CXCR4, CXCL13, PDCD1 | Intermediate, transitioning | Add to state |
| Tex-terminal (Tex-term) | HAVCR2, ENTPD1, CXCL13, TOX, PRDM1 | Terminal, non-responsive to ICB | Validate existing state |

**CD8+ inflammation subtypes**:

| Context | Subtype | Markers | Registry Action |
|---------|---------|---------|-----------------|
| Autoimmune | Teff-autoimmune | TBX21, IFNG, CXCR3, CCL5 | Add to state |
| Autoimmune | Tmem-autoimmune | IL7R, S100A4, CD27, ANXA1 | Add to state |
| irAE | Tex-irAE | PDCD1, CTLA4, ICOS, OX40 | Add to state |
| Cancer | Tex-tumor | PDCD1, HAVCR2, LAG3, TOX | Validate existing |

#### Tier 1b: State (marker_registry)

**TCR clonal expansion-associated states**:
- Expanded clones vs non-expanded show different marker profiles
- Expanded: GZMB, PRF1, CXCL13 (high cytotoxic/exhausted)
- Non-expanded: TCF7, SELL, CCR7 (naive/memory)

These are **state annotations**, not identity markers. Add to `state` section with `use_for_state_annotation = true`.

#### Tier 3: Functional Programs (geneset)

**Exhaustion trajectory gene modules** (from scAtlasVAE):

```json
"CD8_Tex_trajectory": {
  "genes": ["TCF7", "SELL", "IL7R", "GZMK", "CXCR4", "CXCL13", "PDCD1", "TOX", "HAVCR2", "ENTPD1", "PRDM1", "LAG3", "TIGIT", "CTLA4"],
  "description": "CD8+ T cell exhaustion trajectory from progenitor to terminal state"
}
```

**Cytotoxicity program**:
```json
"CD8_cytotoxic": {
  "genes": ["GZMA", "GZMB", "PRF1", "NKG7", "GNLY", "CST7", "IFNG", "TBX21"],
  "description": "CD8+ T cell cytotoxic effector program"
}
```

**Memory formation program**:
```json
"CD8_memory": {
  "genes": ["IL7R", "S100A4", "CD27", "KLRG1", "CD28", "BCL2", "TCF7"],
  "description": "CD8+ T cell memory formation and maintenance program"
}
```

---

### Article 5: TCellSI (iMeta, 2024)

**Role**: T cell state scoring tool. **Geneset-only resource** — not for identity annotation.

**Key feature**: 8 T cell state scoring programs:

| State | Gene Count | Description | Resource |
|-------|-----------|-------------|----------|
| Quiescence | ~15-30 | Resting T cell state | geneset |
| Regulating | ~15-30 | Treg-mediated suppression | geneset |
| Proliferation | ~10-20 | Cell cycling | geneset (overlap with existing) |
| Helper | ~15-30 | Cytokine production | geneset |
| Cytotoxicity | ~15-30 | Killing function | geneset |
| Progenitor exhaustion | ~15-30 | Stem-like exhausted | geneset |
| Terminal exhaustion | ~15-30 | Dysfunctional terminal | geneset |
| Senescence | ~15-30 | Aging dysfunction | geneset |

**Curation decision**: TCellSI is explicitly a **scoring tool**, not a cell identity reference. Its marker gene sets are designed for bulk/pseudo-bulk scoring, not single-cell identity annotation.

**Recommended action**: Add all 8 programs to `genesets_cancer_signatures.json` under a new category `"T_Cell_States"`.

```json
"T_Cell_States": {
  "TCellSI_quiescence": {
    "genes": ["TCF7", "SELL", "CCR7", "LEF1", "MAL", "CD27", "IL7R"],
    "description": "TCellSI quiescence state score (resting/naive-like)"
  },
  "TCellSI_regulating": {
    "genes": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF9", "RTKN2", "IL10", "TGFB1"],
    "description": "TCellSI regulating state score (Treg-mediated suppression)"
  },
  "TCellSI_proliferation": {
    "genes": ["MKI67", "TOP2A", "PCNA", "CCNB1", "CCNB2", "AURKB", "CDK1"],
    "description": "TCellSI proliferation state score"
  },
  "TCellSI_helper": {
    "genes": ["IL2", "IL4", "IL5", "IL13", "IL17A", "IL21", "IFNG", "TNF"],
    "description": "TCellSI helper function state score (cytokine production)"
  },
  "TCellSI_cytotoxicity": {
    "genes": ["GZMA", "GZMB", "PRF1", "NKG7", "GNLY", "CST7", "IFNG", "TBX21"],
    "description": "TCellSI cytotoxicity state score"
  },
  "TCellSI_progenitor_exhaustion": {
    "genes": ["TCF7", "SELL", "IL7R", "PDCD1", "TOX", "LAG3", "CXCR5"],
    "description": "TCellSI progenitor exhaustion state score (stem-like exhausted, ICB-responsive)"
  },
  "TCellSI_terminal_exhaustion": {
    "genes": ["HAVCR2", "ENTPD1", "CXCL13", "TOX", "PRDM1", "LAG3", "TIGIT", "CTLA4"],
    "description": "TCellSI terminal exhaustion state score (dysfunctional, ICB-resistant)"
  },
  "TCellSI_senescence": {
    "genes": ["CDKN1A", "CDKN2A", "GLB1", "TP53", "LMNB1", "HLA-DRA", "CD69"],
    "description": "TCellSI senescence state score"
  }
}
```

**Important**: These genesets are for **module scoring** (e.g., via `scanpy.tl.score_genes` or `scLucid.utils.Manager` program_scoring view). They are **NOT** for direct cell type annotation. The gene lists above are representative; actual TCellSI gene sets should be extracted from the publication or R package.

---

### Article 6: TSTR stress response state (Nat Med, 2023)

**Role**: Discovery of TSTR (T cell stress response) state in tumor-infiltrating T cells. Links stress response to immunotherapy resistance.

**Key findings**:
- TSTR: unique stress response state characterized by **heat shock gene expression**
- Detectable in situ across cancer types
- Located in tertiary lymphoid structures (TLS) or lymphocyte aggregates
- Upregulated after ICB in nonresponsive tumors
- Correlates with immunotherapy resistance

**Resource tier mapping**:

#### Tier 1b: State Marker (marker_registry)

**TSTR state** — add to `marker_registry` state section:

```toml
[[state.minor]]
name = "T cell stress response (TSTR)"
color = "#FF6347"
markers = ["HSPA1A", "HSPB1", "HSPA6", "HSP90AA1", "DNAJB1", "HSPD1", "HSPE1"]

[state.minor.metadata]
kind = "state"
category = "stress_response"
scope = "lineage_restricted"
applies_to = ["T cells", "CD4+ T", "CD8+ T"]
alias_of = "TSTR"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

**Curation note**: TSTR is a **state**, not an identity. It can occur in any T cell subtype under stress conditions (e.g., ICB treatment, tumor microenvironment). `use_for_global_annotation = false` is correct.

**Overlap with existing**: The existing `Stress-high` artifact entry uses HSPA1A, HSPB1, HSPA6, HSP90AA1. TSTR overlaps but is biologically meaningful in tumor context. Recommend:
- Keep `Stress-high` as **artifact/QC** (dissociation stress)
- Add TSTR as **biological state** (tumor microenvironment stress)
- Document distinction in notes

#### Tier 2: Tumor Context (marker_tumor_human.toml)

**TSTR in tumor microenvironment**:

```toml
[[cancer_state]]
name = "TSTR in tumor"
color = "#FF4500"
markers = ["HSPA1A", "HSPB1", "HSPA6", "HSP90AA1", "DNAJB1", "HSPD1"]
negative_markers = ["GZMB", "PRF1", "IFNG"]

[cancer_state.metadata]
kind = "cancer_state"
granularity = "cancer_state"
compartment = "immune"
lineage = "T cells"
cancer_type = ["all"]
scope = "tumor_context"
applies_to = ["T cells", "CD4+ T", "CD8+ T"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Pan-cancer T cell atlas links a cellular stress response state to immunotherapy resistance", year = "2023", doi = "10.1038/s41591-023-02371-y" }
notes = "TSTR cells are found in TLS/lymphocyte aggregates in tumor beds. Associated with ICB nonresponse. Negative cytotoxic markers indicate functional suppression."
```

#### Tier 3: Functional Program (geneset)

**TSTR heat shock program** (for module scoring):

```json
"TSTR_heat_shock": {
  "genes": ["HSPA1A", "HSPB1", "HSPA6", "HSP90AA1", "DNAJB1", "HSPD1", "HSPE1", "HSPA8", "HSPA1B", "HSPH1", "HSP90AB1", "HSPA4", "HSPA5"],
  "description": "T cell stress response (TSTR) heat shock program. Enriched in ICB nonresponsive tumors."
}
```

**ICB response program**:

```json
"T_cell_ICB_response": {
  "genes": ["TCF7", "SELL", "IL7R", "GZMK", "IFNG", "CXCL13", "PRF1", "GZMB"],
  "description": "T cell program associated with immune checkpoint blockade response"
}
```

```json
"T_cell_ICB_resistance": {
  "genes": ["HSPA1A", "HSPB1", "HSPA6", "HSP90AA1", "DNAJB1", "HSPD1", "TOX", "HAVCR2", "ENTPD1"],
  "description": "T cell program associated with immune checkpoint blockade resistance (TSTR + terminal exhaustion)"
}
```

---

### Article 7: Pan-cancer T cell landscape (Science, 2021)

**Role**: Foundational pan-cancer T cell atlas (316 donors, 21 cancer types). Primary source for **exhaustion dynamics** and **tumor-type-specific T cell composition**.

**Key findings**:
- Multiple exhaustion state-transition paths
- Different paths preferred among tumor types
- Predysfunctional → dysfunctional progression
- T cell composition correlates with mutation burden
- T cell composition alone can classify cancer patients

**Resource tier mapping**:

#### Tier 1a: Subtype Identity (marker_registry)

**Predysfunctional T cells** (pre-exhaustion):

```toml
[[subtype]]
name = "CD8+ predysfunctional"
color = "#DAA520"
markers = ["GZMK", "CXCR4", "PDCD1", "TOX", "LAG3"]
negative_markers = ["HAVCR2", "ENTPD1", "CXCL13"]

[subtype.metadata]
kind = "cell_type"
granularity = "subtype"
compartment = "immune"
lineage = "lymphoid"
scope = "lineage_restricted"
applies_to = ["CD8+ T"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = true
use_for_state_annotation = false
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Pan-cancer single-cell landscape of tumor-infiltrating T cells", year = "2021", doi = "10.1126/science.abe6474" }
notes = "Pre-exhaustion CD8+ T cells showing early PDCD1/TOX upregulation but retaining GZMK expression. Negative TIM3/ENTPD1/CXCL13 distinguishes from terminal exhaustion. Intermediate between effector and exhausted."
```

#### Tier 1b: State (marker_registry)

**Exhaustion trajectory states**:

The article describes multiple exhaustion paths. These are **state transitions**, not stable identities:

| State | Markers | Description |
|-------|---------|-------------|
| Predysfunctional | GZMK, CXCR4, PDCD1 (low), TOX (low) | Early dysfunction |
| Transitional exhaustion | GZMK, CXCR4, PDCD1, TOX, LAG3 | Intermediate |
| Terminal exhaustion | HAVCR2, ENTPD1, CXCL13, TOX, PRDM1 | Late dysfunction |

Add to `state` section:

```toml
[[state.minor]]
name = "Predysfunctional T cell"
color = "#DAA520"
markers = ["GZMK", "CXCR4", "PDCD1", "TOX", "LAG3"]
negative_markers = ["HAVCR2", "ENTPD1", "CXCL13"]

[state.minor.metadata]
kind = "state"
category = "immune_dysfunction"
scope = "lineage_restricted"
applies_to = ["CD8+ T"]
alias_of = "Pre-exhaustion"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

```toml
[[state.minor]]
name = "Transitional exhaustion"
color = "#CD853F"
markers = ["GZMK", "CXCR4", "PDCD1", "TOX", "LAG3", "TIGIT"]
negative_markers = ["HAVCR2", "ENTPD1", "PRDM1"]

[state.minor.metadata]
kind = "state"
category = "immune_dysfunction"
scope = "lineage_restricted"
applies_to = ["CD8+ T"]
alias_of = "Tex-transitory"
granularity = "state"
species = "human"
use_for_global_annotation = false
use_for_state_annotation = true
```

#### Tier 2: Tumor Context (marker_tumor_human.toml)

**Tumor-type-specific T cell composition hints**:

The article found different exhaustion paths preferred in different cancer types:
- Liver/colon: higher exhausted T cell fractions
- Lung: two pre-exhaustion states
- Melanoma: linear predysfunctional → dysfunctional progression
- Multiple myeloma: no notable exhausted populations

These are **tumor-type hints**, not universal markers:

```toml
[[tumor_type_hint]]
name = "High exhaustion TME"
color = "#B22222"
markers = ["PDCD1", "HAVCR2", "LAG3", "TOX", "CXCL13", "ENTPD1"]

[tumor_type_hint.metadata]
kind = "tumor_evidence"
granularity = "tumor_type_hint"
cancer_type = ["Liver Cancer", "Colon Cancer"]
scope = "cancer_type_specific"
applies_to = ["Liver", "Colon"]
evidence_tier = "atlas_supported"
source_type = "pan_cancer_atlas"
review_status = "needs_review"
use_for_global_annotation = false
use_for_state_annotation = true
use_for_malignancy_interpretation = false
species = "human"
source = { title = "Pan-cancer single-cell landscape of tumor-infiltrating T cells", year = "2021", doi = "10.1126/science.abe6474" }
notes = "Liver and colon cancers show higher fractions of exhausted T cells compared to lung cancer. Predysfunctional → dysfunctional progression is more linear in melanoma."
```

#### Tier 3: Functional Programs (geneset)

**Exhaustion trajectory scoring** (Science 2021):

```json
"T_exhaustion_trajectory": {
  "genes": ["TCF7", "SELL", "IL7R", "GZMK", "CXCR4", "PDCD1", "TOX", "LAG3", "TIGIT", "HAVCR2", "ENTPD1", "CXCL13", "PRDM1"],
  "description": "Full CD8+ T cell exhaustion trajectory from naive through predysfunctional to terminal exhaustion (Science 2021 pan-cancer atlas)"
}
```

**Predysfunctional signature**:

```json
"T_predysfunctional": {
  "genes": ["GZMK", "CXCR4", "EOMES", "PDCD1", "TOX", "LAG3", "TNFRSF9"],
  "description": "Predysfunctional T cell signature (early exhaustion, GZMK-retaining)"
}
```

**T cell composition classifier**:

```json
"TME_T_cell_enriched": {
  "genes": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "TBX21", "EOMES", "PDCD1", "TOX"],
  "description": "T cell-inflamed tumor microenvironment signature"
}
```

```json
"TME_T_cell_desert": {
  "genes": ["FOXP3", "IL10", "TGFB1", "VEGFA", "ARG1", "IDO1"],
  "description": "T cell-excluded/desert tumor microenvironment signature"
}
```

---

## Comprehensive Registry Update Plan

### A. `marker_registry_human.toml` Updates

#### New Subtype Entries

1. **Tfr (Follicular Regulatory T)** — Article 3
2. **CD8+ predysfunctional** — Article 7

#### New State Entries

1. **T cell stress response (TSTR)** — Article 6
2. **Treg cytotoxic** — Article 3
3. **CD8+ Tex-progenitor** — Articles 3, 4
4. **CD8+ Tex-transitory** — Article 4
5. **Predysfunctional T cell** — Article 7
6. **Transitional exhaustion** — Article 7

#### Existing Entry Updates

| Existing Entry | Update | Source |
|---------------|--------|--------|
| CD4+ Treg | Add tumor-enrichment note | Article 3 |
| CD4+ Th17 | Add late-stage lung cancer note | Article 3 |
| MAIT | Add COVID-19 enrichment note | Article 3 |
| CD8+ Tex | Refine to 3 subtypes (progenitor/transitory/terminal) | Articles 3, 4, 7 |
| CD8+ Tn | Add healthy-abundance note | Article 3 |
| Stress-high (artifact) | Clarify distinction from TSTR state | Article 6 |

### B. `marker_tumor_human.toml` Updates

#### New Cancer State Entries

1. **TSTR in tumor** — Article 6
2. **High exhaustion TME** (liver/colon) — Article 7

### C. `genesets_cancer_signatures.json` Updates

#### New Categories and Programs

**New category: `T_Cell_States`**

| Program | Genes (approx.) | Source | Description |
|---------|----------------|--------|-------------|
| TCellSI_quiescence | ~10-15 | Article 5 | Resting T cell state |
| TCellSI_regulating | ~10-15 | Article 5 | Treg-mediated suppression |
| TCellSI_proliferation | ~10-15 | Article 5 | Cell cycling |
| TCellSI_helper | ~10-15 | Article 5 | Cytokine production |
| TCellSI_cytotoxicity | ~10-15 | Article 5 | Killing function |
| TCellSI_progenitor_exhaustion | ~10-15 | Article 5 | Stem-like exhausted |
| TCellSI_terminal_exhaustion | ~10-15 | Article 5 | Terminal dysfunctional |
| TCellSI_senescence | ~10-15 | Article 5 | Aging dysfunction |

**New category: `T_Cell_Function`**

| Program | Genes (approx.) | Source | Description |
|---------|----------------|--------|-------------|
| CD8_Tex_trajectory | ~15-20 | Articles 4, 7 | Full exhaustion trajectory |
| CD8_cytotoxic | ~8-10 | Articles 4, 7 | Cytotoxic effector |
| CD8_memory | ~8-10 | Article 4 | Memory formation |
| T_predysfunctional | ~8-10 | Article 7 | Early exhaustion |
| TSTR_heat_shock | ~12-15 | Article 6 | Stress response |
| T_cell_ICB_response | ~8-10 | Article 6 | ICB responsive |
| T_cell_ICB_resistance | ~10-12 | Article 6 | ICB resistant |
| TME_T_cell_enriched | ~10-12 | Article 7 | T cell-inflamed TME |
| TME_T_cell_desert | ~8-10 | Article 7 | T cell-excluded TME |

---

## Summary Statistics

| Category | Count | Target Resource |
|----------|-------|-----------------|
| New subtype entries | 2 | marker_registry |
| New state entries | 6 | marker_registry |
| Existing entry updates | 7 | marker_registry |
| New cancer state entries | 2 | marker_tumor |
| New geneset programs | 17 | genesets_cancer_signatures.json |
| Naming alignment notes | 10+ | marker_registry (naming) |
| Validation sources | 3 | marker_registry (evidence tier promotion) |

---

## Curation Notes and Conflicts

### 1. HSPA1A Ambiguity
**Issue**: HSPA1A appears in:
- `artifact.Stress-high` (QC/dissociation artifact)
- `state.TSTR` (biological stress response in tumors)
- CD16hi NK (from Batch 01)
- CD8 Trm (from Batch 01)

**Resolution**:
- Keep `Stress-high` as **artifact** (context: dissociation stress, broad across cell types)
- Add TSTR as **state** (context: tumor microenvironment, T cell-specific)
- Remove HSPA1A from NK/Trm **identity** markers (state marker, not identity)
- Document in notes: "HSPA1A is a stress response marker, not a stable identity marker"

### 2. Exhaustion Nomenclature
**Issue**: Multiple overlapping exhaustion states across papers:
- Article 3: Tex-progenitor, Tex-terminal
- Article 4: Tex-p, Tex-t, Tex-term
- Article 7: Predysfunctional, Transitional, Terminal
- Existing registry: Tex (PDCD1, TOX, LAG3, HAVCR2)

**Resolution**:
- Use **modular nomenclature** (Article 1): describe by marker expression rather than fixed labels
- Registry hierarchy:
  - `state.T cell exhaustion-like` (existing, broad)
  - `state.Tex-progenitor` (TCF7+ SELL+ PDCD1+)
  - `state.Tex-transitory` (GZMK+ PDCD1+ TOX+)
  - `state.Tex-terminal` (HAVCR2+ ENTPD1+ CXCL13+)
  - `state.Predysfunctional` (GZMK+ PDCD1(low) TOX(low))

### 3. TCellSI vs Registry
**Issue**: TCellSI (Article 5) provides 8 state scoring programs. These are **scoring signatures**, not identity markers.

**Resolution**:
- **Do NOT** add TCellSI gene sets to `marker_registry`
- **DO** add to `genesets_cancer_signatures.json` as functional scoring programs
- Cross-reference in notes: "For T cell state scoring, use TCellSI programs in geneset"

### 4. Population Specificity (Article 2)
**Issue**: Asian diversity atlas may reveal population-specific markers.

**Resolution**:
- If population-specific subtypes found: add to `subtype` with `scope = "population_specific"`
- If markers are stable: use as validation source for existing entries
- Abstract not available; recommend checking paper for specific findings

---

## Recommended Curation Workflow for Remaining Papers

For subsequent batches (B cells, Myeloid, etc.), apply the same tier classification:

1. **Read abstract** → identify paper type (atlas/review/tool)
2. **Classify to tier**:
   - Identity markers (3-12 genes) → `marker_registry`
   - State markers (3-15 genes) → `marker_registry` (state section)
   - Tumor-context markers → `marker_tumor`
   - Scoring programs (>5 genes) → `geneset`
3. **Extract markers** → verify official gene symbols
4. **Check conflicts** → compare with existing entries
5. **Assign metadata** → evidence_tier, source_type, review_status, routing flags
6. **Document notes** → biological caveats, source figures/tables

---

## Next Steps

1. **Manual review**: Verify each proposed entry against original paper figures/tables
2. **Gene symbol check**: Confirm all markers are official HGNC symbols
3. **Manager testing**: Test new entries with `get_marker_manager()` views
4. **TCellSI gene sets**: Extract exact gene lists from publication or R package
5. **Article 2 detail**: Check Asian diversity paper for population-specific findings
