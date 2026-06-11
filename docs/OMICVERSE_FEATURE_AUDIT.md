# OmicVerse Feature Audit for scLucid Integration

**Date**: 2026-06-11  
**Source repository**: https://github.com/omicverse/omicverse  
**Source license**: GNU General Public License v3.0 (GPL-3.0)  
**Target project**: scLucid (MIT licensed)  
**Audit purpose**: Identify OmicVerse bulk RNA-seq and spatial transcriptomics capabilities that can inspire scLucid functionality, while avoiding direct GPL source-code incorporation.

## License Constraint

OmicVerse is GPL-3.0 licensed. scLucid is MIT licensed. Directly copying OmicVerse source code, translated code, or derivative implementation structure into scLucid would likely require scLucid to adopt a GPL-compatible distribution model. Therefore, this audit classifies features by how they can be added to scLucid without violating that constraint:

- `clean_room_reimplement`: Statistical formulae or algorithms that are publicly documented and can be implemented independently from first principles.
- `optional_wrapper`: Mature third-party packages with compatible licenses that scLucid can call through optional dependencies.
- `defer`: Useful but complex features (often deep-learning based) that should wait until core APIs stabilize.
- `do_not_include`: Features whose value is low for scLucid's tumor-scRNA-seq focus, or whose only viable implementation path would require copying GPL code.

## Bulk RNA-seq Audit

| Feature | OmicVerse Location | Classification | Rationale |
|---------|-------------------|----------------|-----------|
| CPM/TPM/FPKM/RPKM normalization | `omicverse/bulk/_Deseq2.py` | `clean_room_reimplement` | Standard textbook formulae; no creative expression. |
| DESeq2-style size factor estimation | `omicverse/bulk/_Deseq2.py` | `clean_room_reimplement` | Published algorithm (Anders & Huber 2010). |
| Bulk t-test / Welch DE | `omicverse/bulk/_Deseq2.py` | `clean_room_reimplement` | Standard scipy/statsmodels tests. |
| Bulk DE via pydeseq2 | `omicverse/bulk/_Deseq2.py` | `optional_wrapper` | Wrap the MIT-compatible `pydeseq2` package directly. |
| Bulk DE via edgeR/limma | `omicverse/bulk/_Deseq2.py` | `optional_wrapper` / `defer` | OmicVerse vendors pure-Python ports; scLucid can wrap real edgeR/limma via rpy2 or defer. |
| Continuous/trait-associated DE | `omicverse/bulk/_Deseq2.py` | `clean_room_reimplement` | Pearson/correlation/regression on sample-level statistics. |
| Time-course DE | `omicverse/bulk/_Deseq2.py` | `clean_room_reimplement` | Spline + moderated F-test; publicly documented. |
| Mfuzz-like temporal clustering | `omicverse/bulk/_Deseq2.py` | `optional_wrapper` | Wrap `pymfuzz` or similar; fuzzy c-means is generic. |
| ORA enrichment | `omicverse/bulk/_Enrichment.py` | `optional_wrapper` | Wrap `gseapy` (MIT). |
| GSEA | `omicverse/bulk/_Enrichment.py` | `optional_wrapper` | Wrap `gseapy` (MIT). |
| Enrichment dot/clustermap plots | `omicverse/bulk/_Enrichment.py` | `clean_room_reimplement` | Visualization logic is generic; use scLucid plotting style. |
| WGCNA | `omicverse/bulk/_wgcna.py`, `_Gene_module.py` | `optional_wrapper` | Wrap `PyWGCNA` or similar; avoid vendoring OmicVerse's vendored copy. |
| STRING PPI query | `omicverse/bulk/_network.py` | `clean_room_reimplement` | Calls public STRING REST API; API client can be rewritten. |
| PPI network analysis | `omicverse/bulk/_network.py` | `clean_room_reimplement` | NetworkX-based analysis of STRING results. |
| TCGA data ingestion | `omicverse/bulk/_tcga.py` | `defer` | Useful for tumor validation, but tied to GDC file conventions and large data handling; design separately. |
| TCGA survival analysis | `omicverse/bulk/_tcga.py` | `defer` | Related to TCGA ingestion. |
| ComBat batch correction | `omicverse/bulk/_combat.py` | `optional_wrapper` | Wrap `combat` or `pymc`/`scanpy` combat; standard EB algorithm. |
| Bulk deconvolution (TAPE, Scaden, BayesPrism, OmicsTweezer) | `omicverse/bulk/_decov.py` | `optional_wrapper` / `defer` | scLucid already has BayesPrism/DWLS native code; other backends can be optional wrappers. |
| CHM13 reference utilities | `omicverse/bulk/_chm13.py` | `do_not_include` | Out of scope for tumor scRNA-seq focus. |
| Alignment utilities | `omicverse/bulk/_alignment/` | `do_not_include` | Out of scope. |
| Dynamic tree cut | `omicverse/bulk/_dynamicTree.py` | `optional_wrapper` | Available in `scikit-learn`/`scipy` dendrogram tooling. |
| Gene ID mapping | `omicverse/bulk/_Deseq2.py` | `clean_room_reimplement` | Use public gene annotation resources; implementation is data lookup. |

## Spatial Transcriptomics Audit

| Feature | OmicVerse Location | Classification | Rationale |
|---------|-------------------|----------------|-----------|
| Spatial neighbor graph (KNN / radius) | `omicverse/space/_svg.py` | `clean_room_reimplement` | sklearn.neighbors / scipy.spatial on `obsm['spatial']`. |
| Moran's I spatial autocorrelation | `omicverse/space/_svg.py` | `clean_room_reimplement` | Publicly defined statistic. |
| Geary's C spatial autocorrelation | `omicverse/space/_svg.py` | `clean_room_reimplement` | Publicly defined statistic. |
| Spatially variable genes (Moran/Geary) | `omicverse/space/_svg.py` | `clean_room_reimplement` | Publicly defined statistic + permutation test. |
| Spatially variable genes (PROST) | `omicverse/space/_svg.py` | `optional_wrapper` | Wrap the external PROST package. |
| Spatially variable genes (SpatialDE/SOMDE) | `omicverse/space/_svg.py` | `optional_wrapper` | Wrap external packages or vendored code under compatible license. |
| Visium IO | `omicverse/space/_tools.py` | `clean_room_reimplement` | Scanpy/squidpy already provide readers; scLucid can wrap with audit. |
| Visium crop / rotate / subset window | `omicverse/space/_tools.py` | `clean_room_reimplement` | Coordinate transformation and array slicing. |
| Tissue zones via NMF | `omicverse/space/_tissue_zones.py` | `clean_room_reimplement` | sklearn.decomposition.NMF on features. |
| Spatial clustering domain detection (Squidpy workflow) | `omicverse/space/_cluster.py` | `optional_wrapper` | Wrap `squidpy`/`scanpy`. |
| STAGATE | `omicverse/space/_cluster.py` | `optional_wrapper` / `defer` | External GNN package; optional `spatial-deep` extra. |
| GraphST / CAST / STAligner | `omicverse/space/` | `optional_wrapper` / `defer` | External GNN packages; optional `spatial-deep` extra. |
| Spatial deconvolution (Tangram) | `omicverse/space/_deconvolution.py` | `optional_wrapper` | Wrap `tangram-sc`. |
| Spatial deconvolution (cell2location) | `omicverse/space/_deconvolution.py` | `optional_wrapper` | Wrap `cell2location`. |
| Spatial deconvolution (RCTD) | `omicverse/space/_deconvolution.py` | `optional_wrapper` | Wrap `rctd-py`. |
| Spatial cell-cell communication (COMMOT) | `omicverse/space/_commot.py` | `optional_wrapper` | Wrap `commot`. |
| SpaceFlow embedding | `omicverse/space/_spaceflow.py` | `optional_wrapper` / `defer` | External deep-learning package. |
| GASTON isodepth | `omicverse/space/_gaston.py` | `optional_wrapper` / `defer` | External package. |
| Spatial trajectory tensor (STT) | `omicverse/space/_stt.py` | `optional_wrapper` / `defer` | External package. |
| SpatRio sc→spatial mapping | `omicverse/space/_spatrio.py` | `optional_wrapper` / `defer` | External package. |
| CellCharter | `omicverse/space/_cellcharter.py` | `optional_wrapper` / `defer` | External package. |
| Histology/pathology FM integration | `omicverse/space/histo/` | `defer` | Large scope; depends on pathology foundation models and `wsidata`/`lazyslide`. |

## Priority Recommendation

### Immediate (Phase 1)
Focus only on `clean_room_reimplement` items with high tumor/TME value:

1. Bulk normalization (CPM/TPM/FPKM/RPKM) and median-ratio size factors.
2. Bulk Welch/t-test DE with BH correction.
3. Continuous trait association for deconvolved proportions.
4. Spatial neighbor graph construction.
5. Moran's I / Geary's C spatial autocorrelation.
6. Spatially variable gene detection via Moran's I.
7. Spatial window subsetting and Visium crop/rotate helpers.

### Short-term (Phase 2-3)
1. Optional wrappers for `pydeseq2`, `gseapy`, `combat`.
2. Optional wrappers for `squidpy` spatial clustering and `tangram` deconvolution.
3. Bulk/pseudobulk concordance validator.

### Long-term / Defer
1. Deep-learning spatial models (STAGATE, GraphST, SpaceFlow, GASTON).
2. Histology/pathology foundation model integration.
3. Full TCGA ingestion pipeline.

## References

- OmicVerse repository: https://github.com/omicverse/omicverse
- OmicVerse license: GPL-3.0 (verified 2026-06-11)
- OmicVerse publication: *Nature Communications* 15, 6058 (2024). doi:10.1038/s41467-024-50194-3
- scLucid license: MIT
