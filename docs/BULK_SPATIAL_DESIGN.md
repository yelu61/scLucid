# Bulk and Spatial Design for scLucid

**Status**: Design document for Phase 0  
**Scope**: Define the API shape, storage contracts, diagnostics, and inference semantics for the new `scLucid.analysis.bulk` and `scLucid.tools.spatial` modules.

## Design Principles

All new bulk and spatial functionality must follow scLucid's existing philosophy:

1. **Diagnostic-first**: run `diagnose_*` before applying tests or corrections.
2. **Audit trail**: store parameters, warnings, assumptions, and evidence in `adata.uns["sclucid"]`.
3. **Explicit inference semantics**: every result table must declare its inference level and publication readiness.
4. **Pydantic configs**: public workflows accept config objects.
5. **Tumor-microenvironment focus**: bulk and spatial features should strengthen tumor, TME, malignancy, deconvolution, clinical association, and spatial niche interpretation.
6. **Optional dependency hygiene**: core imports must succeed without Squidpy, PyTorch, Tangram, PyDESeq2, or other heavy packages.

## Namespace and Storage Contract

### Bulk Analysis

Primary storage: `adata.uns["sclucid"]["analysis"]["bulk"]`

| Output | Storage Path | Type |
|--------|-------------|------|
| Diagnostics | `adata.uns["sclucid"]["analysis"]["bulk"]["diagnostics"]` | dict |
| Normalization | `adata.uns["sclucid"]["analysis"]["bulk"]["normalization"]` | dict |
| DE results | `adata.uns["sclucid"]["analysis"]["bulk"]["de"]` or returned DataFrame | DataFrame / dict |
| Deconvolution | `adata.uns["sclucid"]["analysis"]["bulk"]["deconvolution"]` | dict |
| Abundance tests | `adata.uns["sclucid"]["analysis"]["bulk"]["abundance"]` | dict |
| Clinical association | `adata.uns["sclucid"]["analysis"]["bulk"]["clinical"]` | dict |
| Trace/review | `adata.uns["sclucid"]["analysis"]["bulk"]["review_summary"]` | dict |

Backward compatibility: legacy calls through `scLucid.tools.bulk.deconvolve_bulk` continue to write to `adata.uns["sclucid"]["tools"][key_added]`, and additionally mirror a normalized trace into `adata.uns["sclucid"]["analysis"]["bulk"]["deconvolution"]`.

### Spatial Tools

Primary storage: `adata.uns["sclucid"]["tools"]["spatial"]` for utilities; `adata.uns["sclucid"]["analysis"]["spatial"]` for analysis-level statistics.

| Output | Storage Path | Type |
|--------|-------------|------|
| Diagnostics | `adata.uns["sclucid"]["tools"]["spatial"]["diagnostics"]` | dict |
| Neighbor graph | `adata.obsp["spatial_connectivities"]` | sparse matrix |
| Neighbor distances | `adata.obsp["spatial_distances"]` | sparse matrix |
| Moran's I per gene | `adata.uns["sclucid"]["tools"]["spatial"]["moran_i"]` | DataFrame |
| SVG results | `adata.var["spatially_variable"]` + `adata.uns["sclucid"]["tools"]["spatial"]["svg"]` | bool / DataFrame |
| Tissue zones | `adata.obsm["X_tissue_zones"]` + `adata.uns["sclucid"]["tools"]["spatial"]["tissue_zones"]` | matrix / dict |
| Workflow summary | `adata.uns["sclucid"]["tools"]["spatial"]["workflow"]` | dict |

## Inference-Level Tags

Every result DataFrame must contain these columns where applicable:

| Column | Meaning | Example Values |
|--------|---------|----------------|
| `inference_level` | Kind of inference performed | `"sample_level"`, `"exploratory_trait_association"`, `"exploratory_timecourse"`, `"descriptive_sample_level"`, `"exploratory_spatial"` |
| `valid_for_publication_inference` | Whether the comparison can be treated as formal biological inference | `True` / `False` |
| `replicate_requirement_met` | Whether biological replicate counts satisfy the method's assumptions | `True` / `False` |
| `diagnostic_status` | Outcome of the diagnostic gate | `"passed"`, `"warning"`, `"failed"` |
| `result_warning` | Human-readable caveat when inference is limited | str or `None` |

### Decision Rules for `valid_for_publication_inference`

- **Bulk DE**: `True` only when `n_samples_per_condition >= 2` for both groups and the method is sample-level (ttest/welch/pydeseq2/limma). `False` for descriptive single-sample mode.
- **Bulk trait association**: `False` by default (observational association, not causal).
- **Bulk time-course**: `False` by default unless longitudinal replicates are available and modeled.
- **Deconvolved abundance tests**: `True` when the underlying bulk samples have replicates; `False` otherwise.
- **Spatial statistics**: `False` by default (exploratory spatial patterns, not sample-level inference).

## Config Classes

### Bulk Configs

```python
class BulkDiagnosticsConfig(SclucidBaseConfig):
    min_samples_total: int = Field(default=2, ge=2)
    min_samples_per_condition: int = Field(default=2, ge=1)
    require_replicates: bool = Field(default=True)
    max_library_size_cv: Optional[float] = Field(default=None, ge=0)
    max_zero_gene_fraction: Optional[float] = Field(default=None, ge=0, le=1)

class BulkNormalizationConfig(SclucidBaseConfig):
    method: Literal["CPM", "TPM", "FPKM", "RPKM", "deseq2_size_factors"] = "CPM"
    gene_length_col: Optional[str] = None  # for TPM/FPKM/RPKM
    target_sum: float = Field(default=1e6, gt=0)
    pseudocount: float = Field(default=0.0, ge=0)

class BulkDEConfig(SclucidBaseConfig):
    method: Literal["ttest", "welch", "pydeseq2", "limma"] = "welch"
    condition_col: str
    condition1: str
    condition2: str
    covariates: List[str] = Field(default_factory=list)
    min_counts_per_gene: int = Field(default=10, ge=0)
    min_samples_expressing: int = Field(default=2, ge=1)
    p_adjust_method: Literal["fdr_bh", "bonferroni"] = "fdr_bh"
    alpha: float = Field(default=0.05, gt=0, lt=1)
    fallback_to_descriptive: bool = Field(default=False)

class BulkTraitAssociationConfig(SclucidBaseConfig):
    method: Literal["pearson", "spearman", "ols"] = "spearman"
    trait_col: str
    covariates: List[str] = Field(default_factory=list)
    min_samples: int = Field(default=10, ge=5)

class BulkDeconvolutionConfig(SclucidBaseConfig):
    method: Literal["BayesPrism", "DWLS"] = "BayesPrism"
    cell_type_key: str
    sample_key: str = "sampleID"
    min_common_genes: int = Field(default=100, ge=10)
    # method-specific kwargs passed through

class BulkAbundanceConfig(SclucidBaseConfig):
    method: Literal["ttest", "wilcoxon"] = "wilcoxon"
    group_col: str
    group1: str
    group2: str
    p_adjust_method: Literal["fdr_bh", "bonferroni"] = "fdr_bh"

class BulkClinicalAssociationConfig(SclucidBaseConfig):
    method: Literal["pearson", "spearman"] = "spearman"
    clinical_variable: str
    min_samples: int = Field(default=10, ge=5)
```

### Spatial Configs

```python
class SpatialDiagnosticsConfig(SclucidBaseConfig):
    spatial_key: str = "spatial"
    require_image: bool = False
    check_duplicate_coords: bool = True
    min_spots: int = Field(default=10, ge=2)

class SpatialNeighborsConfig(SclucidBaseConfig):
    spatial_key: str = "spatial"
    method: Literal["knn", "radius"] = "knn"
    n_neigh: int = Field(default=6, ge=1)
    radius: Optional[float] = None
    key_added: str = "spatial_neighbors"

class SpatialAutocorrConfig(SclucidBaseConfig):
    spatial_key: str = "spatial"
    mode: Literal["moran", "geary"] = "moran"
    n_permutations: int = Field(default=0, ge=0)
    key_added: str = "moran_i"

class SVGConfig(SclucidBaseConfig):
    spatial_key: str = "spatial"
    method: Literal["moran_i", "pearsonr", "prost"] = "moran_i"
    n_permutations: int = Field(default=100, ge=0)
    alpha: float = Field(default=0.05, gt=0, lt=1)
    layer: Optional[str] = None
    key_added: str = "spatially_variable"

class TissueZonesConfig(SclucidBaseConfig):
    n_components: int = Field(default=5, ge=2)
    method: Literal["nmf"] = "nmf"
    input: Literal["expression", "deconvolution"] = "deconvolution"
    key_added: str = "tissue_zones"

class VisiumIOConfig(SclucidBaseConfig):
    # Reader/wrapper configuration for 10x Visium data
    library_id: Optional[str] = None
    load_images: bool = True
```

## Public API Surface

### Bulk

```python
from scLucid.analysis.bulk import (
    diagnose_bulk_data_quality,
    normalize_bulk_counts,
    estimate_size_factors_median_ratio,
    run_bulk_de,
    run_bulk_trait_association,
    run_bulk_timecourse,
    run_bulk_ora,
    run_bulk_gsea,
    deconvolve_bulk,
    run_bulk_abundance_test,
    correlate_abundance_with_clinical,
)
```

### Spatial

```python
from scLucid.tools.spatial import (
    diagnose_spatial_data_quality,
    build_spatial_neighbors,
    compute_moran_i,
    compute_geary_c,
    compute_spatial_autocorr,
    find_spatially_variable_genes,
    find_tissue_zones,
    subset_spatial_window,
    crop_visium,
    rotate_visium,
    read_visium_10x,
    # legacy high-level workflow
    run_spatial_analysis,
    plot_spatial,
    run_spatial_batch,
    export_spatial_report,
)
```

## Diagnostic Gates

### `diagnose_bulk_data_quality(adata, config)`

Checks:
- `n_obs` >= `min_samples_total`
- Per-condition sample counts
- Zero-inflation fraction
- Library size distribution and CV
- Gene overlap if reference provided
- Condition label completeness
- Likely normalization state (warn if negative values or all-integer counts)

Returns:
```python
{
    "passed": bool,
    "warnings": List[str],
    "replicate_requirement_met": bool,
    "n_samples": int,
    "n_conditions": int,
    "min_replicates": int,
    "max_replicates": int,
    "library_size_cv": float,
    "zero_gene_fraction": float,
    "recommended_method": str,
}
```

### `diagnose_spatial_data_quality(adata, config)`

Checks:
- `obsm[spatial_key]` exists and has shape `(n_obs, 2)`
- No NaN coordinates
- Duplicate coordinate count
- Library ID presence if multi-library
- Image key presence if `require_image=True`
- Spot/cell count
- Platform hint if available

Returns:
```python
{
    "passed": bool,
    "warnings": List[str],
    "n_spots": int,
    "n_duplicate_coords": int,
    "spatial_extent": Tuple[float, float, float, float],
    "platform_hint": Optional[str],
    "image_key_present": bool,
}
```

## Optional Dependency Wrappers

All optional backends follow this pattern:

```python
def _require_optional(package: str, extra: str) -> Any:
    import importlib.util
    spec = importlib.util.find_spec(package)
    if spec is None:
        raise ImportError(
            f"'{package}' is required. Install with: pip install scLucid[{extra}]"
        )
    return importlib.import_module(package)
```

No optional package is imported at module load time.

## Backward Compatibility

- `src/scLucid/tools/bulk.py` remains importable and keeps all public names as thin re-exports or wrappers around `scLucid.analysis.bulk`.
- `src/scLucid/tools/spatial.py` becomes a package; `scLucid.tools.spatial.run_spatial_analysis` and other public names continue to work.
- Legacy storage paths continue to be populated for old API calls, with normalized traces mirrored into the new namespaces.

## Verification

After implementation, the following must succeed:

```bash
# Core import without optional extras
python -c "import scLucid; import scLucid.analysis.bulk; import scLucid.tools.spatial"

# Unit tests
pytest tests/analysis/bulk -v --tb=short
pytest tests/tools/spatial -v --tb=short
pytest tests/analysis/differential_expression/test_de_validation.py -v --tb=short

# Full suite
pytest tests/ -x --tb=short
```
