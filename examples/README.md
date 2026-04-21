# scLucid Examples

This directory contains production-ready Python scripts demonstrating common scLucid workflows. These examples use **Pydantic configuration** for automatic validation and reproducibility.

## Available Examples

### 1. quickstart.py
**Complete pipeline in 5 minutes**

End-to-end workflow demonstrating:
- Quality control with adaptive thresholds
- Preprocessing with HVG selection and batch correction
- Clustering and cell type annotation
- Results visualization and saving

**Use case**: New users wanting to get started quickly
**Data requirements**: PBMC-like dataset (e.g., `data/pbmc_raw.h5ad`)
**Runtime**: ~5 minutes for 10K cells

**Key features**:
- Uses `QCWorkflowConfig` and `WorkflowConfig` for reproducibility
- Demonstrates marker-based annotation with `MarkerManager`
- Saves results at each step for inspection

### 2. qc_pipeline.py
**Quality control workflow**

Comprehensive QC pipeline demonstrating:
- QC metrics calculation (genes, counts, mitochondrial)
- Doublet detection with Scrublet
- Cell filtering with configurable thresholds
- Before/after visualization

**Use case**: Quality control on new datasets
**Data requirements**: Raw count matrix (e.g., `data/raw_counts.h5ad`)
**Runtime**: ~3 minutes for 10K cells

**Key features**:
- Uses `QCThresholds` and `DoubletConfig` for precise control
- Generates diagnostic plots (violin, scatter)
- Reports filtering statistics

### 3. preprocessing.py
**Preprocessing pipeline**

Complete preprocessing workflow demonstrating:
- Normalization with target sum
- Highly variable gene (HVG) selection
- Data scaling and PCA
- Batch correction with Harmony

**Use case**: Preprocessing after QC
**Data requirements**: QC-filtered data (e.g., `results/qc_filtered.h5ad`)
**Runtime**: ~5 minutes for 10K cells

**Key features**:
- Modular approach with individual functions
- Configurable HVG parameters (method, n_top_genes, min/max mean)
- Optional batch correction with multiple methods (Harmony, Scanorama, BBKNN)

### 4. annotation_report.py
**Annotation review report export**

End-to-end example focused on annotation audit and publication-grade review output:
- Hybrid annotation with marker evidence and CellTypist
- Automatic review report export in both PNG and PDF
- Sidecar JSON and Markdown summaries for manual review and pipeline logging

**Use case**: Final annotation review before downstream tumor interpretation
**Data requirements**: Clustered dataset with UMAP-ready embedding or enough preprocessing to compute one
**Runtime**: ~5-10 minutes depending on annotation backend

**Key features**:
- Uses `AnnotationConfig(report=True)` for automatic report export
- Demonstrates `export_annotation_report()` with multi-format output
- Writes reviewer-facing `.md` and machine-readable `.json` sidecars

## Running the Examples

### Prerequisites

```bash
# Install scLucid with all dependencies
pip install sclucid[all]

# Or install from source
cd /path/to/scLucid
pip install -e ".[all]"
```

### Basic Usage

```bash
# Navigate to examples directory
cd examples/

# Run the quickstart example
python quickstart.py

# Run QC pipeline only
python qc_pipeline.py

# Run preprocessing only
python preprocessing.py

# Run annotation report export
python annotation_report.py
```

### Using Your Own Data

Modify the data path in each script:

```python
# In quickstart.py, line 17
adata = sc.read_h5ad("/path/to/your/data.h5ad")
```

## Data Requirements

### Input Format

All examples expect **AnnData objects** (`.h5ad` files) with:
- **Raw UMI counts** in `adata.X`
- **Basic metadata** in `adata.obs` (optional batch, sample columns)

### Recommended Data Structure

```
your_data/
├── raw_counts.h5ad       # Raw UMI matrix
├── metadata.csv          # Optional sample information
└── results/              # Output directory
    ├── qc/
    ├── preprocess/
    └── analysis/
```

## Configuration System

All examples use **Pydantic configuration** for:

✅ **Automatic validation** - Invalid parameters raise clear errors
✅ **Type safety** - Automatic type conversion
✅ **JSON serialization** - Save/load configs from files
✅ **Self-documenting** - All fields have descriptions

### Example: Saving and Loading Configs

```python
from scLucid.qc import QCWorkflowConfig

# Create config
config = QCWorkflowConfig(
    thresholds={"min_genes": 200, "pc_mt": 20.0},
    doublet={"method": "scrublet", "expected_doublet_rate": 0.06}
)

# Save to file
config.to_json_file("my_qc_config.json")

# Load from file
loaded_config = QCWorkflowConfig.from_json_file("my_qc_config.json")

# Use loaded config
adata = run_standard_qc(adata, config=loaded_config)
```

## Customization Tips

### Adjusting QC Thresholds

```python
# In qc_pipeline.py, modify line 43-48
thresholds = QCThresholds(
    min_genes=500,      # Increase for stricter filtering
    max_genes=7000,     # Adjust based on your data
    pc_mt=15.0,         # Lower for stricter MT filtering
    nmads=5.0           # Median absolute deviations
)
```

### Changing HVG Parameters

```python
# In preprocessing.py, modify line 38-43
hvg_config = HVGConfig(
    method="seurat",         # or "cell_ranger"
    n_top_genes=3000,        # Increase for more genes
    min_mean=0.0125,
    max_mean=3,
    span=0.3
)
```

### Batch Correction Methods

```python
# In preprocessing.py, modify line 79-82
integration_config = IntegrationConfig(
    method="harmony",        # Options: "harmony", "scanorama", "bbknn"
    batch_key="batch"        # Must match column in adata.obs
)
```

## Troubleshooting

### Common Issues

**Issue**: `ModuleNotFoundError: No module named 'pydantic'`
```bash
pip install pydantic>=2.0
```

**Issue**: `KeyError: 'batch'` during batch correction
- Ensure `adata.obs["batch"]` column exists
- Or skip batch correction for single-sample data

**Issue**: Memory errors with large datasets
```python
# In config, enable low-memory mode
from scLucid.config import set_config
set_config(low_memory_mode=True, chunk_size=1000)
```

**Issue**: Doublet detection fails
- Ensure scrublet is installed: `pip install scrublet`
- Or use alternative method: `method="doubletdetection"`

## Next Steps

After running these examples:

1. **Explore interactive tutorials**: See `../docs/notebooks/` for Jupyter notebooks
2. **Read API documentation**: See `../docs/source/api/` for detailed function references
3. **Check best practices**: See `../docs/source/best_practices.rst` for recommendations
4. **Review OpenSpec specs**: See `../openspec/specs/` for module requirements

## Additional Resources

- **Documentation**: https://sclucid.readthedocs.io/
- **GitHub Issues**: https://github.com/yelu61/scLucid/issues
- **Examples Data**: Download PBMC datasets from 10x Genomics website

## Contributing Examples

To contribute a new example:

1. Use Pydantic configs for all parameters
2. Include clear comments explaining each step
3. Add data requirements and runtime estimates
4. Test with both synthetic and real data
5. Update this README with your example

---

**Last Updated**: 2026-02-08
**scLucid Version**: 0.1.0
