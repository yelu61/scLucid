"""
Enhanced HVG selection and diagnostics for single-cell RNA-seq data.

Provides config-driven, reproducible, and fully traceable workflows
for highly variable gene (HVG) selection, with batch/sample awareness,
gene-type exclusion, automatic reporting, and large data support.
"""

import logging
from importlib.util import find_spec
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import scipy.sparse
from anndata import AnnData

HAS_VENN = find_spec("matplotlib_venn") is not None

from scLucid.utils.helpers import _show_or_close

from ...runtime import run_joblib_or_sequential
from ...utils import use_layer_as_X
from ..config import HVGConfig
from .plotting import plot_hvg_metrics

log = logging.getLogger(__name__)

__all__ = [
    "find_hvgs",
]


PROTECTED_GENE_PRESETS: Dict[str, List[str]] = {
    "immune_receptor": [
        "TRAC",
        "TRBC1",
        "TRBC2",
        "TRDC",
        "TRGC1",
        "TRGC2",
        "IGHM",
        "IGHD",
        "IGHG1",
        "IGHG2",
        "IGHG3",
        "IGHG4",
        "IGHA1",
        "IGHA2",
        "IGKC",
        "IGLC1",
        "IGLC2",
        "IGLC3",
    ],
    "cytokine": [
        "IL2",
        "IL4",
        "IL6",
        "IL7",
        "IL10",
        "IL12A",
        "IL12B",
        "IL15",
        "IL17A",
        "IL18",
        "IL21",
        "IFNG",
        "TNF",
        "TGFB1",
        "CCL2",
        "CCL3",
        "CCL4",
        "CCL5",
        "CXCL8",
        "CXCL9",
        "CXCL10",
        "CXCL13",
    ],
    "transcription_factor": [
        "TBX21",
        "GATA3",
        "RORC",
        "FOXP3",
        "TCF7",
        "LEF1",
        "PRDM1",
        "BCL6",
        "IRF4",
        "IRF8",
        "SPI1",
        "CEBPA",
        "CEBPB",
        "NFKB1",
        "STAT1",
        "STAT3",
        "STAT5A",
        "STAT5B",
        "MYC",
        "JUN",
        "FOS",
    ],
    "pathway": [
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "MX1",
        "OAS1",
        "OAS2",
        "CXCL10",
        "HLA-DRA",
        "HLA-DRB1",
        "CD274",
        "PDCD1",
        "CTLA4",
        "LAG3",
        "TIGIT",
        "HAVCR2",
    ],
    "tumor_heterogeneity": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "VIM",
        "FN1",
        "CDH1",
        "CDH2",
        "ITGA6",
        "CD44",
        "ALDH1A1",
        "PROM1",
        "MKI67",
        "TOP2A",
        "PCNA",
        "HIF1A",
        "VEGFA",
        "MYC",
        "SOX2",
        "NANOG",
        "POU5F1",
    ],
}


# --- Helper Functions ---#
def _get_hvg_input_matrix(adata: AnnData, input_layer: str):
    """Resolve HVG input matrix from the requested layer."""
    if input_layer == "X":
        return adata.X
    if input_layer in adata.layers:
        return adata.layers[input_layer]
    available = list(adata.layers.keys())
    raise KeyError(
        f"Layer '{input_layer}' not found in adata.layers. Available layers: {available or '[]'}"
    )


def _validate_hvg_input_matrix(X, input_layer: str, method: str) -> None:
    """Validate that HVG input is compatible with the chosen method."""
    if X.shape[0] == 0 or X.shape[1] == 0:
        raise ValueError(f"HVG input layer '{input_layer}' is empty with shape {X.shape}.")

    values = X.data if scipy.sparse.issparse(X) else np.asarray(X)
    if not np.all(np.isfinite(values)):
        raise ValueError(f"HVG input layer '{input_layer}' contains NaN or Inf values.")

    min_val = X.min() if scipy.sparse.issparse(X) else np.min(values)
    if min_val < 0:
        raise ValueError(
            f"HVG input layer '{input_layer}' contains negative values. "
            f"Method '{method}' expects non-negative normalized expression, not regressed/scaled residuals."
        )


def _diagnose_input_for_hvg(X, max_n=10000):
    """Compute and log basic statistics for the input matrix (dense or sparse)."""
    n_cells, n_genes = X.shape
    arr = X
    if hasattr(arr, "toarray"):
        arr = arr.toarray()
    if arr.shape[0] > max_n:
        arr = arr[:max_n]
    if arr.shape[1] > max_n:
        arr = arr[:, :max_n]
    mean = np.mean(arr)
    std = np.std(arr)
    min_val = np.min(arr)
    max_val = np.max(arr)
    zero_frac = np.mean(arr == 0)
    log.info(
        f"[HVG Input] shape={n_cells}x{n_genes}, mean={mean:.2f}, std={std:.2f}, min={min_val:.2f}, max={max_val:.2f}, zero_frac={zero_frac:.2%}"
    )
    return dict(
        mean=mean,
        std=std,
        min=min_val,
        max=max_val,
        zero_frac=zero_frac,
        n_cells=n_cells,
        n_genes=n_genes,
    )


def _matrix_looks_like_counts(X) -> bool:
    """Lightweight check for raw-count-like non-negative integer data."""
    values = X.data if scipy.sparse.issparse(X) else np.asarray(X)
    if values.size == 0:
        return True
    sample = values[: min(values.size, 10000)]
    return bool(np.all(sample >= 0) and np.allclose(sample, np.round(sample)))


def _resolve_hvg_flavor(flavor: str, input_layer: str, X) -> tuple[str, list[str]]:
    """Resolve HVG flavor while keeping the default path dependency-light."""
    notes: list[str] = []
    if flavor != "auto":
        return flavor, notes

    if input_layer in {"counts", "raw_counts"} or _matrix_looks_like_counts(X):
        if find_spec("skmisc") is not None:
            notes.append(
                "flavor='auto' resolved to 'seurat_v3' because the input looks like raw counts "
                "and scikit-misc is available."
            )
            return "seurat_v3", notes
        notes.append(
            "flavor='auto' resolved to 'seurat' because scikit-misc is unavailable; "
            "install scikit-misc to enable seurat_v3 HVG selection on raw counts."
        )
        return "seurat", notes

    notes.append(
        "flavor='auto' resolved to 'seurat' because the HVG input is log-normalized-like; "
        "seurat_v3 expects raw counts."
    )
    return "seurat", notes


def _exclude_genes(
    adata: AnnData,
    hvg_mask: np.ndarray,
    exclude_types: List[str] = ["mitochondrial", "ribosomal"],
    gene_types: Optional[Dict[str, np.ndarray]] = None,
    species: Optional[str] = None,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Exclude specific gene types from HVG selection.

    Enhanced version with:
    - Species-aware detection
    - Better pattern matching
    - Detailed logging

    Args:
        adata: AnnData object
        hvg_mask: Initial HVG mask (boolean array)
        exclude_types: List of gene types to exclude
        gene_types: Pre-computed gene type masks (optional)
        species: Species for detection (optional, auto-detected if None)

    Returns:
        Tuple of (updated HVG mask, exclusion counts dict)
    """
    # Use enhanced detection if gene_types not provided
    if gene_types is None:
        gene_types = _gene_type_detection(adata.var_names, species=species)

    excluded_counts = {}
    updated_mask = hvg_mask.copy()

    log.info(f"\nExcluding gene types: {', '.join(exclude_types)}")

    for gene_type in exclude_types:
        if gene_type in gene_types:
            type_mask = gene_types[gene_type]

            # Count how many HVGs are of this type
            n_hvg_of_type = (hvg_mask & type_mask).sum()
            excluded_counts[gene_type] = int(n_hvg_of_type)

            # Exclude from HVG set
            updated_mask = updated_mask & ~type_mask

            if n_hvg_of_type > 0:
                # Show examples
                excluded_genes = adata.var_names[hvg_mask & type_mask][:5]
                log.info(
                    f"  {gene_type}: excluded {n_hvg_of_type} genes "
                    f"(e.g., {', '.join(excluded_genes)})"
                )
            else:
                log.info(f"  {gene_type}: 0 genes to exclude")
        else:
            log.warning(
                f"  Gene type '{gene_type}' not recognized. "
                f"Available types: {list(gene_types.keys())}"
            )
            excluded_counts[gene_type] = 0

    total_excluded = hvg_mask.sum() - updated_mask.sum()
    log.info(f"\nTotal genes excluded: {total_excluded}")

    return updated_mask, excluded_counts


def _resolve_protected_gene_sets(
    adata: AnnData,
    *,
    presets: Optional[List[str]] = None,
    custom_sets: Optional[Dict[str, List[str]]] = None,
) -> tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """Resolve protected gene presets/custom sets against ``adata.var_names``."""
    requested: Dict[str, List[str]] = {}
    for preset in presets or []:
        if preset not in PROTECTED_GENE_PRESETS:
            log.warning("Unknown protected gene preset '%s'. Skipping.", preset)
            continue
        requested[preset] = list(PROTECTED_GENE_PRESETS[preset])
    for name, genes in (custom_sets or {}).items():
        requested[str(name)] = [str(g) for g in genes]

    upper_to_gene = {str(g).upper(): str(g) for g in adata.var_names}
    present: Dict[str, List[str]] = {}
    missing: Dict[str, List[str]] = {}
    for name, genes in requested.items():
        present_genes: List[str] = []
        missing_genes: List[str] = []
        seen = set()
        for gene in genes:
            key = str(gene).upper()
            if key in upper_to_gene:
                resolved = upper_to_gene[key]
                if resolved not in seen:
                    present_genes.append(resolved)
                    seen.add(resolved)
            else:
                missing_genes.append(str(gene))
        present[name] = present_genes
        missing[name] = missing_genes
    return present, missing


def _apply_hvg_biological_protection(
    adata: AnnData,
    hvg_mask: np.ndarray,
    *,
    presets: Optional[List[str]] = None,
    custom_sets: Optional[Dict[str, List[str]]] = None,
    max_extra_genes: Optional[int] = None,
) -> tuple[np.ndarray, Dict[str, object]]:
    """Include protected biology genes in the final HVG mask and build an audit report."""
    protected_sets, missing_sets = _resolve_protected_gene_sets(
        adata,
        presets=presets,
        custom_sets=custom_sets,
    )
    protected_genes = []
    seen = set()
    for genes in protected_sets.values():
        for gene in genes:
            if gene not in seen:
                protected_genes.append(gene)
                seen.add(gene)

    protected_mask = adata.var_names.isin(protected_genes)
    already_selected = protected_mask & hvg_mask
    selected_genes = set(adata.var_names[np.asarray(hvg_mask, dtype=bool)].astype(str))
    rescued_genes_all = [gene for gene in protected_genes if gene not in selected_genes]
    rescued_genes = list(rescued_genes_all)
    truncated = False

    if max_extra_genes is not None and len(rescued_genes) > int(max_extra_genes):
        rescued_genes = rescued_genes[: int(max_extra_genes)]
        truncated = True
    rescued_mask = adata.var_names.isin(rescued_genes)

    updated_mask = np.asarray(hvg_mask, dtype=bool).copy()
    updated_mask |= np.asarray(rescued_mask, dtype=bool)
    adata.var["hvg_biologically_protected"] = protected_mask
    adata.var["hvg_protection_rescued"] = rescued_mask

    preset_summary = {
        name: {
            "n_requested": len((PROTECTED_GENE_PRESETS.get(name) or custom_sets.get(name, [])) if custom_sets else PROTECTED_GENE_PRESETS.get(name, [])),
            "n_present": len(genes),
            "present_genes": genes,
            "missing_genes": missing_sets.get(name, []),
        }
        for name, genes in protected_sets.items()
    }
    report: Dict[str, object] = {
        "enabled": True,
        "presets": list(presets or []),
        "custom_set_names": list((custom_sets or {}).keys()),
        "n_protected_present": int(protected_mask.sum()),
        "n_already_selected": int(already_selected.sum()),
        "n_rescued": int(rescued_mask.sum()),
        "n_rescued_before_cap": len(rescued_genes_all),
        "rescued_genes": rescued_genes,
        "max_extra_genes": max_extra_genes,
        "truncated": truncated,
        "truncation_policy": (
            "Protected genes are ordered by protected preset order, custom set order, "
            "then gene order within each set; protection_max_extra_genes keeps the first genes in that order."
        ),
        "sets": preset_summary,
        "rationale": (
            "Protected genes are included to preserve interpretable immune, pathway, "
            "transcriptional, or tumor-heterogeneity signals that may not rank highly by variance alone."
        ),
    }
    return updated_mask, report


def _build_hvg_consensus_report(
    adata: AnnData,
    output_key: str,
    *,
    sample_key: Optional[str] = None,
) -> Dict[str, object]:
    """Summarize batch/sample-specific versus consensus HVG evidence when available."""
    report: Dict[str, object] = {
        "output_key": output_key,
        "n_final_hvgs": int(adata.var[output_key].sum()) if output_key in adata.var else 0,
        "sample_key": sample_key,
        "available": False,
    }
    sample_count_key = f"{output_key}_sample_count"
    nbatches_key = f"{output_key}_highly_variable_nbatches"
    if sample_count_key in adata.var:
        counts = pd.to_numeric(adata.var[sample_count_key], errors="coerce").fillna(0)
        report.update(
            {
                "available": True,
                "evidence_key": sample_count_key,
                "n_consensus_all_samples": int((counts == counts.max()).sum()) if counts.max() > 0 else 0,
                "sample_count_distribution": {
                    str(int(k)): int(v) for k, v in counts.value_counts().sort_index().items()
                },
                "median_sample_hvg_count": float(counts[counts > 0].median()) if (counts > 0).any() else 0.0,
            }
        )
    elif nbatches_key in adata.var:
        counts = pd.to_numeric(adata.var[nbatches_key], errors="coerce").fillna(0)
        report.update(
            {
                "available": True,
                "evidence_key": nbatches_key,
                "n_consensus_all_batches": int((counts == counts.max()).sum()) if counts.max() > 0 else 0,
                "batch_count_distribution": {
                    str(int(k)): int(v) for k, v in counts.value_counts().sort_index().items()
                },
            }
        )
    return report


def _compute_hvg_single_sample(
    sample_id: str,
    X_sample: Union[np.ndarray, scipy.sparse.spmatrix],
    var_names: pd.Index,
    n_top_genes: int,
    flavor: str,
    span: Optional[float] = None,
) -> pd.Series:
    """
    Compute HVGs for a single sample (for parallel execution).

    This function is designed to be called in parallel via joblib.
    It takes pre-extracted data to avoid pickling the entire AnnData object.

    Args:
        sample_id: Sample identifier
        X_sample: Expression matrix for this sample (cells × genes)
        var_names: Gene names
        n_top_genes: Number of top HVGs to select
        flavor: HVG selection flavor
        span: LOWESS span parameter (for seurat flavor)

    Returns:
        Boolean Series indicating HVG status for each gene
    """
    n_cells = X_sample.shape[0]
    n_genes = X_sample.shape[1]

    # Validate input
    if n_cells < 50:
        log.warning(
            f"Sample '{sample_id}' has only {n_cells} cells (< 50). "
            "HVG calculation may be unreliable. Returning all False."
        )
        return pd.Series([False] * n_genes, index=var_names, name=sample_id)

    try:
        # Create minimal AnnData for HVG calculation
        temp_adata = AnnData(X=X_sample)
        temp_adata.var_names = var_names

        # Compute HVGs using scanpy
        kwargs = {"n_top_genes": n_top_genes, "flavor": flavor, "inplace": True}
        if flavor == "seurat" and span is not None:
            kwargs["span"] = span

        sc.pp.highly_variable_genes(temp_adata, **kwargs)

        # Extract result
        hvg_mask = temp_adata.var["highly_variable"].values

        return pd.Series(hvg_mask, index=var_names, name=sample_id)

    except Exception as e:
        log.error(
            f"HVG calculation failed for sample '{sample_id}' "
            f"({n_cells} cells, {n_genes} genes): {e}"
        )
        # Return all False on failure
        return pd.Series([False] * n_genes, index=var_names, name=sample_id)


def _compute_hvg_per_sample_parallel(
    adata: AnnData,
    sample_key: str,
    n_top_genes: int,
    flavor: str = "seurat",
    span: Optional[float] = None,
    n_jobs: int = -1,
    min_cells_per_sample: int = 50,
) -> pd.DataFrame:
    """
    Compute HVGs for each sample in parallel.

    This is a memory-efficient implementation that:
    1. Extracts only necessary data for each sample
    2. Uses joblib for parallel processing
    3. Provides progress tracking
    4. Handles errors gracefully

    Args:
        adata: AnnData object
        sample_key: Column in adata.obs identifying samples
        n_top_genes: Number of top HVGs per sample
        flavor: HVG selection method
        span: LOWESS span (for seurat flavor)
        n_jobs: Number of parallel jobs (-1 = all cores)
        min_cells_per_sample: Minimum cells required for HVG calculation

    Returns:
        DataFrame with genes as rows and samples as columns.
        Values are boolean (True = gene is HVG in that sample).
    """
    if sample_key not in adata.obs.columns:
        raise KeyError(f"Sample key '{sample_key}' not found in adata.obs")

    # Get unique samples
    samples = adata.obs[sample_key].unique()
    n_samples = len(samples)

    log.info(f"Computing HVGs for {n_samples} samples in parallel (n_jobs={n_jobs})...")

    # Pre-filter samples by cell count
    valid_samples = []
    for sample in samples:
        n_cells = (adata.obs[sample_key] == sample).sum()
        if n_cells >= min_cells_per_sample:
            valid_samples.append(sample)
        else:
            log.warning(
                f"Skipping sample '{sample}': only {n_cells} cells "
                f"(< {min_cells_per_sample} minimum)"
            )

    if len(valid_samples) == 0:
        raise ValueError(f"No samples have >= {min_cells_per_sample} cells. Cannot compute HVGs.")

    log.info(f"Processing {len(valid_samples)}/{n_samples} samples with sufficient cells")

    # Prepare data for parallel processing
    # Extract X matrices for each sample BEFORE parallelization
    # to avoid repeated slicing overhead
    sample_data = []
    for sample in valid_samples:
        sample_mask = adata.obs[sample_key] == sample
        X_sample = adata[sample_mask, :].X

        # Convert to memory-efficient format
        if scipy.sparse.issparse(X_sample):
            # Ensure CSR format for efficient row slicing
            X_sample = X_sample.tocsr()
        else:
            # Already dense, keep as-is
            pass

        sample_data.append({"sample_id": sample, "X": X_sample, "n_cells": sample_mask.sum()})

    # Define wrapper for parallel execution
    def _process_sample(sample_dict):
        return _compute_hvg_single_sample(
            sample_id=sample_dict["sample_id"],
            X_sample=sample_dict["X"],
            var_names=adata.var_names,
            n_top_genes=n_top_genes,
            flavor=flavor,
            span=span,
        )

    # Execute in parallel with progress bar
    try:
        from tqdm.auto import tqdm

        # Use threading backend for I/O-bound tasks (better for sparse matrices)
        # Use loky backend for CPU-bound tasks
        backend = "loky" if scipy.sparse.issparse(adata.X) else "threading"

        results = run_joblib_or_sequential(
            _process_sample,
            tqdm(sample_data, desc="Computing HVGs per sample", unit="sample"),
            n_jobs=n_jobs,
            backend=backend,
            description="per-sample HVG selection",
        )

    except ImportError:
        # Fallback without tqdm
        log.info("tqdm not available. Running without progress bar...")
        backend = "loky" if scipy.sparse.issparse(adata.X) else "threading"

        results = run_joblib_or_sequential(
            _process_sample,
            sample_data,
            n_jobs=n_jobs,
            backend=backend,
            description="per-sample HVG selection",
        )

    # Combine results into DataFrame
    hvg_df = pd.DataFrame(results).T

    # Log statistics
    total_hvgs_per_sample = hvg_df.sum(axis=0)
    log.info("\nHVG counts per sample:")
    for sample, count in total_hvgs_per_sample.items():
        log.info(f"  {sample}: {count} HVGs")

    log.info(
        f"\nAverage HVGs per sample: {total_hvgs_per_sample.mean():.0f} "
        f"(±{total_hvgs_per_sample.std():.0f})"
    )

    return hvg_df


def _identify_highly_expressed_genes(
    adata: AnnData,
    n_genes: int = 50,
    layer: Optional[str] = None,
) -> np.ndarray:
    """
    Identify top highly expressed genes (to exclude from HVG selection).

    These genes are often housekeeping genes that dominate the dataset
    but don't contribute to biological variation.

    Args:
        adata: AnnData object
        n_genes: Number of top genes to identify
        layer: Layer to use (if None, use adata.X)

    Returns:
        Boolean array indicating highly expressed genes
    """
    import scipy.sparse

    # Get expression matrix
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers")
        X = adata.layers[layer]
    else:
        X = adata.X

    # Compute total expression per gene
    if scipy.sparse.issparse(X):
        gene_total_expr = np.array(X.sum(axis=0)).flatten()
    else:
        gene_total_expr = X.sum(axis=0)

    # Get top N genes
    top_indices = np.argsort(-gene_total_expr)[:n_genes]

    # Create boolean mask
    highly_expressed_mask = np.zeros(adata.n_vars, dtype=bool)
    highly_expressed_mask[top_indices] = True

    # Log top genes
    top_genes = adata.var_names[top_indices][:10]  # Show top 10
    log.info(f"Top highly expressed genes: {', '.join(top_genes)}")

    return highly_expressed_mask


def _identify_sample_specific_genes_parallel(
    adata: AnnData,
    sample_key: str,
    n_genes_per_group: int = 20,
    layer: Optional[str] = None,
    method: str = "t-test",
    n_jobs: int = 1,  # Default to 1 to avoid nested parallelization
) -> np.ndarray:
    """
    Identify sample-specific marker genes (to exclude from HVG selection).

    Improved version with:
    - Better error handling
    - Optional parallel computation
    - More robust DE method

    Args:
        adata: AnnData object
        sample_key: Column identifying samples
        n_genes_per_group: Top N marker genes per sample
        layer: Layer to use for DE analysis
        method: DE method ('t-test', 'wilcoxon', 'logreg')
        n_jobs: Number of parallel jobs for DE calculation

    Returns:
        Boolean array indicating sample-specific genes
    """
    if sample_key not in adata.obs.columns:
        raise KeyError(f"Sample key '{sample_key}' not found in adata.obs")

    n_samples = adata.obs[sample_key].nunique()

    if n_samples <= 1:
        log.info(
            f"Only {n_samples} sample(s) found. " "Skipping sample-specific gene identification."
        )
        return np.zeros(adata.n_vars, dtype=bool)

    log.info(
        f"Identifying sample-specific genes across {n_samples} groups " f"using {method} test..."
    )

    temp_adata = adata

    try:
        # Work on a copy to avoid modifying the original AnnData
        temp_adata = adata.copy()

        # Switch to specified layer if needed
        if layer is not None:
            if layer not in temp_adata.layers:
                raise KeyError(f"Layer '{layer}' not found")
            temp_adata.X = temp_adata.layers[layer]

        # Validate method
        valid_methods = ["t-test", "wilcoxon", "logreg"]
        if method not in valid_methods:
            log.warning(f"Method '{method}' not in {valid_methods}. Defaulting to 't-test'.")
            method = "t-test"

        # Perform differential expression
        sc.tl.rank_genes_groups(
            temp_adata,
            groupby=sample_key,
            method=method,
            n_genes=n_genes_per_group,
            pts=True,
        )

        # Extract marker genes
        marker_genes_df = pd.DataFrame(temp_adata.uns["rank_genes_groups"]["names"])

        # Get unique marker genes across all samples
        marker_genes = set(marker_genes_df.values.flatten())
        marker_genes.discard(None)
        marker_genes = list(marker_genes)

        # Create boolean mask
        sample_specific_mask = adata.var_names.isin(marker_genes)

        log.info(
            f"Identified {len(marker_genes)} unique sample-specific genes "
            f"({sample_specific_mask.sum()} found in current data)"
        )

        # Log top markers per sample
        log.info("Top sample-specific genes per sample:")
        for sample in temp_adata.obs[sample_key].unique()[:5]:
            top_genes = marker_genes_df[sample].head(3).tolist()
            log.info(f"  {sample}: {', '.join(top_genes)}")

        return sample_specific_mask

    except Exception as e:
        log.error(f"Sample-specific gene identification failed: {e}")
        log.exception("Detailed error:")
        return np.zeros(adata.n_vars, dtype=bool)


def _gene_type_detection(
    var_names: pd.Index,
    species: Optional[str] = None,
) -> Dict[str, np.ndarray]:
    """
    Enhanced gene type detection with species-specific patterns.

    Improvements over original _detect_gene_types:
    - Species-aware patterns
    - More comprehensive pattern matching
    - Better logging

    Args:
        var_names: Gene names
        species: Species ('human', 'mouse', 'rat', or None for auto-detect)

    Returns:
        Dict mapping gene type to boolean mask
    """
    # Auto-detect species if not provided
    if species is None:
        species = _infer_species_from_gene_names(var_names)
        log.info(f"Auto-detected species: {species}")

    gene_names = var_names.astype(str)
    n_genes = len(gene_names)

    # Species-specific patterns
    MT_PATTERNS = {
        "human": [
            r"^MT-",  # Standard mitochondrial
            r"^MTRNR",  # MT-rRNA
            r"^MT[ACNT][DOPT]",  # MTND, MTCO, MTATP, MTCYB
        ],
        "mouse": [
            r"^mt-",  # Mouse lowercase
            r"^Mt-",  # Some datasets use title case
        ],
        "rat": [
            r"^Mt-",
            r"^mt-",
        ],
    }

    RIBO_PATTERNS = {
        "human": [
            r"^RP[SL]\d+[A-Z]*$",  # RPS1, RPL1, RPS27A, etc.
            r"^MRPL\d+",  # Mitochondrial ribosomal (large)
            r"^MRPS\d+",  # Mitochondrial ribosomal (small)
        ],
        "mouse": [
            r"^Rp[sl]\d+[a-z]*$",
            r"^Mrpl\d+",
            r"^Mrps\d+",
        ],
        "rat": [
            r"^Rpl\d+",
            r"^Rps\d+",
            r"^Mrpl\d+",
            r"^Mrps\d+",
        ],
    }

    gene_types = {}

    # === 1. Mitochondrial genes ===
    mt_patterns = MT_PATTERNS.get(species, MT_PATTERNS["human"])
    mt_mask = pd.Series([False] * n_genes, index=var_names)

    for pattern in mt_patterns:
        matches = gene_names.str.match(pattern, case=True)
        mt_mask |= matches

    gene_types["mitochondrial"] = mt_mask.values

    # === 2. Ribosomal genes ===
    ribo_patterns = RIBO_PATTERNS.get(species, RIBO_PATTERNS["human"])
    ribo_mask = pd.Series([False] * n_genes, index=var_names)

    for pattern in ribo_patterns:
        matches = gene_names.str.match(pattern, case=True)
        ribo_mask |= matches

    gene_types["ribosomal"] = ribo_mask.values

    # === 3. Hemoglobin genes ===
    # Species-independent patterns
    hb_patterns = [
        r"^HB[ABGDEZQMTP]\d*$",  # Human: HBA1, HBB, HBG1, etc.
        r"^Hb[ab]-",  # Mouse: Hba-a1, Hbb-bs, etc.
    ]

    hb_mask = pd.Series([False] * n_genes, index=var_names)
    for pattern in hb_patterns:
        hb_mask |= gene_names.str.match(pattern, case=True)

    gene_types["hemoglobin"] = hb_mask.values

    # === 4. Heat shock proteins ===
    hsp_patterns = [
        r"^HSP[ABCDE]?\d+[A-Z]*\d*",  # HSPA1A, HSPB1, etc.
        r"^HSPH\d+",
        r"^DNAJ[ABC]\d+",  # HSP40 family
    ]

    hsp_mask = pd.Series([False] * n_genes, index=var_names)
    for pattern in hsp_patterns:
        hsp_mask |= gene_names.str.match(pattern, case=True)

    gene_types["heat_shock"] = hsp_mask.values

    # === 5. Immediate early genes ===
    # Conserved gene symbols across species
    ieg_symbols = [
        "FOS",
        "FOSB",
        "FOSL1",
        "FOSL2",
        "JUN",
        "JUNB",
        "JUND",
        "EGR1",
        "EGR2",
        "EGR3",
        "EGR4",
        "NR4A1",
        "NR4A2",
        "NR4A3",
        "IER2",
        "IER3",
        "IER5",
        "ZFP36",
        "ZFP36L1",
        "ZFP36L2",
        "ATF3",
        "DUSP1",
        "DUSP2",
    ]

    # Add mouse equivalents
    if species == "mouse":
        ieg_symbols.extend(
            [
                "Fos",
                "Fosb",
                "Fosl1",
                "Fosl2",
                "Jun",
                "Junb",
                "Jund",
                "Egr1",
                "Egr2",
                "Egr3",
                "Nr4a1",
                "Nr4a2",
                "Nr4a3",
                "Ier2",
                "Ier3",
                "Ier5",
            ]
        )

    gene_types["immediate_early"] = gene_names.isin(ieg_symbols)

    # === Log detection results ===
    log.info(f"Gene type detection results (species: {species}):")
    for gene_type, mask in gene_types.items():
        n_detected = mask.sum()
        pct = n_detected / n_genes * 100

        if n_detected > 0:
            log.info(f"  {gene_type}: {n_detected} genes ({pct:.2f}%)")

    return gene_types


def _infer_species_from_gene_names(var_names: pd.Index) -> str:
    """
    Infer species from gene naming conventions.

    Args:
        var_names: Gene names

    Returns:
        'human', 'mouse', or 'rat'
    """
    gene_names = var_names.astype(str)

    # Sample genes for efficiency
    sample_size = min(1000, len(gene_names))
    sample_genes = gene_names[:sample_size]

    # Count capitalization patterns
    all_caps = sum(1 for g in sample_genes if len(g) > 1 and g.isupper())
    title_case = sum(1 for g in sample_genes if len(g) > 1 and g[0].isupper() and not g.isupper())

    # Check for species-specific markers
    has_mt_upper = any(g.startswith("MT-") for g in sample_genes)  # Human
    has_mt_lower = any(
        g.startswith("mt-") or g.startswith("Mt-") for g in sample_genes
    )  # Mouse/Rat

    # Decision logic
    if all_caps > title_case * 1.5 or has_mt_upper:
        return "human"
    elif title_case > all_caps * 1.5 or has_mt_lower:
        # Distinguish mouse vs rat (both use title case)
        # Rats often have more Rpl/Rps vs Mrpl/Mrps
        has_rpl = any(g.startswith("Rpl") for g in sample_genes)
        if has_rpl:
            return "rat"
        return "mouse"
    else:
        log.warning(
            "Could not confidently infer species from gene names. " "Defaulting to 'human'."
        )
        return "human"


def adapt_n_top_genes(n_cells: int, method: str = "linear") -> int:
    """Adapt n_top_genes based on dataset size.

    For complex tissues like tumors with >20 cell types, the standard 2000
    HVGs may be insufficient to capture rare populations. This function
    scales HVG count with dataset size.

    Parameters
    ----------
    n_cells : int
        Number of cells in the dataset.
    method : str
        Scaling method: "linear" (default) or "log".

    Returns:
    -------
    int
        Recommended number of HVGs.
    """
    if method == "linear":
        if n_cells < 5000:
            return 2000
        elif n_cells < 20000:
            return 3000
        elif n_cells < 50000:
            return 4000
        else:
            return 5000
    elif method == "log":
        import math
        return min(8000, max(2000, int(1000 * math.log10(max(n_cells, 100)))))
    else:
        raise ValueError(f"Unknown method: {method}. Use 'linear' or 'log'.")


def _write_hvg_report(
    report_path: Path,
    stats: dict,
    n_hvg: int,
    config: HVGConfig,
    gene_type_counts: dict,
    protection_report: Optional[Dict[str, object]] = None,
    consensus_report: Optional[Dict[str, object]] = None,
) -> None:
    """Write a simple markdown HVG report."""
    with open(report_path, "w") as f:
        f.write("# HVG Selection Report\n\n")
        f.write(f"**Method:** {config.method}\n\n")
        f.write(f"**Input shape:** {stats.get('n_cells')} cells × {stats.get('n_genes')} genes\n\n")
        f.write("## Input Statistics\n")
        for k, v in stats.items():
            f.write(f"- {k}: {v:.3g}\n")
        f.write(f"\n## HVG count: {n_hvg}\n")
        f.write("\n## Parameters\n")
        for k, v in config.__dict__.items():
            f.write(f"- {k}: {v}\n")
        f.write("\n## Excluded gene types\n")
        for k, v in gene_type_counts.items():
            f.write(f"- {k}: {v}\n")
        if protection_report:
            f.write("\n## Biological protection\n")
            f.write(f"- enabled: {protection_report.get('enabled')}\n")
            f.write(f"- protected genes present: {protection_report.get('n_protected_present', 0)}\n")
            f.write(f"- rescued into final HVG: {protection_report.get('n_rescued', 0)}\n")
            rescued = protection_report.get("rescued_genes", [])
            if rescued:
                f.write(f"- rescued genes: {', '.join(map(str, rescued[:50]))}\n")
        if consensus_report:
            f.write("\n## Sample/batch consensus\n")
            f.write(f"- available: {consensus_report.get('available')}\n")
            f.write(f"- evidence key: {consensus_report.get('evidence_key', '')}\n")
            for key in ["sample_count_distribution", "batch_count_distribution"]:
                if consensus_report.get(key):
                    f.write(f"- {key}: {consensus_report.get(key)}\n")


def _to_savable_dict(d: dict) -> dict:
    """Recursively convert a dictionary to be h5ad-savable."""
    savable = {}
    for k, v in d.items():
        if isinstance(v, tuple):
            savable[k] = list(v)  # Convert tuples to lists
        elif isinstance(v, dict):
            savable[k] = _to_savable_dict(v)
        elif isinstance(v, (str, int, float, bool, list)) or v is None:
            savable[k] = v
        else:
            # For other complex objects, convert to string representation
            savable[k] = str(v)
    return savable


# --- Main Functions ---#
def find_hvgs(
    adata: AnnData,
    config: Optional[HVGConfig] = None,
    input_layer: str = "normalized",
    preserve_tumor_heterogeneity: bool = False,
    **kwargs,
) -> AnnData:
    """
    Config-driven, reproducible HVG selection with diagnostics and reporting.
    """
    # --- 1. Establish the final configuration ---
    if config is None:
        active_config = HVGConfig()
    else:
        # Create a copy of config and apply kwargs

        active_config = config.model_copy()

    # Apply overrides from kwargs
    ignored_override_keys = {"force", "n_jobs"}
    for key, value in kwargs.items():
        if hasattr(active_config, key):
            setattr(active_config, key, value)
        elif key not in ignored_override_keys:
            log.warning(f"Ignoring unknown parameter: '{key}'")

    # Extract parameters from the final config for use in the function
    force = kwargs.get("force", False)
    report = (
        active_config.report if hasattr(active_config, "report") else kwargs.get("report", False)
    )
    plot = kwargs.get("plot", active_config.plot)
    save_dir = active_config.save_dir
    n_top_genes = active_config.n_top_genes
    # Auto-adapt n_top_genes based on dataset size if enabled and user didn't
    # explicitly override n_top_genes via kwargs or config constructor.
    auto_n_top = getattr(active_config, "auto_n_top_genes", False)
    user_explicit_n_top = kwargs.get("n_top_genes") is not None
    if config is not None:
        user_explicit_n_top = user_explicit_n_top or (
            "n_top_genes" in getattr(config, "model_fields_set", set())
        )
    if auto_n_top and not user_explicit_n_top:
        adaptive_method = getattr(active_config, "auto_n_top_genes_method", "log")
        n_top_genes_adaptive = adapt_n_top_genes(adata.n_obs, method=adaptive_method)
        if n_top_genes_adaptive != n_top_genes:
            log.info(
                f"[HVG] Auto-adapting n_top_genes: {n_top_genes} -> {n_top_genes_adaptive} "
                f"(based on {adata.n_obs} cells, method='{adaptive_method}')"
            )
            n_top_genes = n_top_genes_adaptive
    method = active_config.method
    batch_key = active_config.batch_key
    requested_flavor = active_config.flavor
    exclude_gene_types = active_config.exclude_gene_types
    span = active_config.span

    log.info(f"[HVG] Diagnosing input data from layer '{input_layer}' ...")
    X = _get_hvg_input_matrix(adata, input_layer)
    _validate_hvg_input_matrix(X, input_layer, method)
    stats = _diagnose_input_for_hvg(X)
    flavor, flavor_notes = _resolve_hvg_flavor(requested_flavor, input_layer, X)
    for note in flavor_notes:
        log.info(f"[HVG] {note}")

    output_key = (
        f"highly_variable_{method}_{flavor}" if method == "scanpy" else f"highly_variable_{method}"
    )

    if output_key in adata.var and not force:
        n_existing = adata.var[output_key].sum()
        log.info(
            f"[HVG] Annotations found in '{output_key}' with {n_existing} genes. Use force=True to overwrite."
        )
        return adata

    with use_layer_as_X(adata, input_layer):
        if method == "scanpy":
            log.info(
                f"[HVG] Selecting with scanpy ({flavor}), n_top_genes={n_top_genes}, batch_key={batch_key}"
            )
            sc.pp.highly_variable_genes(
                adata,
                flavor=flavor,
                n_top_genes=n_top_genes,
                batch_key=batch_key,
                span=span,
                inplace=True,
            )
            adata.var[output_key] = adata.var["highly_variable"].copy()
            for metric in [
                "means",
                "dispersions",
                "dispersions_norm",
                "highly_variable_nbatches",
                "highly_variable_intersection",
            ]:
                if metric in adata.var:
                    adata.var[f"{output_key}_{metric}"] = adata.var[metric].copy()
            if output_key != "highly_variable":
                del adata.var["highly_variable"]

        elif method == "custom":
            log.info(f"[HVG] Custom HVG selection (n_top_genes={n_top_genes})")

            # Extract parameters
            sample_key = active_config.sample_key
            min_n_samples = active_config.min_n_samples
            n_highly_expressed_genes = active_config.n_highly_expressed_genes
            n_specific_genes = active_config.n_specific_genes

            # Validate sample_key
            if sample_key not in adata.obs.columns:
                raise KeyError(
                    f"Sample key '{sample_key}' not found in adata.obs. "
                    f"Available columns: {list(adata.obs.columns)}"
                )

            n_samples = adata.obs[sample_key].nunique()
            log.info(f"Found {n_samples} unique samples in column '{sample_key}'")

            if n_samples < min_n_samples:
                raise ValueError(
                    f"Only {n_samples} sample(s) found, but min_n_samples={min_n_samples}. "
                    "Custom HVG selection requires at least min_n_samples samples."
                )

            # === STEP 1: Compute HVGs per sample (PARALLELIZED) ===
            n_jobs = kwargs.get("n_jobs", -1)  # Get from kwargs or use all cores

            hvg_df = _compute_hvg_per_sample_parallel(
                adata,
                sample_key=sample_key,
                n_top_genes=n_top_genes,
                flavor=flavor,
                span=span,
                n_jobs=n_jobs,
                min_cells_per_sample=50,  # Minimum cells to calculate HVGs
            )

            # === STEP 2: Count sample occurrences ===
            # How many samples each gene appears as HVG in
            sample_counts = hvg_df.sum(axis=1)

            # Genes that are HVG in >= min_n_samples samples
            combined_hvgs = sample_counts >= min_n_samples

            # Store sample counts for diagnostics
            adata.var[f"{output_key}_sample_count"] = sample_counts

            log.info(
                f"\nGenes selected as HVG in >= {min_n_samples} samples: "
                f"{combined_hvgs.sum()} genes"
            )

            # Distribution of sample counts
            count_distribution = sample_counts.value_counts().sort_index()
            log.info("\nDistribution of HVG occurrences across samples:")
            for count, n_genes in count_distribution.items():
                if count > 0:
                    log.info(f"  {count} sample(s): {n_genes} genes")

            # === STEP 3: Identify genes to exclude ===
            genes_to_exclude = np.zeros(adata.n_vars, dtype=bool)

            # 3a. Highly expressed genes
            if n_highly_expressed_genes > 0:
                log.info(
                    f"\nIdentifying top {n_highly_expressed_genes} highly expressed genes "
                    "to exclude..."
                )

                highly_expressed_mask = _identify_highly_expressed_genes(
                    adata,
                    n_genes=n_highly_expressed_genes,
                    layer=input_layer,
                )

                genes_to_exclude |= highly_expressed_mask

                # Store for diagnostics
                adata.var[f"{output_key}_highly_expressed"] = highly_expressed_mask

                log.info(f"Marked {highly_expressed_mask.sum()} genes as highly expressed")

            # 3b. Sample-specific marker genes
            sample_specific_mask = np.zeros(adata.n_vars, dtype=bool)
            if n_specific_genes > 0:
                if preserve_tumor_heterogeneity:
                    log.info(
                        "preserve_tumor_heterogeneity=True: Skipping sample-specific gene exclusion to retain inter-tumor heterogeneity."
                    )
                else:
                    log.info(
                        "Identifying sample-specific genes to exclude (Batch Effect removal)..."
                    )

                    sample_specific_mask = _identify_sample_specific_genes_parallel(
                        adata,
                        sample_key=sample_key,
                        n_genes_per_group=n_specific_genes,
                        layer=input_layer,
                        method="t-test",
                        n_jobs=1,  # Avoid nested parallelization
                    )

                    genes_to_exclude |= sample_specific_mask

                # Store for diagnostics
                adata.var[f"{output_key}_sample_specific"] = sample_specific_mask

                log.info(f"Marked {sample_specific_mask.sum()} genes as sample-specific")

            # === STEP 4: Apply exclusions ===
            # Final HVG mask = (selected in >= min_n_samples) AND (not excluded)
            final_hvg_mask = combined_hvgs & ~genes_to_exclude

            # Store result
            adata.var[output_key] = final_hvg_mask

            # Log exclusion statistics
            log.info("\n" + "=" * 60)
            log.info("Custom HVG Selection Summary:")
            log.info("=" * 60)
            log.info(f"Initial HVG candidates: {combined_hvgs.sum()}")
            log.info(
                f"Excluded (highly expressed): {highly_expressed_mask.sum() if n_highly_expressed_genes > 0 else 0}"
            )
            log.info(
                f"Excluded (sample-specific): {sample_specific_mask.sum() if n_specific_genes > 0 else 0}"
            )
            log.info(
                f"Excluded (overlap): "
                f"{(highly_expressed_mask & sample_specific_mask).sum() if (n_highly_expressed_genes > 0 and n_specific_genes > 0) else 0}"
            )
            log.info(f"Final HVGs selected: {final_hvg_mask.sum()}")
            log.info("=" * 60 + "\n")

            # Optional: Show some example HVGs
            if final_hvg_mask.sum() > 0:
                example_hvgs = adata.var_names[final_hvg_mask][:20]
                log.info(f"Example HVGs: {', '.join(example_hvgs)}")

        elif method == "triku":
            try:
                import triku
            except ImportError:
                raise ImportError("Please install triku: pip install triku")
            log.info("[HVG] Running triku method")
            result = triku.tl.triku(adata, return_all=True)
            adata.var[output_key] = result["highly_variable"]
            adata.var[f"{output_key}_score"] = result["score"]
        else:
            raise ValueError(f"Unknown method '{method}'.")

    # Exclude gene types
    gene_type_counts = {}
    if exclude_gene_types:
        log.info(f"[HVG] Excluding gene types: {exclude_gene_types}")

        # Auto-detect species for enhanced pattern matching
        species = kwargs.get("species")

        current_mask = adata.var[output_key]

        # Use enhanced exclusion with species-aware detection
        updated_mask, excluded_counts = _exclude_genes(
            adata,
            current_mask,
            exclude_gene_types,
            gene_types=None,  # Will auto-detect
            species=species,
        )

        adata.var[output_key] = updated_mask
        gene_type_counts.update(excluded_counts)

    protection_report = {"enabled": False}
    if getattr(active_config, "protect_genes", True):
        protected_presets = list(getattr(active_config, "protected_gene_presets", []) or [])
        protected_sets = dict(getattr(active_config, "protected_gene_sets", {}) or {})
        if preserve_tumor_heterogeneity and "tumor_heterogeneity" not in protected_presets:
            protected_presets.append("tumor_heterogeneity")
        if protected_presets or protected_sets:
            protected_mask, protection_report = _apply_hvg_biological_protection(
                adata,
                np.asarray(adata.var[output_key], dtype=bool),
                presets=protected_presets,
                custom_sets=protected_sets,
                max_extra_genes=getattr(active_config, "protection_max_extra_genes", None),
            )
            adata.var[output_key] = protected_mask
            log.info(
                "[HVG] Biological protection added %s genes to '%s'.",
                protection_report.get("n_rescued", 0),
                output_key,
            )

    n_hvg = int(adata.var[output_key].sum())
    log.info(f"[HVG] Final number of highly variable genes: {n_hvg}")
    consensus_report = _build_hvg_consensus_report(
        adata,
        output_key,
        sample_key=getattr(active_config, "sample_key", None),
    )

    # --- Store metadata in .uns ---
    # Use a helper function to ensure the dictionary is savable
    savable_params = _to_savable_dict(active_config.to_dict())  # Pydantic's built-in serialization
    adata.uns.setdefault("sclucid", {}).setdefault("preprocess", {})["hvg"] = {
        "output_key": output_key,
        "method": method,
        "requested_flavor": requested_flavor,
        "resolved_flavor": flavor,
        "flavor_notes": flavor_notes,
        "input_layer": input_layer,
        "params": savable_params,
        "n_hvg": n_hvg,
        "input_stats": stats,
        "excluded_gene_types": gene_type_counts,
        "biological_protection": protection_report,
        "consensus_report": consensus_report,
    }

    if plot:
        plot_save_path = Path(save_dir) / "hvg_metrics.png" if save_dir else None
        if plot_save_path is not None:
            plot_save_path.parent.mkdir(parents=True, exist_ok=True)
        plot_hvg_metrics(
            adata,
            output_key,
            save_path=plot_save_path,
        )
        _show_or_close(plt.gcf())

    if report and save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)
        report_path = Path(save_dir) / "hvg_report.md"
        _write_hvg_report(
            report_path,
            stats,
            n_hvg,
            active_config,
            gene_type_counts,
            protection_report=protection_report,
            consensus_report=consensus_report,
        )
        log.info(f"[HVG] Report written to {report_path}")

    return adata
