"""
Utility functions for single-cell RNA-seq data analysis.

This module provides common helper functions that can be used
across different parts of the analysis pipeline.
"""

import gc
import logging
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union

import anndata
import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import io

from .contracts import LayerKeys, SCLUCID_ROOT, UnsKeys, ensure_sclucid_namespace

log = logging.getLogger(__name__)


def _import_scanpy():
    """Import scanpy lazily so lightweight scLucid imports stay robust."""
    try:
        import scanpy as sc
    except Exception as exc:
        raise ImportError(
            "scanpy is required for this helper. Install/repair scanpy or use "
            "a function that does not depend on 10x loading."
        ) from exc
    return sc


__all__ = [
    "read_10x",
    "load_10x_data",
    "use_layer_as_X",
    "sanitize_for_hdf5",
    "subset_adata",
    "subset_from_annotations",
    "merge_obs_metadata",
]


def read_10x(
    path: Optional[Union[str, Path]] = None,
    *,
    samples: Optional[List[str]] = None,
    base_dir: Optional[Union[str, Path]] = None,
    path_dict: Optional[Dict[str, str]] = None,
    metadata_dicts: Optional[Dict[str, Dict[str, Any]]] = None,
    possible_subpaths: Optional[List[str]] = None,
    var_names: Literal["gene_symbols", "gene_ids"] = "gene_symbols",
    make_unique: bool = True,
    cache: bool = True,
    sample_id: Optional[str] = None,
    species: Optional[str] = None,
    tissue: Optional[str] = None,
    tissue_type: Optional[str] = None,
    cancer_type: Optional[str] = None,
    output_file: Optional[Union[str, Path]] = None,
    compression: Optional[str] = "gzip",
    backup_existing: bool = True,
) -> AnnData:
    """
    Unified 10x Genomics loader for single-sample and multi-sample workflows.

    The function operates in one of two modes, chosen by which arguments are
    provided:

    **Single-sample mode** — pass ``path`` and (optionally) the
    sample-level metadata arguments. Auto-detects whether ``path`` is a
    Cell Ranger ``filtered_feature_bc_matrix`` directory or a ``.h5`` file.

    **Multi-sample mode** — pass ``samples`` together with either
    ``base_dir`` (auto-search the standard Cell Ranger subpaths under
    ``base_dir/<sample>/...``) or ``path_dict`` (explicit
    ``{sample_id: path}`` mapping). Per-sample metadata can be supplied
    via ``metadata_dicts`` (``{column: {sample_id: value}}``). The function
    loads each sample, attaches ``sampleID`` plus any metadata columns,
    falls back to a hand-written mtx parser when scanpy's reader fails, and
    concatenates everything into one AnnData.

    Both modes guarantee that the returned AnnData has ``layers["counts"]``
    populated (when ``X`` looks like raw integer counts) and that user
    metadata is also lifted into ``adata.uns["sclucid"]["analysis_context"]``
    so downstream stages can pick it up without extra arguments.

    Parameters
    ----------
    path : str or Path, optional
        Single-sample input. Either a Cell Ranger output directory or a
        ``.h5`` file produced by Cell Ranger.
    samples : list of str, optional
        Multi-sample sample IDs.
    base_dir : str or Path, optional
        Root directory used to locate ``<base_dir>/<sample>/<subpath>`` for
        each entry in ``samples``.
    path_dict : dict, optional
        Explicit ``{sample_id: path}`` mapping (multi-sample mode).
    metadata_dicts : dict of dicts, optional
        ``{column_name: {sample_id: value}}``. For each sample, the matching
        value is broadcast onto ``adata.obs[column_name]`` before concat.
    possible_subpaths : list of str, optional
        Candidate Cell Ranger subdirectory layouts. Defaults to common
        ``outs/filtered_feature_bc_matrix`` / ``filtered_feature_bc_matrix``
        layouts.
    var_names : {"gene_symbols", "gene_ids"}, default="gene_symbols"
        Whether to use gene symbols or Ensembl IDs as ``var_names``.
    make_unique : bool, default=True
        Make ``var_names`` unique by appending suffixes when symbols collide.
    cache : bool, default=True
        Cache the parsed matrix on disk (single-sample directory mode only).
    sample_id, species, tissue, tissue_type, cancer_type : str, optional
        Single-sample metadata stamped onto every cell of the result and
        lifted into ``analysis_context``.
    output_file : str or Path, optional
        Multi-sample only: if provided, write the concatenated AnnData here.
    compression : str, default="gzip"
        Compression for ``output_file``.
    backup_existing : bool, default=True
        Multi-sample only: rename any existing ``output_file`` to a
        timestamped backup before writing.

    Returns:
    -------
    AnnData
        Loaded data, ready for ``scl.run_pipeline()``.

    Raises:
    ------
    ValueError
        If neither single-sample nor multi-sample arguments are provided, or
        both are provided.
    FileNotFoundError
        If a required path does not exist.

    Examples:
    --------
    Single-sample (wet-lab one-liner):

    >>> import scLucid as scl
    >>> adata = scl.read_10x(
    ...     "data/pbmc3k/filtered_feature_bc_matrix/",
    ...     species="human",
    ...     tissue="PBMC",
    ... )

    Multi-sample (project-style):

    >>> adata = scl.read_10x(
    ...     samples=["P1_tumor", "P1_normal", "P2_tumor"],
    ...     base_dir="data/cellranger",
    ...     metadata_dicts={
    ...         "patient": {"P1_tumor": "P1", "P1_normal": "P1", "P2_tumor": "P2"},
    ...         "condition": {"P1_tumor": "tumor", "P1_normal": "normal", "P2_tumor": "tumor"},
    ...     },
    ...     output_file="results/merged.h5ad",
    ... )
    """
    multi_mode = samples is not None or path_dict is not None
    single_mode = path is not None

    if multi_mode and single_mode:
        raise ValueError(
            "Provide either `path` (single-sample) OR "
            "`samples`/`path_dict` (multi-sample), not both."
        )
    if not multi_mode and not single_mode:
        raise ValueError(
            "read_10x requires either `path` (single-sample mode) "
            "or `samples` with `base_dir`/`path_dict` (multi-sample mode)."
        )

    if single_mode:
        adata = _read_10x_single(
            path=path,
            var_names=var_names,
            make_unique=make_unique,
            cache=cache,
        )
        _attach_counts_layer(adata)
        _attach_sample_metadata(
            adata,
            sample_id=sample_id,
            species=species,
            tissue=tissue,
            tissue_type=tissue_type,
            cancer_type=cancer_type,
        )
        log.info(
            "Loaded 10x AnnData: %d cells x %d genes (counts layer: %s)",
            adata.n_obs,
            adata.n_vars,
            LayerKeys.COUNTS in adata.layers,
        )
        return adata

    return _read_10x_multi(
        samples=samples or list((path_dict or {}).keys()),
        base_dir=base_dir,
        path_dict=path_dict,
        metadata_dicts=metadata_dicts,
        possible_subpaths=possible_subpaths,
        var_names=var_names,
        make_unique=make_unique,
        cache=cache,
        output_file=output_file,
        compression=compression,
        backup_existing=backup_existing,
        species=species,
        tissue=tissue,
        tissue_type=tissue_type,
        cancer_type=cancer_type,
    )


def load_10x_data(
    samples: List[str],
    base_dir: Optional[str] = None,
    path_dict: Optional[Dict[str, str]] = None,
    metadata_dicts: Optional[Dict[str, Dict[str, Any]]] = None,
    possible_subpaths: Optional[List[str]] = None,
    output_file: Optional[str] = None,
    compression: Optional[str] = "gzip",
    backup_existing: bool = True,
    chunk_size: Optional[int] = None,  # noqa: ARG001 — accepted for back-compat
) -> AnnData:
    """
    Backward-compatible multi-sample 10x loader.

    Thin alias for :func:`read_10x` in multi-sample mode. New code should
    call :func:`read_10x` directly; this function is preserved so existing
    scripts and notebooks keep working unchanged.

    All arguments map 1:1 to :func:`read_10x`. ``chunk_size`` is accepted
    but currently unused; multi-sample loads are not chunked.
    """
    return read_10x(
        samples=samples,
        base_dir=base_dir,
        path_dict=path_dict,
        metadata_dicts=metadata_dicts,
        possible_subpaths=possible_subpaths,
        output_file=output_file,
        compression=compression,
        backup_existing=backup_existing,
    )


def _read_10x_single(
    *,
    path: Union[str, Path],
    var_names: str,
    make_unique: bool,
    cache: bool,
) -> AnnData:
    """Read a single Cell Ranger directory or .h5 file."""
    path_obj = Path(path).expanduser().resolve()
    if not path_obj.exists():
        raise FileNotFoundError(f"10x data path does not exist: {path_obj}")

    sc = _import_scanpy()

    if path_obj.is_dir():
        log.info("Reading 10x Cell Ranger directory: %s", path_obj)
        try:
            return sc.read_10x_mtx(
                path_obj,
                var_names=var_names,
                make_unique=make_unique,
                cache=cache,
            )
        except Exception as exc:
            log.warning(
                "scanpy.read_10x_mtx failed for %s: %s. Falling back to "
                "manual mtx parser.",
                path_obj,
                exc,
            )
            return _read_10x_manually(str(path_obj))
    if path_obj.suffix.lower() in {".h5", ".hdf5"}:
        log.info("Reading 10x HDF5 file: %s", path_obj)
        adata = sc.read_10x_h5(path_obj)
        if make_unique:
            adata.var_names_make_unique()
        return adata

    raise ValueError(
        f"Cannot recognise '{path_obj}' as a Cell Ranger output. Expected a "
        "directory containing matrix.mtx[.gz] or a .h5 file."
    )


def _read_10x_multi(
    *,
    samples: List[str],
    base_dir: Optional[Union[str, Path]],
    path_dict: Optional[Dict[str, str]],
    metadata_dicts: Optional[Dict[str, Dict[str, Any]]],
    possible_subpaths: Optional[List[str]],
    var_names: str,
    make_unique: bool,
    cache: bool,
    output_file: Optional[Union[str, Path]],
    compression: Optional[str],
    backup_existing: bool,
    species: Optional[str],
    tissue: Optional[str],
    tissue_type: Optional[str],
    cancer_type: Optional[str],
) -> AnnData:
    """Load multiple 10x samples, attach metadata, concat, optionally write."""
    if path_dict is None:
        if base_dir is None:
            raise ValueError("Either base_dir or path_dict must be provided")
        log.info("Searching for sample paths in %s", base_dir)
        path_dict = _find_sample_paths(str(base_dir), samples, possible_subpaths)

    sample_metadata: Dict[str, Dict[str, Any]] = {}
    if metadata_dicts:
        for sample in samples:
            sample_metadata[sample] = {}
            for metadata_name, metadata_dict in metadata_dicts.items():
                if sample in metadata_dict:
                    sample_metadata[sample][metadata_name] = metadata_dict[sample]

    valid_samples = [s for s in samples if s in path_dict]
    if len(valid_samples) < len(samples):
        log.warning(
            "Found valid paths for %d/%d samples", len(valid_samples), len(samples)
        )

    adata_list: List[AnnData] = []
    sc = _import_scanpy()
    for sample in valid_samples:
        sample_path = path_dict[sample]
        adata: Optional[AnnData] = None

        try:
            log.info("Loading %s with standard method from %s", sample, sample_path)
            adata = sc.read_10x_mtx(
                sample_path,
                var_names=var_names,
                cache=cache,
                make_unique=make_unique,
            )
        except Exception as e:
            log.warning("Standard method failed for %s: %s", sample, e)
            log.info("Attempting robust fallback method for %s...", sample)
            try:
                adata = _read_10x_manually(sample_path)
            except Exception as e2:
                log.error("Robust fallback method also failed for %s: %s", sample, e2)
                continue

        if adata is None:
            continue

        adata.obs["sampleID"] = sample
        if sample in sample_metadata:
            for meta_key, meta_value in sample_metadata[sample].items():
                adata.obs[meta_key] = meta_value
        log.info(
            "Successfully loaded %s: %d cells, %d genes",
            sample,
            adata.n_obs,
            adata.n_vars,
        )
        adata_list.append(adata)
        gc.collect()

    if not adata_list:
        log.error("No samples were loaded successfully.")
        return AnnData()

    log.info("Merging %d samples...", len(adata_list))
    combined = anndata.concat(
        adata_list, join="outer", keys=valid_samples, label="batch", index_unique="_"
    )

    log.info(
        "Combined dataset: %d cells, %d genes", combined.n_obs, combined.n_vars
    )

    if output_file:
        out_path = Path(output_file).expanduser()
        if out_path.exists() and backup_existing:
            backup_path = out_path.with_name(
                f"{out_path.name}.bak.{int(time.time())}"
            )
            log.info("Backing up existing %s to %s", out_path, backup_path)
            out_path.rename(backup_path)
        log.info("Saving combined data to %s", out_path)
        combined.write(out_path, compression=compression)

    _attach_counts_layer(combined)
    _attach_sample_metadata(
        combined,
        sample_id=None,  # Multi-sample mode uses per-sample sampleID column
        species=species,
        tissue=tissue,
        tissue_type=tissue_type,
        cancer_type=cancer_type,
    )

    return combined


def _attach_counts_layer(adata: AnnData) -> None:
    """Populate ``layers["counts"]`` from ``X`` when X looks like raw counts.

    scLucid workflows treat ``layers["counts"]`` as the canonical raw-count
    slot. Missing this layer is a common foot-gun, so we auto-fill it when:

    - the layer is not already present, **and**
    - ``X`` is non-negative integer-valued (heuristic check on a sample).

    Otherwise we log a warning but leave ``X`` alone — the caller is
    responsible for providing real counts.
    """
    if LayerKeys.COUNTS in adata.layers:
        return
    if _looks_like_counts(adata.X):
        adata.layers[LayerKeys.COUNTS] = adata.X.copy()
        log.info(
            "Copied AnnData.X to layers['%s'] (detected integer counts).",
            LayerKeys.COUNTS,
        )
    else:
        log.warning(
            "AnnData has no 'counts' layer and X does not look like raw counts. "
            "Some scLucid workflows require raw counts; consider supplying them "
            "via adata.layers['counts'] before calling run_pipeline()."
        )


def _looks_like_counts(matrix) -> bool:
    """Return True if ``matrix`` appears to hold non-negative integer counts."""
    try:
        if hasattr(matrix, "dtype") and np.issubdtype(matrix.dtype, np.integer):
            return True
        sample = matrix[:64] if hasattr(matrix, "__getitem__") else matrix
        if hasattr(sample, "toarray"):
            sample = sample.toarray()
        sample = np.asarray(sample)
        if sample.size == 0:
            return False
        if np.any(sample < 0):
            return False
        return bool(np.all(sample == sample.astype(int)))
    except Exception:
        return False


def _attach_sample_metadata(
    adata: AnnData,
    *,
    sample_id: Optional[str],
    species: Optional[str],
    tissue: Optional[str],
    tissue_type: Optional[str],
    cancer_type: Optional[str],
) -> None:
    """Stamp sample-level metadata onto ``.obs`` and the analysis context."""
    obs_columns = {
        "sample_id": sample_id,
        "species": species,
        "tissue": tissue,
        "tissue_type": tissue_type,
        "cancer_type": cancer_type,
    }
    for column, value in obs_columns.items():
        if value is None:
            continue
        adata.obs[column] = value

    context_payload = {key: value for key, value in obs_columns.items() if value is not None}
    if not context_payload:
        return

    root = ensure_sclucid_namespace(adata)
    existing = root.get(UnsKeys.ANALYSIS_CONTEXT, {})
    if not isinstance(existing, dict):
        existing = {}
    merged = {**existing, **context_payload}
    root[UnsKeys.ANALYSIS_CONTEXT] = merged
    adata.uns[SCLUCID_ROOT] = root


def _find_sample_paths(
    base_dir: str, samples: List[str], possible_subpaths: Optional[List[str]] = None
) -> Dict[str, str]:
    """
    Find valid 10x Genomics data paths for a list of samples.

    Args:
        base_dir: Base directory containing the sample folders.
        samples: List of sample IDs to search for.
        possible_subpaths: List of possible subdirectory structures to check.
                          Default paths checked are:
                          - "outs/filtered_feature_bc_matrix"
                          - "filtered_feature_bc_matrix"
                          - "" (sample directory itself)

    Returns:
        Dictionary mapping sample IDs to their valid data paths.
    """
    if possible_subpaths is None:
        possible_subpaths = [
            os.path.join("outs", "filtered_feature_bc_matrix"),
            "filtered_feature_bc_matrix",
            "",  # Use sample directory directly
        ]

    found_paths = {}

    for sample in samples:
        found = False
        for subpath in possible_subpaths:
            full_path = os.path.join(base_dir, sample, subpath)
            if os.path.isdir(full_path):
                # Check if directory contains 10x data files
                mtx_files = ["matrix.mtx", "matrix.mtx.gz"]
                feature_files = [
                    "features.tsv",
                    "features.tsv.gz",
                    "genes.tsv",
                    "genes.tsv.gz",
                ]

                has_mtx = any(os.path.exists(os.path.join(full_path, f)) for f in mtx_files)
                has_features = any(
                    os.path.exists(os.path.join(full_path, f)) for f in feature_files
                )

                if has_mtx and has_features:
                    found_paths[sample] = full_path
                    found = True
                    break

        if not found:
            log.warning(f"No valid 10x data path found for sample {sample}")

    return found_paths


def _read_10x_manually(sample_path: str) -> AnnData:
    """
    Manually reads 10x data files as a robust fallback method.

    Args:
        sample_path: Path to the directory containing matrix.mtx.gz, features.tsv.gz, etc.

    Returns:
        An AnnData object.
    """
    log.info(f"Attempting robust manual read from: {sample_path}")

    # --- Find Files (with fallback for different names/compressions) ---
    matrix_file = os.path.join(sample_path, "matrix.mtx.gz")
    if not os.path.exists(matrix_file):
        matrix_file = os.path.join(sample_path, "matrix.mtx")

    features_file = os.path.join(sample_path, "features.tsv.gz")
    if not os.path.exists(features_file):
        features_file = os.path.join(sample_path, "genes.tsv.gz")
    if not os.path.exists(features_file):
        features_file = os.path.join(sample_path, "features.tsv")
    if not os.path.exists(features_file):
        features_file = os.path.join(sample_path, "genes.tsv")

    barcodes_file = os.path.join(sample_path, "barcodes.tsv.gz")
    if not os.path.exists(barcodes_file):
        barcodes_file = os.path.join(sample_path, "barcodes.tsv")

    if not all(os.path.exists(f) for f in [matrix_file, features_file, barcodes_file]):
        raise FileNotFoundError(f"Could not find all required 10x files in {sample_path}")

    # --- Read Files with Explicit Type Control ---
    X = io.mmread(matrix_file).T.tocsr()

    features_df = pd.read_csv(
        features_file,
        sep="\t",
        header=None,
        compression="gzip" if features_file.endswith(".gz") else None,
        dtype=str,  # Crucial: ensure all columns are read as strings
    )
    gene_names = features_df[1] if features_df.shape[1] >= 2 else features_df[0]

    barcodes_df = pd.read_csv(
        barcodes_file,
        sep="\t",
        header=None,
        compression="gzip" if barcodes_file.endswith(".gz") else None,
        dtype=str,  # Crucial: ensure barcodes are strings
    )
    barcodes = barcodes_df[0]

    # --- Create and Sanitize AnnData Object ---
    adata = anndata.AnnData(
        X=X,
        obs=pd.DataFrame(index=barcodes.values),
        var=pd.DataFrame(index=gene_names.values),
    )

    adata.var_names_make_unique()  # Ensure gene names are unique
    adata.layers["counts"] = adata.X.copy()

    if X.shape[0] != len(barcodes):
        raise ValueError(
            f"Mismatch: {X.shape[0]} cells in matrix, " f"but {len(barcodes)} barcodes"
        )

    if X.shape[1] != len(gene_names):
        raise ValueError(
            f"Mismatch: {X.shape[1]} genes in matrix, " f"but {len(gene_names)} gene names"
        )

    # Check for empty matrix
    if X.sum() == 0:
        raise ValueError("Matrix contains no data (all zeros)")

    # Check for genes with zero expression across all cells
    cells_per_gene = (X > 0).sum(axis=0).A1
    if (cells_per_gene == 0).sum() > 0.5 * X.shape[1]:
        log.warning(
            "Over 50% of genes have zero expression. " "Check if matrix is correctly oriented."
        )

    return adata


@contextmanager
def use_layer_as_X(adata: AnnData, layer: Optional[str]):
    """Context manager to temporarily use a layer as adata.X."""
    if layer is None:
        yield
        return

    if layer not in adata.layers:
        log.warning(f"Layer '{layer}' not found in adata.layers. Using adata.X.")
        yield
        return

    X_backup = adata.X.copy()
    adata.X = adata.layers[layer].copy()
    try:
        yield
    finally:
        # Always restore the original .X
        adata.X = X_backup


def sanitize_for_hdf5(obj):
    """
    Make objects HDF5-compatible by:
    1. Converting tuples to lists
    2. Converting integer keys to strings in dictionaries
    3. Handling other non-HDF5 compatible types
    """
    if isinstance(obj, tuple) or isinstance(obj, list):
        return [sanitize_for_hdf5(item) for item in obj]
    elif isinstance(obj, dict):
        return {str(k): sanitize_for_hdf5(v) for k, v in obj.items()}
    elif isinstance(obj, (int, float, str, bool, np.number, np.bool_)) or obj is None:
        return obj
    else:
        # Try to convert other types to string representation
        try:
            return str(obj)
        except:
            return "Unconvertible object"


def subset_adata(
    adata: AnnData,
    filters: Dict[str, Union[Any, List[Any]]],
    keep_raw_genes: bool = True,
    raise_on_empty: bool = True,
) -> AnnData:
    """
    Subset an AnnData object based on metadata criteria, retaining raw gene data.

    This utility function is designed for flexible data exploration, allowing you
    to create a new AnnData object for a specific subset of cells (e.g., one
    sample, or only T cells) that can be re-analyzed from scratch.

    Args:
        adata: The AnnData object to subset.
        filters: Dictionary of metadata filters. Keys are column names in `adata.obs`,
                 and values are the desired value or a list of desired values.
                 Example: `{"sampleID": "sample1", "cell_type": ["T cells", "B cells"]}`
        keep_raw_genes: If True and `adata.raw` exists, the returned object's `.raw`
                        attribute will contain the subset of cells but the full original
                        set of genes, enabling re-running of HVG selection.

    Returns:
        A new, subsetted AnnData object.
    """
    if not isinstance(filters, dict):
        raise TypeError("filters must be a dictionary.")

    initial_cells = adata.n_obs
    combined_mask = pd.Series(True, index=adata.obs_names)

    for key, value in filters.items():
        if key not in adata.obs.columns:
            log.warning(f"Metadata column '{key}' not found in adata.obs. Skipping filter.")
            continue

        if isinstance(value, list):
            mask = adata.obs[key].isin(value)
        else:
            mask = adata.obs[key] == value

        combined_mask &= mask

    final_cells = combined_mask.sum()
    log.info("Subsetting data based on provided filters:")
    log.info(f"  - Initial cells: {initial_cells}")
    log.info(f"  - Final cells after filtering: {final_cells}")

    if final_cells == 0:
        msg = f"No cells remaining after applying filters: {filters}"
        if raise_on_empty:
            raise ValueError(msg)
        else:
            log.warning(msg)
            return AnnData()

    if final_cells < 10:
        log.warning(f"Only {final_cells} cells remaining. " "Results may be unreliable.")

    # The core slicing operation
    adata_subset = adata[combined_mask, :].copy()

    # AnnData slicing automatically handles .raw correctly. If we want to ensure
    # the .raw attribute uses the original var, we can explicitly re-assign it.
    if keep_raw_genes and adata.raw is not None:
        # Create a new raw object from the original raw data, but with subsetted cells
        adata_subset.raw = adata.raw[adata_subset.obs_names, :].copy()
        log.info(f"Subset .raw created, retaining all {adata.raw.n_vars} original genes.")

    return adata_subset


def subset_from_annotations(
    adata_target: AnnData,
    adata_source: AnnData,
    filters: Dict[str, Union[Any, List[Any]]],
    columns_to_merge: Union[str, List[str]],
) -> AnnData:
    """
    Subsets a target AnnData object based on annotations from a source object.

    This is a convenience wrapper for the common sub-clustering workflow where
    annotations (e.g., cell types) are generated on a processed object but the
    subsetting needs to be done on an unprocessed object (e.g., containing
    raw counts for all genes).

    Args:
        adata_target: The AnnData object to be subsetted (e.g., the QC'd object).
        adata_source: The AnnData object containing the annotations in its .obs.
                      Must share the same cell indices as adata_target.
        filters: Dictionary of metadata filters to apply. The keys must be present
                 in the `columns_to_merge`.
        columns_to_merge: A column name or list of column names from `adata_source.obs`
                          to merge into `adata_target.obs` before filtering.

    Returns:
        A new, subsetted AnnData object.
    """
    if isinstance(columns_to_merge, str):
        columns_to_merge = [columns_to_merge]

    # --- Step 1: Merge Annotations ---
    log.info(f"Merging annotations for columns: {columns_to_merge} from source object.")

    # Check if columns exist in the source
    missing_cols = [col for col in columns_to_merge if col not in adata_source.obs]
    if missing_cols:
        raise ValueError(f"Columns {missing_cols} not found in the source AnnData object's .obs")

    annotations = adata_source.obs[columns_to_merge]

    # Use a temporary DataFrame to avoid modifying the original adata_target.obs in case of error
    obs_merged = adata_target.obs.join(annotations)

    # Validate that all cells were matched
    if obs_merged[columns_to_merge[0]].isnull().any():
        unmatched_count = obs_merged[columns_to_merge[0]].isnull().sum()
        log.warning(
            f"Found {unmatched_count} cells in the target object that were not present "
            "in the source object's annotations. These will not be selected."
        )

    # Create a temporary AnnData object with the merged obs for filtering
    temp_adata = adata_target.copy()
    temp_adata.obs = obs_merged

    # --- Step 2: Subset ---
    log.info(f"Subsetting target object based on new annotations with filters: {filters}")

    # Now we can call the original, simple subset_adata function
    adata_subset = subset_adata(temp_adata, filters=filters)

    return adata_subset


def merge_obs_metadata(
    adata: AnnData,
    metadata_path: str,
    left_on: Optional[str] = None,  # If None, uses adata.obs.index
    right_on: Optional[str] = None,  # If None, uses metadata_df.index
    how: str = "left",
    handle_duplicates: str = "warn",  # 'warn', 'error', 'overwrite'
) -> AnnData:
    """
    Merges metadata from an external file into the AnnData object's .obs DataFrame.

    Args:
        adata: The AnnData object to modify.
        metadata_path: Path to the metadata file (.csv, .tsv, or .xlsx).
        left_on: Column in adata.obs to join on. If None, uses the index (cell barcodes).
        right_on: Column in the external file to join on. If None, uses the index.
        how: How to perform the merge (e.g., 'left', 'inner'). Defaults to 'left'.

    Returns:
        The AnnData object with merged metadata (modified in place).
    """
    log.info(f"Loading metadata from {metadata_path}")
    if metadata_path.endswith(".csv"):
        meta_df = pd.read_csv(metadata_path)
    elif metadata_path.endswith((".xlsx", ".xls")):
        meta_df = pd.read_excel(metadata_path)
    elif metadata_path.endswith(".tsv"):
        meta_df = pd.read_csv(metadata_path, sep="\t")
    else:
        raise ValueError("Unsupported file format. Please use .csv, .tsv, or .xlsx.")

    initial_cols = set(adata.obs.columns)
    meta_cols = set(meta_df.columns if right_on is None else meta_df.columns.drop(right_on))

    overlapping = initial_cols & meta_cols

    if overlapping:
        if handle_duplicates == "error":
            raise ValueError(f"Columns already exist in adata.obs: {overlapping}")
        elif handle_duplicates == "warn":
            log.warning(
                f"Columns {overlapping} already exist. " f"New columns will be suffixed with '_new'"
            )
            suffixes = ("", "_new")
        elif handle_duplicates == "overwrite":
            log.info(f"Overwriting columns: {overlapping}")
            # Drop existing columns before merge
            adata.obs.drop(columns=overlapping, inplace=True)
            suffixes = ("", "")
    else:
        suffixes = ("", "")

    # Perform merge
    if left_on is None:
        adata.obs = adata.obs.join(
            meta_df.set_index(right_on) if right_on else meta_df,
            how=how,
            rsuffix="_new" if handle_duplicates == "warn" else "",
        )
    else:
        adata.obs = adata.obs.merge(
            meta_df, left_on=left_on, right_on=right_on, how=how, suffixes=suffixes
        )

    return adata
