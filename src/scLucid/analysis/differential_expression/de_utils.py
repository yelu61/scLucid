"""
Utility functions and result management for differential expression analysis.

This module provides:
- Helper functions for data standardization
- Result storage utilities
- ResultManager class for unified result management
"""

import logging
import re
from pathlib import Path
from typing import Dict, List, Literal, Union

import pandas as pd
from anndata import AnnData

from ...base_config import SclucidBaseConfig
from ...utils import sanitize_for_hdf5

# Canonical implementations live in scanpy_compat.py; re-exported here for
# backward compatibility.

log = logging.getLogger(__name__)


# ==================== Helper Functions ====================


def _safe_filename(s: str) -> str:
    """Convert string to filesystem-safe filename."""
    return re.sub(r"[^0-9a-zA-Z\u4e00-\u9fa5\-_\.]", "_", s)


def _store_results(
    adata: AnnData,
    key: str,
    results: Union[pd.DataFrame, Dict],
    config: SclucidBaseConfig,
    result_type: str = "de",
) -> None:
    """
    Store analysis results to adata.uns with full provenance.

    Creates structure:
    adata.uns['sclucid']['analysis'][result_type][key] = results
    adata.uns['sclucid']['analysis'][result_type][f'{key}_config'] = config_dict

    Args:
        adata: AnnData object
        key: Storage key for results
        results: Results DataFrame or dictionary
        config: Configuration object used
        result_type: Type of analysis ('de', 'enrichment', etc.)
    """
    # Initialize nested structure
    if "sclucid" not in adata.uns:
        adata.uns["sclucid"] = {}
    if "analysis" not in adata.uns["sclucid"]:
        adata.uns["sclucid"]["analysis"] = {}
    if result_type not in adata.uns["sclucid"]["analysis"]:
        adata.uns["sclucid"]["analysis"][result_type] = {}

    # Store results
    adata.uns["sclucid"]["analysis"][result_type][key] = results

    # Store configuration with HDF5 compatibility
    config_dict = {
        k: v for k, v in config.__dict__.items() if not k.startswith("_") and v is not None
    }
    adata.uns["sclucid"]["analysis"][result_type][f"{key}_config"] = sanitize_for_hdf5(config_dict)

    if config.verbose:
        log.info(f"Results stored: adata.uns['sclucid']['analysis']['{result_type}']['{key}']")


# ==================== Result Management ====================


class ResultManager:
    """Unified result saving and loading manager."""

    SUPPORTED_FORMATS = {
        "csv": ".csv",
        "tsv": ".tsv",
        "excel": ".xlsx",
        "pickle": ".pkl",
        "parquet": ".parquet",
    }

    def __init__(self, base_dir: Union[str, Path]):
        """
        Initialize result manager.

        Args:
            base_dir: Base directory for storing results
        """
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

        self.de_dir = self.base_dir / "differential_expression"
        self.enrichment_dir = self.base_dir / "enrichment"
        self.plots_dir = self.base_dir / "plots"

        for dir_path in [self.de_dir, self.enrichment_dir, self.plots_dir]:
            dir_path.mkdir(exist_ok=True)

    def save_deg_results(
        self,
        results: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
        name: str,
        format: str = "csv",
    ) -> List[Path]:
        """
        Save differential expression results to disk.

        Args:
            results: Results DataFrame or dict of DataFrames
            name: Base filename (without extension)
            format: Output format ('csv', 'tsv', 'excel', 'pickle', 'parquet')

        Returns:
            List of saved file paths
        """
        if format not in self.SUPPORTED_FORMATS:
            raise ValueError(
                f"Unsupported format: {format}. "
                f"Supported: {list(self.SUPPORTED_FORMATS.keys())}"
            )

        ext = self.SUPPORTED_FORMATS[format]
        saved_paths = []

        if isinstance(results, pd.DataFrame):
            results = {"all": results}

        for key, df in results.items():
            filename = _safe_filename(f"{name}_{key}{ext}")
            filepath = self.de_dir / filename

            if format == "csv":
                df.to_csv(filepath, index=False)
            elif format == "tsv":
                df.to_csv(filepath, sep="\t", index=False)
            elif format == "excel":
                df.to_excel(filepath, index=False)
            elif format == "pickle":
                df.to_pickle(filepath)
            elif format == "parquet":
                df.to_parquet(filepath, index=False)

            saved_paths.append(filepath)
            log.info(f"Saved DE results to: {filepath}")

        return saved_paths

    def save_enrichment_results(
        self,
        results: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
        name: str,
        format: str = "csv",
    ) -> List[Path]:
        """
        Save enrichment results to disk.

        Args:
            results: Results DataFrame or dict of DataFrames
            name: Base filename
            format: Output format

        Returns:
            List of saved file paths
        """
        if format not in self.SUPPORTED_FORMATS:
            raise ValueError(
                f"Unsupported format: {format}. "
                f"Supported: {list(self.SUPPORTED_FORMATS.keys())}"
            )

        ext = self.SUPPORTED_FORMATS[format]
        saved_paths = []

        if isinstance(results, pd.DataFrame):
            results = {"all": results}

        for key, df in results.items():
            filename = _safe_filename(f"{name}_{key}{ext}")
            filepath = self.enrichment_dir / filename

            if format == "csv":
                df.to_csv(filepath, index=False)
            elif format == "tsv":
                df.to_csv(filepath, sep="\t", index=False)
            elif format == "excel":
                df.to_excel(filepath, index=False)
            elif format == "pickle":
                df.to_pickle(filepath)
            elif format == "parquet":
                df.to_parquet(filepath, index=False)

            saved_paths.append(filepath)
            log.info(f"Saved enrichment results to: {filepath}")

        return saved_paths


def save_results(
    results: Union[pd.DataFrame, Dict[str, pd.DataFrame]],
    name: str,
    outdir: Union[str, Path],
    result_type: Literal["de", "enrichment"] = "de",
    format: str = "csv",
) -> List[Path]:
    """
    Save analysis results to disk.

    Convenience function that creates a ResultManager and saves results.

    Args:
        results: Results DataFrame or dict of DataFrames
        name: Base filename (without extension)
        outdir: Output directory path
        result_type: Type of results ('de' or 'enrichment')
        format: Output format ('csv', 'tsv', 'excel', 'pickle', 'parquet')

    Returns:
        List of saved file paths

    Example:
        >>> markers = find_markers(adata, config)
        >>> save_results(markers, "cluster_markers", "./results", "de")
    """
    manager = ResultManager(outdir)

    if result_type == "de":
        return manager.save_deg_results(results, name, format)
    elif result_type == "enrichment":
        return manager.save_enrichment_results(results, name, format)
    else:
        raise ValueError(f"Unknown result_type: {result_type}")


def load_results(name: str, outdir: str, format: str = "csv") -> pd.DataFrame:
    """
    Load analysis results from disk.

    Args:
        name: Base filename (without extension or _key suffix)
        outdir: Output directory path
        format: File format ('csv', 'tsv', 'excel', 'pickle', 'parquet')

    Returns:
        Loaded DataFrame

    Example:
        >>> results = load_results("cluster_markers_all", "./results", "csv")
    """
    base_dir = Path(outdir)
    ext = ResultManager.SUPPORTED_FORMATS.get(format, f".{format}")
    filepath = base_dir / "differential_expression" / f"{name}{ext}"

    if not filepath.exists():
        raise FileNotFoundError(f"Results file not found: {filepath}")

    if format == "csv":
        return pd.read_csv(filepath)
    elif format == "tsv":
        return pd.read_csv(filepath, sep="\t")
    elif format == "excel":
        return pd.read_excel(filepath)
    elif format == "pickle":
        return pd.read_pickle(filepath)
    elif format == "parquet":
        return pd.read_parquet(filepath)
    else:
        raise ValueError(f"Unsupported format: {format}")
