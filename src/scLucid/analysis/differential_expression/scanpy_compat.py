"""Scanpy version compatibility layer for differential expression analysis.

Centralises column-name normalisation and percentage-scale handling so that
code consuming DE results does not need to know which Scanpy version
produced them.  This is the single source of truth for DE compatibility
glue code.
"""

import logging

import pandas as pd

log = logging.getLogger(__name__)


def _is_0_to_1(series: pd.Series) -> bool:
    """Check if a pandas Series is on a 0-1 scale."""
    s = series.dropna()
    return s.empty or ((s.min() >= 0.0) and (s.max() <= 1.0))


def _to_frac(series: pd.Series) -> pd.Series:
    """Convert Series from 0-100 or 0-1 range to 0-1 fraction scale.

    Handles:
    - None/empty Series
    - Already scaled 0-1 data
    - Percentage scale (0-100)
    - Invalid values (converted to NaN)
    """
    if series is None:
        return pd.Series(dtype=float)

    s_numeric = pd.to_numeric(series, errors="coerce")

    if _is_0_to_1(s_numeric):
        return s_numeric

    log.debug("Detected percentage scale (0-100). Converting to fraction (0-1).")
    return s_numeric.clip(lower=0, upper=100) / 100.0


def standardize_pct_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize percentage column names across different Scanpy versions.

    Converts:
    - pct.1 / pct_nz_group -> pct_nz_group (0-1 scale)
    - pct.2 / pct_nz_reference -> pct_nz_reference (0-1 scale)

    Args:
        df: DataFrame with percentage columns

    Returns:
        DataFrame with standardized column names and values
    """
    df = df.copy()

    # Handle group column
    if "pct_nz_group" in df.columns:
        df["pct_nz_group"] = _to_frac(df["pct_nz_group"])
    elif "pct.1" in df.columns:
        df["pct_nz_group"] = _to_frac(df["pct.1"])
        df.drop(columns=["pct.1"], inplace=True)

    # Handle reference column
    if "pct_nz_reference" in df.columns:
        df["pct_nz_reference"] = _to_frac(df["pct_nz_reference"])
    elif "pct.2" in df.columns:
        df["pct_nz_reference"] = _to_frac(df["pct.2"])
        df.drop(columns=["pct.2"], inplace=True)

    return df


def standardize_enrichment_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Standardize enrichment result column names across GSEApy/Enrichr outputs.

    Args:
        df: Enrichment results DataFrame

    Returns:
        DataFrame with standardized column names
    """
    if df is None or df.empty:
        return pd.DataFrame()

    df = df.copy()

    rename_map = {
        "Adjusted P-value": "pval_adj",
        "P-value": "pval",
        "Overlap": "overlap",
        "Genes": "genes",
        "Gene_set": "gene_set",
        "Term": "term",
        "NES": "nes",
        "FDR q-val": "fdr",
        "FWER p-val": "fwer_pval",
    }

    df.rename(columns=rename_map, inplace=True)

    # Ensure pval_adj exists
    if "pval_adj" not in df.columns and "fdr" in df.columns:
        df["pval_adj"] = df["fdr"]

    return df


__all__ = [
    "standardize_pct_columns",
    "standardize_enrichment_cols",
    "_to_frac",
    "_is_0_to_1",
]
