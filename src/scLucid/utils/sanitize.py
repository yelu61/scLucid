"""HDF5 serialization helpers for scLucid."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


def sanitize_for_hdf5(obj: Any) -> Any:
    """Make objects HDF5-compatible by cleaning numpy/tuple/None types."""
    if isinstance(obj, pd.DataFrame):
        df = obj.copy()
        for col in df.columns:
            if df[col].dtype == object:
                df[col] = df[col].apply(lambda x: "" if x is None or pd.isna(x) else x)
        return df
    elif isinstance(obj, pd.Series):
        if obj.dtype == object:
            return obj.apply(lambda x: "" if x is None or pd.isna(x) else x).to_dict()
        return obj.to_dict()
    elif isinstance(obj, np.ndarray):
        return sanitize_for_hdf5(obj.tolist())
    elif isinstance(obj, tuple) or isinstance(obj, list):
        if obj and all(isinstance(item, dict) for item in obj):
            return {
                str(i): sanitize_for_hdf5(item)
                for i, item in enumerate(obj)
            }
        return [sanitize_for_hdf5(item) for item in obj]
    elif isinstance(obj, dict):
        return {str(k): sanitize_for_hdf5(v) for k, v in obj.items()}
    elif obj is None:
        return ""
    elif isinstance(obj, np.generic):
        return obj.item()
    elif isinstance(obj, (int, float, str, bool)):
        return obj
    else:
        try:
            return str(obj)
        except Exception:
            return "Unconvertible object"
