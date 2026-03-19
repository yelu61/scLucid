"""
Test fixtures for scLucid.

This module provides utilities and data loaders for testing.
"""

# Keep heavy scanpy-backed loaders lazy to avoid import-time dependency failures.
def load_test_data(*args, **kwargs):
    from .data_loader import load_test_data as _load_test_data
    return _load_test_data(*args, **kwargs)


def get_test_config(*args, **kwargs):
    from .data_loader import get_test_config as _get_test_config
    return _get_test_config(*args, **kwargs)


def list_available_datasets(*args, **kwargs):
    from .data_loader import list_available_datasets as _list_available_datasets
    return _list_available_datasets(*args, **kwargs)

__all__ = [
    "load_test_data",
    "get_test_config",
    "list_available_datasets",
]
