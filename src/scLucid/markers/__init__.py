"""
Cell Type Marker Management System for scLucid.

This module provides a robust system for loading, managing, and querying hierarchical
cell type markers from TOML or JSON files. It is central to the biology-aware
features of the scLucid toolkit.

Key exposed components:
- Manager: The main class for handling marker hierarchies.
- CellType: A dataclass representing a single cell type entry.
- get_marker_manager: A factory function to easily build a combined manager.
- KNOWN_SPECIES: A list of built-in species supported.
- MARKER_FORMATS: A list of supported marker file formats.
"""

# Import and expose key functions and classes from the submodule
from .manager import (
    CellType,
    Manager,
    get_marker_manager,
    KNOWN_SPECIES,
    MARKER_FORMATS,
)

# Define what should be accessible when a user does `from scLucid.markers import *`
__all__ = [
    "CellType",
    "Manager",
    "get_marker_manager",
    "KNOWN_SPECIES",
    "MARKER_FORMATS",
]