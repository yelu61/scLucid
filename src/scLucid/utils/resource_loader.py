"""
Unified resource loader for scLucid.
Handles loading of internal data (TOML, JSON, GMT, CSV) robustly.
"""

import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict

# Compat for python < 3.9
if sys.version_info < (3, 9):
    import importlib_resources as resources
else:
    from importlib import resources

# Optional dependencies
try:
    import tomllib
except ImportError:
    import tomli as tomllib  # Fallback

log = logging.getLogger(__name__)

RESOURCE_PACKAGE = "scLucid.resources"


def _resolve_resource(filename: str):
    """Resolve a package resource from a slash-delimited relative path."""
    relative = filename.replace("\\", "/").strip("/")
    try:
        traversable = resources.files(RESOURCE_PACKAGE)
        for part in relative.split("/"):
            traversable = traversable.joinpath(part)
        return traversable
    except (AttributeError, ModuleNotFoundError):
        return None


def get_resource_path(filename: str) -> Path:
    """
    Get the absolute path to a resource file.
    Note: In zipped environments, this might return a temporary path.
    """
    try:
        traversable = _resolve_resource(filename)
        if traversable is None or not traversable.is_file():
            raise FileNotFoundError(filename)
        with resources.as_file(traversable) as path:
            return path
    except (ImportError, FileNotFoundError):
        # Fallback for dev mode (running from source without install)
        # Assuming struct: src/scLucid/utils/resource_loader.py -> src/scLucid/resources/
        dev_path = Path(__file__).parent.parent / "resources" / filename.replace("\\", "/")
        if dev_path.exists():
            return dev_path
        raise FileNotFoundError(f"Resource '{filename}' not found in package '{RESOURCE_PACKAGE}'")


def resource_exists(filename: str) -> bool:
    """Return True when a package/dev resource exists."""
    try:
        get_resource_path(filename)
        return True
    except FileNotFoundError:
        return False


def load_toml(filename: str) -> Dict[str, Any]:
    """Load a TOML file from resources."""
    path = get_resource_path(filename)
    with open(path, "rb") as f:
        return tomllib.load(f)


def load_json(filename: str) -> Dict[str, Any]:
    """Load a JSON file from resources."""
    path = get_resource_path(filename)
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def load_gmt(filename: str) -> Dict[str, list]:
    """Load a GMT file into a dictionary."""
    path = get_resource_path(filename)
    gene_sets = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            # parts[1] is usually description
            genes = [g for g in parts[2:] if g]
            gene_sets[name] = genes
    return gene_sets
