"""Smoke tests for package import stability and public surfaces."""

import importlib
import warnings

import pytest


CORE_MODULES = [
    "scLucid",
    "scLucid.qc",
    "scLucid.preprocess",
    "scLucid.analysis",
    "scLucid.tools",
    "scLucid.config",
]


@pytest.mark.smoke
@pytest.mark.parametrize("module_name", CORE_MODULES)
def test_core_module_imports(module_name):
    """Core package modules should remain importable."""
    module = importlib.import_module(module_name)
    assert module is not None


@pytest.mark.smoke
@pytest.mark.parametrize(
    "module_name",
    ["scLucid", "scLucid.qc", "scLucid.preprocess", "scLucid.analysis", "scLucid.tools", "scLucid.utils"],
)
def test_all_exports_resolve(module_name):
    """Every symbol listed in __all__ must resolve on the module."""
    module = importlib.import_module(module_name)
    exported = getattr(module, "__all__", [])
    assert isinstance(exported, list)
    for symbol in exported:
        assert hasattr(module, symbol), f"{module_name} exports missing symbol: {symbol}"


@pytest.mark.smoke
def test_top_level_aliases_are_consistent():
    """Top-level convenience aliases should mirror full module attributes."""
    import scLucid as scl

    assert scl.pp is scl.preprocess
    assert scl.al is scl.analysis
    assert scl.tl is scl.tools
    assert scl.ut is scl.utils
    assert scl.pl is scl.plotting


@pytest.mark.smoke
def test_optional_import_helper_gracefully_handles_missing_module():
    """Internal optional import helper should return None and warn on missing modules."""
    import scLucid as scl

    if not hasattr(scl, "_import_optional"):
        pytest.skip("No optional import helper exposed")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = scl._import_optional("__module_that_does_not_exist__")

    assert result is None
    assert any("Could not import module" in str(item.message) for item in caught)
