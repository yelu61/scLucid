"""Public API checks for scLucid.tools with optional backend boundaries."""

import pytest

import scLucid.tools as tools


@pytest.mark.unit
@pytest.mark.optional
def test_tools_exports_resolve():
    for symbol in tools.__all__:
        assert hasattr(tools, symbol), f"scLucid.tools missing exported symbol: {symbol}"


@pytest.mark.unit
@pytest.mark.optional
def test_missing_optional_symbol_is_not_exported():
    # pyDWLS currently has incomplete internals in some environments; if it fails to import,
    # the top-level tools module should not claim DWLS in __all__.
    if hasattr(tools, "DWLS"):
        assert "DWLS" in tools.__all__
    else:
        assert "DWLS" not in tools.__all__
