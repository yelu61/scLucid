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
def test_dwls_is_exported():
    # pyDWLS is a first-class scLucid port; the DWLS class must be on the
    # top-level tools namespace and appear in __all__.
    assert hasattr(tools, "DWLS"), "scLucid.tools.DWLS is missing"
    assert "DWLS" in tools.__all__
