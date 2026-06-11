"""Backward-compatible high-level entry points for bulk tools."""

from __future__ import annotations

from .abundance import run_bulk_abundance_test
from .deconvolution import deconvolve_bulk

# Preserve legacy aliases from the original scLucid.tools.bulk module.
differential_abundance = run_bulk_abundance_test
run_deconvolution = deconvolve_bulk

__all__ = [
    "differential_abundance",
    "run_deconvolution",
]
