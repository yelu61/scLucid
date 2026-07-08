#!/usr/bin/env python3
"""Compatibility wrapper for the renamed QC evidence package builder."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.qc.build_qc_evidence_package import (
    DEFAULT_INPUT_ROOT,
    DEFAULT_OUTPUT_DIR,
    _source_rows_from_figure_table,
    build_package,
    main,
)

__all__ = [
    "DEFAULT_INPUT_ROOT",
    "DEFAULT_OUTPUT_DIR",
    "_source_rows_from_figure_table",
    "build_package",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
