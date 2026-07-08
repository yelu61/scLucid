#!/usr/bin/env python3
"""Build a combined QC+preprocess evidence package from executed outputs."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from validation.qc_preprocess.build_review_summary_evidence_tables import (
    _collect_input_paths,
    build_tables,
    main,
    write_outputs,
)

__all__ = ["_collect_input_paths", "build_tables", "write_outputs", "main"]


if __name__ == "__main__":
    main()
