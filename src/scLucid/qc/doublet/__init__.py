"""Doublet detection subpackage."""

from .core import *  # noqa: F401, F403
from .ensemble import (  # noqa: F401
    DoubletEvidenceProfiler,
    predict_doublets,
    predict_doublets_with_profiling,
)

# Re-export internal symbols used by tests and external consumers.
from .algorithms import (  # noqa: F401
    _run_doubletdetection,
    _run_scrublet,
    _run_solo,
)
from .core import (  # noqa: F401
    FINAL_PRED_COL,
    HEURISTIC_PRED_COL,
    HEURISTIC_SCORE_COL,
    LINEAGE_SCORES_KEY,
    _create_doublet_marker_config_from_manager,
)
from .ensemble import _export_doublet_stats  # noqa: F401
from .heuristic import _run_heuristic  # noqa: F401
