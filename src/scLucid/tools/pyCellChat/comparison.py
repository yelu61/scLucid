"""
Comparison analysis between CellChat objects (R-free)
"""

import logging
from typing import Dict, List

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def compare_cellchat_objects(cellchat1, cellchat2, comparison_type: str = "functional") -> Dict:
    """Compare two CellChat objects"""
    from .analysis import compute_network_similarity

    similarity = compute_network_similarity(cellchat1, cellchat2, type=comparison_type)
    diff_pathways = identify_differential_pathways([cellchat1, cellchat2])
    conserved = identify_conserved_pathways([cellchat1, cellchat2])

    return {
        "similarity": similarity,
        "differential_pathways": diff_pathways,
        "conserved_pathways": conserved,
    }


def identify_differential_pathways(cellchat_list: List, thresh: float = 0.05) -> pd.DataFrame:
    """Identify differential pathways across conditions"""
    all_pathways = set(cellchat_list[0].netP["prob"].keys())
    for cellchat in cellchat_list[1:]:
        all_pathways &= set(cellchat.netP["prob"].keys())

    results = []
    for pathway in all_pathways:
        strengths = [cc.netP["prob"][pathway].sum() for cc in cellchat_list]

        if len(strengths) == 2:
            fc = strengths[1] / (strengths[0] + 1e-10)
            log_fc = np.log2(fc)
            results.append(
                {
                    "pathway": pathway,
                    "condition1_strength": strengths[0],
                    "condition2_strength": strengths[1],
                    "log2_fc": log_fc,
                    "significant": abs(log_fc) > 1,
                }
            )

    return pd.DataFrame(results).sort_values("log2_fc", ascending=False)


def identify_conserved_pathways(cellchat_list: List, correlation_thresh: float = 0.7) -> List[str]:
    """Identify conserved pathways across conditions"""
    all_pathways = set(cellchat_list[0].netP["prob"].keys())
    for cellchat in cellchat_list[1:]:
        all_pathways &= set(cellchat.netP["prob"].keys())

    conserved = []
    for pathway in all_pathways:
        matrices = [cc.netP["prob"][pathway].flatten() for cc in cellchat_list]
        correlations = [
            np.corrcoef(matrices[i], matrices[j])[0, 1]
            for i in range(len(matrices))
            for j in range(i + 1, len(matrices))
        ]

        if all(c > correlation_thresh for c in correlations):
            conserved.append(pathway)

    return conserved
