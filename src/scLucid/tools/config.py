"""
Configuration classes for the analysis module of scLucid.
"""

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

__all__ = [
    
]


@dataclass
class InferCNVConfig:
    ref_obs: str = "cell_type"
    ref_keys: Union[str, List[str]] = "Immune"
    window_size: int = 250
    plot_heatmap: bool = True
    find_tumor_cells: bool = True