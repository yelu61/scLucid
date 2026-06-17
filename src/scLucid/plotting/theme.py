"""Publication-ready plotting themes."""

from typing import Dict, Iterable, Optional, Sequence

import matplotlib.pyplot as plt
import seaborn as sns

FONT_FAMILY = ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"]

# Nature style
NATURE_THEME = {
    "font.family": "sans-serif",
    "font.sans-serif": FONT_FAMILY,
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.linewidth": 1.0,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 300,
    "pdf.fonttype": 42,  # 保证导出PDF文字可编辑
    "ps.fonttype": 42,
}

NATURE_COLORS = {
    "palette": [
        "#E64B35",
        "#4DBBD5",
        "#00A087",
        "#3C5488",
        "#F39B7F",
        "#8491B4",
        "#91D1C2",
        "#DC0000",
        "#7E6148",
    ],
    "cmap": "RdYlBu_r",
}

# Science style
SCIENCE_THEME = {
    "font.family": "sans-serif",
    "font.sans-serif": FONT_FAMILY,
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "axes.linewidth": 1.0,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "figure.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

SCIENCE_COLORS = {
    "palette": [
        "#3C5488",
        "#F39B7F",
        "#8491B4",
        "#91D1C2",
        "#DC0000",
        "#00A087",
        "#E64B35",
        "#4DBBD5",
    ],
    "cmap": "viridis",
}


def build_color_palette(
    keys,
    color_maps: Optional[Dict[str, Dict[str, str]]] = None,
    fallback_palette: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """Build a combined color palette by looking up keys across color maps.

    Searches each key in the provided color maps in order and falls back to
    ``fallback_palette`` if none of the maps define the key.

    Parameters
    ----------
    keys
        Iterable of keys to resolve (e.g. sample IDs, group names).
    color_maps
        Ordered mapping ``{"map_name": {key: color}}``. Later maps override
        earlier maps when a key appears in multiple maps. If None, defaults to
        empty dict.
    fallback_palette
        Optional final fallback dict used when a key is not found in any
        ``color_maps``.

    Returns:
    -------
    Dict[str, str]
        Resolved palette with one color per key.

    Examples:
    --------
    >>> palette = build_color_palette(
    ...     ["sample_A", "group_1"],
    ...     color_maps={
    ...         "samples": {"sample_A": "#1E688D"},
    ...         "groups": {"group_1": "#55E08A"},
    ...     },
    ... )
    """
    color_maps = color_maps or {}
    palette: Dict[str, str] = {}
    for key in keys:
        for color_map in color_maps.values():
            if key in color_map:
                palette[key] = color_map[key]
                break
        else:
            if fallback_palette and key in fallback_palette:
                palette[key] = fallback_palette[key]
    return palette


def build_obs_palette(
    adata,
    color_keys: Sequence[str],
    color_maps: Optional[Dict[str, Dict[str, str]]] = None,
    fallback_palette: Optional[Iterable[str]] = None,
) -> Dict[str, str]:
    """Build a Scanpy palette for categorical values in ``adata.obs`` columns.

    ``scanpy.pl.embedding(..., color=[...], palette=...)`` expects colors keyed by
    category values, not by obs column names. This helper collects the displayed
    categories from each requested obs column and resolves them against project
    color maps, filling any missing values with a deterministic fallback palette.
    """
    color_maps = color_maps or {}
    fallback_colors = list(fallback_palette or sns.color_palette("tab20", 20).as_hex())
    palette: Dict[str, str] = {}
    fallback_i = 0

    categories = []
    for key in color_keys:
        if key not in adata.obs.columns:
            continue
        values = adata.obs[key].dropna().astype(str).unique().tolist()
        categories.extend(values)

    for category in dict.fromkeys(categories):
        for color_map in color_maps.values():
            if category in color_map:
                palette[category] = color_map[category]
                break
        else:
            palette[category] = fallback_colors[fallback_i % len(fallback_colors)]
            fallback_i += 1

    return palette


def apply_theme(theme: str = "nature") -> Dict:
    """
    Apply a publication-ready theme.

    Parameters:
        theme: 'nature' or 'science'

    Returns:
        Dict containing color configurations
    """
    if theme.lower() == "nature":
        plt.rcParams.update(NATURE_THEME)
        sns.set_palette(NATURE_COLORS["palette"])
        return NATURE_COLORS
    elif theme.lower() == "science":
        plt.rcParams.update(SCIENCE_THEME)
        sns.set_palette(SCIENCE_COLORS["palette"])
        return SCIENCE_COLORS
    else:
        # Fallback to defaults but try to optimize DPI
        plt.rcParams["figure.dpi"] = 300
        return {"palette": "tab10", "cmap": "viridis"}
