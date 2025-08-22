"""
Global settings for the scLucid toolkit.
Handles logging and plotting defaults.
"""

import logging
from typing import Any, Dict, Optional
from copy import deepcopy

import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc

log = logging.getLogger(__name__)

# --- Store original matplotlib settings for restoration ---
_original_rc_params = mpl.rcParams.copy()

# --- Plotting Style Definitions (Module-level constants) ---

_THEMES = {
    "default": {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "text.color": "black",
        "axes.labelcolor": "black",
        "axes.edgecolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
        "grid.color": "lightgray",
    },
    "dark": {
        "figure.facecolor": "#303030",
        "axes.facecolor": "#303030",
        "savefig.facecolor": "#303030",
        "text.color": "white",
        "axes.labelcolor": "white",
        "axes.edgecolor": "white",
        "xtick.color": "white",
        "ytick.color": "white",
        "grid.color": "#505050",
    },
    "seaborn": {
        "figure.facecolor": "white",
        "axes.facecolor": "#F0F0F0",
        "savefig.facecolor": "white",
        "text.color": "#283747",
        "axes.labelcolor": "#283747",
        "axes.edgecolor": "#283747",
        "xtick.color": "#283747",
        "ytick.color": "#283747",
        "grid.color": "white",
    },
}

_DEFAULT_STYLES = {
    # Figure properties
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.transparent": False,
    # Fonts and text
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "legend.title_fontsize": 12,
    # Axes properties
    "axes.linewidth": 1.0,
    "axes.grid": False,
    "axes.axisbelow": True,
    "axes.spines.top": True,
    "axes.spines.right": True,
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    # Ticks
    "xtick.major.size": 4,
    "xtick.minor.size": 2,
    "xtick.major.width": 1.0,
    "xtick.direction": "out",
    "ytick.major.size": 4,
    "ytick.minor.size": 2,
    "ytick.major.width": 1.0,
    "ytick.direction": "out",
    # Legend properties
    "legend.frameon": False,
    # Line properties
    "lines.linewidth": 1.5,
    "lines.markersize": 6,
}


__all__ = ["setup_logging", "set_figure_params", "reset_figure_params"]


def setup_logging(
    level: str = "INFO",
    file_path: Optional[str] = None,
    log_format: str = "%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
):
    """
    Configures the root logger for the scLucid toolkit.

    Args:
        level: The logging level (e.g., 'DEBUG', 'INFO', 'WARNING').
        file_path: Optional path to a file to save logs.
        log_format: The format string for the log messages.
    """
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {level}")

    handlers = [logging.StreamHandler()]
    if file_path:
        handlers.append(logging.FileHandler(file_path))

    logging.basicConfig(
        level=numeric_level,
        format=log_format,
        handlers=handlers,
        force=True,  # Override any existing configuration
    )
    log.info(f"scLucid logging configured to level {level}.")


def set_figure_params(
    dpi: int = 100,
    dpi_save: int = 300,
    figsize: tuple = (6, 6),
    style: str = None,
    style_dict: Optional[Dict[str, Any]] = None,
    color_theme: str = "default",
):
    """
    Set global plotting parameters for matplotlib and scanpy.

    Args:
        dpi (int): Resolution for figure display.
        dpi_save (int): Resolution for figure saving.
        figsize (tuple): Default figure size in inches.
        style (str): Matplotlib style to use (e.g., 'default', 'classic', 'seaborn-v0_8').
        style_dict (dict): Custom style parameters to override defaults.
        color_theme (str): Color theme to use. Options: 'default', 'dark', 'seaborn'.
    """
    print("Applying global plotting settings...")

    # Set scanpy settings first
    sc.settings.verbosity = 3
    sc.settings.autoshow = True
    sc.settings.autosave = False
    sc.settings.set_figure_params(
        scanpy=True,
        dpi=dpi,
        dpi_save=dpi_save,
        figsize=figsize,
        facecolor=_THEMES.get(color_theme, _THEMES["default"]).get("axes.facecolor"),
        # Let scanpy handle its formats, or enforce them here
        # format="svg",
    )

    # Start with a clean copy of the base styles
    applied_styles = deepcopy(_DEFAULT_STYLES)

    # Apply a base matplotlib style if specified
    if style:
        try:
            plt.style.use(style)
        except Exception as e:
            log.warning(f"Could not apply style '{style}'. Error: {e}. Falling back to default.")
            plt.style.use("default")

    # Update with the selected color theme
    theme_params = _THEMES.get(color_theme, _THEMES["default"])
    applied_styles.update(theme_params)
    
    # Manually set figsize and dpi as they are direct arguments
    applied_styles["figure.figsize"] = figsize
    applied_styles["figure.dpi"] = dpi
    applied_styles["savefig.dpi"] = dpi_save

    # Update with any user-provided custom styles
    if style_dict:
        applied_styles.update(style_dict)

    # Apply the final settings to matplotlib
    plt.rcParams.update(applied_styles)

    # Configure IPython display format
    try:
        from IPython import get_ipython

        ipython = get_ipython()
        if ipython is not None:
            ipython.run_line_magic("matplotlib", "inline")
            # This is a common way to set formats for inline backend
            ipython.run_line_magic("config", "InlineBackend.figure_format = 'retina'")

    except Exception as e:
        log.debug(f"Could not set IPython display format automatically: {e}")

    print(f"Global plotting settings applied with '{color_theme}' color theme.")


def reset_figure_params():
    """Resets matplotlib and scanpy settings to their original default state."""
    print("Resetting plotting settings to matplotlib defaults.")
    mpl.rcParams.update(_original_rc_params)
    # You might also want to reset scanpy settings to its default
    # sc.settings.set_figure_params() # Re-run with defaults if needed