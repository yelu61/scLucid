"""
Global settings for the scLucid toolkit.
Handles logging and plotting defaults.
"""

import logging
from copy import deepcopy
from typing import Any, Dict, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as _scanpy

log = logging.getLogger(__name__)

# --- Backup original matplotlib settings for restoration ---
_original_rc_params = mpl.rcParams.copy()


# --- Backup original scanpy settings for restoration ---
def _get_initial_scanpy_params():
    # For newer scanpy versions
    return {
        "dpi": getattr(_scanpy.settings, "figdir_dpi", 80),
        "dpi_save": getattr(_scanpy.settings, "figdir_dpi_save", 150),
        "figsize": getattr(_scanpy.settings, "figsize", (8, 6)),
        "facecolor": getattr(_scanpy.settings, "facecolor", "white"),
        "autoshow": getattr(_scanpy.settings, "autoshow", True),
        "autosave": getattr(_scanpy.settings, "autosave", False),
        "verbosity": getattr(_scanpy.settings, "verbosity", 1),
    }


_original_scanpy_params = _get_initial_scanpy_params()

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

# --- Academic Font Style Definitions ---
# 三种学术期刊常用的字体风格
_ACADEMIC_FONT_STYLES = {
    "nature": {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica"],
        "font.style": "normal",
        "font.variant": "normal",
        "font.weight": "normal",
        "font.stretch": "normal",
        "description": "Nature/Science style: Clean sans-serif Arial font",
    },
    "cell": {
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.style": "normal",
        "font.variant": "normal",
        "font.weight": "medium",
        "font.stretch": "normal",
        "description": "Cell Press style: Classic Helvetica font",
    },
    "traditional": {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
        "font.style": "normal",
        "font.variant": "normal",
        "font.weight": "normal",
        "font.stretch": "normal",
        "description": "Traditional academic: Classic Times New Roman serif font",
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
    modules: Optional[list] = None,
    log_format: str = "%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
) -> None:
    """
    Configure the root logger for the scLucid toolkit.
    """
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {level}")

    root_logger = logging.getLogger()

    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    handlers = [logging.StreamHandler()]
    if file_path:
        handlers.append(logging.FileHandler(file_path))

    from .config import _config

    _config.verbosity = {"WARNING": 0, "INFO": 1, "DEBUG": 2}[level]

    logging.basicConfig(
        level=numeric_level,
        format=log_format,
        handlers=handlers,
        force=True,
    )

    logging.getLogger("harmonypy").setLevel(logging.ERROR)
    logging.getLogger("scvi").setLevel(logging.WARNING)
    logging.getLogger("numba").setLevel(logging.WARNING)

    log.info(f"scLucid logging configured to level {level}.")


def _configure_scanpy_directly(dpi, dpi_save, figsize, facecolor):
    """
    Configure scanpy settings directly without using set_figure_params method.
    This avoids recursion issues in newer scanpy versions.
    """
    # Set scanpy settings directly
    _scanpy.settings.verbosity = 3
    _scanpy.settings.autoshow = True
    _scanpy.settings.autosave = False

    # Set figure parameters directly
    _scanpy.settings.figdir_dpi = dpi
    _scanpy.settings.figdir_dpi_save = dpi_save
    _scanpy.settings.figsize = figsize

    # Attempt to set facecolor if the attribute exists
    if hasattr(_scanpy.settings, "facecolor"):
        _scanpy.settings.facecolor = facecolor

    # Handle scanpy plotting RC params directly if available
    if hasattr(_scanpy, "plotting") and hasattr(_scanpy.plotting, "rcParams"):
        _scanpy.plotting.rcParams["figure.figsize"] = figsize
        _scanpy.plotting.rcParams["figure.dpi"] = dpi
        _scanpy.plotting.rcParams["savefig.dpi"] = dpi_save
        if "figure.facecolor" in _scanpy.plotting.rcParams:
            _scanpy.plotting.rcParams["figure.facecolor"] = facecolor


def set_figure_params(
    dpi: int = 100,
    dpi_save: int = 300,
    figsize: Tuple[float, float] = (8, 6),
    style: Optional[str] = None,
    style_dict: Optional[Dict[str, Any]] = None,
    color_theme: str = "default",
    font_style: Optional[str] = None,
) -> None:
    """
    Set global plotting parameters for matplotlib and scanpy.

    Args:
        dpi (int): Resolution for figure display.
        dpi_save (int): Resolution for figure saving.
        figsize (tuple): Default figure size in inches.
        style (str): Matplotlib style to use (e.g., 'default', 'classic', 'seaborn-v0_8').
        style_dict (dict): Custom style parameters to override defaults.
        color_theme (str): Color theme to use. Options: 'default', 'dark', 'seaborn'.
        font_style (str): Academic font style to use. Options:
            - 'nature': Nature/Science style (Arial, clean sans-serif)
            - 'cell': Cell Press style (Helvetica, classic sans-serif)
            - 'traditional': Traditional academic (Times New Roman, serif)
            - None: Use default font settings
    """
    log.info("Applying global plotting settings...")

    # Validate color_theme
    if color_theme not in _THEMES:
        log.warning(
            f"Color theme '{color_theme}' not recognized. Using 'default' theme. "
            f"Available: {list(_THEMES.keys())}"
        )
        color_theme = "default"

    # Validate and apply font_style
    if font_style is not None:
        if font_style not in _ACADEMIC_FONT_STYLES:
            log.warning(
                f"Font style '{font_style}' not recognized. Using default font. "
                f"Available academic styles: {list(_ACADEMIC_FONT_STYLES.keys())}"
            )
            font_style = None
        else:
            log.info(
                f"Applying academic font style: '{font_style}' - "
                f"{_ACADEMIC_FONT_STYLES[font_style]['description']}"
            )

    # Configure scanpy directly to avoid recursion issues
    _configure_scanpy_directly(
        dpi=dpi,
        dpi_save=dpi_save,
        figsize=figsize,
        facecolor=_THEMES[color_theme]["axes.facecolor"],
    )

    # Prepare base matplotlib styles
    applied_styles = deepcopy(_DEFAULT_STYLES)
    if style:
        try:
            plt.style.use(style)
            log.info(f"Applied matplotlib style: {style}")
        except Exception as e:
            log.warning(f"Could not apply style '{style}'. Error: {e}. Falling back to default.")
            plt.style.use("default")

    # Apply academic font style if specified
    if font_style:
        font_config = _ACADEMIC_FONT_STYLES[font_style]
        # Remove description field (not a matplotlib parameter)
        font_params = {k: v for k, v in font_config.items() if k != "description"}
        applied_styles.update(font_params)

    # Update with selected color theme
    applied_styles.update(_THEMES[color_theme])
    applied_styles["figure.figsize"] = figsize
    applied_styles["figure.dpi"] = dpi
    applied_styles["savefig.dpi"] = dpi_save

    # Update with any user-provided custom styles
    if style_dict:
        applied_styles.update(style_dict)
        log.info(f"Applied custom style dict: {style_dict}")

    plt.rcParams.update(applied_styles)

    # Set font type for PDF/EPS output to embed fonts (avoid Type 3 fonts)
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42

    # Configure IPython display format if in Jupyter
    try:
        from IPython import get_ipython

        ipython = get_ipython()
        if ipython is not None:
            ipython.run_line_magic("matplotlib", "inline")
            ipython.run_line_magic("config", "InlineBackend.figure_format = 'retina'")
            log.info("IPython inline backend set to retina.")
    except Exception as e:
        log.debug(f"Could not set IPython display format automatically: {e}")

    font_info = f" '{font_style}' font style" if font_style else " default font"
    log.info(f"Global plotting settings applied with '{color_theme}' color theme and{font_info}.")


def reset_figure_params() -> None:
    """
    Reset matplotlib and scanpy settings to their original default state.
    """
    log.info("Resetting plotting settings to matplotlib and scanpy defaults.")
    mpl.rcParams.update(_original_rc_params)

    # Reset scanpy settings directly
    _configure_scanpy_directly(
        dpi=_original_scanpy_params["dpi"],
        dpi_save=_original_scanpy_params["dpi_save"],
        figsize=_original_scanpy_params["figsize"],
        facecolor=_original_scanpy_params["facecolor"],
    )

    _scanpy.settings.autoshow = _original_scanpy_params["autoshow"]
    _scanpy.settings.autosave = _original_scanpy_params["autosave"]
    _scanpy.settings.verbosity = _original_scanpy_params["verbosity"]
    log.info("Plotting settings reset.")


# --- Extra: Warn user if scanpy import path is suspicious, to avoid local conflicts ---
def _check_scanpy_path():
    scanpy_path = getattr(_scanpy, "__file__", "UNKNOWN")
    if "scLucid" in scanpy_path:
        log.warning(
            f"WARNING: Detected scanpy is imported from '{scanpy_path}', "
            "which may indicate a package name conflict. Please check your sys.path "
            "and avoid naming your local modules or packages as 'scanpy'."
        )


_check_scanpy_path()
