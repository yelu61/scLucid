"""
Global settings for the scLucid toolkit.
Handles logging and plotting defaults.
"""

import logging
from typing import Any, Dict, Optional

import matplotlib.pyplot as plt
import scanpy as sc

log = logging.getLogger(__name__)

__all__ = ["setup_logging", "set_figure_params"]


def setup_logging(level: str = "INFO", file_path: Optional[str] = None):
    """
    Configures the root logger for the scLucid toolkit.

    Args:
        level: The logging level (e.g., 'DEBUG', 'INFO', 'WARNING').
        file_path: Optional path to a file to save logs.
    """
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f"Invalid log level: {level}")

    handlers = [logging.StreamHandler()]
    if file_path:
        handlers.append(logging.FileHandler(file_path))

    logging.basicConfig(
        level=numeric_level,
        format="%(asctime)s - %(name)s - [%(levelname)s] - %(message)s",
        handlers=handlers,
        force=True,  # Override any existing configuration
    )
    log.info(f"scLucid logging configured to level {level}.")


def set_figure_params(
    dpi: int = 300,
    figsize: tuple = (6, 6),
    style: str = None,
    style_dict: Optional[Dict[str, Any]] = None,
    color_theme: str = "default",  # 可选: "default", "dark", "seaborn"
):
    """
    Set global plotting parameters for matplotlib and scanpy.

    Args:
        dpi (int): Resolution for figure saving and display.
        figsize (tuple): Default figure size in inches.
        style (str): Matplotlib style to use. Options include: 'default',
                     'classic', 'seaborn', 'ggplot', etc.
        style_dict (dict): Custom style parameters to override defaults.
        color_theme (str): Color theme to use. Options: 'default', 'dark', 'seaborn'.
    """
    print("Applying global plotting settings...")

    # Set scanpy settings
    sc.settings.verbosity = 3
    sc.settings.autoshow = True
    sc.settings.autosave = False
    sc.settings.dpi = dpi
    sc.settings.dpi_save = dpi
    sc.settings.figsize = figsize
    sc.settings.figure_formats = ["png", "svg"]

    # Apply base style if specified
    if style:
        try:
            plt.style.use(style)
        except Exception as e:
            print(f"Warning: Could not apply style '{style}'. Error: {str(e)}")
            print("Falling back to default style.")
            plt.style.use("default")

    # Define color themes
    themes = {
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

    # Apply selected color theme
    theme_params = themes.get(color_theme, themes["default"])

    # Define comprehensive default settings
    default_styles = {
        # Figure properties
        "figure.figsize": figsize,
        "figure.dpi": dpi,
        "savefig.dpi": dpi,
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
        "axes.axisbelow": True,  # Put grid lines behind plot elements
        "axes.spines.top": True,
        "axes.spines.right": True,
        "axes.spines.left": True,
        "axes.spines.bottom": True,
        # Ticks
        "xtick.major.size": 4,
        "xtick.minor.size": 2,
        "xtick.major.width": 1.0,
        "xtick.minor.width": 0.5,
        "xtick.major.pad": 4,
        "xtick.minor.pad": 4,
        "xtick.direction": "out",
        "ytick.major.size": 4,
        "ytick.minor.size": 2,
        "ytick.major.width": 1.0,
        "ytick.minor.width": 0.5,
        "ytick.major.pad": 4,
        "ytick.minor.pad": 4,
        "ytick.direction": "out",
        # Legend properties
        "legend.frameon": False,
        "legend.framealpha": 0.8,
        "legend.edgecolor": "0.8",
        "legend.numpoints": 1,
        "legend.scatterpoints": 3,
        # Grid properties
        "grid.linewidth": 0.8,
        "grid.alpha": 0.5,
        # Line properties
        "lines.linewidth": 1.5,
        "lines.markersize": 6,
        "lines.markeredgewidth": 1.0,
        # Patches properties (for boxplots, etc)
        "patch.linewidth": 1.0,
        "patch.edgecolor": "black",
        # Scatter properties
        "scatter.marker": "o",
        # Violin plots
        "boxplot.showcaps": True,
        "boxplot.showbox": True,
        "boxplot.showfliers": True,
        "boxplot.showmeans": False,
    }

    # Update with the selected theme
    default_styles.update(theme_params)

    # Update with user-provided style dictionary if any
    if style_dict:
        default_styles.update(style_dict)

    # Apply the settings
    plt.rcParams.update(default_styles)

    # Try to set IPython display format safely
    try:
        from IPython import get_ipython

        ipython = get_ipython()
        if ipython is not None:
            # Try both methods to set matplotlib display formats
            try:
                # Newer IPython versions
                ipython.run_line_magic("matplotlib", "inline")
                for fmt in ["png", "svg"]:
                    if fmt in sc.settings.figure_formats:
                        ipython.run_line_magic(
                            "config", f"InlineBackend.figure_formats = ['{fmt}']"
                        )
                        break
            except:
                # Try another approach for older versions
                try:
                    import IPython.display

                    # Older IPython might have this method
                    if hasattr(IPython.display, "set_matplotlib_formats"):
                        IPython.display.set_matplotlib_formats(
                            *sc.settings.figure_formats
                        )
                    else:
                        print(
                            "Note: IPython display format could not be set automatically."
                        )
                except:
                    print(
                        "Note: IPython display format could not be set automatically."
                    )
    except:
        # Not running in IPython or error occurred
        pass

    print(f"Global plotting settings applied with '{color_theme}' color theme.")

    # Return the applied settings dict for reference
    # return default_styles


# --- Initialize with default settings on import ---
#setup_logging()
#set_figure_params()
