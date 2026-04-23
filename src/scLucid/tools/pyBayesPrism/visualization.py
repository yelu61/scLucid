"""
Visualization functions for BayesPrism (R-free)
"""

import logging
from typing import Optional, Tuple

import pandas as pd

log = logging.getLogger(__name__)


def plot_fraction(
    fraction_df: pd.DataFrame,
    figsize: Tuple[int, int] = (12, 6),
    cmap: str = "YlOrRd",
    title: str = "Cell Type Fractions",
    save_path: Optional[str] = None,
):
    """
    Plot cell type fraction heatmap

    Parameters
    ----------
    fraction_df : pd.DataFrame
        Cell type fractions (samples x cell_types)
    figsize : Tuple[int, int]
        Figure size
    cmap : str
        Colormap name
    title : str
        Plot title
    save_path : str, optional
        Path to save figure
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        raise ImportError("matplotlib and seaborn are required for visualization")

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        fraction_df.T,
        cmap=cmap,
        cbar_kws={"label": "Fraction"},
        xticklabels=True,
        yticklabels=True,
        ax=ax,
    )

    ax.set_title(title)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Cell Types")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_correlation(
    fraction_df: pd.DataFrame,
    figsize: Tuple[int, int] = (10, 8),
    title: str = "Cell Type Correlation",
    save_path: Optional[str] = None,
):
    """
    Plot correlation between cell types

    Parameters
    ----------
    fraction_df : pd.DataFrame
        Cell type fractions (samples x cell_types)
    figsize : Tuple[int, int]
        Figure size
    title : str
        Plot title
    save_path : str, optional
        Path to save figure
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        raise ImportError("matplotlib and seaborn are required for visualization")

    corr = fraction_df.corr()

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        corr,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        center=0,
        square=True,
        ax=ax,
    )

    ax.set_title(title)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_stacked_bar(
    fraction_df: pd.DataFrame,
    figsize: Tuple[int, int] = (14, 6),
    title: str = "Cell Type Composition",
    save_path: Optional[str] = None,
):
    """
    Plot stacked bar chart of cell type fractions

    Parameters
    ----------
    fraction_df : pd.DataFrame
        Cell type fractions (samples x cell_types)
    figsize : Tuple[int, int]
        Figure size
    title : str
        Plot title
    save_path : str, optional
        Path to save figure
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required for visualization")

    fig, ax = plt.subplots(figsize=figsize)

    fraction_df.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        colormap="tab20",
    )

    ax.set_title(title)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Fraction")
    ax.legend(title="Cell Types", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_gene_programs(
    W_df: pd.DataFrame,
    top_n: int = 20,
    figsize: Tuple[int, int] = (15, 10),
    save_path: Optional[str] = None,
):
    """
    Plot gene program heatmaps

    Parameters
    ----------
    W_df : pd.DataFrame
        Gene program matrix (genes x programs)
    top_n : int
        Number of top genes to show per program
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save figure
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required for visualization")

    n_programs = W_df.shape[1]
    n_cols = 2
    n_rows = (n_programs + 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten() if n_programs > 1 else [axes]

    for i, program in enumerate(W_df.columns):
        program_values = W_df[program]
        top_genes = program_values.nlargest(top_n)

        axes[i].barh(range(top_n), top_genes.values[::-1])
        axes[i].set_yticks(range(top_n))
        axes[i].set_yticklabels(top_genes.index[::-1], fontsize=8)
        axes[i].set_xlabel("Weight")
        axes[i].set_title(f"{program}")

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, axes


def plot_program_usage(
    H_df: pd.DataFrame,
    figsize: Tuple[int, int] = (12, 6),
    cmap: str = "viridis",
    title: str = "Malignant Gene Program Usage",
    save_path: Optional[str] = None,
):
    """
    Plot program usage heatmap

    Parameters
    ----------
    H_df : pd.DataFrame
        Program usage matrix (programs x samples)
    figsize : Tuple[int, int]
        Figure size
    cmap : str
        Colormap name
    title : str
        Plot title
    save_path : str, optional
        Path to save figure
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        raise ImportError("matplotlib and seaborn are required for visualization")

    fig, ax = plt.subplots(figsize=figsize)

    sns.heatmap(
        H_df,
        cmap=cmap,
        cbar_kws={"label": "Usage"},
        xticklabels=True,
        yticklabels=True,
        ax=ax,
    )

    ax.set_title(title)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Programs")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_cv(
    cv_df: pd.DataFrame,
    figsize: Tuple[int, int] = (10, 6),
    save_path: Optional[str] = None,
):
    """
    Plot coefficient of variation

    Parameters
    ----------
    cv_df : pd.DataFrame
        CV dataframe with 'cell_type' and 'CV' columns
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save figure
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        raise ImportError("matplotlib is required for visualization")

    fig, ax = plt.subplots(figsize=figsize)

    ax.bar(cv_df["cell_type"], cv_df["CV"])
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Coefficient of Variation")
    ax.set_title("Uncertainty Estimation (CV)")
    ax.tick_params(axis="x", rotation=45)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, ax


def plot_validation_scatter(
    predicted: pd.DataFrame,
    actual: pd.DataFrame,
    figsize: Tuple[int, int] = (12, 10),
    save_path: Optional[str] = None,
):
    """
    Plot predicted vs actual proportions for validation

    Parameters
    ----------
    predicted : pd.DataFrame
        Predicted fractions
    actual : pd.DataFrame
        Actual (ground truth) fractions
    figsize : Tuple[int, int]
        Figure size
    save_path : str, optional
        Path to save figure
    """
    try:
        import matplotlib.pyplot as plt
        from scipy.stats import pearsonr
    except ImportError:
        raise ImportError("matplotlib and scipy are required")

    cell_types = predicted.columns
    n_types = len(cell_types)
    n_cols = 3
    n_rows = (n_types + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()

    for i, cell_type in enumerate(cell_types):
        ax = axes[i]

        x = actual[cell_type].values
        y = predicted[cell_type].values

        # Pearson correlation
        r, p = pearsonr(x, y)

        ax.scatter(x, y, alpha=0.6)
        ax.plot([0, 1], [0, 1], "r--", lw=2)

        ax.set_xlabel("Actual Fraction")
        ax.set_ylabel("Predicted Fraction")
        ax.set_title(f"{cell_type}\nr={r:.3f}, p={p:.2e}")

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, axes
