"""Publication-ready plotting themes."""

import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List

FONT_FAMILY = ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif']

# Nature style
NATURE_THEME = {
    'font.family': 'sans-serif',
    'font.sans-serif': FONT_FAMILY,
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.linewidth': 1.0,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'legend.frameon': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
    'pdf.fonttype': 42, # 保证导出PDF文字可编辑
    'ps.fonttype': 42,
}

NATURE_COLORS = {
    'palette': ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#DC0000', '#7E6148'],
    'cmap': 'RdYlBu_r'
}

# Science style
SCIENCE_THEME = {
    'font.family': 'sans-serif',
    'font.sans-serif': FONT_FAMILY,
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.linewidth': 1.0,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'figure.dpi': 300,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
}

SCIENCE_COLORS = {
    'palette': ['#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#DC0000', '#00A087', '#E64B35', '#4DBBD5'],
    'cmap': 'viridis'
}

def apply_theme(theme: str = 'nature') -> Dict:
    """
    Apply a publication-ready theme.
    
    Parameters:
        theme: 'nature' or 'science'
    
    Returns:
        Dict containing color configurations
    """
    if theme.lower() == 'nature':
        plt.rcParams.update(NATURE_THEME)
        sns.set_palette(NATURE_COLORS['palette'])
        return NATURE_COLORS
    elif theme.lower() == 'science':
        plt.rcParams.update(SCIENCE_THEME)
        sns.set_palette(SCIENCE_COLORS['palette'])
        return SCIENCE_COLORS
    else:
        # Fallback to defaults but try to optimize DPI
        plt.rcParams['figure.dpi'] = 300
        return {'palette': 'tab10', 'cmap': 'viridis'}