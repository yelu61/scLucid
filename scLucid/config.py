"""Configuration module for single-cell RNA-seq data analysis"""

import logging
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Union
from pathlib import Path
import scanpy as sc
import matplotlib.pyplot as plt

log = logging.getLogger(__name__)


@dataclass
class QCConfig:
    """Quality Control configuration"""
    # Basic filtering thresholds
    min_genes: int = 200
    max_genes: Optional[int] = 8000
    min_counts: Optional[int] = None
    max_counts: Optional[int] = None
    min_cells: int = 3
    
    # Percentage-based thresholds
    pc_mt: float = 20.0
    pc_hb: float = 20.0
    pc_top20_genes: Optional[float] = 65.0
    
    # Statistical outlier detection
    nmads: float = 5.0
    use_fixed_top20_threshold: bool = False
    
    # Gene pattern definitions
    gene_patterns: Dict[str, str] = field(default_factory=lambda: {
        "mt": r"^(MT|Mt|mt)-",
        "ribo": r"^(RP[SL]|Rp[sl])",
        "hb": r"^(HB|hb)[^(P|p)]",
    })
    
    # QC metrics to calculate and plot
    qc_metrics: List[str] = field(default_factory=lambda: [
        "total_counts", "n_genes_by_counts", "pct_counts_mt", 
        "pct_counts_hb", "pct_counts_ribo", "pct_counts_in_top_20_genes"
    ])
    
    # Plotting options
    plot_violin: bool = True
    plot_scatter: bool = True
    plot_top20: bool = True
    plot_outliers: bool = True
    
    def validate(self):
        """Validate QC configuration parameters"""
        if self.min_genes < 0:
            raise ValueError("min_genes must be non-negative")
        if self.max_genes is not None and self.max_genes <= self.min_genes:
            raise ValueError("max_genes must be greater than min_genes")
        if not 0 <= self.pc_mt <= 100:
            raise ValueError("pc_mt must be between 0 and 100")
        if not 0 <= self.pc_hb <= 100:
            raise ValueError("pc_hb must be between 0 and 100")
        if self.pc_top20_genes is not None and not 0 <= self.pc_top20_genes <= 100:
            raise ValueError("pc_top20_genes must be between 0 and 100")


@dataclass
class DoubletConfig:
    """Doublet detection configuration"""
    # Methods and parameters
    method: str = "scrublet"
    rate_per_1000_cells: float = 0.008
    expected_doublet_rate: Union[float, Dict[str, float]] = 0.1
    
    # Scrublet-specific parameters
    n_pcs: int = 30
    use_heuristics: bool = True
    
    # Marker-based doublet detection
    marker_dict: Optional[Dict[str, List[str]]] = None
    use_raw_for_markers: bool = True
    
    # Plotting
    plot_umap: bool = True
    
    def validate(self):
        """Validate doublet configuration parameters"""
        if self.method not in ["scrublet"]:  # 将来可以添加更多方法
            raise ValueError(f"Unsupported doublet detection method: {self.method}")
        if not 0 <= self.rate_per_1000_cells <= 1:
            raise ValueError("rate_per_1000_cells must be between 0 and 1")


@dataclass
class PreprocessConfig:
    """Preprocessing configuration"""
    # Normalization
    target_sum: float = 1e4
    layer_raw: str = "counts"
    layer_norm: str = "log1p_norm"
    layer_scale: str = "scaled"
    
    # Highly variable genes
    n_top_genes: int = 2000
    flavor: str = "seurat_v3"
    
    # Scaling
    max_value: Optional[float] = 10
    zero_center: bool = True
    
    def validate(self):
        """Validate preprocessing configuration"""
        if self.target_sum <= 0:
            raise ValueError("target_sum must be positive")
        if self.n_top_genes <= 0:
            raise ValueError("n_top_genes must be positive")


@dataclass
class AnalysisConfig:
    """Analysis configuration"""
    # PCA
    n_pcs: int = 50
    use_highly_variable: bool = True
    
    # Neighborhood graph
    n_neighbors: int = 30
    metric: str = "euclidean"
    
    # Clustering
    resolution: float = 0.8
    clustering_method: str = "leiden"
    
    # UMAP
    min_dist: float = 0.5
    spread: float = 1.0
    n_components: int = 2
    
    def validate(self):
        """Validate analysis configuration"""
        if self.n_pcs <= 0:
            raise ValueError("n_pcs must be positive")
        if self.n_neighbors <= 0:
            raise ValueError("n_neighbors must be positive")
        if not 0 <= self.resolution <= 5:
            raise ValueError("resolution should typically be between 0 and 5")


@dataclass
class PathConfig:
    """Path and I/O configuration"""
    # Base directories
    data_dir: Optional[Path] = None
    output_dir: Optional[Path] = None
    figure_dir: Optional[Path] = None
    cache_dir: Optional[Path] = None
    
    # File formats
    figure_formats: List[str] = field(default_factory=lambda: ['png', 'pdf'])
    adata_format: str = 'h5ad'
    
    def __post_init__(self):
        """Convert strings to Path objects and create directories"""
        if self.data_dir:
            self.data_dir = Path(self.data_dir)
        if self.output_dir:
            self.output_dir = Path(self.output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)
        if self.figure_dir:
            self.figure_dir = Path(self.figure_dir)
            self.figure_dir.mkdir(parents=True, exist_ok=True)
        if self.cache_dir:
            self.cache_dir = Path(self.cache_dir)
            self.cache_dir.mkdir(parents=True, exist_ok=True)


@dataclass
class LoggingConfig:
    """Logging configuration"""
    level: str = "INFO"
    format: str = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    file_path: Optional[Path] = None
    
    def setup_logging(self):
        """Setup logging configuration"""
        # Convert string level to logging constant
        numeric_level = getattr(logging, self.level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f'Invalid log level: {self.level}')
        
        # Configure logging
        handlers = [logging.StreamHandler()]
        if self.file_path:
            handlers.append(logging.FileHandler(self.file_path))
        
        logging.basicConfig(
            level=numeric_level,
            format=self.format,
            handlers=handlers,
            force=True  # Override existing configuration
        )


class Config:
    """
    Main configuration class that combines all sub-configurations.
    """
    def __init__(self):
        # General parameters
        self.batch_key: str = "sampleID"
        self.random_seed: int = 42
        
        # Sub-configurations
        self.qc = QCConfig()
        self.doublet = DoubletConfig()
        self.preprocess = PreprocessConfig()
        self.analysis = AnalysisConfig()
        self.paths = PathConfig()
        self.logging = LoggingConfig()
        
        # Plotting configuration
        self.figure = PlottingConfig()
        
        # Initialize logging
        self.logging.setup_logging()

    def update(self, user_params: Dict[str, Any]):
        """
        Update configuration with user-provided parameters.
        Supports nested updates for sub-configurations.
        """
        for key, value in user_params.items():
            if hasattr(self, key):
                if isinstance(getattr(self, key), (QCConfig, DoubletConfig, 
                                                 PreprocessConfig, AnalysisConfig, 
                                                 PathConfig, LoggingConfig)):
                    # Handle nested configuration updates
                    if isinstance(value, dict):
                        for subkey, subvalue in value.items():
                            if hasattr(getattr(self, key), subkey):
                                setattr(getattr(self, key), subkey, subvalue)
                                log.info(f"Configuration updated: {key}.{subkey} = {subvalue}")
                            else:
                                log.warning(f"Unknown configuration parameter: {key}.{subkey}")
                    else:
                        setattr(self, key, value)
                        log.info(f"Configuration updated: {key} = {value}")
                else:
                    setattr(self, key, value)
                    log.info(f"Configuration updated: {key} = {value}")
            else:
                log.warning(f"Unknown configuration parameter: {key}")

    def validate_all(self):
        """Validate all configuration sections"""
        try:
            self.qc.validate()
            self.doublet.validate()
            self.preprocess.validate()
            self.analysis.validate()
            log.info("All configuration parameters validated successfully")
        except ValueError as e:
            log.error(f"Configuration validation failed: {e}")
            raise

    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary for serialization"""
        from dataclasses import asdict
        return {
            'batch_key': self.batch_key,
            'random_seed': self.random_seed,
            'qc': asdict(self.qc),
            'doublet': asdict(self.doublet),
            'preprocess': asdict(self.preprocess),
            'analysis': asdict(self.analysis),
            'paths': asdict(self.paths),
            'logging': asdict(self.logging),
        }

    def save(self, filepath: Union[str, Path]):
        """Save configuration to file"""
        import json
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2, default=str)
        log.info(f"Configuration saved to {filepath}")

    @classmethod
    def load(cls, filepath: Union[str, Path]):
        """Load configuration from file"""
        import json
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        config = cls()
        config.update(data)
        return config


@dataclass
class PlottingConfig:
    """Enhanced plotting configuration"""
    dpi: int = 300
    figsize: tuple = (6, 6)
    style: str = "default"
    color_theme: str = "default"
    
    # Color palettes for different purposes
    categorical_palette: str = "tab20"
    continuous_palette: str = "viridis"
    diverging_palette: str = "RdBu_r"
    
    # Font settings
    font_family: str = "sans-serif"
    font_size: int = 10
    title_size: int = 14
    
    # Scanpy specific
    scanpy_dpi: int = 300
    scanpy_figsize: tuple = (4, 4)
    scanpy_frameon: bool = False


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
    sc.settings.figure_formats = ['png', 'svg']
    
    # Apply base style if specified
    if style:
        try:
            plt.style.use(style)
        except Exception as e:
            print(f"Warning: Could not apply style '{style}'. Error: {str(e)}")
            print("Falling back to default style.")
            plt.style.use('default')
    
    # Define color themes
    themes = {
        "default": {
            'figure.facecolor': 'white',
            'axes.facecolor': 'white',
            'savefig.facecolor': 'white',
            'text.color': 'black',
            'axes.labelcolor': 'black',
            'axes.edgecolor': 'black',
            'xtick.color': 'black',
            'ytick.color': 'black',
            'grid.color': 'lightgray',
        },
        "dark": {
            'figure.facecolor': '#303030',
            'axes.facecolor': '#303030',
            'savefig.facecolor': '#303030',
            'text.color': 'white',
            'axes.labelcolor': 'white',
            'axes.edgecolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'grid.color': '#505050',
        },
        "seaborn": {
            'figure.facecolor': 'white',
            'axes.facecolor': '#F0F0F0',
            'savefig.facecolor': 'white',
            'text.color': '#283747',
            'axes.labelcolor': '#283747',
            'axes.edgecolor': '#283747',
            'xtick.color': '#283747',
            'ytick.color': '#283747',
            'grid.color': 'white',
        }
    }
    
    # Apply selected color theme
    theme_params = themes.get(color_theme, themes["default"])
    
    # Define comprehensive default settings
    default_styles = {
        # Figure properties
        'figure.figsize': figsize,
        'figure.dpi': dpi,
        'savefig.dpi': dpi,
        'savefig.bbox': 'tight',
        'savefig.transparent': False,
        
        # Fonts and text
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'axes.titlesize': 14,
        'axes.titleweight': 'bold',
        'axes.labelsize': 12,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'legend.title_fontsize': 12,
        
        # Axes properties
        'axes.linewidth': 1.0,
        'axes.grid': False,
        'axes.axisbelow': True,  # Put grid lines behind plot elements
        'axes.spines.top': True,
        'axes.spines.right': True,
        'axes.spines.left': True,
        'axes.spines.bottom': True,
        
        # Ticks
        'xtick.major.size': 4,
        'xtick.minor.size': 2,
        'xtick.major.width': 1.0,
        'xtick.minor.width': 0.5,
        'xtick.major.pad': 4,
        'xtick.minor.pad': 4,
        'xtick.direction': 'out',
        
        'ytick.major.size': 4,
        'ytick.minor.size': 2,
        'ytick.major.width': 1.0,
        'ytick.minor.width': 0.5,
        'ytick.major.pad': 4,
        'ytick.minor.pad': 4,
        'ytick.direction': 'out',
        
        # Legend properties
        'legend.frameon': False,
        'legend.framealpha': 0.8,
        'legend.edgecolor': '0.8',
        'legend.numpoints': 1,
        'legend.scatterpoints': 3,
        
        # Grid properties
        'grid.linewidth': 0.8,
        'grid.alpha': 0.5,
        
        # Line properties
        'lines.linewidth': 1.5,
        'lines.markersize': 6,
        'lines.markeredgewidth': 1.0,
        
        # Patches properties (for boxplots, etc)
        'patch.linewidth': 1.0,
        'patch.edgecolor': 'black',
        
        # Scatter properties
        'scatter.marker': 'o',
        
        # Violin plots
        'boxplot.showcaps': True,
        'boxplot.showbox': True,
        'boxplot.showfliers': True,
        'boxplot.showmeans': False
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
                ipython.run_line_magic('matplotlib', 'inline')
                for fmt in ['png', 'svg']:
                    if fmt in sc.settings.figure_formats:
                        ipython.run_line_magic('config', f"InlineBackend.figure_formats = ['{fmt}']")
                        break
            except:
                # Try another approach for older versions
                try:
                    import IPython.display
                    # Older IPython might have this method
                    if hasattr(IPython.display, 'set_matplotlib_formats'):
                        IPython.display.set_matplotlib_formats(*sc.settings.figure_formats)
                    else:
                        print("Note: IPython display format could not be set automatically.")
                except:
                    print("Note: IPython display format could not be set automatically.")
    except:
        # Not running in IPython or error occurred
        pass
    
    print(f"Global plotting settings applied with '{color_theme}' color theme.")
    
    # Return the applied settings dict for reference
    #return default_styles

# Create global configuration instance
settings = Config()

# Convenience function for quick setup
def setup_environment(
    output_dir: Optional[str] = None,
    random_seed: int = 42,
    log_level: str = "INFO",
    **kwargs
):
    """Quick environment setup with common parameters"""
    if output_dir:
        settings.paths.output_dir = Path(output_dir)
        settings.paths.figure_dir = Path(output_dir) / "figures"
    
    settings.random_seed = random_seed
    settings.logging.level = log_level
    settings.logging.setup_logging()
    
    # Update any additional parameters
    if kwargs:
        settings.update(kwargs)
    
    # Validate configuration
    settings.validate_all()
    
    print(f"Environment setup complete. Output directory: {settings.paths.output_dir}")