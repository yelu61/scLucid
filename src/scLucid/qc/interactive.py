"""
Interactive Jupyter widgets for QC exploration.

This module provides interactive components for exploring QC parameters
and visualizing their effects in real-time within Jupyter notebooks.
"""

import logging
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from anndata import AnnData

log = logging.getLogger(__name__)


class InteractiveQCExplorer:
    """
    Interactive QC parameter explorer with live preview.

    Provides sliders and dropdowns to adjust QC parameters and see
    the effects on filtering in real-time.
    """

    def __init__(
        self,
        adata: AnnData,
        sample_key: Optional[str] = None,
    ):
        """
        Initialize the interactive QC explorer.

        Args:
            adata: AnnData object with QC metrics
            sample_key: Key in obs identifying samples
        """
        self.adata = adata
        self.sample_key = sample_key
        self._widgets = {}
        self._output_handler = None

    def create_threshold_sliders(
        self,
        initial_thresholds: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, 'ipywidgets.Widget']:
        """
        Create interactive sliders for QC thresholds.

        Args:
            initial_thresholds: Initial threshold values

        Returns:
            Dictionary of widget objects
        """
        try:
            import ipywidgets as widgets
            from IPython.display import display
        except ImportError:
            raise ImportError(
                "ipywidgets is required for interactive exploration. "
                "Install with: pip install ipywidgets"
            )

        if initial_thresholds is None:
            initial_thresholds = {
                'min_genes': 200,
                'max_pct_mt': 20,
                'min_counts': 1000,
                'max_pct_top_genes': 80,
            }

        # Gene count slider
        min_genes_slider = widgets.IntSlider(
            value=initial_thresholds.get('min_genes', 200),
            min=50,
            max=2000,
            step=50,
            description='Min Genes:',
            continuous_update=False,
        )

        # Mitochondrial percentage slider
        max_mt_slider = widgets.FloatSlider(
            value=initial_thresholds.get('max_pct_mt', 20),
            min=0,
            max=50,
            step=1,
            description='Max MT%:',
            continuous_update=False,
        )

        # Total count slider
        min_counts_slider = widgets.IntSlider(
            value=initial_thresholds.get('min_counts', 1000),
            min=100,
            max=20000,
            step=500,
            description='Min Counts:',
            continuous_update=False,
        )

        # Top genes percentage slider
        max_top_genes_slider = widgets.FloatSlider(
            value=initial_thresholds.get('max_pct_top_genes', 80),
            min=50,
            max=100,
            step=5,
            description='Max Top Genes%:',
            continuous_update=False,
        )

        self._widgets = {
            'min_genes': min_genes_slider,
            'max_pct_mt': max_mt_slider,
            'min_counts': min_counts_slider,
            'max_pct_top_genes': max_top_genes_slider,
        }

        return self._widgets

    def create_sample_selector(self) -> 'ipywidgets.Widget':
        """
        Create sample selector dropdown.

        Returns:
            Dropdown widget
        """
        try:
            import ipywidgets as widgets
        except ImportError:
            raise ImportError("ipywidgets is required")

        if self.sample_key and self.sample_key in self.adata.obs:
            samples = self.adata.obs[self.sample_key].unique().tolist()
        else:
            samples = ['All']

        sample_dropdown = widgets.Dropdown(
            options=samples,
            value=samples[0] if samples else None,
            description='Sample:',
        )

        self._widgets['sample'] = sample_dropdown

        return sample_dropdown

    def create_output_widget(self) -> 'ipywidgets.Widget':
        """
        Create output widget for displaying results.

        Returns:
            Output widget
        """
        try:
            import ipywidgets as widgets
        except ImportError:
            raise ImportError("ipywidgets is required")

        self._output_handler = widgets.Output()

        return self._output_handler

    def on_threshold_change(
        self,
        change: Dict[str, Any],
        preview_func: Optional[Callable[[Dict[str, Any]], None]] = None,
    ):
        """
        Handle threshold change events.

        Args:
            change: Change event data
            preview_func: Function to call with new thresholds
        """
        if self._output_handler is None:
            return

        thresholds = self.get_current_thresholds()

        with self._output_handler:
            self._output_handler.clear_output(wait=True)

            # Display current thresholds
            print("Current Thresholds:")
            for key, value in thresholds.items():
                print(f"  {key}: {value}")

            # Apply filtering and show preview
            if preview_func:
                preview_func(thresholds)
            else:
                self._default_preview(thresholds)

    def _default_preview(self, thresholds: Dict[str, Any]):
        """
        Default preview showing filtering statistics.

        Args:
            thresholds: Current threshold values
        """
        # Apply filters
        mask = np.ones(self.adata.n_obs, dtype=bool)

        if 'min_genes' in thresholds and 'n_genes_by_counts' in self.adata.obs:
            mask &= self.adata.obs['n_genes_by_counts'] >= thresholds['min_genes']

        if 'max_pct_mt' in thresholds and 'pct_counts_mt' in self.adata.obs:
            mask &= self.adata.obs['pct_counts_mt'] <= thresholds['max_pct_mt']

        if 'min_counts' in thresholds and 'total_counts' in self.adata.obs:
            mask &= self.adata.obs['total_counts'] >= thresholds['min_counts']

        if 'max_pct_top_genes' in thresholds and 'pct_counts_in_top_20_genes' in self.adata.obs:
            mask &= self.adata.obs['pct_counts_in_top_20_genes'] <= thresholds['max_pct_top_genes']

        n_retained = mask.sum()
        n_filtered = self.adata.n_obs - n_retained

        print(f"\nFiltering Results:")
        print(f"  Total cells: {self.adata.n_obs}")
        print(f"  Retained: {n_retained} ({n_retained/self.adata.n_obs:.1%})")
        print(f"  Filtered: {n_filtered} ({n_filtered/self.adata.n_obs:.1%})")

    def get_current_thresholds(self) -> Dict[str, Any]:
        """
        Get current threshold values from widgets.

        Returns:
            Dictionary of threshold values
        """
        thresholds = {}
        for key, widget in self._widgets.items():
            if hasattr(widget, 'value'):
                thresholds[key] = widget.value

        return thresholds

    def display(
        self,
        preview_func: Optional[Callable[[Dict[str, Any]], None]] = None,
    ):
        """
        Display the interactive explorer.

        Args:
            preview_func: Custom function for previewing results
        """
        try:
            import ipywidgets as widgets
            from IPython.display import display
        except ImportError:
            raise ImportError("ipywidgets is required")

        # Create widgets if not already created
        if not self._widgets:
            self.create_threshold_sliders()
            self.create_output_widget()

        # Create layout
        sliders = [
            self._widgets['min_genes'],
            self._widgets['max_pct_mt'],
            self._widgets['min_counts'],
            self._widgets['max_pct_top_genes'],
        ]

        # Attach observers
        for widget in sliders:
            widget.observe(
                lambda change: self.on_threshold_change(change, preview_func),
                names='value',
            )

        # Display
        slider_box = widgets.VBox(sliders)
        display(widgets.VBox([slider_box, self._output_handler]))

        # Initial preview
        self.on_threshold_change(None, preview_func)


class InteractiveQCPlotter:
    """
    Interactive plotting for QC metrics with parameter adjustments.
    """

    def __init__(self, adata: AnnData):
        """
        Initialize interactive plotter.

        Args:
            adata: AnnData object with QC metrics
        """
        self.adata = adata
        self._widgets = {}
        self._current_plot_type = 'violin'

    def create_plot_controls(self) -> 'ipywidgets.Widget':
        """
        Create controls for plot customization.

        Returns:
            VBox widget with controls
        """
        try:
            import ipywidgets as widgets
        except ImportError:
            raise ImportError("ipywidgets is required")

        # Plot type selector
        plot_type_dropdown = widgets.Dropdown(
            options=['violin', 'scatter', 'histogram', 'heatmap'],
            value='violin',
            description='Plot Type:',
        )

        # Metric selector
        available_metrics = [
            m for m in self.adata.obs.columns
            if m.startswith(('log1p_', 'pct_', 'n_', 'total_')) or m in ['phase', 'S_score', 'G2M_score']
        ]

        metric_dropdown = widgets.Dropdown(
            options=available_metrics,
            value=available_metrics[0] if available_metrics else None,
            description='Metric:',
        )

        # Figure size sliders
        width_slider = widgets.IntSlider(
            value=10,
            min=5,
            max=20,
            step=1,
            description='Width:',
        )

        height_slider = widgets.IntSlider(
            value=6,
            min=4,
            max=15,
            step=1,
            description='Height:',
        )

        # Color selector
        color_dropdown = widgets.Dropdown(
            options=['None'] + [m for m in self.adata.obs.columns if m != 'metric'],
            value='None',
            description='Color by:',
        )

        # Update button
        update_button = widgets.Button(
            description='Update Plot',
            button_style='primary',
        )

        self._widgets = {
            'plot_type': plot_type_dropdown,
            'metric': metric_dropdown,
            'width': width_slider,
            'height': height_slider,
            'color': color_dropdown,
            'update': update_button,
        }

        # Attach update handler
        update_button.on_click(self._on_update_click)

        return widgets.VBox(list(self._widgets.values()))

    def _on_update_click(self, button):
        """Handle update button click."""
        self.update_plot()

    def update_plot(self):
        """Update the plot based on current widget values."""
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("matplotlib is required")

        plot_type = self._widgets['plot_type'].value
        metric = self._widgets['metric'].value
        width = self._widgets['width'].value
        height = self._widgets['height'].value
        color_by = self._widgets['color'].value

        if color_by == 'None':
            color_by = None

        # Create figure
        fig, ax = plt.subplots(figsize=(width, height))

        if plot_type == 'violin':
            self._plot_violin(ax, metric, color_by)
        elif plot_type == 'scatter':
            self._plot_scatter(ax, metric, color_by)
        elif plot_type == 'histogram':
            self._plot_histogram(ax, metric, color_by)
        elif plot_type == 'heatmap':
            self._plot_heatmap(ax, metric)

        plt.tight_layout()
        plt.show()

    def _plot_violin(self, ax, metric, color_by):
        """Plot violin plot."""
        try:
            import seaborn as sns
        except ImportError:
            raise ImportError("seaborn is required")

        if color_by and color_by in self.adata.obs:
            sns.violinplot(
                data=self.adata.obs,
                x=color_by,
                y=metric,
                ax=ax,
            )
        else:
            sns.violinplot(
                y=self.adata.obs[metric],
                ax=ax,
            )

        ax.set_title(f'{metric} Distribution')
        ax.set_xlabel(color_by if color_by else '')
        ax.set_ylabel(metric)

    def _plot_scatter(self, ax, metric, color_by):
        """Plot scatter plot."""
        x_data = range(len(self.adata.obs))

        if color_by and color_by in self.adata.obs:
            scatter = ax.scatter(
                x_data,
                self.adata.obs[metric],
                c=self.adata.obs[color_by],
                cmap='viridis',
                alpha=0.6,
            )
            plt.colorbar(scatter, ax=ax, label=color_by)
        else:
            ax.scatter(
                x_data,
                self.adata.obs[metric],
                alpha=0.6,
            )

        ax.set_title(f'{metric} Scatter Plot')
        ax.set_xlabel('Cell Index')
        ax.set_ylabel(metric)

    def _plot_histogram(self, ax, metric, color_by):
        """Plot histogram."""
        if color_by and color_by in self.adata.obs:
            # Group by color variable
            for group in self.adata.obs[color_by].unique():
                group_data = self.adata.obs[self.adata.obs[color_by] == group][metric]
                ax.hist(group_data, alpha=0.5, label=group, bins=50)
            ax.legend()
        else:
            ax.hist(self.adata.obs[metric], bins=50, edgecolor='black')

        ax.set_title(f'{metric} Histogram')
        ax.set_xlabel(metric)
        ax.set_ylabel('Frequency')

    def _plot_heatmap(self, ax, metric):
        """Plot metric heatmap (for multiple samples)."""
        # This is a placeholder - actual implementation depends on data structure
        ax.text(0.5, 0.5, 'Heatmap view requires sample grouping',
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'{metric} Heatmap')

    def display(self):
        """Display the interactive plotter."""
        try:
            from IPython.display import display
        except ImportError:
            raise ImportError("IPython is required")

        controls = self.create_plot_controls()
        display(controls)

        # Initial plot
        self.update_plot()


def create_interactive_dashboard(
    adata: AnnData,
    sample_key: Optional[str] = None,
):
    """
    Create a comprehensive interactive dashboard for QC exploration.

    Args:
        adata: AnnData object with QC metrics
        sample_key: Key in obs identifying samples
    """
    try:
        from IPython.display import display
        import ipywidgets as widgets
    except ImportError:
        raise ImportError("ipywidgets and IPython are required")

    # Create tabbed interface
    tab = widgets.Tab()

    # Tab 1: Threshold explorer
    explorer = InteractiveQCExplorer(adata, sample_key)
    tab.children = [widgets.VBox([])]

    tab.set_title(0, 'Threshold Explorer')

    # Tab 2: Interactive plots
    plotter = InteractiveQCPlotter(adata)
    tab.children = [widgets.VBox([]), widgets.VBox([])]
    tab.set_title(1, 'Interactive Plots')

    # Add content to tabs
    with tab.children[0]:
        explorer.display()

    with tab.children[1]:
        plotter.display()

    display(tab)


def interactive_filter_preview(
    adata: AnnData,
    thresholds: Dict[str, Any],
    sample_key: Optional[str] = None,
) -> AnnData:
    """
    Preview filtering with given thresholds.

    Args:
        adata: AnnData object
        thresholds: Threshold dictionary
        sample_key: Sample key for grouping

    Returns:
        Filtered AnnData (preview only, not saved)
    """
    mask = np.ones(adata.n_obs, dtype=bool)

    # Apply each threshold
    if 'min_genes' in thresholds and 'n_genes_by_counts' in adata.obs:
        mask &= adata.obs['n_genes_by_counts'] >= thresholds['min_genes']

    if 'max_genes' in thresholds and 'n_genes_by_counts' in adata.obs:
        mask &= adata.obs['n_genes_by_counts'] <= thresholds['max_genes']

    if 'min_counts' in thresholds and 'total_counts' in adata.obs:
        mask &= adata.obs['total_counts'] >= thresholds['min_counts']

    if 'max_counts' in thresholds and 'total_counts' in adata.obs:
        mask &= adata.obs['total_counts'] <= thresholds['max_counts']

    if 'max_pct_mt' in thresholds and 'pct_counts_mt' in adata.obs:
        mask &= adata.obs['pct_counts_mt'] <= thresholds['max_pct_mt']

    if 'max_pct_hb' in thresholds and 'pct_counts_hb' in adata.obs:
        mask &= adata.obs['pct_counts_hb'] <= thresholds['max_pct_hb']

    if 'max_pct_top_genes' in thresholds and 'pct_counts_in_top_20_genes' in adata.obs:
        mask &= adata.obs['pct_counts_in_top_20_genes'] <= thresholds['max_pct_top_genes']

    return adata[mask].copy()
