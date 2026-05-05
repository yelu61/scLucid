"""
Dash-based web application for scLucid.

This module provides an interactive web interface using Dash framework
with Plotly.js for visualizations.
"""

import plotly.graph_objects as go

# API client for backend communication
import requests
from dash import Dash, Input, Output, State, dcc, html


class QCDashApp:
    """
    Dash application for QC analysis.

    Features:
    - Interactive QC parameter adjustment
    - Real-time Plotly.js visualizations
    - Live filtering preview
    - QC metrics exploration
    """

    def __init__(
        self,
        api_base_url: str = "http://localhost:8000",
        port: int = 8050,
        debug: bool = False,
    ):
        """
        Initialize Dash application.

        Args:
            api_base_url: Base URL for the FastAPI backend
            port: Port to run Dash server on
            debug: Enable debug mode
        """
        self.api_base_url = api_base_url
        self.port = port
        self.debug = debug

        # Initialize Dash app
        self.app = Dash(
            __name__,
            external_stylesheets=[
                "https://codepen.io/chriddyp/pen/bWLwgP.css",
                "https://cdn.jsdelivr.net/npm/dash-bootstrap@1.5.3/dist/dash-bootstrap.min.css",
            ],
        )

        self.app.title = "scLucid QC Dashboard"

        # Setup layout and callbacks
        self._setup_layout()
        self._setup_callbacks()

    def _setup_layout(self):
        """Setup the dashboard layout."""
        self.app.layout = html.Div(
            [
                # Header
                html.Div(
                    [
                        html.H1("scLucid QC Dashboard", className="display-4"),
                        html.P("Interactive Quality Control Analysis", className="lead"),
                    ],
                    className="jumbotron bg-primary text-white p-4 mb-4",
                ),
                # Project selection
                html.Div(
                    [
                        html.Label("Project ID:"),
                        dcc.Input(
                            id="project-id-input",
                            type="text",
                            value="demo_project",
                            className="form-control",
                        ),
                        html.Button(
                            "Load Project",
                            id="load-project-btn",
                            className="btn btn-primary mt-2",
                        ),
                    ],
                    className="card card-body mb-4",
                ),
                # QC Parameters Control Panel
                html.Div(
                    [
                        html.H4("QC Parameters"),
                        html.Div(
                            [
                                html.Label("Min Genes:"),
                                dcc.Slider(
                                    id="min-genes-slider",
                                    min=50,
                                    max=2000,
                                    value=200,
                                    marks={i: str(i) for i in range(100, 2001, 500)},
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.Div(id="min-genes-value"),
                            ],
                            className="mb-3",
                        ),
                        html.Div(
                            [
                                html.Label("Max MT%:"),
                                dcc.Slider(
                                    id="max-mt-slider",
                                    min=0,
                                    max=50,
                                    value=20,
                                    step=1,
                                    marks={i: str(i) for i in range(0, 51, 10)},
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.Div(id="max-mt-value"),
                            ],
                            className="mb-3",
                        ),
                        html.Div(
                            [
                                html.Label("Min Counts:"),
                                dcc.Slider(
                                    id="min-counts-slider",
                                    min=100,
                                    max=20000,
                                    value=1000,
                                    marks={i: str(i) for i in range(1000, 20001, 5000)},
                                    tooltip={"placement": "bottom", "always_visible": True},
                                ),
                                html.Div(id="min-counts-value"),
                            ],
                            className="mb-3",
                        ),
                        html.Button(
                            "Apply Filters",
                            id="apply-filters-btn",
                            className="btn btn-success",
                        ),
                    ],
                    className="card card-body mb-4",
                ),
                # Summary Statistics
                html.Div(
                    [
                        html.H4("Summary Statistics"),
                        html.Div(id="summary-stats"),
                    ],
                    className="card card-body mb-4",
                ),
                # Plotly.js Visualizations Row
                html.Div(
                    [
                        # Scatter Plot
                        html.Div(
                            [
                                html.H5("Genes vs MT% (Plotly.js)"),
                                dcc.Graph(id="scatter-plot"),
                            ],
                            className="card card-body",
                        ),
                        # Violin Plot
                        html.Div(
                            [
                                html.H5("QC Metrics Distribution"),
                                html.Label("Select Metric:"),
                                dcc.Dropdown(id="violin-metric-dropdown"),
                                dcc.Graph(id="violin-plot"),
                            ],
                            className="card card-body",
                        ),
                    ],
                    className="row",
                ),
                # Histogram
                html.Div(
                    [
                        html.H5("Histogram"),
                        html.Label("Select Metric:"),
                        dcc.Dropdown(id="histogram-metric-dropdown"),
                        dcc.Graph(id="histogram-plot"),
                    ],
                    className="card card-body mt-4",
                ),
                # Footer
                html.Div(
                    [
                        html.P(
                            "Generated by scLucid - Interactive QC Analysis", className="text-muted"
                        ),
                    ],
                    className="text-center mt-4",
                ),
                # Store components
                dcc.Store(id="project-data-store"),
                dcc.Store(id="filtered-data-store"),
            ]
        )

    def _setup_callbacks(self):
        """Setup interactive callbacks."""

        @self.app.callback(
            [
                Output("min-genes-value", "children"),
                Output("max-mt-value", "children"),
                Output("min-counts-value", "children"),
            ],
            [
                Input("min-genes-slider", "value"),
                Input("max-mt-slider", "value"),
                Input("min-counts-slider", "value"),
            ],
        )
        def update_slider_values(min_genes, max_mt, min_counts):
            """Update slider value displays."""
            return (f"{min_genes} genes", f"{max_mt}%", f"{min_counts} counts")

        @self.app.callback(
            Output("project-data-store", "data"),
            [Input("load-project-btn", "n_clicks")],
            [State("project-id-input", "value")],
        )
        def load_project(n_clicks, project_id):
            """Load project data from backend."""
            if n_clicks is None or not project_id:
                return None

            try:
                response = requests.post(
                    f"{self.api_base_url}/api/qc/metrics",
                    params={"project_id": project_id},
                    timeout=10,
                )
                response.raise_for_status()
                return response.json()
            except Exception as e:
                print(f"Error loading project: {e}")
                return None

        @self.app.callback(
            [Output("scatter-plot", "figure"), Output("summary-stats", "children")],
            [Input("apply-filters-btn", "n_clicks"), Input("project-data-store", "data")],
            [
                State("min-genes-slider", "value"),
                State("max-mt-slider", "value"),
                State("min-counts-slider", "value"),
                State("project-id-input", "value"),
            ],
        )
        def update_filter_preview(
            n_clicks, project_data, min_genes, max_mt, min_counts, project_id
        ):
            """Update filter preview with Plotly.js scatter plot."""
            if not project_id:
                return go.Figure(), "No project loaded"

            try:
                response = requests.post(
                    f"{self.api_base_url}/api/qc/filter-preview",
                    params={
                        "project_id": project_id,
                        "min_genes": min_genes,
                        "max_mt_percent": max_mt,
                        "min_counts": min_counts,
                    },
                    timeout=10,
                )
                response.raise_for_status()
                data = response.json()

                # Create Plotly.js scatter figure
                scatter_data = data["scatter_data"]
                summary = data["summary"]

                fig = go.Figure()

                # Add retained cells
                retained_indices = [i for i, r in enumerate(scatter_data["retained"]) if r]
                removed_indices = [i for i, r in enumerate(scatter_data["retained"]) if not r]

                if retained_indices:
                    fig.add_trace(
                        go.Scatter(
                            x=[scatter_data["x"][i] for i in retained_indices],
                            y=[scatter_data["y"][i] for i in retained_indices],
                            mode="markers",
                            name="Retained",
                            marker=dict(size=5, color="green", opacity=0.6),
                        )
                    )

                if removed_indices:
                    fig.add_trace(
                        go.Scatter(
                            x=[scatter_data["x"][i] for i in removed_indices],
                            y=[scatter_data["y"][i] for i in removed_indices],
                            mode="markers",
                            name="Removed",
                            marker=dict(size=5, color="red", opacity=0.6),
                        )
                    )

                fig.update_layout(
                    title="QC Filter Preview",
                    xaxis_title="Number of Genes",
                    yaxis_title="MT%",
                    hovermode="closest",
                    height=400,
                )

                # Summary stats
                stats_html = html.Div(
                    [
                        html.P(f"Total Cells: {summary['total']:,}"),
                        html.P(
                            f"Retained: {summary['retained']:,} ({summary['retained_pct']:.1f}%)"
                        ),
                        html.P(f"Removed: {summary['removed']:,}", className="text-danger"),
                    ]
                )

                return fig, stats_html

            except Exception as e:
                return go.Figure(), f"Error: {str(e)}"

        @self.app.callback(
            [
                Output("violin-metric-dropdown", "options"),
                Output("violin-metric-dropdown", "value"),
            ],
            [Input("project-data-store", "data")],
        )
        def update_metric_options(project_data):
            """Update available metric options."""
            if not project_data:
                return [], None

            metrics = [
                "n_genes_by_counts",
                "total_counts",
                "pct_counts_mt",
                "log1p_total_counts",
            ]

            options = [{"label": m, "value": m} for m in metrics]
            return options, metrics[0]

        @self.app.callback(
            Output("violin-plot", "figure"),
            [Input("violin-metric-dropdown", "value"), Input("project-id-input", "value")],
        )
        def update_violin_plot(metric, project_id):
            """Update Plotly.js violin plot."""
            if not metric or not project_id:
                return go.Figure()

            try:
                response = requests.get(
                    f"{self.api_base_url}/api/qc/violin-data/{project_id}",
                    params={"metric": metric},
                    timeout=10,
                )
                response.raise_for_status()
                data = response.json()

                fig = go.Figure()
                fig.add_trace(
                    go.Violin(
                        y=data["data"],
                        name=metric,
                        box_visible=True,
                        meanline_visible=True,
                    )
                )

                fig.update_layout(
                    title=f"{metric} Distribution",
                    yaxis_title=metric,
                    height=400,
                )

                return fig

            except Exception:
                return go.Figure()

        @self.app.callback(
            Output("histogram-plot", "figure"),
            [Input("histogram-metric-dropdown", "value"), Input("project-id-input", "value")],
        )
        def update_histogram(metric, project_id):
            """Update Plotly.js histogram."""
            if not metric or not project_id:
                return go.Figure()

            try:
                response = requests.get(
                    f"{self.api_base_url}/api/qc/histogram-data/{project_id}",
                    params={"metric": metric},
                    timeout=10,
                )
                response.raise_for_status()
                data = response.json()

                fig = go.Figure()
                fig.add_trace(go.Bar(x=data["x"], y=data["y"], name=metric))

                fig.update_layout(
                    title=f"{metric} Histogram",
                    xaxis_title=metric,
                    yaxis_title="Count",
                    height=400,
                )

                return fig

            except Exception:
                return go.Figure()

    def run(self, host: str = "127.0.0.1"):
        """
        Run the Dash server.

        Args:
            host: Host to bind to
        """
        print(f"Starting scLucid QC Dashboard at http://{host}:{self.port}")
        print(f"Backend API at {self.api_base_url}")
        self.app.run_server(host=host, port=self.port, debug=self.debug)


def launch_dashboard(
    api_base_url: str = "http://localhost:8000",
    dash_port: int = 8050,
    host: str = "127.0.0.1",
    debug: bool = False,
):
    """
    Convenience function to launch the QC dashboard.

    Args:
        api_base_url: Base URL for FastAPI backend
        dash_port: Port for Dash frontend
        host: Host to bind to
        debug: Debug mode
    """
    app = QCDashApp(api_base_url=api_base_url, port=dash_port, debug=debug)
    app.run(host=host)
