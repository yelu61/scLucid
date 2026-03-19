"""
Enhanced HTML report generation for QC results.

This module provides comprehensive, publication-ready HTML reports
with embedded visualizations and recommendations.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from anndata import AnnData

log = logging.getLogger(__name__)


class EnhancedQCReport:
    """
    Generate comprehensive HTML reports for QC analysis.

    Features:
    - Publication-ready formatting
    - Embedded interactive visualizations
    - Automatic recommendations
    - Downloadable tables
    - Statistical summaries
    """

    def __init__(
        self,
        adata: AnnData,
        adata_before: Optional[AnnData] = None,
    ):
        """
        Initialize the report generator.

        Args:
            adata: AnnData after QC filtering
            adata_before: AnnData before QC (for comparison)
        """
        self.adata = adata
        self.adata_before = adata_before or adata

    def generate_html_report(
        self,
        output_path: str,
        title: str = "Quality Control Report",
        author: Optional[str] = None,
        include_plots: bool = True,
        include_recommendations: bool = True,
    ):
        """
        Generate a comprehensive HTML report.

        Args:
            output_path: Path to save the HTML report
            title: Report title
            author: Report author
            include_plots: Whether to include plots
            include_recommendations: Whether to include recommendations
        """
        # Gather report data
        report_data = self._gather_report_data(include_plots=include_plots)

        # Generate HTML
        html_content = self._generate_html(
            title=title,
            author=author,
            data=report_data,
            include_plots=include_plots,
            include_recommendations=include_recommendations,
        )

        # Write to file
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        log.info(f"QC report saved to {output_path}")

    def _gather_report_data(self, include_plots: bool = True) -> Dict[str, Any]:
        """Gather all data for the report."""
        data = {
            'metadata': self._get_metadata(),
            'summary': self._get_summary_statistics(),
            'metrics': self._get_metrics_summary(),
            'filtering': self._get_filtering_summary(),
            'plots': self._get_plot_data() if include_plots else {},
            'recommendations': self._get_recommendations(),
            'tables': self._get_table_data(),
        }
        return data

    def _get_metadata(self) -> Dict[str, Any]:
        """Get report metadata."""
        return {
            'date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'n_cells_after': self.adata.n_obs,
            'n_cells_before': self.adata_before.n_obs,
            'n_genes': self.adata.n_vars,
            'retention_rate': self.adata.n_obs / self.adata_before.n_obs,
        }

    def _get_summary_statistics(self) -> Dict[str, Any]:
        """Get summary statistics."""
        stats = {
            'total_cells': self.adata.n_obs,
            'total_genes': self.adata.n_vars,
        }

        # QC metric statistics
        qc_metrics = [
            'n_genes_by_counts',
            'total_counts',
            'log1p_total_counts',
            'pct_counts_mt',
        ]

        for metric in qc_metrics:
            if metric in self.adata.obs:
                values = self.adata.obs[metric].values
                stats[metric] = {
                    'mean': float(np.mean(values)),
                    'median': float(np.median(values)),
                    'std': float(np.std(values)),
                    'min': float(np.min(values)),
                    'max': float(np.max(values)),
                }

        return stats

    def _get_metrics_summary(self) -> List[Dict[str, Any]]:
        """Get per-metric summary."""
        summary = []

        for metric in self.adata.obs.columns:
            if not metric.startswith(('log1p_', 'pct_', 'n_', 'total_')):
                continue

            values = self.adata.obs[metric].values
            summary.append({
                'name': metric,
                'mean': float(np.mean(values)),
                'median': float(np.median(values)),
                'std': float(np.std(values)),
            })

        return summary

    def _get_filtering_summary(self) -> Dict[str, Any]:
        """Get filtering summary."""
        n_before = self.adata_before.n_obs
        n_after = self.adata.n_obs
        n_removed = n_before - n_after

        return {
            'n_before': n_before,
            'n_after': n_after,
            'n_removed': n_removed,
            'retention_rate': n_after / n_before,
            'removal_rate': n_removed / n_before,
        }

    def _get_recommendations(self) -> List[Dict[str, str]]:
        """Generate QC recommendations."""
        recommendations = []

        # Check filtering rate
        retention_rate = self.adata.n_obs / self.adata_before.n_obs

        if retention_rate > 0.95:
            recommendations.append({
                'category': 'Filtering',
                'severity': 'warning',
                'message': 'Very few cells were filtered. Consider if QC thresholds are appropriate.',
            })
        elif retention_rate < 0.5:
            recommendations.append({
                'category': 'Filtering',
                'severity': 'warning',
                'message': 'More than 50% of cells were removed. Review QC parameters.',
            })

        # Check mitochondrial percentage
        if 'pct_counts_mt' in self.adata.obs:
            mt_values = self.adata.obs['pct_counts_mt'].values
            mean_mt = np.mean(mt_values)

            if mean_mt > 20:
                recommendations.append({
                    'category': 'Mitochondrial',
                    'severity': 'info',
                    'message': f'High mean MT% ({mean_mt:.1f}%). Possible cell stress or dying cells.',
                })

        # Check gene counts
        if 'n_genes_by_counts' in self.adata.obs:
            gene_counts = self.adata.obs['n_genes_by_counts'].values

            if np.median(gene_counts) < 500:
                recommendations.append({
                    'category': 'Gene Detection',
                    'severity': 'warning',
                    'message': f'Low median gene count ({np.median(gene_counts):.0}). Possible low-quality data.',
                })

        return recommendations

    def _get_plot_data(self) -> Dict[str, Any]:
        """Get data for plots."""
        plots = {}

        # Prepare data for violin plots
        for metric in ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']:
            if metric in self.adata.obs:
                plots[metric] = {
                    'values': self.adata.obs[metric].tolist(),
                    'name': metric,
                }

        return plots

    def _get_table_data(self) -> Dict[str, Any]:
        """Get data for tables."""
        tables = {}

        # Summary statistics table
        summary_data = []
        for metric, stats in self._get_summary_statistics().items():
            if isinstance(stats, dict):
                summary_data.append({
                    'Metric': metric,
                    'Mean': f"{stats.get('mean', 0):.2f}",
                    'Median': f"{stats.get('median', 0):.2f}",
                    'Std': f"{stats.get('std', 0):.2f}",
                })

        tables['summary'] = summary_data

        # QC metrics table
        qc_metrics_data = []
        for item in self._get_metrics_summary():
            qc_metrics_data.append({
                'Metric': item['name'],
                'Mean': f"{item['mean']:.2f}",
                'Median': f"{item['median']:.2f}",
                'Std': f"{item['std']:.2f}",
            })

        tables['metrics'] = qc_metrics_data

        return tables

    def _generate_html(
        self,
        title: str,
        author: Optional[str],
        data: Dict[str, Any],
        include_plots: bool,
        include_recommendations: bool,
    ) -> str:
        """Generate complete HTML report."""
        metadata = data['metadata']
        summary = data['summary']
        filtering = data['filtering']
        recommendations = data.get('recommendations', [])
        tables = data.get('tables', {})

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f8f9fa;
        }}

        .container {{
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }}

        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 10px;
            margin-bottom: 30px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }}

        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}

        .header .meta {{
            font-size: 1.1em;
            opacity: 0.9;
        }}

        .section {{
            background: white;
            padding: 30px;
            margin-bottom: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}

        .section h2 {{
            color: #667eea;
            border-bottom: 3px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}

        .section h3 {{
            color: #764ba2;
            margin-top: 25px;
            margin-bottom: 15px;
        }}

        .stat-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}

        .stat-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
        }}

        .stat-card .value {{
            font-size: 2em;
            font-weight: bold;
            margin-bottom: 5px;
        }}

        .stat-card .label {{
            font-size: 0.9em;
            opacity: 0.9;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}

        table th {{
            background: #667eea;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: 600;
        }}

        table td {{
            padding: 12px;
            border-bottom: 1px solid #ddd;
        }}

        table tr:hover {{
            background: #f5f5f5;
        }}

        .recommendation {{
            padding: 15px;
            margin: 10px 0;
            border-radius: 5px;
            border-left: 4px solid #ccc;
        }}

        .recommendation.warning {{
            background: #fff3cd;
            border-left-color: #ffc107;
        }}

        .recommendation.info {{
            background: #d1ecf1;
            border-left-color: #17a2b8;
        }}

        .recommendation.success {{
            background: #d4edda;
            border-left-color: #28a745;
        }}

        .badge {{
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: 600;
        }}

        .badge.warning {{
            background: #ffc107;
            color: #000;
        }}

        .badge.info {{
            background: #17a2b8;
            color: #fff;
        }}

        .badge.success {{
            background: #28a745;
            color: #fff;
        }}

        .footer {{
            text-align: center;
            padding: 20px;
            color: #6c757d;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{title}</h1>
            <div class="meta">
                Generated: {metadata['date']}<br>
                {f"Author: {author}" if author else ""}
            </div>
        </div>

        <div class="section">
            <h2>Summary Statistics</h2>
            <div class="stat-grid">
                <div class="stat-card">
                    <div class="value">{metadata['n_cells_after']:,}</div>
                    <div class="label">Cells After QC</div>
                </div>
                <div class="stat-card">
                    <div class="value">{metadata['n_genes']:,}</div>
                    <div class="label">Total Genes</div>
                </div>
                <div class="stat-card">
                    <div class="value">{filtering['retention_rate']:.1%}</div>
                    <div class="label">Retention Rate</div>
                </div>
                <div class="stat-card">
                    <div class="value">{filtering['n_removed']:,}</div>
                    <div class="label">Cells Removed</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h2>Detailed Metrics</h2>
"""

        # Add summary table
        if 'summary' in tables:
            html += self._generate_table('Summary Statistics', tables['summary'])

        # Add metrics table
        if 'metrics' in tables:
            html += self._generate_table('QC Metrics', tables['metrics'])

        # Add recommendations
        html += """
        </div>

        <div class="section">
            <h2>Recommendations</h2>
"""

        if recommendations:
            for rec in recommendations:
                severity_class = rec['severity']
                html += f"""
            <div class="recommendation {severity_class}">
                <strong>{rec['category']}:</strong> <span class="badge {severity_class}">{rec['severity'].upper()}</span><br>
                {rec['message']}
            </div>
"""
        else:
            html += "<p>No specific recommendations. Your data looks good!</p>"

        html += f"""
        </div>

        <div class="footer">
            <p>Generated by scLucid - A Comprehensive System for Single-Cell Analysis</p>
            <p>Report generated on {metadata['date']}</p>
        </div>
    </div>

    <script>
        // Add any interactive JavaScript here
        document.addEventListener('DOMContentLoaded', function() {{
            console.log('QC Report loaded');
        }});
    </script>
</body>
</html>
"""

        return html

    def _generate_table(self, title: str, data: List[Dict[str, str]]) -> str:
        """Generate HTML table from data."""
        if not data:
            return ""

        html = f"<h3>{title}</h3>\n<table>\n"

        # Header
        html += "<thead><tr>"
        for key in data[0].keys():
            html += f"<th>{key}</th>"
        html += "</tr></thead>\n"

        # Body
        html += "<tbody>"
        for row in data:
            html += "<tr>"
            for value in row.values():
                html += f"<td>{value}</td>"
            html += "</tr>\n"
        html += "</tbody>\n"

        html += "</table>\n"

        return html


def generate_qc_html_report(
    adata: AnnData,
    output_path: str,
    adata_before: Optional[AnnData] = None,
    title: str = "Quality Control Report",
    author: Optional[str] = None,
) -> str:
    """
    Convenience function to generate QC HTML report.

    Args:
        adata: AnnData after QC
        output_path: Path to save report
        adata_before: AnnData before QC
        title: Report title
        author: Report author

    Returns:
        Path to generated report
    """
    reporter = EnhancedQCReport(adata, adata_before)
    reporter.generate_html_report(
        output_path=output_path,
        title=title,
        author=author,
    )

    return output_path


class InteractiveReportGenerator:
    """
    Generate interactive HTML reports with embedded JavaScript.

    Creates reports with interactive plots using Plotly.js.
    """

    def __init__(self, adata: AnnData):
        """
        Initialize interactive report generator.

        Args:
            adata: AnnData object with QC results
        """
        self.adata = adata

    def generate_interactive_html(
        self,
        output_path: str,
        title: str = "Interactive QC Report",
    ):
        """
        Generate interactive HTML report with embedded Plotly charts.

        Args:
            output_path: Path to save report
            title: Report title
        """
        try:
            import plotly.graph_objects as go
            import plotly.io as pio
        except ImportError:
            raise ImportError("Plotly is required for interactive reports")

        # Create plots
        plots = self._create_interactive_plots()

        # Generate HTML with embedded plots
        html = self._generate_interactive_html(title, plots)

        # Write to file
        with open(output_path, 'w') as f:
            f.write(html)

        log.info(f"Interactive report saved to {output_path}")

    def _create_interactive_plots(self) -> Dict[str, 'go.Figure']:
        """Create interactive Plotly figures."""
        plots = {}

        # Violin plot for key metrics
        metrics_to_plot = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']
        available_metrics = [m for m in metrics_to_plot if m in self.adata.obs]

        if available_metrics:
            from plotly.subplots import make_subplots

            fig = make_subplots(
                rows=1, cols=len(available_metrics),
                subplot_titles=available_metrics,
            )

            for i, metric in enumerate(available_metrics):
                fig.add_trace(
                    go.Violin(
                        y=self.adata.obs[metric],
                        name=metric,
                        box_visible=True,
                        meanline_visible=True,
                    ),
                    row=1, col=i+1,
                )

            fig.update_layout(
                title_text="QC Metrics Distribution",
                showlegend=False,
                height=400,
            )

            plots['violin'] = fig

        return plots

    def _generate_interactive_html(
        self,
        title: str,
        plots: Dict[str, 'go.Figure'],
    ) -> str:
        """Generate HTML with embedded Plotly charts."""
        import plotly.io as pio

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f8f9fa;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
        }}
        .plot-container {{
            background: white;
            padding: 20px;
            margin-bottom: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>{title}</h1>
            <p>Interactive Quality Control Report</p>
        </div>
"""

        # Add plots
        for plot_name, fig in plots.items():
            plot_html = pio.to_html(fig, full_html=False, include_plotlyjs=False)
            html += f"""
        <div class="plot-container">
            <h3>{plot_name.replace('_', ' ').title()}</h3>
            {plot_html}
        </div>
"""

        html += """
    </div>

    <script>
        // Initialize plots
        console.log('Interactive QC Report loaded');
    </script>
</body>
</html>
"""

        return html
