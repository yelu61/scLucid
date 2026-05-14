"""
Self-contained HTML audit report for a scLucid analysis.

Walks ``adata.uns["sclucid"]`` and renders the recorded review summaries,
workflow configs, config lineage, and contract validation outcomes into one
self-contained HTML file that can be shared with reviewers or collaborators.

The output is pure HTML + inline CSS — no JavaScript, no external assets —
so the file is portable, prints cleanly, and can be opened from any browser
or attached to an email.
"""

from __future__ import annotations

import html as _html
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Optional, Union

from anndata import AnnData

from .contracts import Modules, SCLUCID_ROOT, UnsKeys

log = logging.getLogger(__name__)

# Render order for module sections; unknown modules are appended afterwards.
_DEFAULT_MODULE_ORDER = (
    Modules.QC,
    Modules.PREPROCESS,
    Modules.ANALYSIS,
    Modules.TUMOR,
    Modules.TOOLS,
)

# Reserved top-level keys under ``adata.uns["sclucid"]`` that are not
# per-module sections and should be rendered separately.
_RESERVED_ROOT_KEYS = {
    UnsKeys.PIPELINE_CONTEXT,
    UnsKeys.ANALYSIS_CONTEXT,
    UnsKeys.NAMESPACE_METADATA,
}


def export_audit_report(
    adata: AnnData,
    out: Union[str, Path],
    *,
    title: Optional[str] = None,
    include_full_config: bool = True,
) -> Path:
    """
    Render a self-contained HTML audit report from an scLucid-annotated AnnData.

    Parameters
    ----------
    adata : AnnData
        Data that has been processed by one or more scLucid workflow stages,
        i.e. ``adata.uns["sclucid"]`` is populated.
    out : str or Path
        Output HTML file path. Parent directories are created if needed.
    title : str, optional
        Title for the HTML page; defaults to ``"scLucid Analysis Audit"``.
    include_full_config : bool, default=True
        If ``True``, expand the full ``workflow_config`` dictionary for each
        stage inside collapsible ``<details>`` blocks. Set to ``False`` for a
        more compact report when configs are very large.

    Returns:
    -------
    Path
        The absolute path of the written HTML file.

    Examples:
    --------
    >>> import scLucid as scl
    >>> adata = scl.run_pipeline(adata, stages=["qc", "preprocess", "analysis"])
    >>> scl.export_audit_report(adata, "results/audit_report.html")
    """
    out_path = Path(out).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    sclucid_data = adata.uns.get(SCLUCID_ROOT, {}) or {}
    if not isinstance(sclucid_data, Mapping):
        log.warning(
            "adata.uns[%r] is not a mapping (type=%s); rendering empty report.",
            SCLUCID_ROOT,
            type(sclucid_data).__name__,
        )
        sclucid_data = {}

    html_doc = _render_html(
        adata=adata,
        sclucid_data=sclucid_data,
        title=title or "scLucid Analysis Audit",
        include_full_config=include_full_config,
    )
    out_path.write_text(html_doc, encoding="utf-8")
    log.info("scLucid audit report written to %s", out_path)
    return out_path


# ---------------------------------------------------------------------------
# Internal renderers
# ---------------------------------------------------------------------------


def _render_html(
    *,
    adata: AnnData,
    sclucid_data: Mapping[str, Any],
    title: str,
    include_full_config: bool,
) -> str:
    parts: list[str] = []
    parts.append(_render_head(title))
    parts.append("<body>")
    parts.append(_render_header(adata, title))
    parts.append(_render_dataset_panel(adata, sclucid_data))
    parts.append(_render_pipeline_context_panel(sclucid_data))

    module_keys = _ordered_module_keys(sclucid_data)
    if module_keys:
        parts.append('<section class="modules">')
        parts.append("<h2>Module Reports</h2>")
        for module in module_keys:
            section = sclucid_data.get(module)
            if isinstance(section, Mapping):
                parts.append(
                    _render_module_section(
                        module=module,
                        section=section,
                        include_full_config=include_full_config,
                    )
                )
        parts.append("</section>")
    else:
        parts.append(
            '<section class="empty"><p>No module results were recorded on this '
            "AnnData. Run a scLucid workflow (for example "
            "<code>scl.run_pipeline(adata)</code>) before exporting the report.</p>"
            "</section>"
        )

    parts.append(_render_footer())
    parts.append("</body></html>")
    return "\n".join(parts)


def _render_head(title: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{_html.escape(title)}</title>
<style>
:root {{
  --fg: #1d1f23;
  --muted: #5f6368;
  --bg: #fbfbfd;
  --panel: #ffffff;
  --border: #e1e4e8;
  --accent: #1d6f8a;
  --warn-bg: #fff5d6;
  --warn-fg: #6b4f00;
  --err-bg: #ffe0e0;
  --err-fg: #8b0000;
  --ok-bg: #e8f6ea;
  --ok-fg: #1e7a32;
  --mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  padding: 0;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", sans-serif;
  color: var(--fg);
  background: var(--bg);
  line-height: 1.5;
}}
header.report-header {{
  background: linear-gradient(135deg, #1d6f8a 0%, #154d62 100%);
  color: #fff;
  padding: 2rem 2.5rem 1.5rem;
}}
header.report-header h1 {{ margin: 0 0 0.25rem; font-size: 1.6rem; font-weight: 600; }}
header.report-header .meta {{ font-size: 0.85rem; opacity: 0.9; }}
main, section, .empty {{
  max-width: 1080px;
  margin: 1.5rem auto;
  padding: 0 2rem;
}}
.panel {{
  background: var(--panel);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 1rem 1.25rem;
  margin-bottom: 1rem;
}}
.panel h2, .panel h3 {{ margin-top: 0; color: var(--accent); }}
h2 {{ font-size: 1.2rem; margin-bottom: 0.75rem; }}
h3 {{ font-size: 1.05rem; margin-bottom: 0.5rem; }}
table.kv {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
table.kv td {{
  padding: 0.3rem 0.5rem;
  border-bottom: 1px solid var(--border);
  vertical-align: top;
}}
table.kv td.label {{
  font-weight: 600;
  width: 30%;
  color: var(--muted);
}}
ul.steps {{ margin: 0; padding-left: 1.25rem; }}
ul.steps li {{ font-family: var(--mono); font-size: 0.85rem; padding: 0.1rem 0; }}
.badge {{
  display: inline-block;
  padding: 0.15rem 0.55rem;
  border-radius: 999px;
  font-size: 0.75rem;
  font-weight: 600;
  margin-right: 0.4rem;
}}
.badge.ok {{ background: var(--ok-bg); color: var(--ok-fg); }}
.badge.warn {{ background: var(--warn-bg); color: var(--warn-fg); }}
.badge.err {{ background: var(--err-bg); color: var(--err-fg); }}
.badge.muted {{ background: #eef0f3; color: var(--muted); }}
details {{
  margin: 0.5rem 0;
  padding: 0.4rem 0.75rem;
  border: 1px solid var(--border);
  border-radius: 4px;
  background: #fafbfc;
}}
details > summary {{
  cursor: pointer;
  font-weight: 600;
  color: var(--accent);
}}
pre.config {{
  margin: 0.5rem 0 0;
  padding: 0.75rem;
  background: #0f1419;
  color: #e7eaee;
  border-radius: 4px;
  font-family: var(--mono);
  font-size: 0.8rem;
  overflow-x: auto;
  white-space: pre-wrap;
  word-break: break-word;
}}
.warning-list {{
  margin: 0.5rem 0;
  padding: 0.6rem 0.9rem;
  background: var(--warn-bg);
  color: var(--warn-fg);
  border-left: 3px solid #c89200;
  border-radius: 3px;
}}
.warning-list ul {{ margin: 0; padding-left: 1.2rem; }}
.module-header {{
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  flex-wrap: wrap;
  gap: 0.5rem;
}}
.module-header h3 {{ margin: 0; }}
.module-header .stats {{ font-size: 0.8rem; color: var(--muted); }}
footer.report-footer {{
  text-align: center;
  padding: 1.25rem;
  font-size: 0.8rem;
  color: var(--muted);
}}
@media print {{
  body {{ background: #fff; }}
  header.report-header {{ background: #1d6f8a !important; print-color-adjust: exact; }}
  details {{ background: transparent; }}
  details > summary {{ list-style: none; }}
  details[open] > summary::after {{ content: ""; }}
}}
</style>
</head>
"""


def _render_header(adata: AnnData, title: str) -> str:
    try:
        from .. import __version__ as scl_version
    except Exception:
        scl_version = "unknown"

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    shape = f"{adata.n_obs:,} cells × {adata.n_vars:,} genes"
    return f"""<header class="report-header">
  <h1>{_html.escape(title)}</h1>
  <div class="meta">
    Generated {timestamp} &middot; scLucid {_html.escape(str(scl_version))} &middot; {shape}
  </div>
</header>
"""


def _render_dataset_panel(adata: AnnData, sclucid_data: Mapping[str, Any]) -> str:
    rows = [
        ("Cells (n_obs)", f"{adata.n_obs:,}"),
        ("Genes (n_vars)", f"{adata.n_vars:,}"),
    ]
    if adata.layers:
        rows.append(("Layers", ", ".join(sorted(adata.layers.keys()))))
    if adata.obsm:
        rows.append(("Embeddings (obsm)", ", ".join(sorted(adata.obsm.keys()))))
    if adata.obs.columns.size:
        rows.append(("obs columns", f"{adata.obs.columns.size}"))
    if adata.var.columns.size:
        rows.append(("var columns", f"{adata.var.columns.size}"))

    analysis_ctx = sclucid_data.get(UnsKeys.ANALYSIS_CONTEXT)
    if isinstance(analysis_ctx, Mapping):
        for label in ("species", "dataset_type", "tissue", "cancer_type"):
            if label in analysis_ctx and analysis_ctx[label] is not None:
                rows.append((label, str(analysis_ctx[label])))

    return _render_panel("Dataset", _render_kv_table(rows))


def _render_pipeline_context_panel(sclucid_data: Mapping[str, Any]) -> str:
    pipeline_ctx = sclucid_data.get(UnsKeys.PIPELINE_CONTEXT)
    if not isinstance(pipeline_ctx, Mapping):
        return ""

    rows = []
    for key, value in pipeline_ctx.items():
        if key == "config_lineage":
            continue
        rows.append((key, _scalar_html(value)))

    body = _render_kv_table(rows)
    lineage = pipeline_ctx.get("config_lineage")
    if isinstance(lineage, Mapping):
        body += _render_collapsible("Pipeline config lineage", _pretty_json_html(lineage))
    return _render_panel("Pipeline Context", body)


def _render_module_section(
    *,
    module: str,
    section: Mapping[str, Any],
    include_full_config: bool,
) -> str:
    review = section.get(UnsKeys.REVIEW_SUMMARY)
    workflow_config = section.get(UnsKeys.WORKFLOW_CONFIG)
    steps_executed = section.get(UnsKeys.STEPS_EXECUTED)
    contract = section.get(UnsKeys.CONTRACT)
    lineage = section.get(UnsKeys.CONFIG_LINEAGE)
    artifacts = section.get(UnsKeys.ARTIFACTS)
    errors = section.get(UnsKeys.ERRORS)

    if isinstance(review, Mapping):
        review_steps = review.get("steps_executed")
        if steps_executed is None and review_steps:
            steps_executed = review_steps
        review_warnings = review.get("warnings", [])
        if workflow_config is None and isinstance(review.get("config"), Mapping):
            workflow_config = review["config"]
        if contract is None and isinstance(review.get("contract"), Mapping):
            contract = review["contract"]
        if lineage is None and isinstance(review.get("config_lineage"), Mapping):
            lineage = review["config_lineage"]
        if artifacts is None and isinstance(review.get("artifacts"), Mapping):
            artifacts = review["artifacts"]
    else:
        review_warnings = []

    title_html = _html.escape(module.capitalize())
    badge_html = _module_badges(review=review, contract=contract, errors=errors)
    stats_html = _module_stats(steps_executed=steps_executed, warnings=review_warnings)

    parts = [
        '<div class="panel module">',
        f'<div class="module-header"><h3>{title_html}</h3>'
        f'<div class="stats">{badge_html}{stats_html}</div></div>',
    ]

    if isinstance(review, Mapping):
        meta_rows = []
        for key in ("workflow_name", "schema_version", "generated_at"):
            if key in review:
                meta_rows.append((key, _scalar_html(review[key])))
        if meta_rows:
            parts.append(_render_kv_table(meta_rows))

    if steps_executed:
        parts.append("<h4>Steps executed</h4>")
        items = "".join(
            f"<li>{_html.escape(str(step))}</li>" for step in steps_executed
        )
        parts.append(f"<ul class='steps'>{items}</ul>")

    if review_warnings:
        items = "".join(
            f"<li>{_html.escape(str(w))}</li>" for w in review_warnings
        )
        parts.append(
            f"<div class='warning-list'><strong>{len(review_warnings)} warning(s):</strong>"
            f"<ul>{items}</ul></div>"
        )

    if errors:
        if isinstance(errors, (list, tuple)):
            items = "".join(
                f"<li>{_html.escape(str(e))}</li>" for e in errors
            )
        else:
            items = f"<li>{_html.escape(str(errors))}</li>"
        parts.append(
            f"<div class='warning-list' style='background:var(--err-bg);"
            f"color:var(--err-fg);border-left-color:#8b0000'>"
            f"<strong>{(len(errors) if hasattr(errors, '__len__') else 1)} error(s):</strong>"
            f"<ul>{items}</ul></div>"
        )

    if include_full_config and isinstance(workflow_config, Mapping):
        parts.append(
            _render_collapsible("Effective workflow_config", _pretty_json_html(workflow_config))
        )

    if isinstance(lineage, Mapping) and lineage:
        parts.append(
            _render_collapsible("Config lineage", _pretty_json_html(lineage))
        )

    if isinstance(contract, Mapping) and contract:
        parts.append(
            _render_collapsible("Contract validation", _pretty_json_html(contract))
        )

    if isinstance(artifacts, Mapping) and artifacts:
        rows = []
        for label, value in artifacts.items():
            rows.append((label, _scalar_html(value)))
        parts.append("<h4>Artifacts</h4>")
        parts.append(_render_kv_table(rows))

    # Render any additional non-standard keys (e.g. execution_trace,
    # recommended_params snapshots) so module-specific evidence does not
    # disappear from the audit just because it doesn't fit the canonical
    # envelope.
    rendered_keys = {
        UnsKeys.REVIEW_SUMMARY,
        UnsKeys.WORKFLOW_CONFIG,
        UnsKeys.STEPS_EXECUTED,
        UnsKeys.CONTRACT,
        UnsKeys.CONFIG_LINEAGE,
        UnsKeys.ARTIFACTS,
        UnsKeys.ERRORS,
        UnsKeys.NAMESPACE_METADATA,
    }
    for key, value in section.items():
        if key in rendered_keys:
            continue
        if not isinstance(value, (Mapping, list, tuple)):
            continue
        if not value:
            continue
        parts.append(_render_collapsible(str(key), _pretty_json_html(value)))

    parts.append("</div>")
    return "\n".join(parts)


def _render_footer() -> str:
    return (
        '<footer class="report-footer">Generated by '
        '<code>scLucid.export_audit_report</code></footer>'
    )


# ---------------------------------------------------------------------------
# Small HTML helpers
# ---------------------------------------------------------------------------


def _render_panel(title: str, body_html: str) -> str:
    return f"""<section class="panel">
  <h2>{_html.escape(title)}</h2>
  {body_html}
</section>"""


def _render_kv_table(rows: list[tuple[str, str]]) -> str:
    if not rows:
        return "<p><em>No values.</em></p>"
    cells = "".join(
        f'<tr><td class="label">{_html.escape(str(label))}</td>'
        f'<td>{value}</td></tr>'
        for label, value in rows
    )
    return f'<table class="kv">{cells}</table>'


def _render_collapsible(label: str, inner_html: str) -> str:
    return (
        f'<details><summary>{_html.escape(label)}</summary>'
        f'{inner_html}</details>'
    )


def _pretty_json_html(value: Any) -> str:
    encoded = json.dumps(value, indent=2, default=str, sort_keys=True)
    return f'<pre class="config">{_html.escape(encoded)}</pre>'


def _scalar_html(value: Any) -> str:
    if isinstance(value, (Mapping, list, tuple, set)):
        return _pretty_json_html(value)
    return _html.escape(str(value))


def _module_badges(
    *,
    review: Any,
    contract: Any,
    errors: Any,
) -> str:
    badges: list[str] = []
    if errors:
        badges.append('<span class="badge err">errors</span>')
    if isinstance(contract, Mapping):
        valid = contract.get("valid")
        if valid is True:
            badges.append('<span class="badge ok">contract ✓</span>')
        elif valid is False:
            badges.append('<span class="badge err">contract ✗</span>')
    if isinstance(review, Mapping):
        warnings = review.get("warnings") or []
        if warnings:
            badges.append(
                f'<span class="badge warn">{len(warnings)} warning(s)</span>'
            )
    if not badges:
        badges.append('<span class="badge muted">recorded</span>')
    return "".join(badges)


def _module_stats(*, steps_executed: Any, warnings: Any) -> str:
    bits: list[str] = []
    if steps_executed:
        bits.append(f"{len(steps_executed)} step(s)")
    if warnings:
        bits.append(f"{len(warnings)} warning(s)")
    return " &middot; ".join(bits)


def _ordered_module_keys(sclucid_data: Mapping[str, Any]) -> list[str]:
    keys: list[str] = []
    seen: set[str] = set()
    for module in _DEFAULT_MODULE_ORDER:
        if module in sclucid_data and isinstance(sclucid_data[module], Mapping):
            keys.append(module)
            seen.add(module)
    for key in sclucid_data:
        if key in seen or key in _RESERVED_ROOT_KEYS:
            continue
        if isinstance(sclucid_data[key], Mapping):
            keys.append(key)
            seen.add(key)
    return keys


__all__ = ["export_audit_report"]
