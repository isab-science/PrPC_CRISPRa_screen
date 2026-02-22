from __future__ import annotations

import os
import sys
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from prpcscreen.visualization.volcano_and_flashlight_plots import (
    _build_sublibrary_series,
    _candidate_mask,
    flashlight_plot,
    write_interactive_volcano_html,
)
from prpcscreen.visualization.plotly_exports import write_plotly_interactive_html

DEBUG_ENV_DEFAULT = os.environ.get("PRPCSCREEN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}


def debug_log(message: str, enabled: bool) -> None:
    """Emit debug lines for candidate landscape plotting."""
    if enabled:
        print(f"[08_plot_candidate_landscape] {message}", file=sys.stderr)


def _write_interactive_flashlight_html(
    df: pd.DataFrame,
    output_png: str,
    rank_col: str = "Mean_log2",
    genomics_excel: str | None = None,
) -> Path:
    keep = _candidate_mask(df)
    plot_df = df.loc[keep].copy()
    plot_df[rank_col] = pd.to_numeric(plot_df[rank_col], errors="coerce")
    plot_df = plot_df[np.isfinite(plot_df[rank_col])].copy()

    pos_mask = plot_df.get("Is_pos_ctrl", pd.Series(False, index=plot_df.index)).astype(bool)
    nt_mask = plot_df.get("Is_NT_ctrl", pd.Series(False, index=plot_df.index)).astype(bool) & ~pos_mask
    gene_mask = ~(pos_mask | nt_mask)

    if "Gene_symbol" in plot_df.columns:
        labels = plot_df["Gene_symbol"].fillna("").astype(str).str.strip()
    else:
        labels = pd.Series("", index=plot_df.index, dtype=object)
    labels = labels.mask(labels.str.lower().eq("nan"), "")
    labels = labels.mask(labels.eq("") & pos_mask, "PRNP")
    labels = labels.mask(labels.eq("") & nt_mask, "NT_CTRL")
    if "Entrez_ID" in plot_df.columns:
        entrez = pd.to_numeric(plot_df["Entrez_ID"], errors="coerce")
        fallback = entrez.map(lambda v: f"Entrez {int(v)}" if pd.notna(v) else "")
        labels = labels.mask(labels.eq(""), fallback)
    labels = labels.mask(labels.eq(""), "UNLABELED")

    if genomics_excel:
        sublibrary_series, sublibrary_options = _build_sublibrary_series(plot_df, genomics_excel)
    else:
        sublibrary_series = pd.Series([""] * len(plot_df), index=plot_df.index, dtype=object)
        sublibrary_options = []
    sublibrary_series = sublibrary_series.fillna("").astype(str).str.strip()
    sublibrary_series = sublibrary_series.where(sublibrary_series != "", "Unknown")

    plot_df = plot_df.assign(
        _rank_value=plot_df[rank_col].astype(float),
        _label=labels.astype(str),
        _is_gene=gene_mask.astype(bool),
        _sublibrary=sublibrary_series,
    )
    plot_df = plot_df.sort_values("_rank_value", kind="mergesort").reset_index(drop=True)
    plot_df["_rank"] = np.arange(1, len(plot_df) + 1, dtype=int)

    customdata = [
        [label, str(label).upper(), bool(is_gene), sublibrary]
        for label, is_gene, sublibrary in zip(
            plot_df["_label"].tolist(),
            plot_df["_is_gene"].tolist(),
            plot_df["_sublibrary"].tolist(),
        )
    ]

    traces = [
        {
            "x": plot_df["_rank"].astype(float).tolist(),
            "y": plot_df["_rank_value"].astype(float).tolist(),
            "customdata": customdata,
            "mode": "lines",
            "type": "scatter",
            "line": {"color": "#1f77b4", "width": 1.2},
            "name": rank_col,
            "hovertemplate": (
                "Gene: %{customdata[0]}<br>Sublibrary: %{customdata[3]}<br>Rank: %{x}<br>Mean log2: %{y:.4f}<extra></extra>"
            ),
        },
        {
            "x": [],
            "y": [],
            "text": [],
            "mode": "markers",
            "type": "scatter",
            "marker": {"size": 8, "color": "#d50000", "line": {"color": "#ffffff", "width": 0.8}},
            "showlegend": False,
            "visible": False,
            "hovertemplate": "Gene: %{text}<br>Rank: %{x}<br>Mean log2: %{y:.4f}<extra>Labeled</extra>",
        },
        {
            "x": [],
            "y": [],
            "text": [],
            "mode": "text",
            "type": "scatter",
            "textposition": "top center",
            "textfont": {"size": 11, "color": "#111111"},
            "showlegend": False,
            "visible": False,
            "hoverinfo": "skip",
        },
    ]
    layout = {
        "paper_bgcolor": "#f6f6f6",
        "plot_bgcolor": "#ffffff",
        "margin": {"l": 70, "r": 30, "t": 35, "b": 60},
        "xaxis": {"title": "Rank"},
        "yaxis": {"title": rank_col},
    }
    extra_controls_html = (
        "    <div class=\"controls\">\n"
        "      <label for=\"flashlight-sublibrary-filter\">Sublibrary</label>\n"
        "      <select id=\"flashlight-sublibrary-filter\"></select>\n"
        "      <span id=\"flashlight-sublibrary-status\" style=\"color: #475569; font-size: 13px; min-height: 1em;\"></span>\n"
        "    </div>\n"
        "    <div class=\"controls\">\n"
        "      <label for=\"label-mode\">Label points</label>\n"
        "      <select id=\"label-mode\">\n"
        "        <option value=\"top10\">Strongest modulators: Top 10</option>\n"
        "        <option value=\"top20\">Strongest modulators: Top 20</option>\n"
        "        <option value=\"top50\">Strongest modulators: Top 50</option>\n"
        "        <option value=\"genes\">Gene list (textbox)</option>\n"
        "      </select>\n"
        "      <input id=\"label-gene-input\" type=\"text\" placeholder=\"GENE1, GENE2\" style=\"min-width: 320px; border: 1px solid #d0d7de; border-radius: 6px; padding: 6px 8px; font-size: 13px;\" />\n"
        "      <button id=\"label-apply\" type=\"button\">Apply</button>\n"
        "      <button id=\"label-clear\" type=\"button\">Clear</button>\n"
        "      <span id=\"label-status\" style=\"color: #475569; font-size: 13px; min-height: 1em;\"></span>\n"
        "    </div>\n"
    )
    extra_script = (
        "const TRACE_LINE = 0;\n"
        "const TRACE_LABEL_MARKERS = 1;\n"
        "const TRACE_LABEL_TEXT = 2;\n"
        "const SUBLIBRARY_OPTIONS = "
        + json.dumps([str(opt) for opt in sublibrary_options], separators=(",", ":"))
        + ";\n"
        "const FLASHLIGHT_POINT_SOURCE = {\n"
        "  x: traces[TRACE_LINE].x || [],\n"
        "  y: traces[TRACE_LINE].y || [],\n"
        "  custom: traces[TRACE_LINE].customdata || []\n"
        "};\n"
        "function normalizedSub(value) {\n"
        "  return (value || '').toString().trim().toLowerCase();\n"
        "}\n"
        "function activeSublibrary() {\n"
        "  const sel = document.getElementById('flashlight-sublibrary-filter');\n"
        "  const value = sel ? sel.value : '__all__';\n"
        "  return value && value !== '__all__' ? value : '';\n"
        "}\n"
        "function populateSublibraryFilter() {\n"
        "  const sel = document.getElementById('flashlight-sublibrary-filter');\n"
        "  sel.innerHTML = '';\n"
        "  const all = document.createElement('option');\n"
        "  all.value = '__all__';\n"
        "  all.textContent = 'All sublibraries';\n"
        "  sel.appendChild(all);\n"
        "  for (const sub of SUBLIBRARY_OPTIONS) {\n"
        "    if (!sub) continue;\n"
        "    const opt = document.createElement('option');\n"
        "    opt.value = sub;\n"
        "    opt.textContent = sub;\n"
        "    sel.appendChild(opt);\n"
        "  }\n"
        "}\n"
        "function filteredPoints() {\n"
        "  const selected = activeSublibrary();\n"
        "  const key = normalizedSub(selected);\n"
        "  const out = [];\n"
        "  const xs = FLASHLIGHT_POINT_SOURCE.x || [];\n"
        "  const ys = FLASHLIGHT_POINT_SOURCE.y || [];\n"
        "  const cs = FLASHLIGHT_POINT_SOURCE.custom || [];\n"
        "  const n = Math.min(xs.length, ys.length, cs.length);\n"
        "  for (let i = 0; i < n; i += 1) {\n"
        "    const c = cs[i];\n"
        "    const sub = Array.isArray(c) ? String(c[3] || '').trim() : '';\n"
        "    if (selected && normalizedSub(sub) !== key) continue;\n"
        "    out.push({x: Number(xs[i]), y: Number(ys[i]), custom: c});\n"
        "  }\n"
        "  return out.filter((p) => Number.isFinite(p.x) && Number.isFinite(p.y));\n"
        "}\n"
        "function applySublibraryFilter() {\n"
        "  const points = filteredPoints();\n"
        "  const xVals = points.map((p) => p.x);\n"
        "  const yVals = points.map((p) => p.y);\n"
        "  const cVals = points.map((p) => p.custom);\n"
        "  Plotly.restyle(PLOT_DIV_ID, {x: [xVals], y: [yVals], customdata: [cVals]}, [TRACE_LINE]);\n"
        "  const sub = activeSublibrary();\n"
        "  const label = sub ? `Sublibrary: ${sub}` : 'All sublibraries';\n"
        "  document.getElementById('flashlight-sublibrary-status').textContent = `${label} (${points.length.toLocaleString()} points)`;\n"
        "  clearLabels('');\n"
        "}\n"
        "function parseGeneInput(raw) {\n"
        "  if (!raw) return [];\n"
        "  const tokens = raw\n"
        "    .split(/[\\s,;]+/)\n"
        "    .map((t) => t.trim())\n"
        "    .filter(Boolean)\n"
        "    .map((t) => t.toUpperCase());\n"
        "  return [...new Set(tokens)];\n"
        "}\n"
        "function clearLabels(statusMessage = '') {\n"
        "  Plotly.restyle(PLOT_DIV_ID, {x: [[]], y: [[]], text: [[]], visible: [false]}, [TRACE_LABEL_MARKERS]);\n"
        "  Plotly.restyle(PLOT_DIV_ID, {x: [[]], y: [[]], text: [[]], visible: [false]}, [TRACE_LABEL_TEXT]);\n"
        "  document.getElementById('label-status').textContent = statusMessage;\n"
        "}\n"
        "function renderLabels(points, statusMessage) {\n"
        "  const xVals = points.map((p) => p.x);\n"
        "  const yVals = points.map((p) => p.y);\n"
        "  const labels = points.map((p) => p.gene || p.geneUpper || 'UNLABELED');\n"
        "  const hasPoints = points.length > 0;\n"
        "  Plotly.restyle(PLOT_DIV_ID, {x: [xVals], y: [yVals], text: [labels], visible: [hasPoints]}, [TRACE_LABEL_MARKERS]);\n"
        "  Plotly.restyle(PLOT_DIV_ID, {x: [xVals], y: [yVals], text: [labels], visible: [hasPoints]}, [TRACE_LABEL_TEXT]);\n"
        "  document.getElementById('label-status').textContent = statusMessage;\n"
        "}\n"
        "function buildBestGeneMap() {\n"
        "  const tr = traces[TRACE_LINE] || {};\n"
        "  const xs = tr.x || [];\n"
        "  const ys = tr.y || [];\n"
        "  const custom = tr.customdata || [];\n"
        "  const best = new Map();\n"
        "  const n = Math.min(xs.length, ys.length, custom.length);\n"
        "  for (let i = 0; i < n; i += 1) {\n"
        "    const c = custom[i];\n"
        "    const gene = Array.isArray(c) ? String(c[0] || '').trim() : '';\n"
        "    const geneUpper = Array.isArray(c) ? String(c[1] || '').trim().toUpperCase() : gene.toUpperCase();\n"
        "    const isGene = Array.isArray(c) ? Boolean(c[2]) : true;\n"
        "    if (!isGene || !geneUpper) continue;\n"
        "    const x = Number(xs[i]);\n"
        "    const y = Number(ys[i]);\n"
        "    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;\n"
        "    const rec = {x, y, gene, geneUpper};\n"
        "    const prev = best.get(geneUpper);\n"
        "    if (!prev || Math.abs(y) > Math.abs(prev.y)) best.set(geneUpper, rec);\n"
        "  }\n"
        "  return best;\n"
        "}\n"
        "function selectTopStrongest(n) {\n"
        "  const points = Array.from(buildBestGeneMap().values());\n"
        "  points.sort((a, b) => (Math.abs(b.y) - Math.abs(a.y)) || (b.y - a.y));\n"
        "  return points.slice(0, n);\n"
        "}\n"
        "function applyLabels() {\n"
        "  const mode = document.getElementById('label-mode').value;\n"
        "  if (mode === 'genes') {\n"
        "    const requested = parseGeneInput(document.getElementById('label-gene-input').value);\n"
        "    if (requested.length === 0) {\n"
        "      clearLabels('Enter at least one gene symbol.');\n"
        "      return;\n"
        "    }\n"
        "    const bestByGene = buildBestGeneMap();\n"
        "    const found = [];\n"
        "    const missing = [];\n"
        "    for (const g of requested) {\n"
        "      if (bestByGene.has(g)) found.push(bestByGene.get(g));\n"
        "      else missing.push(g);\n"
        "    }\n"
        "    let msg = `Labeled ${found.length} gene(s) from textbox.`;\n"
        "    if (missing.length) msg += ` Not found: ${missing.join(', ')}`;\n"
        "    if (found.length === 0) {\n"
        "      clearLabels(msg);\n"
        "      return;\n"
        "    }\n"
        "    renderLabels(found, msg);\n"
        "    return;\n"
        "  }\n"
        "  const n = mode === 'top20' ? 20 : (mode === 'top50' ? 50 : 10);\n"
        "  const selected = selectTopStrongest(n);\n"
        "  renderLabels(selected, `Labeled ${selected.length} strongest modulator(s) (top ${n}).`);\n"
        "}\n"
        "function syncLabelInputState() {\n"
        "  const isGeneMode = document.getElementById('label-mode').value === 'genes';\n"
        "  const input = document.getElementById('label-gene-input');\n"
        "  input.disabled = !isGeneMode;\n"
        "  if (isGeneMode) {\n"
        "    input.placeholder = 'GENE1, GENE2';\n"
        "  } else {\n"
        "    input.placeholder = 'Switch mode to Gene list to type genes';\n"
        "  }\n"
        "}\n"
        "document.getElementById('label-mode').addEventListener('change', syncLabelInputState);\n"
        "document.getElementById('flashlight-sublibrary-filter').addEventListener('change', applySublibraryFilter);\n"
        "document.getElementById('label-apply').addEventListener('click', applyLabels);\n"
        "document.getElementById('label-clear').addEventListener('click', () => {\n"
        "  document.getElementById('label-gene-input').value = '';\n"
        "  clearLabels('');\n"
        "});\n"
        "document.getElementById('label-gene-input').addEventListener('keydown', (ev) => {\n"
        "  if (ev.key === 'Enter') {\n"
        "    ev.preventDefault();\n"
        "    applyLabels();\n"
        "  }\n"
        "});\n"
        "syncLabelInputState();\n"
        "populateSublibraryFilter();\n"
        "applySublibraryFilter();\n"
    )
    output_html = Path(output_png).with_name(Path(output_png).stem + "_interactive.html")
    return write_plotly_interactive_html(
        output_html=output_html,
        traces=traces,
        layout=layout,
        title="Candidate ranking flashlight",
        filename_base="candidate_flashlight_ranked_meanlog2",
        extra_controls_html=extra_controls_html,
        extra_script=extra_script,
    )


def run_landscape_cli() -> None:
    # Parse analyzed input plus output destinations.
    parser = argparse.ArgumentParser(description="Generate interactive volcano and flashlight plots.")
    parser.add_argument("input_csv")
    parser.add_argument("flashlight_png")
    parser.add_argument(
        "--volcano_html",
        required=True,
        help="Interactive volcano HTML output with on/off toggles by category.",
    )
    parser.add_argument(
        "--genomics_excel",
        default=None,
        help="Optional genomics workbook used to preload sublibrary metadata for volcano filtering.",
    )
    parser.add_argument("--debug", action="store_true", help="Enable verbose debug logging.")
    args = parser.parse_args()
    debug_enabled = DEBUG_ENV_DEFAULT or args.debug

    # Load analyzed table once and generate outputs.
    df = pd.read_csv(args.input_csv)
    debug_log(f"Loaded analyzed data: {args.input_csv} ({len(df)} rows)", debug_enabled)

    html_path = write_interactive_volcano_html(df, args.volcano_html, genomics_excel=args.genomics_excel)
    debug_log(f"Wrote interactive volcano HTML: {html_path}", debug_enabled)

    fig_f, _ = flashlight_plot(df)
    Path(args.flashlight_png).parent.mkdir(parents=True, exist_ok=True)
    fig_f.savefig(args.flashlight_png, dpi=600, bbox_inches="tight")
    plt.close(fig_f)
    debug_log(f"Wrote flashlight figure: {args.flashlight_png}", debug_enabled)

    flashlight_html = _write_interactive_flashlight_html(
        df,
        args.flashlight_png,
        rank_col="Mean_log2",
        genomics_excel=args.genomics_excel,
    )
    debug_log(f"Wrote interactive flashlight HTML: {flashlight_html}", debug_enabled)


if __name__ == "__main__":
    run_landscape_cli()

