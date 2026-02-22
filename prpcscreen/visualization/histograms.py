from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from prpcscreen.visualization.volcano_and_flashlight_plots import _build_sublibrary_series


def _histogram_inputs(df: pd.DataFrame, use_column: str) -> tuple[pd.Series, pd.Series, pd.Series, np.ndarray]:
    vals = pd.to_numeric(df[use_column], errors="coerce")
    is_nt = df.get("Is_NT_ctrl", pd.Series(False, index=df.index)).astype(bool)
    is_pos = df.get("Is_pos_ctrl", pd.Series(False, index=df.index)).astype(bool)
    is_gene = pd.to_numeric(df.get("Entrez_ID"), errors="coerce").notna() if "Entrez_ID" in df.columns else ~(is_nt | is_pos)

    valid = vals.notna() & (is_gene | is_nt | is_pos)
    vals_all = vals[valid]
    if vals_all.empty:
        bins = np.linspace(-1.0, 1.0, 70)
    else:
        lo = float(vals_all.min())
        hi = float(vals_all.max())
        if np.isclose(lo, hi):
            lo -= 0.5
            hi += 0.5
        bins = np.linspace(lo, hi, 70)

    gene_vals = vals[is_gene & vals.notna()]
    nt_vals = vals[is_nt & vals.notna()]
    pos_vals = vals[is_pos & vals.notna()]
    return gene_vals, nt_vals, pos_vals, bins


def three_histograms(df: pd.DataFrame, use_column: str):
    gene_vals, nt_vals, pos_vals, bins = _histogram_inputs(df, use_column)

    fig, ax = plt.subplots(figsize=(7, 5))

    # Back-to-front draw order to preserve visibility.
    ax.hist(pos_vals, bins=bins, color="#de2d26", alpha=0.75, edgecolor="#a50f15", linewidth=0.5, zorder=1)
    ax.hist(gene_vals, bins=bins, color="black", alpha=0.35, edgecolor="#8c8c8c", linewidth=0.5, zorder=2)
    ax.hist(nt_vals, bins=bins, color="#2171b5", alpha=0.75, edgecolor="#08519c", linewidth=0.5, zorder=3)

    ax.set_xlabel(use_column)
    ax.set_ylabel("Count")
    ax.tick_params(labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=(0, 0, 0, 0.35), edgecolor="#8c8c8c", linewidth=0.8, label="Genes in CRISPRa library"),
        plt.Rectangle((0, 0), 1, 1, facecolor=(33 / 255, 113 / 255, 181 / 255, 0.75), edgecolor="#08519c", linewidth=0.8, label="Non-targeting controls"),
        plt.Rectangle((0, 0), 1, 1, facecolor=(222 / 255, 45 / 255, 38 / 255, 0.75), edgecolor="#a50f15", linewidth=0.8, label="Positive controls"),
    ]
    ax.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=7.5)

    fig.tight_layout()
    return fig, ax


def write_interactive_histogram_html(
    df: pd.DataFrame,
    output_html: str | Path,
    use_column: str,
    genomics_excel: str | None = None,
    write_gzip_sidecar: bool = True,
) -> Path:
    gene_vals, nt_vals, pos_vals, bins = _histogram_inputs(df, use_column)
    vals = pd.to_numeric(df[use_column], errors="coerce")
    is_nt = df.get("Is_NT_ctrl", pd.Series(False, index=df.index)).astype(bool)
    is_pos = df.get("Is_pos_ctrl", pd.Series(False, index=df.index)).astype(bool)
    is_gene = pd.to_numeric(df.get("Entrez_ID"), errors="coerce").notna() if "Entrez_ID" in df.columns else ~(is_nt | is_pos)
    gene_mask = is_gene & vals.notna()
    if genomics_excel:
        sublibrary_series, sublibrary_options = _build_sublibrary_series(df, genomics_excel)
    else:
        sublibrary_series = pd.Series([""] * len(df), index=df.index, dtype=object)
        sublibrary_options = []
    sublibrary_options = [str(opt) for opt in sublibrary_options if str(opt).strip().lower() != "unknown"]
    sublibrary_series = sublibrary_series.fillna("").astype(str).str.strip()
    sublibrary_series = sublibrary_series.where(sublibrary_series != "", "Unknown")
    gene_points = pd.DataFrame(
        {
            "x": vals[gene_mask].astype(float),
            "sublibrary": sublibrary_series[gene_mask].astype(str),
        }
    ).reset_index(drop=True)
    start = float(bins[0])
    end = float(bins[-1])
    size = float((end - start) / max(1, (len(bins) - 1)))

    traces = [
        {
            "x": gene_points["x"].astype(float).tolist(),
            "customdata": gene_points["sublibrary"].tolist(),
            "type": "histogram",
            "name": "Genes in CRISPRa library",
            "marker": {"color": "rgba(0,0,0,0.35)", "line": {"color": "#8c8c8c", "width": 1}},
            "xbins": {"start": start, "end": end, "size": size},
            "hovertemplate": "Sublibrary: %{customdata}<br>Value: %{x:.4f}<extra></extra>",
        },
        {
            "x": nt_vals.astype(float).tolist(),
            "type": "histogram",
            "name": "Non-targeting controls",
            "marker": {"color": "rgba(33,113,181,0.75)", "line": {"color": "#08519c", "width": 1}},
            "xbins": {"start": start, "end": end, "size": size},
        },
        {
            "x": pos_vals.astype(float).tolist(),
            "type": "histogram",
            "name": "Positive controls",
            "marker": {"color": "rgba(222,45,38,0.75)", "line": {"color": "#a50f15", "width": 1}},
            "xbins": {"start": start, "end": end, "size": size},
        },
    ]
    layout = {
        "title": {"text": f"Histogram ({use_column})"},
        "paper_bgcolor": "#f6f6f6",
        "plot_bgcolor": "#ffffff",
        "barmode": "overlay",
        "margin": {"l": 70, "r": 40, "t": 60, "b": 60},
        "xaxis": {"title": use_column},
        "yaxis": {"title": "Count"},
        "showlegend": True,
        "legend": {"orientation": "v", "x": 1.02, "y": 0.5, "xanchor": "left", "yanchor": "middle"},
    }

    traces_json = json.dumps(traces, separators=(",", ":"))
    layout_json = json.dumps(layout, separators=(",", ":"))
    column_json = json.dumps(use_column, separators=(",", ":"))
    sublibrary_options_json = json.dumps([str(opt) for opt in sublibrary_options], separators=(",", ":"))

    html = (
        "<!doctype html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "  <title>Interactive Histogram</title>\n"
        "  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>\n"
        "  <style>\n"
        "    body { margin: 0; font-family: Segoe UI, Arial, sans-serif; background: #f6f6f6; color: #1f1f1f; }\n"
        "    .wrap { max-width: 1400px; margin: 0 auto; padding: 16px; }\n"
        "    .controls { display: flex; gap: 18px; flex-wrap: wrap; align-items: center; margin-bottom: 10px; font-size: 14px; }\n"
        "    .controls label { display: inline-flex; align-items: center; gap: 6px; user-select: none; }\n"
        "    .controls select { border: 1px solid #d0d7de; border-radius: 6px; padding: 6px 8px; font-size: 13px; background: #ffffff; }\n"
        "    .controls button { border: 1px solid #c9ced6; border-radius: 6px; background: #ffffff; padding: 6px 10px; cursor: pointer; font-size: 13px; }\n"
        "    .controls button:hover { background: #f4f6f8; }\n"
        "    #export-status { color: #475569; font-size: 13px; min-height: 1em; }\n"
        "    #hist { width: 100%; height: 78vh; min-height: 620px; border: 1px solid #d8d8d8; background: #ffffff; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n"
        "  <div class=\"wrap\">\n"
        "    <div class=\"controls\">\n"
        "      <label><input id=\"toggle-genes\" type=\"checkbox\" checked> Genes in CRISPRa library</label>\n"
        "      <label><input id=\"toggle-nt\" type=\"checkbox\" checked> Non-targeting controls</label>\n"
        "      <label><input id=\"toggle-pos\" type=\"checkbox\" checked> Positive controls</label>\n"
        "      <label for=\"hist-sublibrary-filter\">Gene sublibrary</label>\n"
        "      <select id=\"hist-sublibrary-filter\"></select>\n"
        "      <span id=\"hist-sublibrary-status\" style=\"color: #475569; font-size: 13px; min-height: 1em;\"></span>\n"
        "    </div>\n"
        "    <div class=\"controls\">\n"
        "      <label for=\"export-quality\">Publication-quality figure</label>\n"
        "      <select id=\"export-quality\">\n"
        "        <option value=\"low\">Low (1200x800)</option>\n"
        "        <option value=\"medium\" selected>Medium (1800x1200)</option>\n"
        "        <option value=\"high\">High (2400x1600)</option>\n"
        "      </select>\n"
        "      <button id=\"export-png\" type=\"button\">Export PNG</button>\n"
        "      <span id=\"export-status\"></span>\n"
        "    </div>\n"
        "    <div id=\"hist\"></div>\n"
        "  </div>\n"
        "  <script>\n"
        "    const COLUMN_NAME = "
        + column_json
        + ";\n"
        "    const traces = "
        + traces_json
        + ";\n"
        "    const HIST_GENE_SOURCE = { x: traces[0].x || [], sublibrary: traces[0].customdata || [] };\n"
        "    const HIST_SUBLIBRARY_OPTIONS = "
        + sublibrary_options_json
        + ";\n"
        "    const layout = "
        + layout_json
        + ";\n"
        "    Plotly.newPlot('hist', traces, layout, {responsive: true});\n"
        "    function normalizedSub(value) {\n"
        "      return (value || '').toString().trim().toLowerCase();\n"
        "    }\n"
        "    function activeSublibrary() {\n"
        "      const sel = document.getElementById('hist-sublibrary-filter');\n"
        "      const value = sel ? sel.value : '__all__';\n"
        "      return value && value !== '__all__' ? value : '';\n"
        "    }\n"
        "    function populateSublibraryFilter() {\n"
        "      const sel = document.getElementById('hist-sublibrary-filter');\n"
        "      sel.innerHTML = '';\n"
        "      const all = document.createElement('option');\n"
        "      all.value = '__all__';\n"
        "      all.textContent = 'All sublibraries';\n"
        "      sel.appendChild(all);\n"
        "      for (const sub of HIST_SUBLIBRARY_OPTIONS) {\n"
        "        if (!sub) continue;\n"
        "        const opt = document.createElement('option');\n"
        "        opt.value = sub;\n"
        "        opt.textContent = sub;\n"
        "        sel.appendChild(opt);\n"
        "      }\n"
        "    }\n"
        "    function applySublibraryFilter() {\n"
        "      const selected = activeSublibrary();\n"
        "      const key = normalizedSub(selected);\n"
        "      const xs = HIST_GENE_SOURCE.x || [];\n"
        "      const ss = HIST_GENE_SOURCE.sublibrary || [];\n"
        "      const n = Math.min(xs.length, ss.length);\n"
        "      const xVals = [];\n"
        "      const sVals = [];\n"
        "      for (let i = 0; i < n; i += 1) {\n"
        "        const sub = (ss[i] || 'Unknown').toString().trim() || 'Unknown';\n"
        "        if (selected && normalizedSub(sub) !== key) continue;\n"
        "        const xVal = Number(xs[i]);\n"
        "        if (!Number.isFinite(xVal)) continue;\n"
        "        xVals.push(xVal);\n"
        "        sVals.push(sub);\n"
        "      }\n"
        "      Plotly.restyle('hist', {x: [xVals], customdata: [sVals]}, [0]);\n"
        "      const label = selected ? `Sublibrary: ${selected}` : 'All sublibraries';\n"
        "      document.getElementById('hist-sublibrary-status').textContent = `${label} (${xVals.length.toLocaleString()} genes)`;\n"
        "    }\n"
        "    function applyVisibility() {\n"
        "      const genesOn = document.getElementById('toggle-genes').checked;\n"
        "      const ntOn = document.getElementById('toggle-nt').checked;\n"
        "      const posOn = document.getElementById('toggle-pos').checked;\n"
        "      Plotly.restyle('hist', {visible: genesOn ? true : 'legendonly'}, [0]);\n"
        "      Plotly.restyle('hist', {visible: ntOn ? true : 'legendonly'}, [1]);\n"
        "      Plotly.restyle('hist', {visible: posOn ? true : 'legendonly'}, [2]);\n"
        "    }\n"
        "    function getExportPreset() {\n"
        "      const mode = document.getElementById('export-quality').value;\n"
        "      if (mode === 'low') return { label: 'low', width: 1200, height: 800 };\n"
        "      if (mode === 'high') return { label: 'high', width: 2400, height: 1600 };\n"
        "      return { label: 'medium', width: 1800, height: 1200 };\n"
        "    }\n"
        "    function exportPublicationPng() {\n"
        "      const preset = getExportPreset();\n"
        "      const col = COLUMN_NAME.replace(/[^a-zA-Z0-9_]+/g, '_');\n"
        "      const filename = `distribution_${col}_${preset.label}`;\n"
        "      document.getElementById('export-status').textContent = `Exporting ${preset.label} PNG...`;\n"
        "      Plotly.downloadImage('hist', {\n"
        "        format: 'png',\n"
        "        filename,\n"
        "        width: preset.width,\n"
        "        height: preset.height,\n"
        "        scale: 1,\n"
        "      }).then(() => {\n"
        "        document.getElementById('export-status').textContent = `Saved ${filename}.png (${preset.width}x${preset.height}).`;\n"
        "      }).catch((err) => {\n"
        "        const msg = (err && err.message) ? err.message : 'Export failed.';\n"
        "        document.getElementById('export-status').textContent = msg;\n"
        "      });\n"
        "    }\n"
        "    document.getElementById('toggle-genes').addEventListener('change', applyVisibility);\n"
        "    document.getElementById('toggle-nt').addEventListener('change', applyVisibility);\n"
        "    document.getElementById('toggle-pos').addEventListener('change', applyVisibility);\n"
        "    document.getElementById('hist-sublibrary-filter').addEventListener('change', applySublibraryFilter);\n"
        "    document.getElementById('export-png').addEventListener('click', exportPublicationPng);\n"
        "    populateSublibraryFilter();\n"
        "    applySublibraryFilter();\n"
        "    applyVisibility();\n"
        "  </script>\n"
        "</body>\n"
        "</html>\n"
    )

    output_path = Path(output_html)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    if write_gzip_sidecar:
        gzip_path = output_path.with_suffix(output_path.suffix + ".gz")
        payload = html.encode("utf-8")
        with gzip_path.open("wb") as raw_fh:
            with gzip.GzipFile(fileobj=raw_fh, mode="wb", compresslevel=9, mtime=0) as gz_fh:
                gz_fh.write(payload)
    return output_path
