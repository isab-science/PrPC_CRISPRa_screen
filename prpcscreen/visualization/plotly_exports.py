from __future__ import annotations

import gzip
import json
from pathlib import Path


def write_plotly_interactive_html(
    output_html: str | Path,
    traces: list[dict],
    layout: dict,
    title: str,
    filename_base: str,
    *,
    plot_div_id: str = "plot",
    extra_controls_html: str = "",
    extra_script: str = "",
    write_gzip_sidecar: bool = True,
) -> Path:
    traces_json = json.dumps(traces, separators=(",", ":"))
    layout_json = json.dumps(layout, separators=(",", ":"))
    title_json = json.dumps(title, separators=(",", ":"))
    filename_json = json.dumps(filename_base, separators=(",", ":"))
    div_id_json = json.dumps(plot_div_id, separators=(",", ":"))

    html = (
        "<!doctype html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "  <title>"
        + title
        + "</title>\n"
        "  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>\n"
        "  <style>\n"
        "    body { margin: 0; font-family: Segoe UI, Arial, sans-serif; background: #f6f6f6; color: #1f1f1f; }\n"
        "    .wrap { max-width: 1500px; margin: 0 auto; padding: 16px; }\n"
        "    h1 { font-size: 18px; margin: 0 0 10px; font-weight: 700; }\n"
        "    .controls { display: flex; gap: 14px; flex-wrap: wrap; align-items: center; margin-bottom: 10px; font-size: 14px; }\n"
        "    .controls label { display: inline-flex; align-items: center; gap: 6px; user-select: none; }\n"
        "    .controls select { border: 1px solid #d0d7de; border-radius: 6px; padding: 6px 8px; font-size: 13px; background: #ffffff; }\n"
        "    .controls button { border: 1px solid #c9ced6; border-radius: 6px; background: #ffffff; padding: 6px 10px; cursor: pointer; font-size: 13px; }\n"
        "    .controls button:hover { background: #f4f6f8; }\n"
        "    #export-status { color: #475569; font-size: 13px; min-height: 1em; }\n"
        "    #"
        + plot_div_id
        + " { width: 100%; height: 78vh; min-height: 620px; border: 1px solid #d8d8d8; background: #ffffff; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n"
        "  <div class=\"wrap\">\n"
        "    <h1 id=\"plot-title\"></h1>\n"
        + extra_controls_html
        + "    <div class=\"controls\">\n"
        "      <label for=\"export-quality\">Publication SVG</label>\n"
        "      <select id=\"export-quality\">\n"
        "        <option value=\"low\">Low (1200x800)</option>\n"
        "        <option value=\"medium\" selected>Medium (1800x1200)</option>\n"
        "        <option value=\"high\">High (2400x1600)</option>\n"
        "      </select>\n"
        "      <button id=\"export-svg\" type=\"button\">Export SVG</button>\n"
        "      <span id=\"export-status\"></span>\n"
        "    </div>\n"
        "    <div id=\""
        + plot_div_id
        + "\"></div>\n"
        "  </div>\n"
        "  <script>\n"
        "    const PLOT_DIV_ID = "
        + div_id_json
        + ";\n"
        "    const FIGURE_TITLE = "
        + title_json
        + ";\n"
        "    const FILE_BASE = "
        + filename_json
        + ";\n"
        "    const traces = "
        + traces_json
        + ";\n"
        "    const layout = "
        + layout_json
        + ";\n"
        "    document.getElementById('plot-title').textContent = FIGURE_TITLE;\n"
        "    Plotly.newPlot(PLOT_DIV_ID, traces, layout, {responsive: true});\n"
        "    function getExportPreset() {\n"
        "      const mode = document.getElementById('export-quality').value;\n"
        "      if (mode === 'low') return { label: 'low', width: 1200, height: 800 };\n"
        "      if (mode === 'high') return { label: 'high', width: 2400, height: 1600 };\n"
        "      return { label: 'medium', width: 1800, height: 1200 };\n"
        "    }\n"
        "    function exportPublicationSvg() {\n"
        "      const preset = getExportPreset();\n"
        "      const filename = `${FILE_BASE}_${preset.label}`;\n"
        "      document.getElementById('export-status').textContent = `Exporting ${preset.label} SVG...`;\n"
        "      Plotly.downloadImage(PLOT_DIV_ID, {\n"
        "        format: 'svg',\n"
        "        filename,\n"
        "        width: preset.width,\n"
        "        height: preset.height,\n"
        "        scale: 1,\n"
        "      }).then(() => {\n"
        "        document.getElementById('export-status').textContent = `Saved ${filename}.svg (${preset.width}x${preset.height}).`;\n"
        "      }).catch((err) => {\n"
        "        const msg = (err && err.message) ? err.message : 'Export failed.';\n"
        "        document.getElementById('export-status').textContent = msg;\n"
        "      });\n"
        "    }\n"
        "    document.getElementById('export-svg').addEventListener('click', exportPublicationSvg);\n"
        + extra_script
        + "\n"
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

