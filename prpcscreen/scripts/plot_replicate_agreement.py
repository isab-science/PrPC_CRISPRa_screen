from __future__ import annotations

import os
import sys
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from prpcscreen.visualization.replicate_scatter_plots import replicate_diagnostics
from prpcscreen.visualization.plotly_exports import write_plotly_interactive_html

DEBUG_ENV_DEFAULT = os.environ.get("PRPCSCREEN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}


def debug_log(message: str, enabled: bool) -> None:
    """Emit debug lines for replicate agreement plotting."""
    if enabled:
        print(f"[06_plot_replicate_agreement] {message}", file=sys.stderr)


def _candidate_mask(df: pd.DataFrame) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    nt = df.get("Is_NT_ctrl", pd.Series(False, index=df.index)).astype(bool)
    pos = df.get("Is_pos_ctrl", pd.Series(False, index=df.index)).astype(bool)
    is_control = nt | pos

    if "Gene_symbol" in df.columns:
        gene = df["Gene_symbol"].astype(str).str.strip()
        gene_ok = df["Gene_symbol"].notna() & gene.ne("")
        if "Entrez_ID" in df.columns:
            entrez_ok = pd.to_numeric(df["Entrez_ID"], errors="coerce").notna()
            mask &= is_control | gene_ok | entrez_ok
        else:
            mask &= is_control | gene_ok
    elif "Entrez_ID" in df.columns:
        mask &= is_control | pd.to_numeric(df["Entrez_ID"], errors="coerce").notna()
    return mask


def _write_interactive_replicate_html(df: pd.DataFrame, stem: str, output_png: str) -> Path:
    keep = _candidate_mask(df)
    x = pd.to_numeric(df.get(f"{stem}_rep1"), errors="coerce")
    y = pd.to_numeric(df.get(f"{stem}_rep2"), errors="coerce")
    valid = keep & x.notna() & y.notna()

    if int(valid.sum()) == 0:
        traces = [
            {
                "x": [0.0],
                "y": [0.0],
                "mode": "text",
                "type": "scatter",
                "text": [f"No valid {stem} replicate values found"],
                "hoverinfo": "skip",
            }
        ]
        layout = {
            "paper_bgcolor": "#f6f6f6",
            "plot_bgcolor": "#ffffff",
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
        }
    else:
        is_nt = df.get("Is_NT_ctrl", pd.Series(False, index=df.index)).astype(bool)
        is_pos = df.get("Is_pos_ctrl", pd.Series(False, index=df.index)).astype(bool)
        is_gene = (~is_nt) & (~is_pos)

        xv = x[valid].astype(float)
        yv = y[valid].astype(float)
        mv = (xv + yv) / 2.0
        dv = yv - xv
        av = mv.abs()
        adv = dv.abs()

        gene_vals = valid & is_gene
        nt_vals = valid & is_nt
        pos_vals = valid & is_pos

        x_all = pd.concat([xv, yv])
        lo = float(np.nanpercentile(x_all, 0.5))
        hi = float(np.nanpercentile(x_all, 99.5))
        if (not np.isfinite(lo)) or (not np.isfinite(hi)) or hi <= lo:
            lo, hi = float(np.nanmin(x_all)), float(np.nanmax(x_all))
        if (not np.isfinite(lo)) or (not np.isfinite(hi)) or hi <= lo:
            lo, hi = -1.0, 1.0
        pad = (hi - lo) * 0.05
        x_lo, x_hi = lo - pad, hi + pad

        bias = float(np.nanmean(dv))
        sd = float(np.nanstd(dv))

        trend_x: list[float] = []
        trend_y: list[float] = []
        if av.size >= 100:
            bins = pd.qcut(av, q=min(20, max(5, int(np.sqrt(av.size)))), duplicates="drop")
            trend = (
                pd.DataFrame({"abs_mean": av, "abs_delta": adv, "bin": bins})
                .groupby("bin", observed=True)
                .median(numeric_only=True)
            )
            trend_x = trend["abs_mean"].astype(float).tolist()
            trend_y = trend["abs_delta"].astype(float).tolist()

        corr_x: list[float] = []
        corr_y: list[float] = []
        if av.size >= 60:
            bin_codes = pd.qcut(av, q=min(10, max(5, int(np.sqrt(av.size) / 2))), duplicates="drop")
            tmp = pd.DataFrame({"x": xv, "y": yv, "abs_mean": av, "bin": bin_codes})
            rows: list[tuple[float, float]] = []
            for _, sub in tmp.groupby("bin", observed=True):
                if len(sub) < 30:
                    continue
                sx = float(sub["x"].std())
                sy = float(sub["y"].std())
                if np.isclose(sx, 0) or np.isclose(sy, 0):
                    continue
                rows.append((float(sub["abs_mean"].median()), float(sub["x"].corr(sub["y"]))))
            if rows:
                rows = sorted(rows, key=lambda r: r[0])
                corr_x = [r[0] for r in rows]
                corr_y = [r[1] for r in rows]

        def _subset(series: pd.Series, mask: pd.Series) -> list[float]:
            return series[mask[series.index]].astype(float).tolist()

        traces = [
            {
                "x": _subset(x, gene_vals),
                "y": _subset(y, gene_vals),
                "mode": "markers",
                "type": "scatter",
                "name": "Genes",
                "marker": {"color": "#7a7a7a", "size": 5, "opacity": 0.25},
                "xaxis": "x",
                "yaxis": "y",
                "hovertemplate": "Rep1: %{x:.4f}<br>Rep2: %{y:.4f}<extra>Genes</extra>",
            },
            {
                "x": _subset(x, nt_vals),
                "y": _subset(y, nt_vals),
                "mode": "markers",
                "type": "scatter",
                "name": "Non-targeting",
                "marker": {"color": "#d62728", "size": 6, "opacity": 0.9},
                "xaxis": "x",
                "yaxis": "y",
                "hovertemplate": "Rep1: %{x:.4f}<br>Rep2: %{y:.4f}<extra>Non-targeting</extra>",
            },
            {
                "x": _subset(x, pos_vals),
                "y": _subset(y, pos_vals),
                "mode": "markers",
                "type": "scatter",
                "name": "Positive",
                "marker": {"color": "#1f77b4", "size": 6, "opacity": 0.9},
                "xaxis": "x",
                "yaxis": "y",
                "hovertemplate": "Rep1: %{x:.4f}<br>Rep2: %{y:.4f}<extra>Positive</extra>",
            },
            {
                "x": [x_lo, x_hi],
                "y": [x_lo, x_hi],
                "mode": "lines",
                "type": "scatter",
                "name": "Identity",
                "line": {"color": "#111111", "width": 1, "dash": "dash"},
                "xaxis": "x",
                "yaxis": "y",
                "showlegend": False,
                "hoverinfo": "skip",
            },
            {
                "x": _subset(mv, gene_vals),
                "y": _subset(dv, gene_vals),
                "mode": "markers",
                "type": "scatter",
                "name": "Bland-Altman genes",
                "marker": {"color": "#7a7a7a", "size": 5, "opacity": 0.25},
                "xaxis": "x2",
                "yaxis": "y2",
                "showlegend": False,
                "hovertemplate": "Mean: %{x:.4f}<br>Delta: %{y:.4f}<extra>Genes</extra>",
            },
            {
                "x": _subset(mv, nt_vals),
                "y": _subset(dv, nt_vals),
                "mode": "markers",
                "type": "scatter",
                "name": "Bland-Altman Non-targeting",
                "marker": {"color": "#d62728", "size": 6, "opacity": 0.9},
                "xaxis": "x2",
                "yaxis": "y2",
                "showlegend": False,
                "hovertemplate": "Mean: %{x:.4f}<br>Delta: %{y:.4f}<extra>Non-targeting</extra>",
            },
            {
                "x": _subset(mv, pos_vals),
                "y": _subset(dv, pos_vals),
                "mode": "markers",
                "type": "scatter",
                "name": "Bland-Altman Positive",
                "marker": {"color": "#1f77b4", "size": 6, "opacity": 0.9},
                "xaxis": "x2",
                "yaxis": "y2",
                "showlegend": False,
                "hovertemplate": "Mean: %{x:.4f}<br>Delta: %{y:.4f}<extra>Positive</extra>",
            },
            {
                "x": [float(mv.min()), float(mv.max())],
                "y": [bias, bias],
                "mode": "lines",
                "type": "scatter",
                "line": {"color": "#111111", "width": 1},
                "xaxis": "x2",
                "yaxis": "y2",
                "showlegend": False,
                "hoverinfo": "skip",
            },
            {
                "x": [float(mv.min()), float(mv.max())],
                "y": [bias + 1.96 * sd, bias + 1.96 * sd],
                "mode": "lines",
                "type": "scatter",
                "line": {"color": "#444444", "width": 1, "dash": "dot"},
                "xaxis": "x2",
                "yaxis": "y2",
                "showlegend": False,
                "hoverinfo": "skip",
            },
            {
                "x": [float(mv.min()), float(mv.max())],
                "y": [bias - 1.96 * sd, bias - 1.96 * sd],
                "mode": "lines",
                "type": "scatter",
                "line": {"color": "#444444", "width": 1, "dash": "dot"},
                "xaxis": "x2",
                "yaxis": "y2",
                "showlegend": False,
                "hoverinfo": "skip",
            },
            {
                "x": _subset(av, gene_vals),
                "y": _subset(adv, gene_vals),
                "mode": "markers",
                "type": "scatter",
                "marker": {"color": "#7a7a7a", "size": 5, "opacity": 0.2},
                "xaxis": "x3",
                "yaxis": "y3",
                "showlegend": False,
                "hovertemplate": "|Mean|: %{x:.4f}<br>|Delta|: %{y:.4f}<extra>Genes</extra>",
            },
            {
                "x": _subset(av, nt_vals),
                "y": _subset(adv, nt_vals),
                "mode": "markers",
                "type": "scatter",
                "marker": {"color": "#d62728", "size": 6, "opacity": 0.8},
                "xaxis": "x3",
                "yaxis": "y3",
                "showlegend": False,
                "hovertemplate": "|Mean|: %{x:.4f}<br>|Delta|: %{y:.4f}<extra>Non-targeting</extra>",
            },
            {
                "x": _subset(av, pos_vals),
                "y": _subset(adv, pos_vals),
                "mode": "markers",
                "type": "scatter",
                "marker": {"color": "#1f77b4", "size": 6, "opacity": 0.8},
                "xaxis": "x3",
                "yaxis": "y3",
                "showlegend": False,
                "hovertemplate": "|Mean|: %{x:.4f}<br>|Delta|: %{y:.4f}<extra>Positive</extra>",
            },
            {
                "x": trend_x,
                "y": trend_y,
                "mode": "lines",
                "type": "scatter",
                "line": {"color": "#111111", "width": 1.5},
                "xaxis": "x3",
                "yaxis": "y3",
                "showlegend": False,
                "hovertemplate": "Median |Mean|: %{x:.4f}<br>Median |Delta|: %{y:.4f}<extra>Trend</extra>",
            },
            {
                "x": corr_x,
                "y": corr_y,
                "mode": "lines+markers",
                "type": "scatter",
                "line": {"color": "#2f4b7c", "width": 1.5},
                "marker": {"color": "#2f4b7c", "size": 7},
                "xaxis": "x4",
                "yaxis": "y4",
                "showlegend": False,
                "hovertemplate": "|Mean|: %{x:.4f}<br>corr: %{y:.4f}<extra>Binned corr</extra>",
            },
            {
                "x": [float(av.min()), float(av.max())],
                "y": [0.0, 0.0],
                "mode": "lines",
                "type": "scatter",
                "line": {"color": "#777777", "width": 1, "dash": "dot"},
                "xaxis": "x4",
                "yaxis": "y4",
                "showlegend": False,
                "hoverinfo": "skip",
            },
        ]

        layout = {
            "paper_bgcolor": "#f6f6f6",
            "plot_bgcolor": "#ffffff",
            "margin": {"l": 70, "r": 30, "t": 70, "b": 60},
            "showlegend": True,
            "legend": {"orientation": "h", "x": 0.0, "y": 1.12},
            "xaxis": {"domain": [0.0, 0.45], "title": f"{stem}_rep1", "range": [x_lo, x_hi]},
            "yaxis": {"domain": [0.56, 1.0], "title": f"{stem}_rep2", "range": [x_lo, x_hi]},
            "xaxis2": {"domain": [0.55, 1.0], "title": "Mean effect"},
            "yaxis2": {"domain": [0.56, 1.0], "title": "Rep2 - Rep1"},
            "xaxis3": {"domain": [0.0, 0.45], "title": "|Mean effect|"},
            "yaxis3": {"domain": [0.0, 0.44], "title": "|Rep2 - Rep1|"},
            "xaxis4": {"domain": [0.55, 1.0], "title": "|Mean effect|"},
            "yaxis4": {"domain": [0.0, 0.44], "title": "corr(rep1, rep2)", "range": [-1.0, 1.0]},
            "annotations": [
                {"text": "Replicate Scatter", "xref": "paper", "yref": "paper", "x": 0.225, "y": 1.04, "showarrow": False},
                {"text": "Bland-Altman", "xref": "paper", "yref": "paper", "x": 0.775, "y": 1.04, "showarrow": False},
                {"text": "Error vs Effect Magnitude", "xref": "paper", "yref": "paper", "x": 0.225, "y": 0.48, "showarrow": False},
                {"text": "Binned Correlation vs |Effect|", "xref": "paper", "yref": "paper", "x": 0.775, "y": 0.48, "showarrow": False},
            ],
        }

    interpretation_legend_html = (
        "    <details style=\"margin: 8px 0 12px 0; border: 1px solid #d0d7de; border-radius: 8px; background: #ffffff;\">\n"
        "      <summary style=\"cursor: pointer; padding: 10px 12px; font-weight: 700;\">Interpretation Legend: Replicate Diagnostics</summary>\n"
        "      <div style=\"padding: 10px 12px 12px 12px; font-size: 13px; line-height: 1.45; color: #1f2937;\">\n"
        "        <p style=\"margin: 0 0 8px 0;\"><strong>What this figure measures:</strong> replicate reproducibility across the same perturbations. Each point is one perturbation (gene or control).</p>\n"
        "        <p style=\"margin: 0 0 8px 0;\"><strong>Panel 1 - Replicate Scatter:</strong> x = rep1, y = rep2, with identity line <code>y = x</code>. Points near the line indicate agreement. Pearson <code>r</code> summarizes linear agreement, while outliers show discordant wells/genes. Controls (NT/positive) are highlighted for sanity-checking assay behavior.</p>\n"
        "        <p style=\"margin: 0 0 8px 0;\"><strong>Panel 2 - Bland-Altman:</strong> x = mean effect <code>(rep1 + rep2)/2</code>, y = difference <code>(rep2 - rep1)</code>. Center line is bias <code>mean(rep2 - rep1)</code>. Dotted lines are limits of agreement: <code>bias +/- 1.96 * SD(rep2 - rep1)</code>. This tests systematic shift and spread of replicate disagreement.</p>\n"
        "        <p style=\"margin: 0 0 8px 0;\"><strong>Panel 3 - Error vs Effect Magnitude:</strong> x = <code>|mean effect|</code>, y = <code>|rep2 - rep1|</code>, with a median trend. This quantifies heteroscedasticity (whether disagreement grows at stronger effects).</p>\n"
        "        <p style=\"margin: 0 0 8px 0;\"><strong>Panel 4 - Binned Correlation vs |Effect|:</strong> points are correlation estimates <code>corr(rep1, rep2)</code> computed within bins of <code>|mean effect|</code>. This shows whether reproducibility is stable across weak-to-strong signal regimes.</p>\n"
        "        <p style=\"margin: 0 0 8px 0;\"><strong>Practical quality targets (rules of thumb):</strong> Pearson <code>r &gt;= 0.70</code> acceptable, <code>&gt;= 0.80</code> good, <code>&gt;= 0.90</code> excellent (assay-dependent). Bland-Altman bias near 0 and most points within LoA is desirable. Error-vs-effect trend should be flat or gently rising; steep growth suggests noisy extremes. Binned correlations should stay positive and preferably high across bins.</p>\n"
        "        <p style=\"margin: 0 0 8px 0;\"><strong>What indicates concern:</strong> broad scatter away from identity line, strong non-zero bias, fan-shaped Bland-Altman structure, rapidly increasing <code>|delta|</code> with effect size, or bin correlations collapsing toward 0/negative in high-effect bins.</p>\n"
        "        <p style=\"margin: 0;\"><strong>Interpretation note:</strong> thresholds depend on assay dynamic range, replicate count, and control design. Always interpret with control separation (NT vs positive), hit concordance, and known biology.</p>\n"
        "      </div>\n"
        "    </details>\n"
    )

    output_html = Path(output_png).with_name(Path(output_png).stem + "_interactive.html")
    file_base = f"replicate_agreement_{stem.lower()}"
    return write_plotly_interactive_html(
        output_html=output_html,
        traces=traces,
        layout=layout,
        title=f"Replicate Diagnostics ({stem})",
        filename_base=file_base,
        extra_controls_html=interpretation_legend_html,
    )


def run_concordance_cli() -> None:
    # Parse analyzed input, output figure path, and metric stem selector.
    parser = argparse.ArgumentParser(description="Generate replicate diagnostics plots.")
    parser.add_argument("input_csv")
    parser.add_argument("output_png")
    parser.add_argument("--stem", default="Raw")
    parser.add_argument("--debug", action="store_true", help="Enable verbose debug logging.")
    args = parser.parse_args()
    debug_enabled = DEBUG_ENV_DEFAULT or args.debug

    # Load input and render concordance between replicate columns.
    df = pd.read_csv(args.input_csv)
    debug_log(f"Loaded analyzed data: {args.input_csv} ({len(df)} rows)", debug_enabled)
    debug_log(f"Using replicate stem: {args.stem}", debug_enabled)
    fig, _ = replicate_diagnostics(df, stem=args.stem)
    Path(args.output_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    debug_log(f"Wrote replicate agreement figure: {args.output_png}", debug_enabled)

    html_path = _write_interactive_replicate_html(df, args.stem, args.output_png)
    debug_log(f"Wrote interactive replicate diagnostics HTML: {html_path}", debug_enabled)


if __name__ == "__main__":
    run_concordance_cli()
