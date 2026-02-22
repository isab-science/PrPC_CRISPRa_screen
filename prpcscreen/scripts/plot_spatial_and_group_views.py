from __future__ import annotations

import os
import sys
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from prpcscreen.visualization.box_plots import beebox_plates
from prpcscreen.visualization.heatmaps import heatmap_for_plate
from prpcscreen.visualization.plotly_exports import write_plotly_interactive_html

DEBUG_ENV_DEFAULT = os.environ.get("PRPCSCREEN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}


def debug_log(message: str, enabled: bool) -> None:
    """Emit debug lines for spatial/group summary plotting."""
    if enabled:
        print(f"[plot_spatial_and_group_views] {message}", file=sys.stderr)


def _write_interactive_heatmap_html(df: pd.DataFrame, plate: str, column: str, output_png: Path) -> Path:
    sub = df[df["Plate_number_384"].astype(str) == str(plate)].copy().sort_values("Well_number_384")
    values = pd.to_numeric(sub[column], errors="coerce").to_numpy(dtype=float)
    if values.size != 384:
        raise ValueError(f"Expected 384 wells for plate {plate}, got {values.size}")
    matrix = values.reshape(16, 24)
    y_labels = [chr(ord("A") + i) for i in range(16)]
    x_labels = [str(i + 1) for i in range(24)]
    traces = [
        {
            "type": "heatmap",
            "z": matrix.tolist(),
            "x": x_labels,
            "y": y_labels,
            "colorscale": "Viridis",
            "colorbar": {"title": column},
            "hovertemplate": "Row %{y} / Col %{x}<br>Value: %{z:.4f}<extra></extra>",
        }
    ]
    layout = {
        "paper_bgcolor": "#f6f6f6",
        "plot_bgcolor": "#ffffff",
        "margin": {"l": 70, "r": 30, "t": 35, "b": 70},
        "xaxis": {"title": "Column"},
        "yaxis": {"title": "Row", "autorange": "reversed"},
    }
    output_html = output_png.with_name(output_png.stem + "_interactive.html")
    return write_plotly_interactive_html(
        output_html=output_html,
        traces=traces,
        layout=layout,
        title=f"Plate {plate} heatmap ({column})",
        filename_base=f"plate_heatmap_{column.lower()}_plate{plate}",
    )


def _write_interactive_grouped_html(df: pd.DataFrame, column: str, output_png: str) -> Path:
    plot_df = df.copy()
    if "Target_flag" in plot_df:
        flag = plot_df["Target_flag"].fillna("").astype(str).str.strip().str.lower()
        own_non_targeting = flag.isin({"own non-targeting control", "own non targeting control", "own nt control"})
        plot_df["group"] = "Gene"
        plot_df.loc[plot_df["Is_pos_ctrl"].astype(bool), "group"] = "Positive"
        plot_df.loc[plot_df["Is_NT_ctrl"].astype(bool), "group"] = "Non-targeting"
        plot_df.loc[own_non_targeting, "group"] = "Own Non-targeting"
    else:
        plot_df["group"] = "Gene"
        plot_df.loc[plot_df["Is_pos_ctrl"].astype(bool), "group"] = "Positive"
        plot_df.loc[plot_df["Is_NT_ctrl"].astype(bool), "group"] = "Non-targeting"

    order = [g for g in ["Gene", "Non-targeting", "Own Non-targeting", "Positive"] if g in plot_df["group"].unique()]
    palette = {
        "Gene": "#90caf9",
        "Non-targeting": "#3b82f6",
        "Own Non-targeting": "#0ea5e9",
        "Positive": "#ef4444",
    }
    traces: list[dict] = []
    for group in order:
        vals = pd.to_numeric(plot_df.loc[plot_df["group"] == group, column], errors="coerce").dropna().astype(float)
        traces.append(
            {
                "type": "violin",
                "name": group,
                "y": vals.tolist(),
                "box": {"visible": True},
                "meanline": {"visible": False},
                "points": False,
                "line": {"color": palette.get(group, "#1f1f1f"), "width": 1.0},
                "fillcolor": palette.get(group, "#90caf9"),
                "opacity": 0.6,
                "hovertemplate": f"{group}<br>{column}: %{{y:.4f}}<extra></extra>",
            }
        )

    layout = {
        "paper_bgcolor": "#f6f6f6",
        "plot_bgcolor": "#ffffff",
        "margin": {"l": 70, "r": 30, "t": 35, "b": 70},
        "xaxis": {"title": "Group"},
        "yaxis": {"title": column},
        "violinmode": "group",
    }
    output_html = Path(output_png).with_name(Path(output_png).stem + "_interactive.html")
    return write_plotly_interactive_html(
        output_html=output_html,
        traces=traces,
        layout=layout,
        title=f"Grouped violin/box ({column})",
        filename_base=f"grouped_boxplot_{column.lower()}",
    )


def plate_sort_key(plate: str) -> tuple[int, int | str]:
    # Keep numeric plate labels in natural order, then non-numeric labels.
    try:
        return (0, int(str(plate)))
    except ValueError:
        return (1, str(plate))


def parse_plate_selector(selector: str, available_complete_plates: list[str]) -> list[str]:
    """
    Parse a heatmap plate selector.

    Supported formats:
    - single number: "1"
    - numeric range: "1-4"
    - numeric series: "1,2,6"
    - all plates: "all"
    """
    requested = str(selector).strip().lower()
    if not requested:
        raise ValueError("Heatmap plate selector is required.")

    values: list[str]
    sorted_available = sorted(available_complete_plates, key=plate_sort_key)
    has_numeric_labels = any(re.fullmatch(r"\d+", p) for p in available_complete_plates)
    ordinal_index = {str(i + 1): plate for i, plate in enumerate(sorted_available)}
    if requested == "all":
        values = sorted_available
    else:
        m_range = re.fullmatch(r"(\d+)\s*-\s*(\d+)", requested)
        if m_range:
            start = int(m_range.group(1))
            end = int(m_range.group(2))
            if start > end:
                raise ValueError("Heatmap plate range must be ascending (for example, '1-4').")
            values = [str(v) for v in range(start, end + 1)]
        elif re.fullmatch(r"\d+(?:\s*,\s*\d+)+", requested):
            values = [v.strip() for v in requested.split(",")]
        elif re.fullmatch(r"\d+", requested):
            values = [requested]
        else:
            raise ValueError(
                "Invalid heatmap plate selector. Use one number (for example '1'), "
                "a range ('1-4'), a series ('1,2,6'), or 'all'."
            )

    selected: list[str] = []
    missing: list[str] = []
    available_set = set(available_complete_plates)
    for value in values:
        if value in available_set:
            if value not in selected:
                selected.append(value)
        elif (not has_numeric_labels) and (value in ordinal_index):
            mapped_plate = ordinal_index[value]
            if mapped_plate not in selected:
                selected.append(mapped_plate)
            print(
                f"[plot_spatial_and_group_views] INFO: Interpreting selector '{value}' as plate position -> '{mapped_plate}'.",
                file=sys.stderr,
            )
        else:
            missing.append(value)

    if missing:
        print(
            f"[plot_spatial_and_group_views] WARNING: Requested plate(s) not found as complete 384-well sets: {', '.join(missing)}",
            file=sys.stderr,
        )

    if not selected:
        available_preview = ", ".join(sorted(available_complete_plates, key=plate_sort_key)[:20])
        raise ValueError(
            "No requested plate was available as a complete 384-well set. "
            f"Requested='{selector}'. Available (first 20): {available_preview}"
        )

    return sorted(selected, key=plate_sort_key)


def build_heatmap_outputs(path: str, selected_plates: list[str]) -> list[Path]:
    # Keep the exact output path for one plate; fan out with suffixes for multi-plate runs.
    base = Path(path)
    if len(selected_plates) == 1:
        return [base]
    suffix = base.suffix or ".png"
    stem = base.stem if base.suffix else base.name
    return [base.with_name(f"{stem}_plate{plate}{suffix}") for plate in selected_plates]


def run_spatial_cli() -> None:
    # Parse analyzed input, heatmap target, combined violin/box target, and plate selector.
    parser = argparse.ArgumentParser(description="Generate heatmap and combined violin/box plots.")
    parser.add_argument("input_csv")
    parser.add_argument("heatmap_png")
    parser.add_argument("boxplot_png")
    parser.add_argument("--plate", default="1")
    parser.add_argument("--debug", action="store_true", help="Enable verbose debug logging.")
    args = parser.parse_args()
    debug_enabled = DEBUG_ENV_DEFAULT or args.debug

    # Load analyzed dataset once and create both spatial and grouped summaries.
    df = pd.read_csv(args.input_csv)
    debug_log(f"Loaded analyzed data: {args.input_csv} ({len(df)} rows)", debug_enabled)
    if "Plate_number_384" not in df.columns:
        raise ValueError("Missing required column: Plate_number_384")

    plate_series = df["Plate_number_384"].astype(str)
    counts = plate_series.value_counts()
    complete_plates = [str(p) for p, c in counts.items() if int(c) == 384]
    if not complete_plates:
        available = ", ".join(sorted(counts.index.astype(str).tolist())[:20])
        raise ValueError(f"No complete 384-well plates found. Available (first 20): {available}")

    selected_plates = parse_plate_selector(args.plate, complete_plates)
    debug_log(f"Heatmap plate selector: {args.plate}", debug_enabled)
    debug_log(f"Heatmap plates resolved: {selected_plates}", debug_enabled)

    output_paths = build_heatmap_outputs(args.heatmap_png, selected_plates)
    Path(args.heatmap_png).parent.mkdir(parents=True, exist_ok=True)
    for plate_for_heatmap, output_path in zip(selected_plates, output_paths):
        fig_h, _ = heatmap_for_plate(df, plate_for_heatmap, "Raw_rep1")
        fig_h.savefig(output_path, dpi=180, bbox_inches="tight")
        plt.close(fig_h)
        debug_log(f"Wrote heatmap figure for plate {plate_for_heatmap}: {output_path}", debug_enabled)
        heatmap_html = _write_interactive_heatmap_html(df, plate_for_heatmap, "Raw_rep1", output_path)
        debug_log(f"Wrote interactive heatmap HTML for plate {plate_for_heatmap}: {heatmap_html}", debug_enabled)

    fig_b, _ = beebox_plates(df, "Raw_rep1")
    Path(args.boxplot_png).parent.mkdir(parents=True, exist_ok=True)
    fig_b.savefig(args.boxplot_png, dpi=180, bbox_inches="tight")
    plt.close(fig_b)
    debug_log(f"Wrote violin/box figure: {args.boxplot_png}", debug_enabled)
    grouped_html = _write_interactive_grouped_html(df, "Raw_rep1", args.boxplot_png)
    debug_log(f"Wrote interactive grouped violin/box HTML: {grouped_html}", debug_enabled)


if __name__ == "__main__":
    run_spatial_cli()

