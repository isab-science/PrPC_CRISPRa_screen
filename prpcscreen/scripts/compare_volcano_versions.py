from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from prpcscreen.visualization.volcano_and_flashlight_plots import _prepare_volcano_dataframe, volcano_plot


DEBUG_ENV_DEFAULT = os.environ.get("PRPCSCREEN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}


def debug_log(message: str, enabled: bool) -> None:
    if enabled:
        print(f"[compare_volcano_versions] {message}", file=sys.stderr)


def _original_candidate_mask(df: pd.DataFrame) -> pd.Series:
    """Candidate filtering from historical baseline commit c856def."""
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


def _prepare_original_volcano_dataframe(df: pd.DataFrame, x_col: str, p_col: str) -> tuple[pd.DataFrame, pd.Series]:
    keep = _original_candidate_mask(df)
    x = pd.to_numeric(df.loc[keep, x_col], errors="coerce")
    p = pd.to_numeric(df.loc[keep, p_col], errors="coerce")
    y = -np.log10(np.clip(p, 1e-300, 1.0))

    plot_df = df.loc[keep].copy()
    plot_df["_x"] = x
    plot_df["_y"] = y
    plot_df = plot_df[np.isfinite(plot_df["_x"]) & np.isfinite(plot_df["_y"])].copy()
    return plot_df, keep


def _save_original_volcano_png(df: pd.DataFrame, output_png: str | Path, x_col: str, p_col: str) -> None:
    plot_df, keep = _prepare_original_volcano_dataframe(df, x_col=x_col, p_col=p_col)
    pos_mask = plot_df.get("Is_pos_ctrl", pd.Series(False, index=plot_df.index)).astype(bool)
    nt_mask = plot_df.get("Is_NT_ctrl", pd.Series(False, index=plot_df.index)).astype(bool) & ~pos_mask
    gene_mask = ~(pos_mask | nt_mask)

    fig, ax = plt.subplots(figsize=(6.3, 5.1))
    log2fc_cutoff = 1.0
    p_cutoff = 0.05
    y_cutoff = -np.log10(p_cutoff)
    x_min, x_max = -5.5, 3.0
    y_min, y_max = -0.05, 2.25
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.fill_between([x_min, -log2fc_cutoff], y_cutoff, y_max, color="#EBEBEB", zorder=0)
    ax.fill_between([log2fc_cutoff, x_max], y_cutoff, y_max, color="#EBEBEB", zorder=0)
    ax.axhline(0, color="grey", ls=":", lw=0.6, zorder=1)
    ax.axvline(0, color="grey", ls=":", lw=0.6, zorder=1)
    ax.axhline(y_cutoff, color="#E8E8E8", lw=0.6, zorder=1)
    ax.axvline(-log2fc_cutoff, color="#E8E8E8", lw=0.6, zorder=1)
    ax.axvline(log2fc_cutoff, color="#E8E8E8", lw=0.6, zorder=1)
    ax.scatter(plot_df.loc[gene_mask, "_x"], plot_df.loc[gene_mask, "_y"], s=6, alpha=0.30, c="#000000", label="Genes")
    ax.scatter(plot_df.loc[pos_mask, "_x"], plot_df.loc[pos_mask, "_y"], s=6, alpha=0.50, c="#de2d26", label="Positive controls")
    ax.scatter(plot_df.loc[nt_mask, "_x"], plot_df.loc[nt_mask, "_y"], s=6, alpha=0.50, c="#3182bd", label="Negative controls")
    ax.set_xlabel(x_col)
    ax.set_ylabel(f"-log10({p_col})")
    removed = int((~keep).sum())
    ax.set_title("Volcano plot (baseline)")
    if removed:
        ax.text(0.99, 0.01, f"Filtered rows: {removed}", transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#555555")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, loc="center left", bbox_to_anchor=(1.02, 0.5))
    out = Path(output_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)


def _well_coordinate_384(series: pd.Series) -> pd.Series:
    nums = pd.to_numeric(series, errors="coerce")
    out: list[str] = []
    for value in nums.tolist():
        if pd.isna(value):
            out.append("NA")
            continue
        well = int(value)
        if 1 <= well <= 384:
            row_idx = (well - 1) // 24
            col_idx = (well - 1) % 24 + 1
            out.append(f"{chr(ord('A') + row_idx)}{col_idx:02d}")
        else:
            out.append(str(well))
    return pd.Series(out, index=series.index, dtype=object)


def _row_key(df: pd.DataFrame) -> pd.Series:
    if {"Plate_number_384", "Well_number_384"}.issubset(df.columns):
        plate = df["Plate_number_384"].astype(str).str.strip()
        well = pd.to_numeric(df["Well_number_384"], errors="coerce")
        well_text = well.map(lambda v: str(int(v)) if pd.notna(v) else "")
        has_pw = plate.ne("") & well_text.ne("")
        key = pd.Series(df.index.astype(str), index=df.index, dtype=object)
        key.loc[has_pw] = plate.loc[has_pw] + "::" + well_text.loc[has_pw]
        return key
    return pd.Series(df.index.astype(str), index=df.index, dtype=object)


def _control_and_gene_label(df: pd.DataFrame) -> pd.Series:
    symbol = df.get("Gene_symbol", pd.Series("", index=df.index)).fillna("").astype(str).str.strip()
    pos = df.get("Is_pos_ctrl", pd.Series(False, index=df.index)).astype(bool)
    nt = df.get("Is_NT_ctrl", pd.Series(False, index=df.index)).astype(bool) & ~pos
    symbol = symbol.mask(symbol.eq("") & pos, "PRNP")
    symbol = symbol.mask(symbol.eq("") & nt, "NT_control")
    symbol = symbol.mask(symbol.eq(""), "NA")
    return symbol


def _projection_for_compare(df: pd.DataFrame, mode: str, x_col: str, p_col: str) -> pd.DataFrame:
    if mode == "original":
        plot_df, _ = _prepare_original_volcano_dataframe(df, x_col=x_col, p_col=p_col)
    elif mode == "current":
        plot_df, _ = _prepare_volcano_dataframe(df, x_col=x_col, p_col=p_col)
    else:
        raise ValueError(f"Unknown mode: {mode}")

    keys = _row_key(plot_df)
    plate = plot_df.get("Plate_number_384", pd.Series("NA", index=plot_df.index)).astype(str)
    well = pd.to_numeric(plot_df.get("Well_number_384", pd.Series(np.nan, index=plot_df.index)), errors="coerce")
    coord = _well_coordinate_384(plot_df.get("Well_number_384", pd.Series(np.nan, index=plot_df.index)))
    symbol = _control_and_gene_label(plot_df)
    plasmid = plot_df.get("Plasmid_ID", pd.Series("", index=plot_df.index)).fillna("").astype(str)
    is_pos = plot_df.get("Is_pos_ctrl", pd.Series(False, index=plot_df.index)).astype(bool)
    is_nt = plot_df.get("Is_NT_ctrl", pd.Series(False, index=plot_df.index)).astype(bool) & ~is_pos

    return pd.DataFrame(
        {
            "row_key": keys,
            "x": plot_df["_x"].astype(float).to_numpy(),
            "y": plot_df["_y"].astype(float).to_numpy(),
            "plate": plate.to_numpy(),
            "well_number": well.to_numpy(),
            "well_coordinate": coord.to_numpy(),
            "gene_label": symbol.to_numpy(),
            "plasmid_id": plasmid.to_numpy(),
            "is_pos_ctrl": is_pos.to_numpy(),
            "is_nt_ctrl": is_nt.to_numpy(),
        }
    )


def _discordant_table(df: pd.DataFrame, x_col: str, p_col: str, tol: float) -> pd.DataFrame:
    orig = _projection_for_compare(df, mode="original", x_col=x_col, p_col=p_col).add_prefix("orig_")
    curr = _projection_for_compare(df, mode="current", x_col=x_col, p_col=p_col).add_prefix("curr_")
    merged = orig.merge(curr, left_on="orig_row_key", right_on="curr_row_key", how="outer", indicator=True)
    merged["row_key"] = merged["orig_row_key"].combine_first(merged["curr_row_key"])

    both = merged["_merge"].eq("both")
    changed = both & (
        (np.abs(merged["orig_x"] - merged["curr_x"]) > tol)
        | (np.abs(merged["orig_y"] - merged["curr_y"]) > tol)
    )
    only_original = merged["_merge"].eq("left_only")
    only_current = merged["_merge"].eq("right_only")
    discordant = merged.loc[changed | only_original | only_current].copy()

    discordant["discordance"] = np.select(
        [only_original.loc[discordant.index], only_current.loc[discordant.index], changed.loc[discordant.index]],
        ["only_in_original", "only_in_current", "changed_value"],
        default="unknown",
    )
    discordant["delta_x"] = discordant["curr_x"] - discordant["orig_x"]
    discordant["delta_y"] = discordant["curr_y"] - discordant["orig_y"]

    for col in ["plate", "well_number", "well_coordinate", "gene_label", "plasmid_id", "is_pos_ctrl", "is_nt_ctrl"]:
        discordant[col] = discordant[f"curr_{col}"].combine_first(discordant[f"orig_{col}"])

    cols = [
        "row_key",
        "discordance",
        "plate",
        "well_number",
        "well_coordinate",
        "gene_label",
        "plasmid_id",
        "is_pos_ctrl",
        "is_nt_ctrl",
        "orig_x",
        "orig_y",
        "curr_x",
        "curr_y",
        "delta_x",
        "delta_y",
    ]
    return discordant[cols].sort_values(["discordance", "plate", "well_number"], na_position="last")


def _discordant_interactive_html(discordant: pd.DataFrame, output_html: str | Path, x_col: str, p_col: str, base_ref: str) -> None:
    only_orig = discordant[discordant["discordance"].eq("only_in_original")].copy()
    only_curr = discordant[discordant["discordance"].eq("only_in_current")].copy()
    changed = discordant[discordant["discordance"].eq("changed_value")].copy()

    def _hover_lines(df_sub: pd.DataFrame, side: str) -> list[str]:
        out: list[str] = []
        for _, row in df_sub.iterrows():
            x_val = row["orig_x"] if side == "orig" else row["curr_x"]
            y_val = row["orig_y"] if side == "orig" else row["curr_y"]
            out.append(
                "<br>".join(
                    [
                        f"Row: {row['row_key']}",
                        f"Gene: {row['gene_label']}",
                        f"Plate: {row['plate']}",
                        f"Coordinate: {row['well_coordinate']}",
                        f"Plasmid: {row['plasmid_id']}",
                        f"Status: {row['discordance']}",
                        f"{x_col}: {x_val:.5f}" if pd.notna(x_val) else f"{x_col}: NA",
                        f"-log10({p_col}): {y_val:.5f}" if pd.notna(y_val) else f"-log10({p_col}): NA",
                        f"Delta x (current-original): {row['delta_x']:.5f}" if pd.notna(row["delta_x"]) else "Delta x (current-original): NA",
                        f"Delta y (current-original): {row['delta_y']:.5f}" if pd.notna(row["delta_y"]) else "Delta y (current-original): NA",
                    ]
                )
            )
        return out

    line_x: list[float | None] = []
    line_y: list[float | None] = []
    for _, row in changed.iterrows():
        line_x.extend([float(row["orig_x"]), float(row["curr_x"]), None])
        line_y.extend([float(row["orig_y"]), float(row["curr_y"]), None])

    traces = [
        {
            "x": line_x,
            "y": line_y,
            "mode": "lines",
            "name": "Shift (original -> current)",
            "line": {"color": "#9ca3af", "width": 1},
            "hoverinfo": "skip",
            "showlegend": True,
        },
        {
            "x": only_orig["orig_x"].astype(float).tolist(),
            "y": only_orig["orig_y"].astype(float).tolist(),
            "text": _hover_lines(only_orig, side="orig"),
            "mode": "markers",
            "name": "Only in original",
            "marker": {"size": 6, "color": "#1d4ed8", "opacity": 0.85, "symbol": "square"},
            "hovertemplate": "%{text}<extra></extra>",
        },
        {
            "x": only_curr["curr_x"].astype(float).tolist(),
            "y": only_curr["curr_y"].astype(float).tolist(),
            "text": _hover_lines(only_curr, side="curr"),
            "mode": "markers",
            "name": "Only in current",
            "marker": {"size": 6, "color": "#dc2626", "opacity": 0.85, "symbol": "circle"},
            "hovertemplate": "%{text}<extra></extra>",
        },
        {
            "x": changed["orig_x"].astype(float).tolist(),
            "y": changed["orig_y"].astype(float).tolist(),
            "text": _hover_lines(changed, side="orig"),
            "mode": "markers",
            "name": "Changed (original position)",
            "marker": {"size": 6, "color": "#f59e0b", "opacity": 0.9, "symbol": "circle-open"},
            "hovertemplate": "%{text}<extra></extra>",
        },
        {
            "x": changed["curr_x"].astype(float).tolist(),
            "y": changed["curr_y"].astype(float).tolist(),
            "text": _hover_lines(changed, side="curr"),
            "mode": "markers",
            "name": "Changed (current position)",
            "marker": {"size": 6, "color": "#16a34a", "opacity": 0.9, "symbol": "circle"},
            "hovertemplate": "%{text}<extra></extra>",
        },
    ]

    all_x = pd.concat(
        [
            only_orig["orig_x"],
            only_curr["curr_x"],
            changed["orig_x"],
            changed["curr_x"],
        ],
        ignore_index=True,
    ).dropna()
    all_y = pd.concat(
        [
            only_orig["orig_y"],
            only_curr["curr_y"],
            changed["orig_y"],
            changed["curr_y"],
        ],
        ignore_index=True,
    ).dropna()
    if all_x.empty:
        x_range = [-2.0, 2.0]
    else:
        x_min, x_max = float(all_x.min()), float(all_x.max())
        pad = max(0.1, 0.05 * (x_max - x_min if x_max > x_min else 1.0))
        x_range = [x_min - pad, x_max + pad]
    if all_y.empty:
        y_range = [0.0, 2.5]
    else:
        y_min, y_max = float(all_y.min()), float(all_y.max())
        pad = max(0.05, 0.05 * (y_max - y_min if y_max > y_min else 1.0))
        y_range = [max(-0.05, y_min - pad), y_max + pad]

    layout = {
        "title": {
            "text": (
                f"Discordant volcano datapoints only "
                f"(original={base_ref}, current=HEAD; n={len(discordant)})"
            )
        },
        "paper_bgcolor": "#f6f6f6",
        "plot_bgcolor": "#ffffff",
        "margin": {"l": 70, "r": 40, "t": 70, "b": 60},
        "xaxis": {"title": x_col, "range": x_range, "zeroline": False},
        "yaxis": {"title": f"-log10({p_col})", "range": y_range, "zeroline": False},
        "showlegend": True,
        "legend": {"orientation": "v", "x": 1.02, "y": 0.5, "xanchor": "left", "yanchor": "middle"},
    }
    if len(discordant) == 0:
        layout["annotations"] = [
            {
                "xref": "paper",
                "yref": "paper",
                "x": 0.5,
                "y": 0.5,
                "text": "No discordant datapoints found.",
                "showarrow": False,
                "font": {"size": 18, "color": "#334155"},
            }
        ]

    html = (
        "<!doctype html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "  <title>Discordant Volcano Points</title>\n"
        "  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>\n"
        "  <style>\n"
        "    body { margin: 0; font-family: Segoe UI, Arial, sans-serif; background: #f6f6f6; color: #1f1f1f; }\n"
        "    .wrap { max-width: 1400px; margin: 0 auto; padding: 16px; }\n"
        "    #volcano { width: 100%; height: 80vh; min-height: 640px; border: 1px solid #d8d8d8; background: #ffffff; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n"
        "  <div class=\"wrap\">\n"
        "    <div id=\"volcano\"></div>\n"
        "  </div>\n"
        "  <script>\n"
        "    const traces = "
        + json.dumps(traces)
        + ";\n"
        "    const layout = "
        + json.dumps(layout)
        + ";\n"
        "    Plotly.newPlot('volcano', traces, layout, {responsive: true});\n"
        "  </script>\n"
        "</body>\n"
        "</html>\n"
    )

    out = Path(output_html)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")


def run_compare_cli() -> None:
    parser = argparse.ArgumentParser(
        description="Generate baseline/current volcano plots and an interactive volcano with only discordant datapoints."
    )
    parser.add_argument("input_csv")
    parser.add_argument("original_volcano_png")
    parser.add_argument("current_volcano_png")
    parser.add_argument("discordant_volcano_html")
    parser.add_argument("--discordant_csv", default=None, help="Optional CSV export of discordant datapoints.")
    parser.add_argument("--x_col", default="Mean_log2")
    parser.add_argument("--p_col", default="p_value_log2")
    parser.add_argument("--base_ref", default="c856def", help="Label used in output metadata/title.")
    parser.add_argument("--tol", type=float, default=1e-12, help="Absolute tolerance for x/y equality checks.")
    parser.add_argument("--debug", action="store_true", help="Enable verbose debug logging.")
    args = parser.parse_args()
    debug_enabled = DEBUG_ENV_DEFAULT or args.debug

    df = pd.read_csv(args.input_csv)
    debug_log(f"Loaded analyzed data: {args.input_csv} ({len(df)} rows)", debug_enabled)

    _save_original_volcano_png(df, args.original_volcano_png, x_col=args.x_col, p_col=args.p_col)
    debug_log(f"Wrote original volcano PNG: {args.original_volcano_png}", debug_enabled)

    fig_curr, _ = volcano_plot(df, x_col=args.x_col, p_col=args.p_col)
    out_curr = Path(args.current_volcano_png)
    out_curr.parent.mkdir(parents=True, exist_ok=True)
    fig_curr.savefig(out_curr, dpi=600, bbox_inches="tight")
    plt.close(fig_curr)
    debug_log(f"Wrote current volcano PNG: {args.current_volcano_png}", debug_enabled)

    discordant = _discordant_table(df, x_col=args.x_col, p_col=args.p_col, tol=args.tol)
    debug_log(
        "Discordant counts: "
        f"total={len(discordant)}, "
        f"only_original={(discordant['discordance'] == 'only_in_original').sum()}, "
        f"only_current={(discordant['discordance'] == 'only_in_current').sum()}, "
        f"changed_value={(discordant['discordance'] == 'changed_value').sum()}",
        debug_enabled,
    )

    if args.discordant_csv:
        out_csv = Path(args.discordant_csv)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        discordant.to_csv(out_csv, index=False)
        debug_log(f"Wrote discordant CSV: {out_csv}", debug_enabled)

    _discordant_interactive_html(
        discordant=discordant,
        output_html=args.discordant_volcano_html,
        x_col=args.x_col,
        p_col=args.p_col,
        base_ref=args.base_ref,
    )
    debug_log(f"Wrote discordant interactive HTML: {args.discordant_volcano_html}", debug_enabled)

    print(
        "Discordant summary: "
        f"total={len(discordant)}, "
        f"only_in_original={(discordant['discordance'] == 'only_in_original').sum()}, "
        f"only_in_current={(discordant['discordance'] == 'only_in_current').sum()}, "
        f"changed_value={(discordant['discordance'] == 'changed_value').sum()}"
    )


if __name__ == "__main__":
    run_compare_cli()
