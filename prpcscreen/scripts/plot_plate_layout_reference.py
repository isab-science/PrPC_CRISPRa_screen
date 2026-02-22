from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd


def _well_to_row_col(well_number: int) -> tuple[int, int]:
    idx = int(well_number) - 1
    return idx // 24, (idx % 24) + 1


def _draw_grid(
    ax: plt.Axes,
    color_map: dict[tuple[int, int], str],
    title: str,
    legend_items: list[tuple[str, str]],
) -> None:
    ax.set_xlim(0, 24)
    ax.set_ylim(16, 0)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=11, pad=10)

    for r in range(16):
        for c in range(1, 25):
            color = color_map.get((r, c), "#f3f4f6")
            rect = patches.Rectangle(
                (c - 1, r),
                1,
                1,
                facecolor=color,
                edgecolor="#c7ccd1",
                linewidth=0.6,
            )
            ax.add_patch(rect)

    ax.set_xticks([i - 0.5 for i in range(1, 25)])
    ax.set_xticklabels([str(i) for i in range(1, 25)], fontsize=7)
    ax.set_yticks([i + 0.5 for i in range(16)])
    ax.set_yticklabels([chr(ord("A") + i) for i in range(16)], fontsize=8)
    ax.tick_params(length=0)

    handles = [patches.Patch(facecolor=color, edgecolor="#9aa0a6", label=label) for label, color in legend_items]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=8, frameon=False)


def build_figure(integrated_csv: Path, output_png: Path) -> None:
    df = pd.read_csv(integrated_csv)
    if not {"Well_number_384", "Is_NT_ctrl", "Is_pos_ctrl", "Plate_number_384"}.issubset(df.columns):
        raise ValueError(
            "Integrated CSV must contain Well_number_384, Is_NT_ctrl, Is_pos_ctrl, Plate_number_384."
        )

    df = df.copy()
    df["Is_NT_ctrl"] = df["Is_NT_ctrl"].astype(bool)
    df["Is_pos_ctrl"] = df["Is_pos_ctrl"].astype(bool)

    n_plates = int(df["Plate_number_384"].nunique())
    freq = (
        df.groupby("Well_number_384", dropna=True)
        .agg(nt=("Is_NT_ctrl", "sum"), pos=("Is_pos_ctrl", "sum"))
        .reset_index()
    )
    majority = max(1, n_plates // 2)

    campaign_map: dict[tuple[int, int], str] = {}
    for _, row in freq.iterrows():
        r, c = _well_to_row_col(int(row["Well_number_384"]))
        nt = int(row["nt"])
        pos = int(row["pos"])
        if nt >= majority:
            campaign_map[(r, c)] = "#2f78b7"  # Non-targeting primary
        elif pos >= majority:
            campaign_map[(r, c)] = "#d84a4a"  # Pos primary
        elif nt > 0:
            campaign_map[(r, c)] = "#8dbce3"  # Non-targeting subset
        elif pos > 0:
            campaign_map[(r, c)] = "#f2a3a3"  # Pos subset
        else:
            campaign_map[(r, c)] = "#f3f4f6"  # Library/default

    # TR-FRET correction reference wells used by the algorithm.
    trfret_map: dict[tuple[int, int], str] = {(r, c): "#f3f4f6" for r in range(16) for c in range(1, 25)}
    odd_rows = [0, 2, 4, 6, 8, 10, 12, 14]
    even_rows = [1, 3, 5, 7, 9, 11, 13, 15]
    for r in odd_rows:
        trfret_map[(r, 1)] = "#4b9c8f"  # ch1 baseline set
        trfret_map[(r, 24)] = "#7dcfbf"  # final recenter set
    for r in even_rows:
        trfret_map[(r, 1)] = "#8d78cc"  # k numerator set
        trfret_map[(r, 24)] = "#b6a6e8"  # k denominator/ch2 baseline set

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), constrained_layout=True)

    _draw_grid(
        axes[0],
        campaign_map,
        title=f"Control Layout In Current Dataset (n={n_plates} plates)",
        legend_items=[
            ("Non-targeting control (majority plates)", "#2f78b7"),
            ("Positive control (majority plates)", "#d84a4a"),
            ("Non-targeting control (subset of plates)", "#8dbce3"),
            ("Positive control (subset of plates)", "#f2a3a3"),
        ],
    )

    _draw_grid(
        axes[1],
        trfret_map,
        title="TR-FRET Correction Reference Wells",
        legend_items=[
            ("Odd rows, col 1 (ch1 baseline)", "#4b9c8f"),
            ("Odd rows, col 24 (final recenter)", "#7dcfbf"),
            ("Even rows, col 1 (k numerator)", "#8d78cc"),
            ("Even rows, col 24 (k denominator/ch2 baseline)", "#b6a6e8"),
        ],
    )

    fig.suptitle("384-Well Plate Reference Map (Rows A-P, Columns 1-24)", fontsize=13, y=1.02)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[2]
    integrated = repo_root / "results" / "01_integrated.csv"
    out = repo_root / "docs" / "figures" / "plate_layout_reference.png"
    build_figure(integrated, out)
    print(f"Wrote {out}")
