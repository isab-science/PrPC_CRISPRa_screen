from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd


def beebox_plates(df: pd.DataFrame, column: str, split_nt: bool = True):
    plot_df = df.copy()
    if split_nt and "Target_flag" in plot_df:
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

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    order = [g for g in ["Gene", "Non-targeting", "Own Non-targeting", "Positive"] if g in plot_df["group"].unique()]
    data = [pd.to_numeric(plot_df.loc[plot_df["group"] == g, column], errors="coerce").dropna() for g in order]

    # Violin layer (distribution shape) + box layer (quartiles/median) in one panel.
    violin = ax.violinplot(data, showmeans=False, showmedians=False, showextrema=False)
    for body in violin["bodies"]:
        body.set_facecolor("#90caf9")
        body.set_edgecolor("#4f83cc")
        body.set_alpha(0.55)

    box = ax.boxplot(data, labels=order, patch_artist=True, widths=0.22, showfliers=False)
    for patch in box["boxes"]:
        patch.set_facecolor("#ffffff")
        patch.set_edgecolor("#1f1f1f")
        patch.set_linewidth(1.0)
    for key in ("whiskers", "caps", "medians"):
        for line in box[key]:
            line.set_color("#1f1f1f")
            line.set_linewidth(1.0)

    ax.set_ylabel(column)
    ax.set_title(f"Violin + box plot: {column}")
    return fig, ax
