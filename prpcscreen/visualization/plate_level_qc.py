from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_plate_qualities(
    rep1: np.ndarray,
    rep2: np.ndarray,
    ylabel: str = "Quality metric",
    plate_labels: list[str] | None = None,
):
    rep1_vals = np.asarray(rep1, dtype=float)
    rep2_vals = np.asarray(rep2, dtype=float)
    n = len(rep1_vals)
    labels = list(plate_labels) if plate_labels and len(plate_labels) == n else [str(i + 1) for i in range(n)]

    # Reorder plates by ascending pair mean so weakest-quality pairs appear first.
    pair_mean = np.nanmean(np.vstack([rep1_vals, rep2_vals]), axis=0)
    sort_key = np.where(np.isfinite(pair_mean), pair_mean, np.inf)
    order = np.argsort(sort_key, kind="stable")
    rep1_vals = rep1_vals[order]
    rep2_vals = rep2_vals[order]
    labels = [labels[i] for i in order]

    x = np.arange(1, n + 1)
    fig, ax = plt.subplots(figsize=(13, 3.8))

    all_vals = np.concatenate([rep1_vals, rep2_vals])
    all_vals = all_vals[np.isfinite(all_vals)]
    y_lo = min(-0.2, float(all_vals.min()) - 0.05) if all_vals.size else -0.2
    y_hi = max(1.0, float(all_vals.max()) + 0.05) if all_vals.size else 1.0

    pair_colors = np.array(["#d62728" if i % 2 == 0 else "#1f77b4" for i in range(n)], dtype=object)
    for idx, xi in enumerate(x):
        color = str(pair_colors[idx])
        y1 = rep1_vals[idx]
        y2 = rep2_vals[idx]
        if np.isfinite(y1) and np.isfinite(y2):
            ax.plot([xi, xi], [y1, y2], color=color, lw=1.2, zorder=2)
        if np.isfinite(y1):
            ax.scatter([xi], [y1], s=18, c=color, edgecolors=color, linewidths=0.6, zorder=3)
        if np.isfinite(y2):
            ax.scatter([xi], [y2], s=18, c=color, edgecolors=color, linewidths=0.6, zorder=3)

    ax.set_xlim(0.4, len(x) + 0.6)
    ax.set_ylim(y_lo, y_hi)
    ax.set_xlabel("Plates")
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5)
    for idx, tick in enumerate(ax.get_xticklabels()):
        tick.set_color(str(pair_colors[idx]))

    ax.tick_params(axis="x", labelsize=6)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return fig, ax


def get_quality_metric(df: pd.DataFrame, fn, filter_nt: bool = False) -> tuple[np.ndarray, np.ndarray]:
    mask = np.ones(len(df), dtype=bool)
    if filter_nt and "Target_flag" in df:
        flag = df["Target_flag"].fillna("").astype(str).str.strip().str.lower()
        own_non_targeting = flag.isin({"own non-targeting control", "own non targeting control", "own nt control"})
        mask &= ~own_non_targeting.to_numpy()

    rep1, rep2 = [], []
    for _, plate_df in df[mask].groupby("Plate_number_384"):
        rep1.append(fn(plate_df, "Raw_rep1"))
        rep2.append(fn(plate_df, "Raw_rep2"))
    return np.asarray(rep1, dtype=float), np.asarray(rep2, dtype=float)
