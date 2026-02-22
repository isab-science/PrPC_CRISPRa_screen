from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plate_well_plot(df: pd.DataFrame, use_column: str = "Raw_rep1", by_row: bool = False):
    fig, ax = plt.subplots(figsize=(11, 3.8))
    order_col = "Well_number_384" if not by_row else "x_position"
    if order_col in df.columns:
        sub = df.sort_values(order_col)
        x_label = "Well order"
    else:
        sub = df.reset_index(drop=True)
        x_label = "Feature order"
    vals = pd.to_numeric(sub[use_column], errors="coerce")
    ax.plot(np.arange(len(vals)), vals, lw=0.8)
    ax.set_xlabel(x_label)
    ax.set_ylabel(use_column)
    ax.set_title(f"Series plot: {use_column}")
    return fig, ax
