from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def heatmap_384(values: np.ndarray, title: str = "Heatmap"):
    arr = np.asarray(values, dtype=float).reshape(16, 24)
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    im = ax.imshow(arr, cmap="viridis", aspect="auto")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    return fig, ax


def heatmap_for_plate(df: pd.DataFrame, plate: int | str, column: str):
    sub = df[df["Plate_number_384"].astype(str) == str(plate)].copy()
    sub = sub.sort_values("Well_number_384")
    values = pd.to_numeric(sub[column], errors="coerce").to_numpy()
    if values.size != 384:
        raise ValueError(f"Expected 384 wells for plate {plate}, got {values.size}")
    return heatmap_384(values, title=f"Plate {plate}: {column}")