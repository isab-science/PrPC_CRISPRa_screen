from __future__ import annotations

import numpy as np


def setup_mappings() -> dict[str, np.ndarray]:
    plate_names = ["A1", "A2", "B1", "B2"]
    scheme_2rows = np.array([plate_names[:2] * 12, plate_names[2:] * 12], dtype=object)
    scheme = np.vstack([scheme_2rows for _ in range(8)])

    idx_384 = np.arange(1, 385, dtype=int).reshape(16, 24)
    out: dict[str, np.ndarray] = {}
    for plate in plate_names:
        mask = scheme == plate
        rows = [idx_384[r, :][mask[r, :]] for r in range(16)]
        out[plate] = np.vstack(rows)
    return out


def convert_well_numbers(well_numbers: np.ndarray | list[int]) -> np.ndarray:
    numbers = np.asarray(well_numbers, dtype=int)
    mappings = setup_mappings()
    inv: dict[int, int] = {}
    for arr in mappings.values():
        for i96, i384 in enumerate(arr.ravel(order="C"), start=1):
            inv[int(i384)] = i96
    return np.array([inv.get(int(n), -1) for n in numbers], dtype=int)