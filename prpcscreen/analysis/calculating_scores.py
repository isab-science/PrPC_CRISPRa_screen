from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats

NORMALIZATION_METHODS = (
    "all Non-targeting",
    "own Non-targeting",
    "own Non-targeting, with factor",
    "genes",
    "genes and all Non-targeting",
    "genes and own Non-targeting",
)

_NORM_METHOD_ALIASES = {
    "all nt": "all Non-targeting",
    "own nt": "own Non-targeting",
    "own nt, with factor": "own Non-targeting, with factor",
    "genes and all nt": "genes and all Non-targeting",
    "genes and own nt": "genes and own Non-targeting",
}

P_VALUE_FLOOR = 1e-300


def _two_sided_p_from_z(z: pd.Series | np.ndarray) -> pd.Series:
    # Use log-survival for numerical stability at extreme z-scores.
    z_vals = pd.to_numeric(pd.Series(z), errors="coerce")
    log_p = np.log(2.0) + stats.norm.logsf(np.abs(z_vals))
    p = np.exp(log_p)
    return pd.Series(np.clip(p, P_VALUE_FLOOR, 1.0), index=z_vals.index, dtype=float)


def _plate_numbers(df: pd.DataFrame) -> pd.Series:
    raw = df["Plate_number_384"]
    numeric = pd.to_numeric(raw, errors="coerce")
    if numeric.notna().any():
        return numeric

    # Support non-numeric plate labels (for example Roman numerals) via stable factorization.
    text = raw.astype(str).str.strip()
    text = text.replace({"": np.nan, "nan": np.nan, "None": np.nan})
    codes, _ = pd.factorize(text, sort=True)
    return pd.Series(np.where(codes < 0, np.nan, codes + 1), index=df.index, dtype=float)


def _median_safe(values: pd.Series) -> float:
    vals = pd.to_numeric(values, errors="coerce").dropna()
    if vals.empty:
        return np.nan
    return float(vals.median())


def _canonical_norm_method(norm_method: str) -> str:
    raw = str(norm_method).strip()
    key = raw.lower()
    if key in _NORM_METHOD_ALIASES:
        return _NORM_METHOD_ALIASES[key]
    return raw


def _own_non_targeting_mask(df: pd.DataFrame) -> pd.Series:
    if "Target_flag" not in df.columns:
        return pd.Series(False, index=df.index, dtype=bool)
    flag = df["Target_flag"].fillna("").astype(str).str.strip().str.lower()
    return flag.isin({"own non-targeting control", "own non targeting control", "own nt control"})


def obtain_nt_factors(df: pd.DataFrame, use_column: str, take_log2: bool = False) -> dict[int, float]:
    values = pd.to_numeric(df[use_column], errors="coerce")
    if take_log2:
        values = np.log2(np.clip(values, 1e-12, None))

    plate_numbers = _plate_numbers(df)
    own_nt = _own_non_targeting_mask(df)
    tubingen_nt = df["Is_NT_ctrl"].astype(bool) & ~own_nt

    factors: dict[int, float] = {}
    for plate in sorted(pd.unique(plate_numbers[tubingen_nt].dropna())):
        is_plate = plate_numbers.eq(plate)
        med_own = _median_safe(values[is_plate & own_nt])
        med_tub = _median_safe(values[is_plate & tubingen_nt])
        factors[int(plate)] = med_tub / med_own if med_own and not np.isnan(med_own) else np.nan
    return factors


def norm_plates(
    df: pd.DataFrame,
    use_column: str,
    fold_nt: bool = False,
    percent_activation: bool = False,
    take_log2: bool = False,
    norm_method: str = "all Non-targeting",
) -> pd.Series:
    norm_method = _canonical_norm_method(norm_method)
    if norm_method not in NORMALIZATION_METHODS:
        raise ValueError(f"Unknown normalization method: {norm_method}")

    plate_numbers = _plate_numbers(df)
    values = pd.to_numeric(df[use_column], errors="coerce")
    if take_log2:
        values = np.log2(np.clip(values, 1e-12, None))

    own_nt = _own_non_targeting_mask(df)
    any_nt = df["Is_NT_ctrl"].astype(bool) if "Is_NT_ctrl" in df else pd.Series(False, index=df.index)
    gene = df["Entrez_ID"].notna() if "Entrez_ID" in df else pd.Series(False, index=df.index)
    pos = df["Is_pos_ctrl"].astype(bool) if "Is_pos_ctrl" in df else pd.Series(False, index=df.index)

    if norm_method == "all Non-targeting":
        norm_group = any_nt
    elif norm_method == "own Non-targeting":
        norm_group = own_nt
    elif norm_method == "own Non-targeting, with factor":
        factors = obtain_nt_factors(df, use_column, take_log2=take_log2)
        norm_group = own_nt
    elif norm_method == "genes":
        norm_group = gene
    elif norm_method == "genes and all Non-targeting":
        # R-compatible behavior: use gene-population baseline.
        norm_group = gene
    elif norm_method == "genes and own Non-targeting":
        norm_group = gene
    else:
        norm_group = any_nt

    out = pd.Series(np.nan, index=df.index, dtype=float)
    for plate in sorted(pd.unique(plate_numbers.dropna())):
        is_plate = plate_numbers.eq(plate)
        baseline = _median_safe(values[is_plate & norm_group])
        if norm_method == "own Non-targeting, with factor":
            baseline *= factors.get(int(plate), 1.0)

        if np.isnan(baseline):
            continue

        if fold_nt:
            out[is_plate] = values[is_plate] / baseline
        elif percent_activation:
            pos_med = _median_safe(values[is_plate & pos])
            denom = pos_med - baseline
            out[is_plate] = np.where(np.isclose(denom, 0), np.nan, (values[is_plate] - baseline) / denom * 100)
        else:
            out[is_plate] = values[is_plate] - baseline

    return out


def calculate_ssmd(df: pd.DataFrame, rep1_column: str, rep2_column: str) -> pd.Series:
    rep1 = pd.to_numeric(df[rep1_column], errors="coerce")
    rep2 = pd.to_numeric(df[rep2_column], errors="coerce")
    mean_diff = rep1 - rep2
    pooled_sd = np.sqrt(rep1.var(skipna=True) + rep2.var(skipna=True))
    return mean_diff / pooled_sd if pooled_sd and not np.isnan(pooled_sd) else pd.Series(np.nan, index=df.index)


def calculate_ssmd_moderated(
    df: pd.DataFrame,
    rep1_column: str,
    rep2_column: str,
    plate_column: str = "Plate_number_384",
    nt_column: str = "Is_NT_ctrl",
) -> pd.Series:
    """
    R-compatible moderated SSMD:
    - per row: mean(rep1, rep2) / sqrt(0.5 * var_row + 0.5 * s0_plate)
    - s0_plate = median row-variance among Non-targeting controls on that plate
    """
    rep1 = pd.to_numeric(df[rep1_column], errors="coerce")
    rep2 = pd.to_numeric(df[rep2_column], errors="coerce")
    mean_row = (rep1 + rep2) / 2.0
    var_row = ((rep1 - rep2) ** 2) / 2.0  # sample var of two values

    nt_mask = df.get(nt_column, pd.Series(False, index=df.index)).astype(bool)
    plate_numbers = _plate_numbers(df) if plate_column == "Plate_number_384" else pd.to_numeric(df.get(plate_column), errors="coerce")
    out = pd.Series(np.nan, index=df.index, dtype=float)

    for plate in sorted(pd.unique(plate_numbers.dropna())):
        is_plate = plate_numbers.eq(plate)
        nt_vars = var_row[is_plate & nt_mask].dropna()
        if nt_vars.empty:
            continue
        s0 = float(nt_vars.median())
        denom = np.sqrt(0.5 * var_row[is_plate] + 0.5 * s0)
        out.loc[is_plate] = np.where(np.isclose(denom, 0), np.nan, mean_row[is_plate] / denom)

    return out


def calculate_p_from_ssmd_t(ssmd: pd.Series | np.ndarray, df: int = 1) -> pd.Series:
    vals = pd.to_numeric(pd.Series(ssmd), errors="coerce")
    p = 2.0 * stats.t.sf(np.abs(vals), df=df)
    return pd.Series(np.clip(p, P_VALUE_FLOOR, 1.0), index=vals.index, dtype=float)


def calculate_t(df: pd.DataFrame, rep1_column: str, rep2_column: str) -> pd.Series:
    rep1 = pd.to_numeric(df[rep1_column], errors="coerce")
    rep2 = pd.to_numeric(df[rep2_column], errors="coerce")
    return rep1 - rep2


def calculate_p(df: pd.DataFrame, rep1_column: str, rep2_column: str) -> pd.Series:
    rep1 = pd.to_numeric(df[rep1_column], errors="coerce")
    rep2 = pd.to_numeric(df[rep2_column], errors="coerce")
    diffs = rep1 - rep2
    sd = diffs.std(skipna=True)
    if not sd or np.isnan(sd):
        return pd.Series(np.nan, index=df.index)
    z = (diffs - diffs.mean(skipna=True)) / sd
    # Use survival function for numerical stability and clip to avoid exact zeros.
    return _two_sided_p_from_z(z).reindex(df.index)


def calculate_p_vs_nt(
    df: pd.DataFrame,
    value_column: str,
    plate_column: str = "Plate_number_384",
    nt_column: str = "Is_NT_ctrl",
) -> pd.Series:
    values = pd.to_numeric(df.get(value_column), errors="coerce")
    if values.notna().sum() == 0:
        return pd.Series(np.nan, index=df.index, dtype=float)

    nt_mask = df.get(nt_column, pd.Series(False, index=df.index)).astype(bool)
    if nt_mask.sum() == 0:
        return pd.Series(np.nan, index=df.index, dtype=float)

    # Compute Non-targeting-based z-scores per plate; fill gaps from global Non-targeting stats.
    out = pd.Series(np.nan, index=df.index, dtype=float)
    plate_numbers = _plate_numbers(df) if plate_column == "Plate_number_384" else pd.to_numeric(df.get(plate_column), errors="coerce")

    for plate in sorted(pd.unique(plate_numbers.dropna())):
        is_plate = plate_numbers.eq(plate)
        nt_vals = values[is_plate & nt_mask].dropna()
        if nt_vals.size < 2:
            continue
        nt_sd = nt_vals.std()
        if np.isclose(nt_sd, 0) or np.isnan(nt_sd):
            continue
        z = (values[is_plate] - nt_vals.mean()) / nt_sd
        out.loc[is_plate] = _two_sided_p_from_z(z).to_numpy()

    unresolved = out.isna() & values.notna()
    if unresolved.any():
        nt_vals_global = values[nt_mask].dropna()
        if nt_vals_global.size >= 2:
            nt_sd_global = nt_vals_global.std()
            if not np.isclose(nt_sd_global, 0) and not np.isnan(nt_sd_global):
                z_global = (values[unresolved] - nt_vals_global.mean()) / nt_sd_global
                out.loc[unresolved] = _two_sided_p_from_z(z_global).to_numpy()

    unresolved = out.isna() & values.notna()
    if unresolved.any():
        # Last fallback to avoid fully-empty volcanoes when Non-targeting variance is degenerate.
        vals_sd = values.std(skipna=True)
        if not np.isclose(vals_sd, 0) and not np.isnan(vals_sd):
            z_all = (values[unresolved] - values.mean(skipna=True)) / vals_sd
            out.loc[unresolved] = _two_sided_p_from_z(z_all).to_numpy()

    return out.clip(lower=P_VALUE_FLOOR, upper=1.0)


def calculate_z_prime(df: pd.DataFrame, use_column: str, filter_nt: bool = False) -> float:
    pos = pd.to_numeric(df.loc[df["Is_pos_ctrl"].astype(bool), use_column], errors="coerce").dropna()
    nt_mask = df["Is_NT_ctrl"].astype(bool)
    if filter_nt and "Target_flag" in df:
        nt_mask &= ~_own_non_targeting_mask(df)
    nt = pd.to_numeric(df.loc[nt_mask, use_column], errors="coerce").dropna()
    if pos.empty or nt.empty:
        return np.nan
    denom = abs(pos.mean() - nt.mean())
    if np.isclose(denom, 0):
        return np.nan
    return 1 - (3 * (pos.std() + nt.std()) / denom)


def calculate_ssmd_ctrls(df: pd.DataFrame, use_column: str, filter_nt: bool = False) -> float:
    pos = pd.to_numeric(df.loc[df["Is_pos_ctrl"].astype(bool), use_column], errors="coerce").dropna()
    nt_mask = df["Is_NT_ctrl"].astype(bool)
    if filter_nt and "Target_flag" in df:
        nt_mask &= ~_own_non_targeting_mask(df)
    nt = pd.to_numeric(df.loc[nt_mask, use_column], errors="coerce").dropna()
    if pos.empty or nt.empty:
        return np.nan
    denom = np.sqrt(pos.var() + nt.var())
    if np.isclose(denom, 0):
        return np.nan
    return (pos.mean() - nt.mean()) / denom
