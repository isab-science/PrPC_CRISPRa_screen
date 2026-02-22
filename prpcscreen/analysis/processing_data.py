from __future__ import annotations

import numpy as np
import pandas as pd

from .calculating_scores import (
    calculate_p,
    calculate_p_from_ssmd_t,
    calculate_ssmd,
    calculate_ssmd_moderated,
    norm_plates,
)


def normalize_with_nt_controls(df: pd.DataFrame, norm_method: str = "genes and all Non-targeting") -> pd.DataFrame:
    out = df.copy()
    has_glo = "CellTiterGlo_raw" in out.columns
    if has_glo:
        out["CellTiterGlo_foldNT"] = norm_plates(out, "CellTiterGlo_raw", fold_nt=True, norm_method=norm_method)
    else:
        out["CellTiterGlo_foldNT"] = np.nan

    for ri in (1, 2):
        raw_col = f"Raw_rep{ri}"
        out[f"DeltaNT_rep{ri}"] = norm_plates(out, raw_col, norm_method=norm_method)
        out[f"FoldNT_rep{ri}"] = norm_plates(out, raw_col, fold_nt=True, norm_method=norm_method)
        out[f"PercActivation_rep{ri}"] = norm_plates(out, raw_col, percent_activation=True, norm_method=norm_method)

        raw = pd.to_numeric(out[raw_col], errors="coerce")
        out[f"Raw_log2_rep{ri}"] = np.log2(np.clip(raw, 1e-12, None))
        out[f"Log2FC_rep{ri}"] = norm_plates(out, raw_col, take_log2=True, norm_method=norm_method)

        if has_glo:
            out[f"Raw_Glo_rep{ri}"] = raw / pd.to_numeric(out["CellTiterGlo_foldNT"], errors="coerce")
            out[f"DeltaNT_Glo_rep{ri}"] = norm_plates(out, f"Raw_Glo_rep{ri}", norm_method=norm_method)
            out[f"FoldNT_Glo_rep{ri}"] = norm_plates(out, f"Raw_Glo_rep{ri}", fold_nt=True, norm_method=norm_method)
            out[f"Log2FC_Glo_rep{ri}"] = norm_plates(out, f"Raw_Glo_rep{ri}", take_log2=True, norm_method=norm_method)
        else:
            out[f"Raw_Glo_rep{ri}"] = np.nan
            out[f"DeltaNT_Glo_rep{ri}"] = np.nan
            out[f"FoldNT_Glo_rep{ri}"] = np.nan
            out[f"Log2FC_Glo_rep{ri}"] = np.nan
    return out


def run_ssmd_stats(df: pd.DataFrame, norm_method: str = "genes and all Non-targeting") -> pd.DataFrame:
    out = normalize_with_nt_controls(df, norm_method=norm_method)

    tracked_pairs = [
        ("Raw_rep1", "Raw_rep2", "raw"),
        ("Log2FC_rep1", "Log2FC_rep2", "log2"),
        ("Log2FC_Glo_rep1", "Log2FC_Glo_rep2", "log2_glo"),
    ]
    for rep1, rep2, stem in tracked_pairs:
        if rep1 in out and rep2 in out:
            out[f"SSMD_{stem}"] = calculate_ssmd(out, rep1, rep2)
            out[f"SSMD_mod_{stem}"] = calculate_ssmd_moderated(out, rep1, rep2)
            out[f"Mean_{stem}"] = pd.to_numeric(out[rep1], errors="coerce").add(pd.to_numeric(out[rep2], errors="coerce")).div(2)
            # R-compatible volcano p-values: moderated SSMD interpreted via t(df=1).
            out[f"p_value_{stem}"] = calculate_p_from_ssmd_t(out[f"SSMD_mod_{stem}"], df=1)
            # Keep replicate-difference p-value as a reproducibility diagnostic.
            out[f"p_value_repro_{stem}"] = calculate_p(out, rep1, rep2)
    return out


def create_hit_lists(
    df: pd.DataFrame,
    p_column: str = "p_value_log2",
    effect_column: str = "Mean_log2",
    p_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
) -> dict[str, pd.DataFrame]:
    out = df.copy()
    p_vals = pd.to_numeric(out.get(p_column), errors="coerce")
    effect = pd.to_numeric(out.get(effect_column), errors="coerce")

    hit_mask = p_vals.lt(p_cutoff) & effect.abs().gt(log2fc_cutoff)
    down = out[hit_mask & effect.lt(0)].copy()
    up = out[hit_mask & effect.gt(0)].copy()

    if "Gene_symbol" in out.columns:
        reordered = out.sort_values(by=["Gene_symbol"], na_position="last")
    else:
        reordered = out

    return {
        "original_df": out,
        "reordered_df": reordered,
        "down_hits": down,
        "up_hits": up,
        "hits": pd.concat([down, up], axis=0),
    }
