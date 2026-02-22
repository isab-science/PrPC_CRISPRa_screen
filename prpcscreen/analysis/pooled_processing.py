from __future__ import annotations

from dataclasses import dataclass
import re
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .calculating_scores import calculate_p

P_VALUE_FLOOR = 1e-300

DEFAULT_REFERENCE_REGEX = r"(?i)^Negative_R\d+$"
DEFAULT_TREATMENT_REGEX = r"(?i)^Positive_R\d+$"
DEFAULT_NT_REGEX = r"(?i)(?:^control(?:_|$)|^nt(?:_|$)|non[-_ ]?target|negative[_ ]?control|nt[_ -]?ctrl)"
DEFAULT_POS_REGEX = r"(?i)(?:^prnp$|positive[_ ]?control|pos[_ -]?ctrl)"

_TRUE_TOKENS = {"1", "true", "t", "yes", "y", "on"}
_FALSE_TOKENS = {"0", "false", "f", "no", "n", "off"}

_PREFERRED_SHEETS = (
    "Primary_RawData",
    "Primary_Raw",
    "RawData",
    "davide_analysis",
    "Davide_analysis",
)


@dataclass(frozen=True)
class PooledAnalysisConfig:
    reference_cols: tuple[str, ...]
    treatment_cols: tuple[str, ...]
    pseudocount: float = 1.0
    pvalue_method: str = "welch"
    p_cutoff: float = 0.05
    log2fc_cutoff: float = 0.3
    nt_regex: str = DEFAULT_NT_REGEX
    pos_regex: str = DEFAULT_POS_REGEX


def _replicate_sort_key(name: str) -> tuple[int, int, str]:
    match = re.search(r"(\d+)$", str(name))
    if match:
        return (0, int(match.group(1)), str(name))
    return (1, 0, str(name))


def _to_numeric_frame(df: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    for col in columns:
        out[col] = pd.to_numeric(df[col], errors="coerce")
    return out


def _coerce_bool_series(values: pd.Series | None, index: pd.Index, default: bool = False) -> pd.Series:
    if values is None:
        return pd.Series(default, index=index, dtype=bool)

    if pd.api.types.is_bool_dtype(values):
        return values.fillna(default).astype(bool)

    text = values.fillna("").astype(str).str.strip().str.lower()
    out = pd.Series(default, index=index, dtype=bool)
    out.loc[text.isin(_TRUE_TOKENS)] = True
    out.loc[text.isin(_FALSE_TOKENS)] = False
    return out


def benjamini_hochberg(p_values: pd.Series | np.ndarray) -> pd.Series:
    p = pd.to_numeric(pd.Series(p_values), errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)

    valid = p.dropna().astype(float)
    if valid.empty:
        return out.fillna(1.0)

    n = len(valid)
    order = np.argsort(valid.to_numpy())
    ranked = valid.to_numpy()[order]
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, P_VALUE_FLOOR, 1.0)

    restored = np.empty(n, dtype=float)
    restored[order] = adjusted
    out.loc[valid.index] = restored
    return out.fillna(1.0)


def detect_replicate_columns(df: pd.DataFrame, pattern: str) -> list[str]:
    regex = re.compile(pattern)
    hits = [str(col) for col in df.columns if regex.search(str(col))]
    return sorted(hits, key=_replicate_sort_key)


def resolve_replicate_columns(
    df: pd.DataFrame,
    reference_cols: Sequence[str] | None = None,
    treatment_cols: Sequence[str] | None = None,
    reference_regex: str = DEFAULT_REFERENCE_REGEX,
    treatment_regex: str = DEFAULT_TREATMENT_REGEX,
) -> tuple[list[str], list[str]]:
    if reference_cols:
        ref = [str(col) for col in reference_cols]
    else:
        ref = detect_replicate_columns(df, reference_regex)

    if treatment_cols:
        treat = [str(col) for col in treatment_cols]
    else:
        treat = detect_replicate_columns(df, treatment_regex)

    missing = [col for col in ref + treat if col not in df.columns]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing requested replicate columns: {joined}")

    if len(ref) == 0 or len(treat) == 0:
        raise ValueError(
            "Could not determine pooled replicate columns. "
            f"Detected reference={ref} treatment={treat}. "
            "Pass --reference-cols and --treatment-cols explicitly."
        )

    return ref, treat


def load_pooled_table(path: str | Path, sheet: str | None = None) -> tuple[pd.DataFrame, str | None]:
    input_path = Path(path)
    suffix = input_path.suffix.lower()

    if suffix in {".xlsx", ".xls", ".xlsm"}:
        if sheet:
            return pd.read_excel(input_path, sheet_name=sheet), sheet

        workbook = pd.ExcelFile(input_path)
        names = list(workbook.sheet_names)
        name_lookup = {name.lower(): name for name in names}
        for candidate in _PREFERRED_SHEETS:
            if candidate.lower() in name_lookup:
                chosen = name_lookup[candidate.lower()]
                return pd.read_excel(input_path, sheet_name=chosen), chosen

        chosen = names[0]
        return pd.read_excel(input_path, sheet_name=chosen), chosen

    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(input_path, sep=None, engine="python"), None
    return pd.read_csv(input_path), None


def _first_existing(columns: Sequence[str], candidates: Sequence[str]) -> str | None:
    lookup = {str(col).strip().lower(): str(col) for col in columns}
    for candidate in candidates:
        hit = lookup.get(str(candidate).strip().lower())
        if hit is not None:
            return hit
    return None


def ensure_standard_annotation_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    columns = [str(col) for col in out.columns]

    if "Gene_symbol" not in out.columns:
        gene_col = _first_existing(
            columns,
            ("Gene_symbol", "Gene", "GeneSymbol", "Gene_TSS", "All_targeted_genes", "label"),
        )
        if gene_col is None:
            out["Gene_symbol"] = ""
        else:
            out["Gene_symbol"] = out[gene_col]

    if "Entrez_ID" not in out.columns:
        entrez_col = _first_existing(columns, ("Entrez_ID", "Entrez", "Entrez ID"))
        if entrez_col is None:
            out["Entrez_ID"] = np.nan
        else:
            out["Entrez_ID"] = out[entrez_col]

    if "TSS_ID" not in out.columns:
        tss_col = _first_existing(columns, ("TSS_ID", "TSS", "Gene_TSS"))
        out["TSS_ID"] = out[tss_col] if tss_col is not None else ""

    if "Plasmid_ID" not in out.columns:
        plasmid_col = _first_existing(columns, ("Plasmid_ID", "Guide", "sgRNA_sequence", "sgRNA", "guide"))
        out["Plasmid_ID"] = out[plasmid_col] if plasmid_col is not None else ""

    return out


def derive_control_flags(
    df: pd.DataFrame,
    nt_regex: str = DEFAULT_NT_REGEX,
    pos_regex: str = DEFAULT_POS_REGEX,
) -> tuple[pd.Series, pd.Series]:
    nt = _coerce_bool_series(df.get("Is_NT_ctrl"), index=df.index, default=False)
    pos = _coerce_bool_series(df.get("Is_pos_ctrl"), index=df.index, default=False)

    gene_text = df.get("Gene_symbol", pd.Series("", index=df.index)).fillna("").astype(str)
    if nt_regex:
        nt = nt | gene_text.str.contains(nt_regex, regex=True, na=False)
    if pos_regex:
        pos = pos | gene_text.str.contains(pos_regex, regex=True, na=False)

    if "diffExpress" in df.columns:
        is_control = df["diffExpress"].fillna("").astype(str).str.strip().str.lower().eq("control")
        nt = nt | is_control

    nt = nt & ~pos
    return nt.astype(bool), pos.astype(bool)


def estimate_size_factors(counts: pd.DataFrame) -> pd.Series:
    numeric = counts.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        log_vals = np.log(numeric.replace(0.0, np.nan))
    geometric_means = np.exp(log_vals.mean(axis=1, skipna=True))

    size_factors = {}
    for col in numeric.columns:
        ratios = numeric[col] / geometric_means
        ratios = ratios.replace([np.inf, -np.inf], np.nan).dropna()
        ratios = ratios[ratios > 0]
        size_factors[col] = float(ratios.median()) if not ratios.empty else np.nan

    sf = pd.Series(size_factors, dtype=float)

    if sf.isna().any() or (sf <= 0).any():
        library_sizes = numeric.sum(axis=0)
        library_median = float(library_sizes.median()) if library_sizes.median() > 0 else 1.0
        fallback = library_sizes / library_median
        sf = sf.where(sf > 0, fallback)

    sf = sf.replace([np.inf, -np.inf], np.nan)
    if sf.isna().all():
        sf[:] = 1.0
    else:
        positive = sf[sf > 0]
        fill_value = float(positive.median()) if not positive.empty else 1.0
        sf = sf.fillna(fill_value).clip(lower=1e-12)
        sf = sf / float(np.exp(np.log(sf).mean()))

    return sf.astype(float)


def prepare_pooled_integrated_table(
    df: pd.DataFrame,
    reference_cols: Sequence[str],
    treatment_cols: Sequence[str],
    nt_regex: str = DEFAULT_NT_REGEX,
    pos_regex: str = DEFAULT_POS_REGEX,
) -> pd.DataFrame:
    out = ensure_standard_annotation_columns(df)
    nt, pos = derive_control_flags(out, nt_regex=nt_regex, pos_regex=pos_regex)
    out["Is_NT_ctrl"] = nt
    out["Is_pos_ctrl"] = pos

    for col in list(reference_cols) + list(treatment_cols):
        out[col] = pd.to_numeric(out[col], errors="coerce")

    if "Raw_rep1" not in out.columns and len(treatment_cols) >= 1:
        out["Raw_rep1"] = pd.to_numeric(out[treatment_cols[0]], errors="coerce")
    if "Raw_rep2" not in out.columns and len(treatment_cols) >= 2:
        out["Raw_rep2"] = pd.to_numeric(out[treatment_cols[1]], errors="coerce")
    if "Raw_ref1" not in out.columns and len(reference_cols) >= 1:
        out["Raw_ref1"] = pd.to_numeric(out[reference_cols[0]], errors="coerce")
    if "Raw_ref2" not in out.columns and len(reference_cols) >= 2:
        out["Raw_ref2"] = pd.to_numeric(out[reference_cols[1]], errors="coerce")

    out["Analysis_mode"] = "pooled"
    return out


def compute_pooled_analysis(df: pd.DataFrame, config: PooledAnalysisConfig) -> pd.DataFrame:
    out = prepare_pooled_integrated_table(
        df,
        reference_cols=config.reference_cols,
        treatment_cols=config.treatment_cols,
        nt_regex=config.nt_regex,
        pos_regex=config.pos_regex,
    )

    sample_cols = list(config.reference_cols) + list(config.treatment_cols)
    counts = _to_numeric_frame(out, sample_cols).fillna(0.0).clip(lower=0.0)
    size_factors = estimate_size_factors(counts)
    norm_counts = counts.div(size_factors, axis=1)
    log_norm = np.log2(norm_counts + float(config.pseudocount))

    for col in sample_cols:
        out[f"SizeFactor_{col}"] = float(size_factors[col])
        out[f"{col}_norm"] = norm_counts[col]
        out[f"{col}_log2norm"] = log_norm[col]

    ref_cols = list(config.reference_cols)
    treat_cols = list(config.treatment_cols)
    log_ref = log_norm[ref_cols]
    log_treat = log_norm[treat_cols]

    out["Mean_log2_ref"] = log_ref.mean(axis=1)
    out["Mean_log2_treat"] = log_treat.mean(axis=1)
    out["Mean_log2"] = out["Mean_log2_treat"] - out["Mean_log2_ref"]
    out["log2Ratio"] = out["Mean_log2"]

    pair_count = min(len(ref_cols), len(treat_cols))
    num_log2fc_reps = max(2, pair_count)
    for idx in range(num_log2fc_reps):
        col_name = f"Log2FC_rep{idx + 1}"
        if idx < pair_count:
            out[col_name] = log_treat.iloc[:, idx] - log_ref.iloc[:, idx]
        else:
            out[col_name] = np.nan

    if config.pvalue_method == "paired":
        if pair_count >= 2:
            _, p_vals = stats.ttest_rel(
                log_treat.iloc[:, :pair_count].to_numpy(dtype=float),
                log_ref.iloc[:, :pair_count].to_numpy(dtype=float),
                axis=1,
                nan_policy="omit",
            )
        else:
            p_vals = np.full(len(out), np.nan, dtype=float)
    else:
        equal_var = config.pvalue_method == "student"
        _, p_vals = stats.ttest_ind(
            log_treat.to_numpy(dtype=float),
            log_ref.to_numpy(dtype=float),
            axis=1,
            equal_var=equal_var,
            nan_policy="omit",
        )

    p_vals = np.nan_to_num(p_vals, nan=1.0, posinf=1.0, neginf=P_VALUE_FLOOR)
    p_vals = np.clip(p_vals, P_VALUE_FLOOR, 1.0)
    fdr_vals = benjamini_hochberg(p_vals)

    out["p_value_log2"] = p_vals
    out["pValue"] = p_vals
    out["fdr_log2"] = fdr_vals
    out["fdr"] = fdr_vals
    out["p_value_log2_method"] = f"ttest_{config.pvalue_method}_log2norm"

    if out["Log2FC_rep1"].notna().any() and out["Log2FC_rep2"].notna().any():
        out["p_value_repro_log2"] = calculate_p(out, "Log2FC_rep1", "Log2FC_rep2")
    else:
        out["p_value_repro_log2"] = np.nan

    effect = pd.to_numeric(out["Mean_log2"], errors="coerce")
    control_mask = out["Is_NT_ctrl"].astype(bool) | out["Is_pos_ctrl"].astype(bool)
    sig = pd.to_numeric(out["p_value_log2"], errors="coerce") < float(config.p_cutoff)
    is_up = sig & effect.ge(float(config.log2fc_cutoff))
    is_down = sig & effect.le(-float(config.log2fc_cutoff))

    out["diffExpress"] = np.where(
        control_mask,
        "Control",
        np.where(is_up, "UP", np.where(is_down, "DOWN", "NO")),
    )
    out["label"] = out["Gene_symbol"].fillna("").astype(str).str.strip()
    out["label"] = out["label"].where(out["label"].ne(""), out["Plasmid_ID"].fillna("").astype(str).str.strip())
    out["label"] = out["label"].replace("", "UNLABELED")

    out["Hit_p_cutoff"] = float(config.p_cutoff)
    out["Hit_log2fc_cutoff"] = float(config.log2fc_cutoff)
    out["Mean_log2FC"] = out["Mean_log2"]

    return out
