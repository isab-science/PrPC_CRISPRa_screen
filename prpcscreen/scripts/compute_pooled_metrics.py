from __future__ import annotations

import os
import sys
import argparse
from pathlib import Path

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from prpcscreen.analysis.pooled_processing import (
    DEFAULT_NT_REGEX,
    DEFAULT_POS_REGEX,
    DEFAULT_REFERENCE_REGEX,
    DEFAULT_TREATMENT_REGEX,
    PooledAnalysisConfig,
    compute_pooled_analysis,
    load_pooled_table,
    prepare_pooled_integrated_table,
    resolve_replicate_columns,
)
from prpcscreen.analysis.processing_data import create_hit_lists

DEBUG_ENV_DEFAULT = os.environ.get("PRPCSCREEN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}


def debug_log(message: str, enabled: bool) -> None:
    if enabled:
        print(f"[compute_pooled_metrics] {message}", file=sys.stderr)


def run_metrics_cli() -> None:
    parser = argparse.ArgumentParser(description="Analyze pooled-screen replicate counts and produce figure-compatible metrics.")
    parser.add_argument("input_table", help="CSV/TSV/XLSX input table with replicate columns (for example Negative_R1.. and Positive_R1..).")
    parser.add_argument("output_csv", help="Analyzed output CSV (figure-compatible schema).")
    parser.add_argument("--integrated_csv", default=None, help="Optional standardized integrated CSV output.")
    parser.add_argument("--hits_csv", default=None, help="Optional hit subset CSV output.")
    parser.add_argument("--sheet", default=None, help="Excel sheet to read when input is XLSX.")

    parser.add_argument("--reference-cols", nargs="+", default=None, help="Reference condition replicate columns.")
    parser.add_argument("--treatment-cols", nargs="+", default=None, help="Treatment condition replicate columns.")
    parser.add_argument("--reference-regex", default=DEFAULT_REFERENCE_REGEX, help="Regex used to auto-detect reference columns.")
    parser.add_argument("--treatment-regex", default=DEFAULT_TREATMENT_REGEX, help="Regex used to auto-detect treatment columns.")

    parser.add_argument("--pseudocount", type=float, default=1.0, help="Pseudocount added before log2 transform (default: 1.0).")
    parser.add_argument(
        "--pvalue-method",
        choices=("welch", "student", "paired"),
        default="welch",
        help="Row-wise test on log2-normalized counts: welch, student (equal variance), or paired.",
    )
    parser.add_argument("--p-cutoff", type=float, default=0.05, help="P-value cutoff for hit labeling.")
    parser.add_argument("--log2fc-cutoff", type=float, default=0.3, help="Absolute Mean_log2 cutoff for hit labeling.")
    parser.add_argument("--nt-regex", default=DEFAULT_NT_REGEX, help="Regex for non-targeting control inference from Gene_symbol.")
    parser.add_argument("--pos-regex", default=DEFAULT_POS_REGEX, help="Regex for positive-control inference from Gene_symbol.")
    parser.add_argument("--debug", action="store_true", help="Enable verbose debug logging.")
    args = parser.parse_args()

    debug_enabled = DEBUG_ENV_DEFAULT or args.debug

    input_path = Path(args.input_table)
    debug_log(f"Loading pooled input: {input_path}", debug_enabled)
    pooled_df, used_sheet = load_pooled_table(input_path, sheet=args.sheet)
    if used_sheet:
        debug_log(f"Loaded Excel sheet: {used_sheet}", debug_enabled)
    debug_log(f"Rows loaded: {len(pooled_df)}", debug_enabled)

    reference_cols, treatment_cols = resolve_replicate_columns(
        pooled_df,
        reference_cols=args.reference_cols,
        treatment_cols=args.treatment_cols,
        reference_regex=args.reference_regex,
        treatment_regex=args.treatment_regex,
    )
    debug_log(f"Reference columns: {reference_cols}", debug_enabled)
    debug_log(f"Treatment columns: {treatment_cols}", debug_enabled)

    integrated = prepare_pooled_integrated_table(
        pooled_df,
        reference_cols=reference_cols,
        treatment_cols=treatment_cols,
        nt_regex=args.nt_regex,
        pos_regex=args.pos_regex,
    )
    if args.integrated_csv:
        integrated_path = Path(args.integrated_csv)
        integrated_path.parent.mkdir(parents=True, exist_ok=True)
        integrated.to_csv(integrated_path, index=False)
        debug_log(f"Wrote integrated table: {integrated_path}", debug_enabled)

    config = PooledAnalysisConfig(
        reference_cols=tuple(reference_cols),
        treatment_cols=tuple(treatment_cols),
        pseudocount=float(args.pseudocount),
        pvalue_method=args.pvalue_method,
        p_cutoff=float(args.p_cutoff),
        log2fc_cutoff=float(args.log2fc_cutoff),
        nt_regex=args.nt_regex,
        pos_regex=args.pos_regex,
    )
    analyzed = compute_pooled_analysis(pooled_df, config)
    output_path = Path(args.output_csv)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    analyzed.to_csv(output_path, index=False)
    debug_log(f"Wrote analyzed table: {output_path} ({len(analyzed)} rows)", debug_enabled)

    if args.hits_csv:
        hits = create_hit_lists(
            analyzed,
            p_column="p_value_log2",
            effect_column="Mean_log2",
            p_cutoff=float(args.p_cutoff),
            log2fc_cutoff=float(args.log2fc_cutoff),
        )
        hits_path = Path(args.hits_csv)
        hits_path.parent.mkdir(parents=True, exist_ok=True)
        hits["hits"].to_csv(hits_path, index=False)
        debug_log(f"Wrote hits table: {hits_path} ({len(hits['hits'])} rows)", debug_enabled)


if __name__ == "__main__":
    run_metrics_cli()

