from __future__ import annotations

import os
import sys
import argparse
import pandas as pd

# Ensure package imports work when this file is executed directly.
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from prpcscreen.analysis.processing_data import create_hit_lists, run_ssmd_stats
from prpcscreen.analysis.calculating_scores import NORMALIZATION_METHODS

DEBUG_ENV_DEFAULT = os.environ.get("PRPCSCREEN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}


def debug_log(message: str, enabled: bool) -> None:
    """Emit debug lines for metric computation flow."""
    if enabled:
        print(f"[03_compute_screen_metrics] {message}", file=sys.stderr)


def run_metrics_cli() -> None:
    # Parse input/output paths and optional hit export.
    parser = argparse.ArgumentParser(description="Analyze integrated PrP screen data.")
    parser.add_argument("input_csv")
    parser.add_argument("output_csv")
    parser.add_argument("--hits_csv", default=None)
    parser.add_argument(
        "--norm-method",
        default="genes and all Non-targeting",
        choices=NORMALIZATION_METHODS,
        help="Per-plate normalization baseline strategy (default: genes and all Non-targeting).",
    )
    parser.add_argument("--debug", action="store_true", help="Enable verbose debug logging.")
    args = parser.parse_args()
    debug_enabled = DEBUG_ENV_DEFAULT or args.debug

    # Read integrated screen data and coerce replicate columns to numeric ranges.
    df = pd.read_csv(args.input_csv)
    debug_log(f"Loaded integrated data: {args.input_csv} ({len(df)} rows)", debug_enabled)
    for col in ("Raw_rep1", "Raw_rep2"):
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="coerce").clip(lower=1)
            debug_log(f"Normalized numeric safety floor for {col}", debug_enabled)

    # Run primary analysis routine (normalization + SSMD + p-value derivation).
    debug_log(f"Normalization method: {args.norm_method}", debug_enabled)
    analyzed = run_ssmd_stats(df, norm_method=args.norm_method)
    analyzed.to_csv(args.output_csv, index=False)
    debug_log(f"Wrote analyzed dataset: {args.output_csv} ({len(analyzed)} rows)", debug_enabled)

    # Optionally materialize filtered hit calls for downstream review.
    if args.hits_csv:
        hits = create_hit_lists(analyzed, p_column="p_value_log2", effect_column="Mean_log2")
        hits["hits"].to_csv(args.hits_csv, index=False)
        debug_log(f"Wrote hits dataset: {args.hits_csv} ({len(hits['hits'])} rows)", debug_enabled)


if __name__ == "__main__":
    run_metrics_cli()

