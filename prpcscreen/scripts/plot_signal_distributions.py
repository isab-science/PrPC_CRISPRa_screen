from __future__ import annotations

import os
import sys
import argparse
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from prpcscreen.visualization.histograms import write_interactive_histogram_html

DEBUG_ENV_DEFAULT = os.environ.get("PRPCSCREEN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}


def debug_log(message: str, enabled: bool) -> None:
    """Emit debug lines for histogram plotting."""
    if enabled:
        print(f"[07_plot_signal_distributions] {message}", file=sys.stderr)


def run_distribution_cli() -> None:
    # Parse source and selected metric for interactive histograms.
    parser = argparse.ArgumentParser(description="Generate interactive PrP histograms.")
    parser.add_argument("input_csv")
    parser.add_argument(
        "--output_html",
        required=True,
        help="Interactive histogram HTML with category on/off toggles and publication-quality PNG export.",
    )
    parser.add_argument("--column", default="Log2FC_rep1")
    parser.add_argument(
        "--genomics_excel",
        default=None,
        help="Optional genomics workbook used to preload sublibrary metadata for histogram filtering.",
    )
    parser.add_argument("--debug", action="store_true", help="Enable verbose debug logging.")
    args = parser.parse_args()
    debug_enabled = DEBUG_ENV_DEFAULT or args.debug

    # Generate interactive histogram from analyzed data.
    df = pd.read_csv(args.input_csv)
    debug_log(f"Loaded analyzed data: {args.input_csv} ({len(df)} rows)", debug_enabled)
    debug_log(f"Histogram base column: {args.column}", debug_enabled)
    html_path = write_interactive_histogram_html(
        df,
        args.output_html,
        use_column=args.column,
        genomics_excel=args.genomics_excel,
    )
    debug_log(f"Wrote interactive histogram HTML: {html_path}", debug_enabled)


if __name__ == "__main__":
    run_distribution_cli()
