from __future__ import annotations

import os
import sys
import argparse
import pandas as pd

# Ensure package imports work when this file is executed directly.
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from prpcscreen.misc.converting_plate_layouts import convert_well_numbers

DEBUG_ENV_DEFAULT = os.environ.get("PRPCSCREEN_DEBUG", "").strip().lower() in {"1", "true", "yes", "on"}


def debug_log(message: str, enabled: bool) -> None:
    """Emit a debug line to stderr when debug mode is enabled."""
    if enabled:
        print(f"[01_remap_plate_coordinates] {message}", file=sys.stderr)


def run_remap_cli() -> None:
    # Parse CLI inputs for source and destination layout tables.
    parser = argparse.ArgumentParser(description="Convert 384-well numbers to 96-well numbering.")
    parser.add_argument("input_csv", help="CSV with Well_number_384 column")
    parser.add_argument("output_csv", help="Output CSV path")
    parser.add_argument("--debug", action="store_true", help="Enable verbose debug logging.")
    args = parser.parse_args()
    debug_enabled = DEBUG_ENV_DEFAULT or args.debug

    # Load input layout and verify row volume for quick sanity checks.
    df = pd.read_csv(args.input_csv)
    debug_log(f"Loaded input layout: {args.input_csv} ({len(df)} rows)", debug_enabled)

    # Perform coordinate remapping from 384 well numbering to 96 mapping.
    df["Well_number_96"] = convert_well_numbers(df["Well_number_384"].to_numpy())

    # Persist transformed table for downstream plate-processing workflows.
    df.to_csv(args.output_csv, index=False)
    debug_log(f"Wrote remapped layout: {args.output_csv}", debug_enabled)


if __name__ == "__main__":
    run_remap_cli()

