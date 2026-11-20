#!/usr/bin/env python3
"""
Convert a plain (x, y, cell_type) CSV into Einzell's dataset JSON format.

Usage:
    python convert.py input.csv output.json
"""
import argparse
import json
import sys

import pandas as pd

from einzell_convert.core import ConversionError, compute_dataset_output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", help="CSV with at least x, y, cell_type columns")
    parser.add_argument("output_json", help="path to write the Einzell-compatible dataset JSON")
    args = parser.parse_args()

    df = pd.read_csv(args.input_csv)
    try:
        output_map = compute_dataset_output(df)
    except ConversionError as e:
        print(f"error: {e}", file=sys.stderr)
        sys.exit(1)

    with open(args.output_json, "w") as f:
        # allow_nan=False: fail loudly here rather than silently writing an
        # Infinity/NaN token that only breaks later when something tries to
        # JSON-decode it strictly (exactly the bug just fixed in
        # einzell-web/census_api.py's own dataset pipeline).
        json.dump(output_map, f, allow_nan=False)
    print(f"wrote {len(output_map)} cells to {args.output_json}")


if __name__ == "__main__":
    main()
