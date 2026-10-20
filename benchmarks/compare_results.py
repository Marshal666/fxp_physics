#!/usr/bin/env python3
"""Compare deterministic benchmark outputs, ignoring machine-dependent timings."""

import argparse
import csv
from pathlib import Path

KEY = ("scene", "bodies", "ticks", "warmup", "queries", "sample")
STATE = (
    "peak_pairs",
    "peak_contacts",
    "peak_awake",
    "checksum",
    "query_checksum",
    "query_hits",
)


def read_results(path):
    rows = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        missing = set(KEY + STATE) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for row in reader:
            key = tuple(row[field] for field in KEY)
            if key in rows:
                raise ValueError(f"{path}: duplicate scenario/sample {key}")
            rows[key] = tuple(int(row[field]) for field in STATE)
    if not rows:
        raise ValueError(f"{path}: no measurements")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="directory containing downloaded CSV artifacts")
    parser.add_argument("--expected-files", type=int, default=8)
    args = parser.parse_args()
    paths = sorted(args.directory.rglob("*.csv"))
    if len(paths) != args.expected_files or len(paths) < 2:
        parser.error(f"expected {args.expected_files} CSV files (at least two), found {len(paths)}")
    try:
        baseline = read_results(paths[0])
        for path in paths[1:]:
            actual = read_results(path)
            if actual.keys() != baseline.keys():
                raise ValueError(f"{path}: scenario/sample set differs from {paths[0]}")
            for key, expected in baseline.items():
                for field, want, got in zip(STATE, expected, actual[key]):
                    if want != got:
                        raise ValueError(f"{path}: {key}: {field} is {got}; "
                                         f"{paths[0]} has {want}")
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))
    print(f"{len(paths)} platform/backend results matched {len(baseline)} scenario/sample rows exactly.")


if __name__ == "__main__":
    main()
