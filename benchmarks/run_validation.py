#!/usr/bin/env python3
"""Measure bounded CI workloads with one and four workers, including 20 Hz stacks."""

import argparse
import csv
import os
from pathlib import Path
import subprocess


def workloads():
    yield ["--suite", "scaling", "--hz", "20", "--ticks", "8"]
    yield ["--scene", "playground", "--ticks", "240"]
    for phase in ("settled", "disturbed"):
        yield [
            "--scene", "stack", "--bodies", "100", "--hz", "20",
            "--phase", phase, "--warmup", "600", "--ticks", "40",
        ]
    yield ["--scene", "pile", "--bodies", "128", "--hz", "20", "--ticks", "40"]
    yield [
        "--scene", "pile", "--bodies", "64", "--hz", "20", "--substeps", "2",
        "--iterations", "8", "--ticks", "40", "--queries", "4",
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    executable = args.executable
    if os.name == "nt" and not executable.suffix:
        executable = executable.with_suffix(".exe")
    executable = executable.resolve(strict=True)
    args.output.mkdir(parents=True, exist_ok=True)
    for workers in (1, 4):
        path = args.output / f"workers-{workers}.csv"
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = None
            for workload in workloads():
                command = [
                    str(executable), "--workers", str(workers), "--samples", "1",
                    *workload,
                ]
                print(" ".join(command), flush=True)
                result = subprocess.check_output(command, text=True)
                reader = csv.DictReader(result.splitlines())
                if writer is None:
                    writer = csv.DictWriter(stream, fieldnames=reader.fieldnames)
                    writer.writeheader()
                elif reader.fieldnames != writer.fieldnames:
                    raise ValueError("inconsistent benchmark CSV columns")
                rows = list(reader)
                if not rows:
                    raise ValueError(f"benchmark returned no rows: {command}")
                writer.writerows(rows)
                stream.flush()
        print(f"Wrote {path}", flush=True)


if __name__ == "__main__":
    main()
