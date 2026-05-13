#!/usr/bin/env python3
"""Filter a catalog/manifest to selected phenotype IDs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--ids", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    ids = {line.strip() for line in args.ids.read_text().splitlines() if line.strip()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest.open() as inp, args.out.open("w", newline="") as out:
        reader = csv.DictReader(inp, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        n = 0
        for row in reader:
            if row["phenotype_id"] in ids:
                writer.writerow(row)
                n += 1
    print(f"wrote {n} rows to {args.out}")


if __name__ == "__main__":
    main()

