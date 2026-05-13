#!/usr/bin/env python3
"""Convert one Pan-UKBB flat file stream to a minimal LDSC sumstats file."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sys
import time
from pathlib import Path


NA = {"", "NA", "nan", "NaN"}


def load_snps(path: Path) -> set[str]:
    with path.open() as f:
        return {line.strip() for line in f if line.strip()}


def open_text_input(path: str):
    if path == "-":
        return gzip.open(sys.stdin.buffer, "rt")
    return gzip.open(path, "rt")


def is_missing(value: str) -> bool:
    return value in NA


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="-", help="Input .tsv.bgz path or '-' for stdin")
    parser.add_argument("--ld-snps", required=True, type=Path)
    parser.add_argument("--phenotype-id", required=True)
    parser.add_argument("--qc", required=True)
    parser.add_argument("--n", required=True, type=float)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--stats-out", required=True, type=Path)
    args = parser.parse_args()

    started = time.time()
    ld_snps = load_snps(args.ld_snps)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.stats_out.parent.mkdir(parents=True, exist_ok=True)
    out_tmp = args.out.with_name(f"{args.out.name}.tmp.{os.getpid()}")
    stats_tmp = args.stats_out.with_name(f"{args.stats_out.name}.tmp.{os.getpid()}")

    n_in = 0
    n_ld = 0
    n_written = 0
    n_missing = 0
    n_low_conf = 0
    n_bad_se = 0

    try:
        with open_text_input(args.input) as inp, gzip.open(out_tmp, "wt", compresslevel=6) as out:
            reader = csv.DictReader(inp, delimiter="\t")
            required = [
                "chr",
                "pos",
                "ref",
                "alt",
                "beta_EUR",
                "se_EUR",
                "low_confidence_EUR",
            ]
            missing_cols = [c for c in required if c not in reader.fieldnames]
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")

            out.write("SNP\tA1\tA2\tZ\tN\n")
            for row in reader:
                n_in += 1
                snp = f"{row['chr']}:{row['pos']}:{row['ref']}:{row['alt']}"
                if snp not in ld_snps:
                    continue
                n_ld += 1

                if row["low_confidence_EUR"].lower() == "true":
                    n_low_conf += 1
                    continue

                beta = row["beta_EUR"]
                se = row["se_EUR"]
                if is_missing(beta) or is_missing(se):
                    n_missing += 1
                    continue
                beta_f = float(beta)
                se_f = float(se)
                if not math.isfinite(beta_f) or not math.isfinite(se_f) or se_f <= 0:
                    n_bad_se += 1
                    continue
                z = beta_f / se_f
                if not math.isfinite(z):
                    n_bad_se += 1
                    continue
                out.write(f"{snp}\t{row['alt']}\t{row['ref']}\t{z:.8g}\t{args.n:.0f}\n")
                n_written += 1

        elapsed = time.time() - started
        stats = {
            "phenotype_id": args.phenotype_id,
            "qc": args.qc,
            "n_eur": args.n,
            "input_rows": n_in,
            "ld_snps_seen": n_ld,
            "written_rows": n_written,
            "missing_eur": n_missing,
            "low_confidence_eur": n_low_conf,
            "bad_se": n_bad_se,
            "elapsed_seconds": elapsed,
            "rows_per_second": n_in / elapsed if elapsed else None,
            "out": str(args.out),
        }
        with stats_tmp.open("w") as f:
            json.dump(stats, f, indent=2, sort_keys=True)
            f.write("\n")
        os.replace(out_tmp, args.out)
        os.replace(stats_tmp, args.stats_out)
    finally:
        for tmp in [out_tmp, stats_tmp]:
            if tmp.exists():
                tmp.unlink()
    print(json.dumps(stats, sort_keys=True))


if __name__ == "__main__":
    main()
