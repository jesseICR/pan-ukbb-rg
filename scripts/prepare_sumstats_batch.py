#!/usr/bin/env python3
"""Prepare benchmark LDSC sumstats files in parallel."""

from __future__ import annotations

import argparse
import csv
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def public_url(path: str) -> str:
    prefix = "s3://pan-ukb-us-east-1/"
    if path.startswith(prefix):
        return "https://pan-ukb-us-east-1.s3.amazonaws.com/" + path[len(prefix) :]
    return path


def run_one(row: dict[str, str], args: argparse.Namespace) -> str:
    phenotype_id = row["phenotype_id"]
    out = args.out_dir / f"{phenotype_id}.sumstats.gz"
    stats_out = args.stats_dir / f"{phenotype_id}.json"
    if out.exists() and stats_out.exists() and not args.force:
        return f"skip\t{phenotype_id}"

    url = row.get("public_url") or public_url(row["aws_path"])
    cmd = (
        f"set -o pipefail; "
        f"curl -fL --retry 5 --retry-delay 5 --silent --show-error {shlex.quote(url)} | "
        f"{shlex.quote(sys.executable)} scripts/convert_pan_flat_to_ldsc.py "
        f"--input - "
        f"--ld-snps {shlex.quote(str(args.ld_snps))} "
        f"--phenotype-id {shlex.quote(phenotype_id)} "
        f"--qc {shlex.quote(row['phenotype_qc_EUR'])} "
        f"--n {shlex.quote(row['N_EUR_for_ldsc'])} "
        f"--out {shlex.quote(str(out))} "
        f"--stats-out {shlex.quote(str(stats_out))}"
    )
    log_path = args.log_dir / f"{phenotype_id}.log"
    with log_path.open("w") as log:
        proc = subprocess.run(["bash", "-lc", cmd], stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"{phenotype_id} failed; see {log_path}")
    return f"done\t{phenotype_id}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotypes", required=True, type=Path)
    parser.add_argument("--ld-snps", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--stats-dir", required=True, type=Path)
    parser.add_argument("--log-dir", default=Path("logs/prepare_sumstats"), type=Path)
    parser.add_argument("--jobs", default=16, type=int)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.stats_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    with args.phenotypes.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, row, args) for row in rows]
        for future in as_completed(futures):
            print(future.result(), flush=True)


if __name__ == "__main__":
    main()
