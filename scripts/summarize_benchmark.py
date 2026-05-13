#!/usr/bin/env python3
"""Summarize benchmark outputs and extrapolate to all Pan-UKBB EUR GWAS pairs."""

from __future__ import annotations

import argparse
import csv
import glob
import json
import statistics
from datetime import datetime
from pathlib import Path


FULL_EUR_GWAS = 7160
FULL_EUR_PAIRS = FULL_EUR_GWAS * (FULL_EUR_GWAS - 1) // 2


def read_prepare_stats(path: Path) -> list[dict[str, object]]:
    rows = []
    for fn in sorted(glob.glob(str(path / "*.json"))):
        with open(fn) as f:
            rows.append(json.load(f))
    return rows


def read_driver_summary(path: Path) -> list[dict[str, str]]:
    summary = path / "driver_summary.tsv"
    if not summary.exists():
        return []
    with summary.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def count_result_rows(path: Path) -> int:
    n = 0
    for fn in glob.glob(str(path / "*.r2")):
        with open(fn) as f:
            header = next(f, None)
            for line in f:
                if line.strip():
                    n += 1
    return n


def wall_seconds_from_files(path: Path) -> float | None:
    starts = [Path(fn).stat().st_mtime for fn in glob.glob(str(path / "*.rg-list.txt"))]
    ends = [Path(fn).stat().st_mtime for fn in glob.glob(str(path / "*.r2"))]
    if not starts or not ends:
        return None
    return max(ends) - min(starts)


def fmt_hours(seconds: float) -> str:
    return f"{seconds / 3600:.2f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotypes", required=True, type=Path)
    parser.add_argument("--prepare-stats", required=True, type=Path)
    parser.add_argument("--rg-dir", required=True, type=Path)
    parser.add_argument("--out-md", required=True, type=Path)
    parser.add_argument("--out-tsv", required=True, type=Path)
    args = parser.parse_args()

    with args.phenotypes.open() as f:
        phenotypes = list(csv.DictReader(f, delimiter="\t"))
    n_phenos = len(phenotypes)
    benchmark_pairs = n_phenos * (n_phenos - 1) // 2

    prep = read_prepare_stats(args.prepare_stats)
    driver = read_driver_summary(args.rg_dir)
    done = [r for r in driver if r["status"] in ("done", "skip_done")]
    measured_done = [r for r in driver if r["status"] == "done"]
    elapsed = sum(float(r["elapsed_seconds"]) for r in measured_done)
    pairs_run = sum(int(r["pairs"]) for r in done)
    rg_rows = count_result_rows(args.rg_dir)
    wall_seconds = wall_seconds_from_files(args.rg_dir)
    sec_per_pair = elapsed / sum(int(r["pairs"]) for r in measured_done) if measured_done else None
    full_seconds_one_worker = sec_per_pair * FULL_EUR_PAIRS if sec_per_pair else None

    prep_elapsed = sum(float(r["elapsed_seconds"]) for r in prep)
    prep_by_qc: dict[str, int] = {}
    for row in phenotypes:
        prep_by_qc[row["phenotype_qc_EUR"]] = prep_by_qc.get(row["phenotype_qc_EUR"], 0) + 1

    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w") as f:
        f.write("metric\tvalue\n")
        f.write(f"generated_at\t{datetime.now().isoformat()}\n")
        f.write(f"benchmark_phenotypes\t{n_phenos}\n")
        f.write(f"benchmark_pairs\t{benchmark_pairs}\n")
        f.write(f"ldsc_pairs_completed\t{pairs_run}\n")
        f.write(f"result_rows_written\t{rg_rows}\n")
        f.write(f"ldsc_elapsed_seconds_measured\t{elapsed:.6f}\n")
        if wall_seconds is not None:
            f.write(f"ldsc_wall_seconds\t{wall_seconds:.6f}\n")
            f.write(f"ldsc_wall_hours\t{wall_seconds / 3600:.4f}\n")
        if sec_per_pair is not None:
            f.write(f"ldsc_seconds_per_pair\t{sec_per_pair:.8f}\n")
            f.write(f"full_eur_pairs\t{FULL_EUR_PAIRS}\n")
            f.write(f"full_eur_single_worker_hours\t{full_seconds_one_worker / 3600:.4f}\n")
        f.write(f"prepare_sumstats_elapsed_seconds_sum\t{prep_elapsed:.6f}\n")

    lines = [
        "# Benchmark 90 Summary",
        "",
        f"Generated: {datetime.now().isoformat()}",
        "",
        "## Phenotypes",
        "",
        f"- Benchmark phenotypes: {n_phenos}",
        f"- Benchmark unordered pairs: {benchmark_pairs:,}",
        f"- QC group counts: {prep_by_qc}",
        "",
        "## LDSC Runtime",
        "",
        f"- LDSC pairs completed: {pairs_run:,}",
        f"- `.r2` result rows written: {rg_rows:,}",
        f"- Measured LDSC elapsed seconds across completed lead jobs: {elapsed:,.1f}",
    ]
    if wall_seconds is not None:
        lines.append(f"- Approximate local wall-clock runtime: {fmt_hours(wall_seconds)} hours")
    if sec_per_pair is not None:
        lines.extend(
            [
                f"- Effective seconds per pair: {sec_per_pair:.6f}",
                f"- Full all-EUR pairs: {FULL_EUR_PAIRS:,}",
                f"- Extrapolated single-worker LDSC hours: {fmt_hours(full_seconds_one_worker)}",
                "",
                "Parallel wall time is approximately this single-worker total divided by the number of independent LDSC processes that fit in RAM and IO.",
            ]
        )
    else:
        lines.append("- No measured LDSC jobs found yet.")

    if prep:
        written = [int(r["written_rows"]) for r in prep]
        lines.extend(
            [
                "",
                "## Sumstats Preparation",
                "",
                f"- Sum of per-file preparation seconds: {prep_elapsed:,.1f}",
                f"- Median output SNPs per phenotype: {statistics.median(written):,.0f}",
                f"- Min/max output SNPs per phenotype: {min(written):,} / {max(written):,}",
            ]
        )

    args.out_md.write_text("\n".join(lines) + "\n")
    print(args.out_md.read_text())


if __name__ == "__main__":
    main()
