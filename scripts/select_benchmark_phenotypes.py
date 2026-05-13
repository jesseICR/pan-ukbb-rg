#!/usr/bin/env python3
"""Select a fixed-size stratified benchmark set from the Pan-UKBB manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import random
from pathlib import Path


def parse_config(path: Path) -> tuple[int, dict[str, int]]:
    seed = 20260513
    groups: dict[str, int] = {}
    in_groups = False
    with path.open() as f:
        for raw in f:
            line = raw.rstrip()
            if not line or line.lstrip().startswith("#"):
                continue
            if line.startswith("seed:"):
                seed = int(line.split(":", 1)[1].strip())
            elif line.startswith("qc_groups:"):
                in_groups = True
            elif in_groups and raw.startswith("  "):
                key, value = line.split(":", 1)
                groups[key.strip()] = int(value.strip())
            elif in_groups:
                in_groups = False
    if not groups:
        raise ValueError(f"No qc_groups found in {path}")
    return seed, groups


def as_int(value: str) -> int:
    if value in ("", "NA"):
        return 0
    return int(float(value))


def public_url(path: str) -> str:
    prefix = "s3://pan-ukb-us-east-1/"
    if path.startswith(prefix):
        return "https://pan-ukb-us-east-1.s3.amazonaws.com/" + path[len(prefix) :]
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    seed, groups = parse_config(args.config)
    rng = random.Random(seed)

    by_qc: dict[str, list[dict[str, str]]] = {qc: [] for qc in groups}
    with gzip.open(args.manifest, "rt", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pops = set((row.get("pops") or "").split(","))
            if "EUR" not in pops:
                continue
            qc = row["phenotype_qc_EUR"]
            if qc not in by_qc:
                continue
            cases = as_int(row.get("n_cases_EUR", "0"))
            controls = as_int(row.get("n_controls_EUR", "0"))
            n = cases + controls if controls else cases
            row = dict(row)
            row["N_EUR_for_ldsc"] = str(n)
            row["phenotype_id"] = row["filename"].removesuffix(".tsv.bgz")
            row["public_url"] = public_url(row["aws_path"])
            by_qc[qc].append(row)

    selected: list[dict[str, str]] = []
    for qc, n in groups.items():
        rows = by_qc[qc]
        if len(rows) < n:
            raise ValueError(f"QC group {qc} has only {len(rows)} rows, need {n}")
        rows = sorted(rows, key=lambda r: r["filename"])
        selected.extend(rng.sample(rows, n))

    selected = sorted(selected, key=lambda r: (r["phenotype_qc_EUR"], r["filename"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    keep = [
        "phenotype_id",
        "phenotype_qc_EUR",
        "trait_type",
        "phenocode",
        "pheno_sex",
        "coding",
        "modifier",
        "description",
        "N_EUR_for_ldsc",
        "sldsc_25bin_h2_z_EUR",
        "aws_path",
        "public_url",
        "size_in_bytes",
        "filename",
    ]
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keep, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)

    for qc, n in groups.items():
        got = sum(1 for row in selected if row["phenotype_qc_EUR"] == qc)
        print(f"{qc}\t{got}")
    total_bytes = sum(as_int(row["size_in_bytes"]) for row in selected)
    print(f"selected\t{len(selected)}")
    print(f"compressed_input_gib\t{total_bytes / 1024**3:.2f}")


if __name__ == "__main__":
    main()
