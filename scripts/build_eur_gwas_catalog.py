#!/usr/bin/env python3
"""Build a searchable Pan-UKBB EUR GWAS catalog for this pipeline."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def as_int(value: str) -> int:
    if value in ("", "NA"):
        return 0
    return int(float(value))


def first_present(row: dict[str, str], names: list[str]) -> str:
    for name in names:
        value = row.get(name, "")
        if value not in ("", "NA"):
            return value
    return "NA"


def public_url(path: str) -> str:
    prefix = "s3://pan-ukb-us-east-1/"
    if path.startswith(prefix):
        return "https://pan-ukb-us-east-1.s3.amazonaws.com/" + path[len(prefix) :]
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotype-manifest", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    keep = [
        "phenotype_id",
        "trait_type",
        "phenocode",
        "pheno_sex",
        "coding",
        "modifier",
        "description",
        "description_more",
        "coding_description",
        "category",
        "n_cases_EUR",
        "n_controls_EUR",
        "N_EUR_for_ldsc",
        "phenotype_qc_EUR",
        "h2_observed_EUR",
        "h2_observed_se_EUR",
        "h2_liability_EUR",
        "h2_liability_se_EUR",
        "h2_z_EUR",
        "lambda_gc_EUR",
        "in_max_independent_set",
        "filename",
        "aws_path",
        "public_url",
        "size_in_bytes",
    ]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(args.phenotype_manifest, "rt", newline="") as inp, args.out.open(
        "w", newline=""
    ) as out:
        reader = csv.DictReader(inp, delimiter="\t")
        writer = csv.DictWriter(out, fieldnames=keep, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            pops = set((row.get("pops") or "").split(","))
            if "EUR" not in pops:
                continue
            cases = as_int(row.get("n_cases_EUR", "0"))
            controls = as_int(row.get("n_controls_EUR", "0"))
            n_eur = cases + controls if controls else cases
            out_row = dict(row)
            out_row["phenotype_id"] = row["filename"].removesuffix(".tsv.bgz")
            out_row["N_EUR_for_ldsc"] = str(n_eur)
            out_row["h2_observed_EUR"] = first_present(
                row, ["sldsc_25bin_h2_observed_EUR", "rhemc_25bin_50rv_h2_observed_EUR"]
            )
            out_row["h2_observed_se_EUR"] = first_present(
                row,
                [
                    "sldsc_25bin_h2_observed_se_EUR",
                    "rhemc_25bin_50rv_h2_observed_se_EUR",
                ],
            )
            out_row["h2_liability_EUR"] = first_present(
                row, ["sldsc_25bin_h2_liability_EUR", "rhemc_25bin_50rv_h2_liability_EUR"]
            )
            out_row["h2_liability_se_EUR"] = first_present(
                row,
                [
                    "sldsc_25bin_h2_liability_se_EUR",
                    "rhemc_25bin_50rv_h2_liability_se_EUR",
                ],
            )
            out_row["h2_z_EUR"] = first_present(
                row, ["sldsc_25bin_h2_z_EUR", "rhemc_25bin_50rv_h2_z_EUR"]
            )
            out_row["public_url"] = public_url(row["aws_path"])
            writer.writerow(out_row)
            n += 1
    print(f"wrote {n} EUR GWAS rows to {args.out}")


if __name__ == "__main__":
    main()
