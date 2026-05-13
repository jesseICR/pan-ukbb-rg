#!/usr/bin/env python3
"""Validate the derived EUR catalog against the source Pan phenotype manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def as_int(value: str) -> int:
    if value in ("", "NA"):
        return 0
    return int(float(value))


def n_eur(row: dict[str, str]) -> str:
    cases = as_int(row.get("n_cases_EUR", "0"))
    controls = as_int(row.get("n_controls_EUR", "0"))
    return str(cases + controls if controls else cases)


def phenotype_id(row: dict[str, str]) -> str:
    return row["filename"].removesuffix(".tsv.bgz")


def read_source(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with gzip.open(path, "rt", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if "EUR" not in set((row.get("pops") or "").split(",")):
                continue
            rows[phenotype_id(row)] = row
    return rows


def read_catalog(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as f:
        return {row["phenotype_id"]: row for row in csv.DictReader(f, delimiter="\t")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotype-manifest", required=True, type=Path)
    parser.add_argument("--catalog", required=True, type=Path)
    args = parser.parse_args()

    source = read_source(args.phenotype_manifest)
    catalog = read_catalog(args.catalog)
    errors: list[str] = []

    missing = sorted(set(source) - set(catalog))
    extra = sorted(set(catalog) - set(source))
    if missing:
        errors.append(f"catalog missing {len(missing)} source EUR rows; first={missing[:5]}")
    if extra:
        errors.append(f"catalog has {len(extra)} extra rows; first={extra[:5]}")

    checks = [
        ("phenocode", lambda s: s["phenocode"], lambda c: c["phenocode"]),
        ("trait_type", lambda s: s["trait_type"], lambda c: c["trait_type"]),
        ("pheno_sex", lambda s: s["pheno_sex"], lambda c: c["pheno_sex"]),
        ("coding", lambda s: s["coding"], lambda c: c["coding"]),
        ("modifier", lambda s: s["modifier"], lambda c: c["modifier"]),
        ("description", lambda s: s["description"], lambda c: c["description"]),
        ("phenotype_qc_EUR", lambda s: s["phenotype_qc_EUR"], lambda c: c["phenotype_qc_EUR"]),
        ("aws_path", lambda s: s["aws_path"], lambda c: c["aws_path"]),
        ("filename", lambda s: s["filename"], lambda c: c["filename"]),
        ("N_EUR_for_ldsc", n_eur, lambda c: c["N_EUR_for_ldsc"]),
    ]
    mismatches = 0
    for pid in sorted(set(source) & set(catalog)):
        s = source[pid]
        c = catalog[pid]
        for name, s_get, c_get in checks:
            if s_get(s) != c_get(c):
                mismatches += 1
                if mismatches <= 20:
                    errors.append(
                        f"mismatch {pid} {name}: source={s_get(s)!r} catalog={c_get(c)!r}"
                    )

    print(f"source_eur_rows\t{len(source)}")
    print(f"catalog_rows\t{len(catalog)}")
    print(f"mismatches\t{mismatches}")
    if errors:
        for error in errors:
            print(f"ERROR\t{error}")
        raise SystemExit(1)
    print("validation\tPASS")


if __name__ == "__main__":
    main()

