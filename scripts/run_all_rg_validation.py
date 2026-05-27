#!/usr/bin/env python3
"""Validate rg-batch against the existing Rust ldsc rg command on random pairs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import statistics
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


THREAD_ENV = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


def read_traits(manifest: Path, sumstats_dir: Path) -> list[tuple[int, str]]:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "phenotype_id" not in (reader.fieldnames or []):
            raise SystemExit(f"{manifest} does not contain a phenotype_id column")
        rows = [(idx, row["phenotype_id"]) for idx, row in enumerate(reader) if row.get("phenotype_id")]
    missing = [trait for _, trait in rows if not (sumstats_dir / f"{trait}.sumstats.gz").exists()]
    if missing:
        preview = ", ".join(missing[:10])
        raise SystemExit(f"{len(missing)} traits are missing prepared sumstats; run make setup. First: {preview}")
    return rows


def select_disjoint_pairs(
    traits: list[tuple[int, str]], n_pairs: int, seed: int
) -> list[dict[str, str]]:
    if len(traits) < 2 * n_pairs:
        raise SystemExit(f"Need at least {2 * n_pairs} traits for {n_pairs} disjoint pairs")
    rng = random.Random(seed)
    shuffled = traits[:]
    rng.shuffle(shuffled)
    out = []
    for pair_idx in range(n_pairs):
        (i, p1) = shuffled[2 * pair_idx]
        (j, p2) = shuffled[2 * pair_idx + 1]
        if j < i:
            i, j = j, i
            p1, p2 = p2, p1
        out.append(
            {
                "pair_id": str(pair_idx + 1),
                "i": str(i),
                "j": str(j),
                "p1": p1,
                "p2": p2,
                "lead": p1,
                "target": p2,
            }
        )
    return out


def write_pairs(path: Path, pairs: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["pair_id", "i", "j", "p1", "p2", "lead", "target"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(pairs)
    tmp.replace(path)


def resolve_m(ld_prefix: Path, m_override: str) -> str:
    if m_override != "auto":
        return m_override
    path = Path(str(ld_prefix) + ".l2.M_5_50")
    if not path.exists():
        raise SystemExit(f"M override is auto but M file does not exist: {path}")
    with path.open() as handle:
        return handle.read().split()[0]


def run_hybrid(args: argparse.Namespace, pairs_path: Path, m_snps: str) -> tuple[Path, float]:
    out = args.out_dir / "hybrid_rg.tsv"
    cmd = [
        str(args.ldsc_bin),
        "rg-batch",
        "--rayon-threads",
        str(args.hybrid_threads),
        "--pairs",
        str(pairs_path),
        "--sumstats-dir",
        str(args.sumstats_dir),
        "--ref-ld",
        str(Path(str(args.ld_prefix) + ".l2.ldscore.gz")),
        "--w-ld",
        str(Path(str(args.ld_prefix) + ".l2.ldscore.gz")),
        "--M",
        m_snps,
        "--out",
        str(out),
        "--verbose-timing",
    ]
    log = args.out_dir / "hybrid_rg.log"
    started = time.time()
    with log.open("w") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT)
    elapsed = time.time() - started
    if proc.returncode != 0:
        raise RuntimeError(f"hybrid rg-batch failed; see {log}")
    return out, elapsed


def parse_rust_rg_summary(stdout_path: Path) -> dict[str, str]:
    lines = stdout_path.read_text().splitlines()
    for line in reversed(lines):
        parts = line.split()
        if len(parts) >= 12 and parts[0].endswith(".sumstats.gz") and parts[1].endswith(".sumstats.gz"):
            return {
                "p1": Path(parts[0]).name.removesuffix(".sumstats.gz"),
                "p2": Path(parts[1]).name.removesuffix(".sumstats.gz"),
                "rg": parts[2],
                "se": parts[3],
                "z": parts[4],
                "p": parts[5],
                "h2_obs": parts[6],
                "h2_obs_se": parts[7],
                "h2_int": parts[8],
                "h2_int_se": parts[9],
                "gcov_int": parts[10],
                "gcov_int_se": parts[11],
            }
    raise RuntimeError(f"could not parse Rust rg summary from {stdout_path}")


def run_baseline_one(
    args: argparse.Namespace, pair: dict[str, str], m_snps: str
) -> tuple[dict[str, str], float]:
    pair_id = pair["pair_id"]
    log = args.out_dir / "baseline_logs" / f"pair_{int(pair_id):04d}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    out_prefix = args.out_dir / "baseline_logs" / f"pair_{int(pair_id):04d}"
    p1 = args.sumstats_dir / f"{pair['p1']}.sumstats.gz"
    p2 = args.sumstats_dir / f"{pair['p2']}.sumstats.gz"
    cmd = [
        str(args.ldsc_bin),
        "rg",
        "--rayon-threads",
        str(args.baseline_threads),
        "--rg",
        f"{p1},{p2}",
        "--ref-ld",
        str(Path(str(args.ld_prefix) + ".l2.ldscore.gz")),
        "--w-ld",
        str(Path(str(args.ld_prefix) + ".l2.ldscore.gz")),
        "--M",
        m_snps,
        "--out",
        str(out_prefix),
    ]
    env = os.environ.copy()
    env.update(THREAD_ENV)
    started = time.time()
    with log.open("w") as handle:
        proc = subprocess.run(cmd, stdout=handle, stderr=subprocess.STDOUT, env=env)
    elapsed = time.time() - started
    if proc.returncode != 0:
        raise RuntimeError(f"baseline rg failed for pair {pair_id}; see {log}")
    row = parse_rust_rg_summary(log)
    row.update(
        {
            "pair_id": pair_id,
            "requested_p1": pair["p1"],
            "requested_p2": pair["p2"],
            "elapsed_seconds": f"{elapsed:.6f}",
        }
    )
    return row, elapsed


def run_baseline(
    args: argparse.Namespace, pairs: list[dict[str, str]], m_snps: str
) -> tuple[Path, float]:
    started = time.time()
    rows: list[dict[str, str]] = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_baseline_one, args, pair, m_snps) for pair in pairs]
        for future in as_completed(futures):
            row, _elapsed = future.result()
            rows.append(row)
            print(f"baseline_done\t{row['pair_id']}", flush=True)
    rows.sort(key=lambda row: int(row["pair_id"]))
    out = args.out_dir / "baseline_rg.tsv"
    with out.open("w", newline="") as handle:
        fields = [
            "pair_id",
            "requested_p1",
            "requested_p2",
            "p1",
            "p2",
            "rg",
            "se",
            "z",
            "p",
            "h2_obs",
            "h2_obs_se",
            "h2_int",
            "h2_int_se",
            "gcov_int",
            "gcov_int_se",
            "elapsed_seconds",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return out, time.time() - started


def to_float(value: str) -> float:
    try:
        out = float(value)
    except ValueError:
        return math.nan
    return out if math.isfinite(out) else math.nan


def read_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        return {row["pair_id"]: row for row in csv.DictReader(handle, delimiter="\t")}


def r2(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return math.nan
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    ssx = sum((x - mx) ** 2 for x in xs)
    ssy = sum((y - my) ** 2 for y in ys)
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if ssx <= 0 or ssy <= 0:
        return math.nan
    return (cov / math.sqrt(ssx * ssy)) ** 2


def compare(hybrid_path: Path, baseline_path: Path) -> dict[str, object]:
    hybrid = read_rows(hybrid_path)
    baseline = read_rows(baseline_path)
    out: dict[str, object] = {"matched_rows": sum(pair_id in baseline for pair_id in hybrid)}
    for col in ["rg", "se"]:
        xs = []
        ys = []
        deltas = []
        rel = []
        for pair_id, hrow in hybrid.items():
            brow = baseline.get(pair_id)
            if brow is None:
                continue
            x = to_float(hrow[col])
            y = to_float(brow[col])
            if math.isfinite(x) and math.isfinite(y):
                xs.append(x)
                ys.append(y)
                delta = x - y
                deltas.append(delta)
                if y != 0:
                    rel.append(abs(delta / y))
        out[f"{col}_finite_pairs"] = len(xs)
        out[f"{col}_r2"] = r2(xs, ys)
        out[f"{col}_max_abs_delta"] = max(map(abs, deltas)) if deltas else math.nan
        out[f"{col}_mean_abs_delta"] = statistics.mean(map(abs, deltas)) if deltas else math.nan
        out[f"{col}_max_abs_rel_delta"] = max(rel) if rel else math.nan
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=Path("data/catalog/eur_gwas_manifest.tsv"), type=Path)
    parser.add_argument("--sumstats-dir", default=Path("data/sumstats/eur"), type=Path)
    parser.add_argument("--ld-prefix", default=Path("data/ld/UKBB.EUR"), type=Path)
    parser.add_argument(
        "--ldsc-bin",
        default=Path("external/ldsc-rs-rg-batch-target/release/ldsc"),
        type=Path,
    )
    parser.add_argument("--out-dir", default=Path("results/all_rg_validation"), type=Path)
    parser.add_argument("--n-pairs", default=100, type=int)
    parser.add_argument("--seed", default=20260527, type=int)
    parser.add_argument("--jobs", default=8, type=int)
    parser.add_argument("--hybrid-threads", default=8, type=int)
    parser.add_argument("--baseline-threads", default=1, type=int)
    parser.add_argument("--m-snps", default="auto")
    parser.add_argument("--min-rg-r2", default=0.999, type=float)
    parser.add_argument("--min-se-r2", default=0.999, type=float)
    parser.add_argument("--max-abs-delta", default=0.001, type=float)
    args = parser.parse_args()

    if not args.ldsc_bin.exists():
        raise SystemExit(f"Patched ldsc-rs binary not found: {args.ldsc_bin}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    traits = read_traits(args.manifest, args.sumstats_dir)
    pairs = select_disjoint_pairs(traits, args.n_pairs, args.seed)
    pairs_path = args.out_dir / "random_disjoint_pairs.tsv"
    write_pairs(pairs_path, pairs)
    m_snps = resolve_m(args.ld_prefix, args.m_snps)

    hybrid_path, hybrid_elapsed = run_hybrid(args, pairs_path, m_snps)
    baseline_path, baseline_elapsed = run_baseline(args, pairs, m_snps)
    summary = compare(hybrid_path, baseline_path)
    summary.update(
        {
            "n_pairs_requested": args.n_pairs,
            "seed": args.seed,
            "pair_selection": "random disjoint trait pairs",
            "unique_traits": len({trait for pair in pairs for trait in (pair["p1"], pair["p2"])}),
            "hybrid_elapsed_seconds": hybrid_elapsed,
            "baseline_elapsed_seconds": baseline_elapsed,
            "hybrid_path": str(hybrid_path),
            "baseline_path": str(baseline_path),
            "pairs_path": str(pairs_path),
        }
    )
    summary_path = args.out_dir / "comparison_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    failures = []
    matched_rows = int(summary.get("matched_rows", 0))
    rg_r2 = float(summary.get("rg_r2", math.nan))
    se_r2 = float(summary.get("se_r2", math.nan))
    rg_max_abs_delta = float(summary.get("rg_max_abs_delta", math.inf))
    se_max_abs_delta = float(summary.get("se_max_abs_delta", math.inf))
    rg_finite_pairs = int(summary.get("rg_finite_pairs", 0))
    se_finite_pairs = int(summary.get("se_finite_pairs", 0))

    if matched_rows != args.n_pairs:
        failures.append(f"matched_rows {matched_rows} != {args.n_pairs}")
    if rg_finite_pairs == 0:
        failures.append("rg_finite_pairs is 0")
    if se_finite_pairs == 0:
        failures.append("se_finite_pairs is 0")
    if not math.isfinite(rg_r2) or rg_r2 < args.min_rg_r2:
        failures.append(f"rg_r2 {summary.get('rg_r2')} < {args.min_rg_r2}")
    if not math.isfinite(se_r2) or se_r2 < args.min_se_r2:
        failures.append(f"se_r2 {summary.get('se_r2')} < {args.min_se_r2}")
    if not math.isfinite(rg_max_abs_delta) or rg_max_abs_delta > args.max_abs_delta:
        failures.append(f"rg_max_abs_delta {summary.get('rg_max_abs_delta')} > {args.max_abs_delta}")
    if not math.isfinite(se_max_abs_delta) or se_max_abs_delta > args.max_abs_delta:
        failures.append(f"se_max_abs_delta {summary.get('se_max_abs_delta')} > {args.max_abs_delta}")
    if failures:
        raise SystemExit("validation failed: " + "; ".join(failures))


if __name__ == "__main__":
    main()
