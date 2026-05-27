#!/usr/bin/env python3
"""Run all pairwise EUR LDSC genetic correlations using patched ldsc-rs rg-batch."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import TextIO


HEADER = ["pair_id", "i", "j", "p1", "p2", "lead", "target"]


def open_text(path: Path, mode: str) -> TextIO:
    if path.name.endswith(".gz") or path.name.endswith(".gz.tmp"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def count_data_rows(path: Path) -> int:
    with open_text(path, "r") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def pair_id(n_traits: int, i: int, j: int) -> int:
    # One-based, row-major over the upper triangle, matching the earlier validation files.
    return i * (2 * n_traits - i - 1) // 2 + (j - i)


def read_traits(manifest: Path, sumstats_dir: Path, allow_missing: bool) -> list[str]:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "phenotype_id" not in (reader.fieldnames or []):
            raise SystemExit(f"{manifest} does not contain a phenotype_id column")
        traits = [row["phenotype_id"] for row in reader if row.get("phenotype_id")]

    seen: set[str] = set()
    deduped: list[str] = []
    for trait in traits:
        if trait in seen:
            continue
        seen.add(trait)
        deduped.append(trait)

    missing = [t for t in deduped if not (sumstats_dir / f"{t}.sumstats.gz").exists()]
    if missing and not allow_missing:
        preview = ", ".join(missing[:10])
        raise SystemExit(
            f"{len(missing)} traits are missing prepared sumstats. "
            f"Run `make setup` first. First missing: {preview}"
        )
    if missing:
        missing_set = set(missing)
        deduped = [t for t in deduped if t not in missing_set]
    return deduped


def write_traits(path: Path, traits: list[str], sumstats_dir: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait_index", "trait_id", "sumstats_path"])
        for idx, trait in enumerate(traits):
            writer.writerow([idx, trait, sumstats_dir / f"{trait}.sumstats.gz"])
    tmp.replace(path)


def shard_rows(
    traits: list[str], block_i: int, block_j: int, block_size: int
) -> list[tuple[int, int, int, str, str, str, str]]:
    n = len(traits)
    i0 = block_i * block_size
    i1 = min(i0 + block_size, n)
    j0 = block_j * block_size
    j1 = min(j0 + block_size, n)
    rows = []
    for i in range(i0, i1):
        start_j = max(i + 1, j0)
        for j in range(start_j, j1):
            rows.append((pair_id(n, i, j), i, j, traits[i], traits[j], traits[i], traits[j]))
    return rows


def write_shard(path: Path, rows: list[tuple[int, int, int, str, str, str, str]]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)
    tmp.replace(path)


def prepare_shards(args: argparse.Namespace) -> list[dict[str, str]]:
    traits = read_traits(args.manifest, args.sumstats_dir, args.allow_missing_sumstats)
    if len(traits) < 2:
        raise SystemExit("Need at least two prepared traits")

    meta_dir = args.out_dir / "metadata"
    pair_dir = args.out_dir / "pair_shards"
    meta_dir.mkdir(parents=True, exist_ok=True)
    pair_dir.mkdir(parents=True, exist_ok=True)
    write_traits(meta_dir / "traits.tsv", traits, args.sumstats_dir)

    n = len(traits)
    n_blocks = (n + args.trait_block_size - 1) // args.trait_block_size
    manifest_rows: list[dict[str, str]] = []
    total_pairs = 0
    shard_id = 0

    for block_i in range(n_blocks):
        for block_j in range(block_i, n_blocks):
            rows = shard_rows(traits, block_i, block_j, args.trait_block_size)
            if not rows:
                continue
            shard = pair_dir / f"pairs.block_{block_i:03d}_{block_j:03d}.tsv"
            if args.force_shards or not shard.exists():
                write_shard(shard, rows)
            total_pairs += len(rows)
            manifest_rows.append(
                {
                    "shard_id": str(shard_id),
                    "block_i": str(block_i),
                    "block_j": str(block_j),
                    "n_pairs": str(len(rows)),
                    "path": str(shard.resolve()),
                }
            )
            shard_id += 1

    expected = n * (n - 1) // 2
    if total_pairs != expected:
        raise SystemExit(f"internal error: planned {total_pairs} pairs, expected {expected}")

    manifest = meta_dir / "shards.tsv"
    tmp = manifest.with_suffix(manifest.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["shard_id", "block_i", "block_j", "n_pairs", "path"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)
    tmp.replace(manifest)
    print(
        f"prepared {len(manifest_rows)} shards for {n} traits and {total_pairs} pairs",
        flush=True,
    )
    return manifest_rows


def read_shard_manifest(out_dir: Path) -> list[dict[str, str]]:
    manifest = out_dir / "metadata" / "shards.tsv"
    if not manifest.exists():
        raise SystemExit(f"Missing shard manifest: {manifest}. Run with --prepare-only first.")
    with manifest.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def result_path(args: argparse.Namespace, shard_path: Path) -> Path:
    base = shard_path.stem
    suffix = ".rg.tsv.gz" if args.compress_output else ".rg.tsv"
    return args.out_dir / "rg_shards" / f"{base}{suffix}"


def done_path(args: argparse.Namespace, shard_path: Path) -> Path:
    return args.out_dir / "rg_shards" / f"{shard_path.stem}.done"


def is_complete(args: argparse.Namespace, row: dict[str, str]) -> bool:
    shard = Path(row["path"])
    out = result_path(args, shard)
    done = done_path(args, shard)
    if not out.exists() or not done.exists():
        return False
    return count_data_rows(out) == int(row["n_pairs"])


def resolve_m(ld_prefix: Path, m_override: str | None) -> str:
    if m_override and m_override != "auto":
        return m_override
    m_file = Path(str(ld_prefix) + ".l2.M_5_50")
    if not m_file.exists():
        raise SystemExit(f"M override is auto but M file does not exist: {m_file}")
    with m_file.open() as handle:
        return handle.read().split()[0]


def gzip_move(src: Path, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    with src.open("rb") as inp, gzip.open(tmp, "wb", compresslevel=6) as out:
        shutil.copyfileobj(inp, out)
    tmp.replace(dest)
    src.unlink()


def run_one_shard(args: argparse.Namespace, row: dict[str, str], m_snps: str) -> tuple[str, str]:
    shard = Path(row["path"])
    expected = int(row["n_pairs"])
    out = result_path(args, shard)
    done = done_path(args, shard)
    lock = args.out_dir / "locks" / f"{shard.stem}.lock"
    log = args.out_dir / "logs" / f"{shard.stem}.log"
    out.parent.mkdir(parents=True, exist_ok=True)
    log.parent.mkdir(parents=True, exist_ok=True)
    lock.parent.mkdir(parents=True, exist_ok=True)

    if not args.force and is_complete(args, row):
        return shard.stem, "skip"
    try:
        lock.mkdir()
    except FileExistsError:
        return shard.stem, "locked"

    tmp = out.parent / f"{shard.stem}.rg.tsv.tmp.{os.getpid()}"
    if tmp.exists():
        tmp.unlink()

    cmd = [
        str(args.ldsc_bin),
        "rg-batch",
        "--rayon-threads",
        str(args.rayon_threads),
        "--pairs",
        str(shard),
        "--sumstats-dir",
        str(args.sumstats_dir),
        "--ref-ld",
        str(Path(str(args.ld_prefix) + ".l2.ldscore.gz")),
        "--w-ld",
        str(Path(str(args.ld_prefix) + ".l2.ldscore.gz")),
        "--M",
        m_snps,
        "--out",
        str(tmp),
        "--verbose-timing",
    ]
    if args.no_check_alleles:
        cmd.append("--no-check-alleles")

    started = time.time()
    try:
        with log.open("w") as handle:
            if Path("/usr/bin/time").exists():
                full_cmd = ["/usr/bin/time", "-v", *cmd]
            else:
                full_cmd = cmd
            proc = subprocess.run(full_cmd, stdout=handle, stderr=subprocess.STDOUT)
        if proc.returncode != 0:
            raise RuntimeError(f"{shard.stem} failed; see {log}")
        observed = count_data_rows(tmp)
        if observed != expected:
            raise RuntimeError(f"{shard.stem} produced {observed} rows, expected {expected}")
        if args.compress_output:
            gzip_move(tmp, out)
        else:
            tmp.replace(out)
        done.write_text(
            json.dumps(
                {
                    "shard": shard.stem,
                    "pairs": expected,
                    "threads": args.rayon_threads,
                    "elapsed_seconds": time.time() - started,
                    "output": str(out),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        return shard.stem, "done"
    finally:
        if tmp.exists():
            tmp.unlink()
        try:
            lock.rmdir()
        except OSError:
            pass


def print_progress(args: argparse.Namespace, rows: list[dict[str, str]]) -> None:
    total_shards = len(rows)
    done_shards = 0
    total_pairs = 0
    done_pairs = 0
    bad = []
    for row in rows:
        expected = int(row["n_pairs"])
        total_pairs += expected
        shard = Path(row["path"])
        out = result_path(args, shard)
        done = done_path(args, shard)
        if out.exists() and done.exists():
            observed = count_data_rows(out)
            if observed == expected:
                done_shards += 1
                done_pairs += observed
            else:
                bad.append(f"{shard.stem}: {observed} rows, expected {expected}")
    pct = 100.0 * done_pairs / total_pairs if total_pairs else 0.0
    print(f"shards: {done_shards}/{total_shards}")
    print(f"pairs:  {done_pairs}/{total_pairs} ({pct:.2f}%)")
    if bad:
        print("bad shards:")
        for item in bad[:20]:
            print(f"  {item}")


def collect_results(args: argparse.Namespace, rows: list[dict[str, str]]) -> None:
    out = args.collect_out or (args.out_dir / "rg.tsv.gz")
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    wrote_header = False
    n_rows = 0
    with open_text(tmp, "w") as out_handle:
        for row in rows:
            shard = Path(row["path"])
            expected = int(row["n_pairs"])
            path = result_path(args, shard)
            if not path.exists():
                if args.allow_incomplete_collect:
                    continue
                raise SystemExit(f"missing result: {path}")
            observed = count_data_rows(path)
            if observed != expected and not args.allow_incomplete_collect:
                raise SystemExit(f"{path} has {observed} rows, expected {expected}")
            with open_text(path, "r") as inp:
                header = inp.readline()
                if not wrote_header:
                    out_handle.write(header)
                    wrote_header = True
                for line in inp:
                    out_handle.write(line)
                    n_rows += 1
    tmp.replace(out)
    print(f"wrote {n_rows} rows to {out}")


def run_shards(args: argparse.Namespace, rows: list[dict[str, str]]) -> None:
    if args.max_shards:
        rows = rows[: args.max_shards]
    pending = [row for row in rows if args.force or not is_complete(args, row)]
    print(f"pending shards: {len(pending)}/{len(rows)}")
    if args.dry_run:
        for row in pending[:20]:
            print(f"would run {Path(row['path']).name} ({row['n_pairs']} pairs)")
        if len(pending) > 20:
            print(f"... {len(pending) - 20} more")
        return
    m_snps = resolve_m(args.ld_prefix, args.m_snps)
    failures = []
    with ThreadPoolExecutor(max_workers=args.max_parallel_shards) as pool:
        active = set()
        iterator = iter(pending)
        while True:
            while len(active) < args.max_parallel_shards:
                try:
                    row = next(iterator)
                except StopIteration:
                    break
                active.add(pool.submit(run_one_shard, args, row, m_snps))
            if not active:
                break
            done, active = wait(active, return_when=FIRST_COMPLETED)
            for future in done:
                try:
                    shard, status = future.result()
                    print(f"{status}\t{shard}", flush=True)
                except Exception as exc:  # noqa: BLE001
                    print(f"failed\t{exc}", file=sys.stderr, flush=True)
                    failures.append(str(exc))
    if failures:
        raise SystemExit(f"{len(failures)} shard(s) failed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=Path("data/catalog/eur_gwas_manifest.tsv"), type=Path)
    parser.add_argument("--sumstats-dir", default=Path("data/sumstats/eur"), type=Path)
    parser.add_argument("--ld-prefix", default=Path("data/ld/UKBB.EUR"), type=Path)
    parser.add_argument(
        "--ldsc-bin",
        default=Path("external/ldsc-rs-rg-batch-target/release/ldsc"),
        type=Path,
    )
    parser.add_argument("--out-dir", default=Path("results/all_rg"), type=Path)
    parser.add_argument("--trait-block-size", default=256, type=int)
    parser.add_argument("--rayon-threads", default=50, type=int)
    parser.add_argument("--max-parallel-shards", default=1, type=int)
    parser.add_argument("--m-snps", default="auto")
    parser.add_argument("--allow-missing-sumstats", action="store_true")
    parser.add_argument("--force-shards", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--collect", action="store_true")
    parser.add_argument("--collect-out", type=Path)
    parser.add_argument("--allow-incomplete-collect", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--max-shards", type=int)
    parser.add_argument("--no-compress-output", dest="compress_output", action="store_false")
    parser.add_argument("--no-check-alleles", action="store_true")
    parser.set_defaults(compress_output=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.trait_block_size < 2:
        raise SystemExit("--trait-block-size must be >= 2")
    if args.max_parallel_shards < 1:
        raise SystemExit("--max-parallel-shards must be >= 1")
    if args.rayon_threads < 1:
        raise SystemExit("--rayon-threads must be >= 1")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.progress or args.collect:
        rows = read_shard_manifest(args.out_dir)
    else:
        rows = prepare_shards(args)

    if args.prepare_only:
        return
    if args.progress:
        print_progress(args, rows)
        return
    if args.collect:
        collect_results(args, rows)
        return
    if not args.ldsc_bin.exists():
        raise SystemExit(f"Patched ldsc-rs binary not found: {args.ldsc_bin}")
    run_shards(args, rows)
    print_progress(args, read_shard_manifest(args.out_dir))


if __name__ == "__main__":
    main()
