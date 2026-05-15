#!/usr/bin/env python3
"""Run Neale-patched LDSC for one Pan-UKBB EUR GWAS against all others."""

from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


THREAD_ENV = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


def default_ldsc_python() -> str:
    env = os.environ.get("LDSC_PYTHON")
    if env:
        return env
    local_env = Path(".envs/ldsc-neale/bin/python")
    if local_env.exists():
        return str(local_env)
    return "python"


def norm(value: str | None) -> str:
    return (value or "").strip().lower()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def resolve_lead(rows: list[dict[str, str]], args: argparse.Namespace) -> dict[str, str]:
    matches = rows
    if args.phenotype_id:
        matches = [r for r in matches if r["phenotype_id"] == args.phenotype_id]
    if args.phenocode:
        matches = [r for r in matches if r["phenocode"] == args.phenocode]
    if args.query:
        q = norm(args.query)
        matches = [
            r
            for r in matches
            if q
            in norm(
                " ".join(
                    [
                        r.get("phenotype_id", ""),
                        r.get("phenocode", ""),
                        r.get("description", ""),
                        r.get("description_more", ""),
                        r.get("category", ""),
                    ]
                )
            )
        ]
    for field, value in [
        ("pheno_sex", args.sex),
        ("coding", args.coding),
        ("modifier", args.modifier),
        ("trait_type", args.trait_type),
    ]:
        if value is not None:
            matches = [r for r in matches if r.get(field, "") == value]

    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise SystemExit("No matching EUR GWAS found.")

    print("Ambiguous phenotype selector; matching rows:")
    fields = ["phenotype_id", "phenocode", "description", "pheno_sex", "coding", "modifier", "phenotype_qc_EUR"]
    print("\t".join(fields))
    for row in matches[:50]:
        print("\t".join(row.get(f, "") for f in fields))
    if len(matches) > 50:
        print(f"... {len(matches) - 50} more matches")
    raise SystemExit("Refine with --phenotype-id, --sex, --coding, or --modifier.")


def chunks(items: list[str], n_chunks: int) -> list[list[str]]:
    n_chunks = max(1, min(n_chunks, len(items)))
    size = int(math.ceil(len(items) / n_chunks))
    return [items[i : i + size] for i in range(0, len(items), size)]


def build_cmd(args: argparse.Namespace, list_path: Path, out_prefix: Path) -> list[str]:
    cmd = [
        args.ldsc_python,
        str(args.ldsc_dir / "ldsc.py"),
        "--rg",
        str(list_path),
        "--rg-file",
        "--ref-ld",
        str(args.ld_prefix),
        "--w-ld",
        str(args.ld_prefix),
        "--n-blocks",
        str(args.n_blocks),
        "--out",
        str(out_prefix),
        "--write-rg",
    ]
    if args.mem_gb_per_job and args.mem_gb_per_job > 0:
        prlimit = shutil.which("prlimit")
        if not prlimit:
            return cmd
        bytes_limit = str(int(args.mem_gb_per_job * 1024**3))
        cmd = [prlimit, f"--as={bytes_limit}", "--"] + cmd
    return cmd


def run_chunk(task: tuple[int, list[str], dict[str, str], argparse.Namespace]) -> dict[str, object]:
    idx, target_ids, lead, args = task
    chunk_name = f"chunk_{idx:04d}"
    out_prefix = args.out_dir / "chunks" / chunk_name
    list_path = args.out_dir / "chunks" / f"{chunk_name}.rg-list.txt"
    log_path = args.out_dir / "logs" / f"{chunk_name}.log"
    result_path = Path(str(out_prefix) + ".r2")
    if result_path.exists() and not args.force:
        return {
            "chunk": chunk_name,
            "targets": len(target_ids),
            "status": "skip_done",
            "elapsed_seconds": 0.0,
        }

    lead_file = args.sumstats_dir / f"{lead['phenotype_id']}.sumstats.gz"
    files = [lead_file] + [args.sumstats_dir / f"{pid}.sumstats.gz" for pid in target_ids]
    missing = [str(p) for p in files if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing prepared LDSC sumstats: {missing[:5]}")
    list_path.write_text(",".join(str(p) for p in files))

    env = os.environ.copy()
    env.update(THREAD_ENV)
    cmd = build_cmd(args, list_path, out_prefix)
    started = time.time()
    with log_path.open("w") as log:
        proc = subprocess.run(cmd, cwd=str(args.ldsc_dir), stdout=log, stderr=subprocess.STDOUT, env=env)
    elapsed = time.time() - started
    if proc.returncode != 0:
        raise RuntimeError(f"LDSC failed for {chunk_name}; see {log_path}")
    return {
        "chunk": chunk_name,
        "targets": len(target_ids),
        "status": "done",
        "elapsed_seconds": elapsed,
    }


def combine_results(out_dir: Path) -> int:
    def clean_trait(value: str) -> str:
        name = Path(value).name
        if name.endswith(".sumstats.gz"):
            return name[: -len(".sumstats.gz")]
        return value

    combined = out_dir / "rg.tsv"
    n = 0
    with combined.open("w") as out:
        wrote_header = False
        p1_idx = None
        p2_idx = None
        for path in sorted((out_dir / "chunks").glob("chunk_*.r2")):
            with path.open() as inp:
                header = inp.readline().split()
                if not wrote_header:
                    p1_idx = header.index("p1") if "p1" in header else None
                    p2_idx = header.index("p2") if "p2" in header else None
                    out.write("\t".join(header) + "\n")
                    wrote_header = True
                for line in inp:
                    parts = line.split()
                    if parts:
                        if p1_idx is not None and p1_idx < len(parts):
                            parts[p1_idx] = clean_trait(parts[p1_idx])
                        if p2_idx is not None and p2_idx < len(parts):
                            parts[p2_idx] = clean_trait(parts[p2_idx])
                        out.write("\t".join(parts) + "\n")
                        n += 1
    return n


def write_lead_metadata(out_dir: Path, lead: dict[str, str]) -> None:
    fields = [
        "phenotype_id",
        "phenocode",
        "description",
        "description_more",
        "category",
        "trait_type",
        "pheno_sex",
        "coding",
        "modifier",
        "phenotype_qc_EUR",
        "N_EUR_for_ldsc",
        "h2_observed_EUR",
        "h2_liability_EUR",
        "h2_z_EUR",
        "lambda_gc_EUR",
        "filename",
        "aws_path",
        "public_url",
    ]
    with (out_dir / "lead.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerow(lead)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", default=Path("data/catalog/eur_gwas_manifest.tsv"), type=Path)
    parser.add_argument("--sumstats-dir", default=Path("data/sumstats/eur"), type=Path)
    parser.add_argument("--ldsc-dir", default=Path("external/ldsc-neale"), type=Path)
    parser.add_argument("--ldsc-python", default=default_ldsc_python())
    parser.add_argument("--ld-prefix", default=Path("data/ld/UKBB.EUR"), type=Path)
    parser.add_argument("--out-base", default=Path("results/one_vs_all"), type=Path)
    parser.add_argument("--phenotype-id")
    parser.add_argument("--phenocode")
    parser.add_argument("--query")
    parser.add_argument("--sex")
    parser.add_argument("--coding")
    parser.add_argument("--modifier")
    parser.add_argument("--trait-type")
    parser.add_argument("--jobs", default=16, type=int)
    parser.add_argument("--chunks", type=int, help="Number of target chunks; default is --jobs.")
    parser.add_argument("--n-blocks", default=200, type=int)
    parser.add_argument("--mem-gb-per-job", default=8.0, type=float)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.manifest = args.manifest.resolve()
    args.sumstats_dir = args.sumstats_dir.resolve()
    args.ldsc_dir = args.ldsc_dir.resolve()
    args.ld_prefix = args.ld_prefix.resolve()
    args.out_base = args.out_base.resolve()
    ldsc_python_path = Path(args.ldsc_python)
    if not ldsc_python_path.is_absolute() and ldsc_python_path.exists():
        args.ldsc_python = str(ldsc_python_path.resolve())

    rows = read_manifest(args.manifest)
    lead = resolve_lead(rows, args)
    target_ids = [r["phenotype_id"] for r in rows if r["phenotype_id"] != lead["phenotype_id"]]
    target_chunks = chunks(target_ids, args.chunks or args.jobs)

    safe_id = lead["phenotype_id"]
    args.out_dir = args.out_base / safe_id
    (args.out_dir / "chunks").mkdir(parents=True, exist_ok=True)
    (args.out_dir / "logs").mkdir(parents=True, exist_ok=True)
    write_lead_metadata(args.out_dir, lead)

    expected = [args.sumstats_dir / f"{lead['phenotype_id']}.sumstats.gz"] + [
        args.sumstats_dir / f"{pid}.sumstats.gz" for pid in target_ids
    ]
    missing = [p for p in expected if not p.exists()]
    if args.dry_run:
        print(f"lead\t{lead['phenotype_id']}")
        print(f"phenocode\t{lead['phenocode']}")
        print(f"description\t{lead['description']}")
        print(f"targets\t{len(target_ids)}")
        print(f"chunks\t{len(target_chunks)}")
        print(f"jobs\t{args.jobs}")
        print(f"mem_gb_per_job\t{args.mem_gb_per_job}")
        print(f"prepared_sumstats_present\t{len(expected) - len(missing)}")
        print(f"prepared_sumstats_missing\t{len(missing)}")
        if missing:
            print("first_missing")
            for path in missing[:10]:
                print(path)
        return

    summary_path = args.out_dir / "driver_summary.tsv"
    with summary_path.open("w") as summary:
        summary.write("chunk\ttargets\tstatus\telapsed_seconds\n")
        tasks = [(i, target_chunk, lead, args) for i, target_chunk in enumerate(target_chunks, start=1)]
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_chunk, task) for task in tasks]
            for future in as_completed(futures):
                res = future.result()
                summary.write(
                    f"{res['chunk']}\t{res['targets']}\t{res['status']}\t{res['elapsed_seconds']:.6f}\n"
                )
                summary.flush()
                print(res, flush=True)

    rows_written = combine_results(args.out_dir)
    print(f"lead\t{lead['phenotype_id']}")
    print(f"targets\t{len(target_ids)}")
    print(f"chunks\t{len(target_chunks)}")
    print(f"combined_results\t{args.out_dir / 'rg.tsv'}")
    print(f"result_rows\t{rows_written}")


if __name__ == "__main__":
    main()
