#!/usr/bin/env python3
"""Run Neale-patched LDSC over the upper triangle of a phenotype list."""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def run_lead(task: tuple[int, list[str], argparse.Namespace]) -> dict[str, object]:
    i, phenotype_ids, args = task
    lead = phenotype_ids[i]
    out_prefix = args.out_dir / lead
    list_path = args.out_dir / f"{lead}.rg-list.txt"
    log_path = args.out_dir / f"{lead}.driver.log"
    pair_count = len(phenotype_ids) - i - 1
    if pair_count == 0:
        return {"lead": lead, "pairs": 0, "status": "skip_last", "elapsed_seconds": 0}
    if Path(str(out_prefix) + ".r2").exists() and not args.force:
        return {"lead": lead, "pairs": pair_count, "status": "skip_done", "elapsed_seconds": 0}

    files = [args.sumstats_dir / f"{pid}.sumstats.gz" for pid in phenotype_ids[i:]]
    missing = [str(p) for p in files if not p.exists()]
    if missing:
        raise FileNotFoundError(f"Missing sumstats for {lead}: {missing[:3]}")
    list_path.write_text(",".join(str(p) for p in files))

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
    env = os.environ.copy()
    env.update(
        {
            "OPENBLAS_NUM_THREADS": "1",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
            "VECLIB_MAXIMUM_THREADS": "1",
        }
    )
    started = time.time()
    with log_path.open("w") as log:
        proc = subprocess.run(
            cmd, cwd=str(args.ldsc_dir), stdout=log, stderr=subprocess.STDOUT, env=env
        )
    elapsed = time.time() - started
    if proc.returncode != 0:
        raise RuntimeError(f"LDSC failed for {lead}; see {log_path}")
    return {"lead": lead, "pairs": pair_count, "status": "done", "elapsed_seconds": elapsed}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phenotypes", required=True, type=Path)
    parser.add_argument("--sumstats-dir", required=True, type=Path)
    parser.add_argument("--ldsc-dir", required=True, type=Path)
    parser.add_argument("--ldsc-python", required=True)
    parser.add_argument("--ld-prefix", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--jobs", default=4, type=int)
    parser.add_argument("--n-blocks", default=200, type=int)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.phenotypes = args.phenotypes.resolve()
    args.sumstats_dir = args.sumstats_dir.resolve()
    args.ldsc_dir = args.ldsc_dir.resolve()
    args.ld_prefix = args.ld_prefix.resolve()
    args.out_dir = args.out_dir.resolve()
    ldsc_python_path = Path(args.ldsc_python)
    if not ldsc_python_path.is_absolute() and ldsc_python_path.exists():
        args.ldsc_python = str(ldsc_python_path.resolve())
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with args.phenotypes.open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    phenotype_ids = [row["phenotype_id"] for row in rows]

    tasks = [(i, phenotype_ids, args) for i in range(len(phenotype_ids))]
    summary_path = args.out_dir / "driver_summary.tsv"
    with summary_path.open("w") as summary:
        summary.write("lead\tpairs\tstatus\telapsed_seconds\n")
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_lead, task) for task in tasks]
            for future in as_completed(futures):
                res = future.result()
                summary.write(
                    f"{res['lead']}\t{res['pairs']}\t{res['status']}\t{res['elapsed_seconds']:.6f}\n"
                )
                summary.flush()
                print(res, flush=True)


if __name__ == "__main__":
    main()
