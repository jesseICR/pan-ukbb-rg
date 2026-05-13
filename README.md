# Pan-UKBB EUR One-vs-All Genetic Correlations

This repository prepares public Pan-UKBB EUR GWAS summary statistics for LDSC and computes genetic correlations between one selected GWAS and every other EUR GWAS.

The workflow uses the Neale Lab batching idea: convert Pan-UKBB flat files into compact LDSC `.sumstats.gz` files, then run the patched Neale LDSC fork with `--rg-file` and `--write-rg` so one lead phenotype can be compared against target chunks in parallel.

## Quick Start

Clone the repo, then run the one-time setup:

```bash
git clone https://github.com/jesseICR/pan-ukbb-rg.git
cd pan-ukbb-rg
```

Run setup:

```bash
make setup
```

Then compute all genetic correlations involving one GWAS:

```bash
python3 scripts/run_one_vs_all.py --phenocode 20016
```

`20016` is the Pan-UKBB phenocode for fluid intelligence score. Replace it with any phenocode in `data/catalog/eur_gwas_manifest.tsv`.

You can also select by phenotype ID or text query:

```bash
python3 scripts/run_one_vs_all.py --phenotype-id continuous-20016-both_sexes-irnt
python3 scripts/run_one_vs_all.py --query "fluid intelligence"
```

If a selector matches multiple GWAS, the script prints the matching rows and asks you to use a more specific selector.

## Requirements

You need:

- a Unix-like shell with `bash` and `make`
- `curl`, `gzip`, `awk`, and `git`
- Python 3.10 or newer for this repo's scripts
- `micromamba`, `mamba`, or `conda` to create the isolated LDSC Python 2.7 environment
- public internet access to Pan-UKBB HTTPS URLs
- at least 150 GiB free disk for the prepared data cache
- enough RAM/CPU for parallel LDSC jobs

No AWS account, AWS credentials, or AWS CLI are required. Pan-UKBB data are downloaded from public HTTPS URLs such as:

```text
https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/phenotype_manifest.tsv.bgz
```

By default, `make setup` uses `micromamba`. If you use another environment manager:

```bash
make setup ENV_MANAGER=mamba
make setup ENV_MANAGER=conda
```

The LDSC environment is created under `.envs/ldsc-neale/`, which is local generated state and is not meant to be committed.

## Runtime And Cores

The default is `JOBS=16` for both setup and one-vs-all rg. That may still be too many jobs for a laptop or shared workstation. Use fewer jobs if the machine becomes unresponsive:

```bash
make setup JOBS=8
python3 scripts/run_one_vs_all.py --phenocode 20016 --jobs 8
```

Setup parallelism helps because each job streams, decompresses, filters, and writes a separate GWAS. Scaling is not perfect because network bandwidth, disk writes, gzip decompression, and remote throttling become shared bottlenecks.

From the 90-phenotype benchmark in this repo:

```text
compressed input: 138.9 GiB
wall time at 16 jobs: ~31 min
sum of per-file preparation time: 19,654 sec
observed parallel efficiency at 16 jobs: ~66%
```

Extrapolated setup time for the full ~6.7 TiB input:

```text
8 jobs:  ~34 h ideal, ~40-55 h realistic
16 jobs: ~17 h ideal, ~26 h based on observed benchmark
32 jobs: ~8 h ideal, probably ~13-24 h depending on network/disk saturation
```

GPUs are not useful for this workflow. LDSC is CPU-bound Python, and setup is mostly streaming, decompression, filtering, and writing.

## Is Setup Resumable?

Yes. `make setup` is designed to be rerun.

The large setup step is `make prepare-all-sumstats`. For each phenotype, the converter writes temporary files first and only marks that phenotype complete after both of these final files exist:

```text
data/sumstats/eur/<phenotype_id>.sumstats.gz
results/prepare_all/prepare_stats/<phenotype_id>.json
```

If the computer shuts down mid-run, rerun:

```bash
make setup
```

Completed phenotypes are skipped. Interrupted or incomplete phenotypes are streamed and converted again. The raw Pan-UKBB flat files are not stored locally.

## Important Files

After `make setup`, the most useful files are:

```text
data/catalog/eur_gwas_manifest.tsv
```

The searchable EUR GWAS catalog. Use this to find phenocodes, phenotype IDs, descriptions, QC labels, sample sizes, and source URLs.

```text
data/sumstats/eur/<phenotype_id>.sumstats.gz
```

The compact LDSC-ready EUR sumstats cache. These are the files used by `run_one_vs_all.py`.

```text
results/prepare_all/prepare_stats/<phenotype_id>.json
```

Per-GWAS conversion diagnostics: input rows, LD-score SNPs seen, written rows, skipped low-confidence rows, skipped missing EUR stats, and elapsed time.

```text
data/ld/UKBB.EUR.l2.ldscore.gz
data/ld/UKBB.EUR.l2.M
data/ld/UKBB.EUR.l2.M_5_50
data/ld/UKBB.EUR.snps
```

EUR LD-score reference files used by LDSC. `UKBB.EUR.snps` is the extracted list of SNP IDs used to filter the Pan-UKBB flat files.

```text
external/ldsc-neale/
.envs/ldsc-neale/
```

The patched LDSC code and its isolated Python 2.7 environment.

After `run_one_vs_all.py`, the main output is:

```text
results/one_vs_all/<lead_phenotype_id>/rg.tsv
```

This is the combined table of genetic correlations between the selected lead GWAS and all other EUR GWAS.
The `p1` and `p2` columns are phenotype IDs from `data/catalog/eur_gwas_manifest.tsv`.

Other useful one-vs-all outputs:

```text
results/one_vs_all/<lead_phenotype_id>/lead.tsv
results/one_vs_all/<lead_phenotype_id>/driver_summary.tsv
results/one_vs_all/<lead_phenotype_id>/logs/
results/one_vs_all/<lead_phenotype_id>/chunks/
```

`lead.tsv` records the selected GWAS metadata. `driver_summary.tsv` records per-chunk runtime/status. `logs/` is for troubleshooting. `chunks/` contains the raw Neale LDSC `.r2` chunk outputs that are combined into `rg.tsv`.

Large public inputs and generated outputs under `data/`, `external/`, `results/`, `logs/`, and `.envs/` are gitignored.

## Dry Runs

A dry run is optional. It does not stream GWAS data and does not run LDSC. It only resolves the phenotype and reports whether the prepared sumstats cache exists.

```bash
python3 scripts/run_one_vs_all.py --phenocode 20016 --dry-run
```

This is useful before launching a long one-vs-all run or when checking whether `make setup` completed.

## What Setup Does

`make setup` runs these steps:

```text
make setup-ldsc-env
make validate-catalog
make prepare-ldscores
make setup-ldsc
make prepare-all-sumstats
```

`make setup-ldsc-env` creates `.envs/ldsc-neale/` from `envs/ldsc-neale.yml`.

`make validate-catalog` downloads the public Pan-UKBB manifests, builds `data/catalog/eur_gwas_manifest.tsv`, and validates it against the source manifest.

`make prepare-ldscores` downloads the Pan-UKBB EUR LD-score files and extracts `data/ld/UKBB.EUR.snps`.

`make setup-ldsc` checks out the patched Neale LDSC fork from `https://github.com/astheeggeggs/ldsc.git` at commit `a4ee4c8aa065a1c9a586c3b678e9b3040bbebafc`.

`make prepare-all-sumstats` streams every Pan-UKBB EUR GWAS flat file and writes compact LDSC-ready sumstats.

Expected data volume:

```text
streamed from Pan-UKBB: ~6.7 TiB
stored after setup: ~70-100 GiB in data/sumstats/eur/
recommended free disk: >=150 GiB
```

## Sumstats Conversion

The Pan-UKBB per-phenotype flat files are not pre-munged LDSC files. This repo converts EUR columns to minimal LDSC-compatible files with:

- `SNP = chr:pos:ref:alt`
- `A1 = alt`
- `A2 = ref`
- `Z = beta_EUR / se_EUR`
- `N = n_cases_EUR + n_controls_EUR` for binary traits, or `n_cases_EUR` otherwise

Rows are skipped if they are outside the EUR LD-score SNP list, have `low_confidence_EUR == true`, have missing EUR beta/se, or have nonpositive/nonfinite standard error.

## Optional Benchmark

The benchmark is separate from the main one-vs-all workflow. It estimates local LDSC runtime using 90 selected phenotypes: 30 `PASS`, 30 `h2_z_insignificant`, and 30 `not_EUR_plus_1`.

```bash
make benchmark90
```

Benchmark outputs are written under:

```text
benchmarks/benchmark90/
results/benchmark90/
```
