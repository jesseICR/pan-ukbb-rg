# Design

## Why use the Neale batching pattern?

Stock LDSC can accept multiple traits in `--rg`, but Neale's fork adds `--rg-file` and `--write-rg`, which makes it practical to submit one job per lead phenotype. For phenotype `i`, the input list contains phenotype `i` plus phenotypes `i+1..N`. Across all lead phenotypes, this covers each unordered pair once.

This reduces scheduler overhead from one job per pair to one job per phenotype, while still performing the same number of pairwise LDSC regressions.

## Why benchmark locally?

LDSC is CPU-bound Python and does not use GPUs. Local benchmarking gives a concrete seconds-per-pair estimate for the hardware that runs the workflow. Any cloud estimate should be treated as optional planning context; the main benchmark output is wall-clock runtime and CPU-hours.

## Neale script sources

`Nealelab/UKBB_ldsc_scripts` documents the Neale Lab UKBB SNP-heritability workflow and includes `round1/rg_single.py` for one target phenotype against chunks of significant UKBB phenotypes. That repository is useful background for sumstats assumptions and cloud batching, but it is not the direct all-pairs implementation used here.

For the all-pairs genetic-correlation benchmark, this repo follows the later `astheeggeggs/UKBB_ldsc_r2` pattern and the patched `astheeggeggs/ldsc` fork, because that fork exposes `--rg-file` and `--write-rg` for upper-triangle batching.

## Full run denominator

From the Pan-UKBB phenotype manifest, the current all-QC EUR set has 7,160 GWAS rows:

```text
7,160 * 7,159 / 2 = 25,629,220 unordered pairs
```

The benchmark uses 90 GWAS rows:

```text
90 * 89 / 2 = 4,005 unordered pairs
```
