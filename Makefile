SHELL := /bin/bash

PYTHON ?= python3
LDSC_ENV_PREFIX ?= .envs/ldsc-neale
LDSC_PYTHON ?= $(LDSC_ENV_PREFIX)/bin/python
ENV_MANAGER ?= micromamba
JOBS ?= 16
BENCH_N ?= 90

PAN_BASE := https://pan-ukb-us-east-1.s3.amazonaws.com
MANIFEST := data/manifests/phenotype_manifest.tsv.bgz
H2_MANIFEST := data/manifests/h2_manifest.tsv.bgz
BENCH_PHENOS := data/benchmark90/phenotypes.tsv
CATALOG := data/catalog/eur_gwas_manifest.tsv
SUMSTATS_DIR := data/sumstats/eur
LD_PREFIX := data/ld/UKBB.EUR
LDSC_DIR := external/ldsc-neale

.PHONY: all init setup fetch-manifests catalog validate-catalog select-benchmark prepare-ldscores setup-ldsc setup-ldsc-env prepare-sumstats prepare-all-sumstats one-vs-all one-vs-all-dry-run run-benchmark summarize benchmark90 hardware clean-small

all: setup

init:
	@git init

fetch-manifests:
	@mkdir -p data/manifests
	curl -fL --retry 5 --retry-delay 5 -o $(MANIFEST) $(PAN_BASE)/sumstats_release/phenotype_manifest.tsv.bgz
	curl -fL --retry 5 --retry-delay 5 -o $(H2_MANIFEST) $(PAN_BASE)/sumstats_release/h2_manifest.tsv.bgz

catalog: fetch-manifests
	@mkdir -p data/catalog
	$(PYTHON) scripts/build_eur_gwas_catalog.py \
		--phenotype-manifest $(MANIFEST) \
		--out $(CATALOG)

validate-catalog: catalog
	$(PYTHON) scripts/validate_eur_gwas_catalog.py \
		--phenotype-manifest $(MANIFEST) \
		--catalog $(CATALOG)

select-benchmark: fetch-manifests
	@mkdir -p data/benchmark90
	$(PYTHON) scripts/select_benchmark_phenotypes.py \
		--manifest $(MANIFEST) \
		--config config/benchmark90.yaml \
		--out $(BENCH_PHENOS)

prepare-ldscores:
	@mkdir -p data/ld
	curl -fL --retry 5 --retry-delay 5 -o $(LD_PREFIX).l2.ldscore.gz $(PAN_BASE)/ld_release/UKBB.EUR.l2.ldscore.gz
	curl -fL --retry 5 --retry-delay 5 -o $(LD_PREFIX).l2.M $(PAN_BASE)/ld_release/UKBB.EUR.l2.M
	curl -fL --retry 5 --retry-delay 5 -o $(LD_PREFIX).l2.M_5_50 $(PAN_BASE)/ld_release/UKBB.EUR.l2.M_5_50
	gzip -dc $(LD_PREFIX).l2.ldscore.gz | awk 'NR>1 {print $$2}' > data/ld/UKBB.EUR.snps

setup-ldsc:
	bash scripts/setup_neale_ldsc.sh

setup-ldsc-env:
	@if [ -x "$(LDSC_ENV_PREFIX)/bin/python" ] && bash scripts/check_ldsc_env.sh "$(LDSC_ENV_PREFIX)/bin/python" >/dev/null; then \
		echo "LDSC environment already exists at $(LDSC_ENV_PREFIX)"; \
	else \
		if [ -e "$(LDSC_ENV_PREFIX)" ]; then \
			echo "Removing incomplete or invalid LDSC environment at $(LDSC_ENV_PREFIX)"; \
			rm -rf "$(LDSC_ENV_PREFIX)"; \
		fi; \
		if [ "$(ENV_MANAGER)" = "conda" ] || [ "$(ENV_MANAGER)" = "mamba" ]; then \
			$(ENV_MANAGER) env create -y -p $(LDSC_ENV_PREFIX) -f envs/ldsc-neale.yml; \
		else \
			$(ENV_MANAGER) create -y -p $(LDSC_ENV_PREFIX) -f envs/ldsc-neale.yml; \
		fi; \
		bash scripts/check_ldsc_env.sh "$(LDSC_ENV_PREFIX)/bin/python"; \
	fi

setup: setup-ldsc-env validate-catalog prepare-ldscores setup-ldsc prepare-all-sumstats

prepare-sumstats: $(BENCH_PHENOS) prepare-ldscores
	mkdir -p data/benchmark90/sumstats results/benchmark90/prepare_stats logs/prepare_sumstats
	$(PYTHON) scripts/prepare_sumstats_batch.py \
		--phenotypes $(BENCH_PHENOS) \
		--ld-snps data/ld/UKBB.EUR.snps \
		--out-dir data/benchmark90/sumstats \
		--stats-dir results/benchmark90/prepare_stats \
		--jobs $(JOBS)

prepare-all-sumstats: catalog prepare-ldscores
	mkdir -p $(SUMSTATS_DIR) results/prepare_all/prepare_stats logs/prepare_all_sumstats
	$(PYTHON) scripts/prepare_sumstats_batch.py \
		--phenotypes $(CATALOG) \
		--ld-snps data/ld/UKBB.EUR.snps \
		--out-dir $(SUMSTATS_DIR) \
		--stats-dir results/prepare_all/prepare_stats \
		--log-dir logs/prepare_all_sumstats \
		--jobs $(JOBS)

one-vs-all: setup-ldsc setup-ldsc-env catalog prepare-ldscores
	@if [ -z "$(PHENOCODE)$(PHENOTYPE_ID)$(QUERY)" ]; then \
		echo "Set one selector, e.g. make one-vs-all PHENOCODE=20016 JOBS=16"; \
		exit 2; \
	fi
	$(PYTHON) scripts/run_one_vs_all.py \
		--manifest $(CATALOG) \
		--sumstats-dir $(SUMSTATS_DIR) \
		--ldsc-dir $(LDSC_DIR) \
		--ldsc-python "$(LDSC_PYTHON)" \
		--ld-prefix $(LD_PREFIX) \
		--jobs $(JOBS) \
		$(if $(PHENOCODE),--phenocode $(PHENOCODE),) \
		$(if $(PHENOTYPE_ID),--phenotype-id $(PHENOTYPE_ID),) \
		$(if $(QUERY),--query "$(QUERY)",)

one-vs-all-dry-run: setup-ldsc setup-ldsc-env catalog prepare-ldscores
	@if [ -z "$(PHENOCODE)$(PHENOTYPE_ID)$(QUERY)" ]; then \
		echo "Set one selector, e.g. make one-vs-all-dry-run PHENOCODE=20016"; \
		exit 2; \
	fi
	$(PYTHON) scripts/run_one_vs_all.py \
		--manifest $(CATALOG) \
		--sumstats-dir $(SUMSTATS_DIR) \
		--ldsc-dir $(LDSC_DIR) \
		--ldsc-python "$(LDSC_PYTHON)" \
		--ld-prefix $(LD_PREFIX) \
		--jobs $(JOBS) \
		--dry-run \
		$(if $(PHENOCODE),--phenocode $(PHENOCODE),) \
		$(if $(PHENOTYPE_ID),--phenotype-id $(PHENOTYPE_ID),) \
		$(if $(QUERY),--query "$(QUERY)",)

run-benchmark: setup-ldsc setup-ldsc-env prepare-sumstats
	mkdir -p results/benchmark90/rg logs/ldsc
	$(PYTHON) scripts/run_ldsc_triangle.py \
		--phenotypes $(BENCH_PHENOS) \
		--sumstats-dir data/benchmark90/sumstats \
		--ldsc-dir $(LDSC_DIR) \
		--ldsc-python "$(LDSC_PYTHON)" \
		--ld-prefix $(LD_PREFIX) \
		--out-dir results/benchmark90/rg \
		--jobs $(JOBS)

summarize:
	mkdir -p benchmarks/benchmark90
	$(PYTHON) scripts/summarize_benchmark.py \
		--phenotypes $(BENCH_PHENOS) \
		--prepare-stats results/benchmark90/prepare_stats \
		--rg-dir results/benchmark90/rg \
		--out-md benchmarks/benchmark90/summary.md \
		--out-tsv benchmarks/benchmark90/summary.tsv

hardware:
	mkdir -p results/benchmark90
	bash scripts/hardware_report.sh > results/benchmark90/hardware.txt

benchmark90: fetch-manifests select-benchmark run-benchmark summarize

clean-small:
	rm -rf benchmarks/benchmark90/*.tmp logs/*.tmp
