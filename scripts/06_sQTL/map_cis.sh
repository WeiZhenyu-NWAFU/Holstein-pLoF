#!/usr/bin/env bash
set -euo pipefail
# Run in an assay-specific work directory containing milk/{genotypes,phenotypes,covFile}.
# OmiGA installation used in the source command: v1.1.3-beta.2+250831.
# Resolve `omiga` through PATH in the configured environment.
tis=${1:?Usage: bash map_cis.sh milk 30 20 NEW_OUTPUT_DIR}
peer_num=${2:?Set PEER factor count}
threads=${3:?Set thread count}
output_dir=${4:?Set a new output directory}
cov_file="${tis}/covFile/milk.peer${peer_num}.tsv"
[[ -f "$cov_file" ]] || { echo 'Missing covariate file' >&2; exit 1; }
[[ ! -e "$output_dir" ]] || { echo 'Output directory already exists' >&2; exit 1; }
mkdir -p "$output_dir"
omiga --mode cis --verbose --debug --threads $threads --force-double-precision \
      --covariates "$cov_file" --calcu-variant-threshold \
      --genotype "${tis}/genotypes/${tis}" --phenotype "${tis}/phenotypes/${tis}.leafcutter.sorted.bed.gz" \
      --prefix "${tis}" --output-dir "$output_dir" --pheno-group "${tis}/phenotypes/${tis}.leafcutter.phenotype_groups_omiga.txt"
cis_result_file="${output_dir}/${tis}.cis_qtl.txt.gz"
[[ -f "$cis_result_file" ]] || { echo 'cis output missing' >&2; exit 1; }
omiga --mode cis_independent --verbose --debug --threads $threads --force-double-precision \
      --covariates "$cov_file" --calcu-variant-threshold \
      --genotype "${tis}/genotypes/${tis}" --phenotype "${tis}/phenotypes/${tis}.leafcutter.sorted.bed.gz" \
      --cis-file "$cis_result_file" \
      --prefix "${tis}" --output-dir "$output_dir" --pheno-group "${tis}/phenotypes/${tis}.leafcutter.phenotype_groups_omiga.txt"
