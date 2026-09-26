#!/usr/bin/env bash
set -euo pipefail
: "${VCF:?Set phased/imputed genotype VCF}"
: "${OUTDIR:?Set a new output directory}"
[[ ! -e "$OUTDIR" ]] || { echo 'Output directory already exists' >&2; exit 1; }
mkdir -p "$OUTDIR"
plink --vcf "$VCF" --maf 0.01 --chr-set 29 --keep-allele-order \
  --make-bed --out "$OUTDIR/milk"
