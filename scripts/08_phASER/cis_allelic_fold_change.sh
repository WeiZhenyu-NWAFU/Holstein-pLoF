#!/usr/bin/env bash
set -euo pipefail
# Prepared gene-level GW-phased count matrix and allele-specific variant–gene pairs.
: "${PHASER_DIR:?Set phASER installation directory}"
: "${GW_BED:?Set indexed GW-phased count BED.gz}"
: "${VCF:?Set phased VCF}"
: "${PAIRS:?Set variant-gene pair TSV}"
: "${SAMPLE_MAP:?Set VCF-to-count-matrix sample map}"
: "${OUTPUT:?Set a new output filename}"
[[ ! -e "$OUTPUT" ]] || { echo 'Output already exists' >&2; exit 1; }
python2 "$PHASER_DIR/phaser_pop/phaser_cis_var.py" \
  --bed "$GW_BED" --vcf "$VCF" --pair "$PAIRS" --map "$SAMPLE_MAP" \
  --o "$OUTPUT" --t 20
