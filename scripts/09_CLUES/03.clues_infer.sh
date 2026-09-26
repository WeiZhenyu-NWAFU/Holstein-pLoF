#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${GENEALOGY_DIR:?Set GENEALOGY_DIR; see SOFTWARE.md}"
: "${CLUES_DIR:?Set CLUES_DIR; see SOFTWARE.md}"

chrom=${1}
pos=${2}

python "${CLUES_DIR}"/inference.py \
        --times clues_infer/${chrom}_${pos}/${chrom}_${pos} \
        --coal "${GENEALOGY_DIR}"/Chr${chrom}_Holstein_reinfer.coal \
        --out clues_infer/${chrom}_${pos}/${chrom}_${pos} \
	--burnin 100 \
	--timeBins timeBins.txt
