#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${RELATE_DIR:?Set RELATE_DIR; see SOFTWARE.md}"
: "${GENEALOGY_DIR:?Set GENEALOGY_DIR; see SOFTWARE.md}"

chrom=${1}
pos=${2}

"${RELATE_DIR}"/scripts/SampleBranchLengths/SampleBranchLengths.sh \
        -i "${GENEALOGY_DIR}"/Chr${chrom}_Holstein_reinfer \
        -o clues_infer/${chrom}_${pos}/${chrom}_${pos} \
        -m 1.26e-8 \
        --coal "${GENEALOGY_DIR}"/Chr${chrom}_Holstein_reinfer.coal \
        --format b \
        --first_bp ${pos} \
        --last_bp ${pos} \
        --seed 1 \
        --num_samples 200
