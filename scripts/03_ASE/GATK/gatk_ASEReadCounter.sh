#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA; see SOFTWARE.md}"
sample=$1
gatk ASEReadCounter \
	--reference "${REFERENCE_FASTA}" \
	--input $sample.bam \
	--variant sites.vcf.gz \
	--output $sample.ASE.table \
	--showHidden
