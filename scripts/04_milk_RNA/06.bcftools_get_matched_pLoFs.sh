#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${PROJECT_ROOT:?Set PROJECT_ROOT; see SOFTWARE.md}"
sample=$1

bcftools isec \
	-c none \
	-p ${sample}/matched_pLoF \
	"${PROJECT_ROOT}"/574_sample/02.annotation/04.venn_collect/v10_sLOF_a4/upper_bound.pLoF.vcf.gz \
	${sample}/${sample}.vcf.gz
