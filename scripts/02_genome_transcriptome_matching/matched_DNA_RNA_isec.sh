#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${PROJECT_ROOT:?Set PROJECT_ROOT; see SOFTWARE.md}"
sample=$1

bcftools isec \
	-c none \
	-p lung-${sample}/matched_pLoF \
	"${PROJECT_ROOT}"/533_sample/pLoF_variant/533_indvduals_vcf/533_sample_private_variants_vcf/liuwei-chongcexu-${sample}/liuwei-chongcexu-${sample}.vcf.gz \
	lung-${sample}/lung-${sample}.vcf.gz
