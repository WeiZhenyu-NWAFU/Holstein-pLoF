#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
gatk SelectVariants \
	--variant  upper_bound.pLoF.vcf.gz \
	--output   upper_bound.pLoF.filter.vcf.gz \
	--max-indel-size 50 \
	--exclude-non-variants true \
        --remove-unused-alternates true

