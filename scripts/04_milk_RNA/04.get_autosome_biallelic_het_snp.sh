#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${PROJECT_ROOT:?Set PROJECT_ROOT; see SOFTWARE.md}"
## based GATK) v4.3.0.0

sample=$1


gatk --java-options "-Xmx300G -XX:ParallelGCThreads=1  -Djava.io.tmpdir=tmp/${sample}" SelectVariants \
        --variant  "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.vcf.gz \
        --output  "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.biallenic_snp.vcf.gz \
        --select-type-to-include  SNP \
        --exclude-non-variants true \
        --remove-unused-alternates true \
	--restrict-alleles-to BIALLELIC

bcftools view \
	--regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 \
	-i 'GT="RA"' \
	-o "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.autosome_biallelic_het_snp.vcf.gz \
	-Oz \
	"${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.biallenic_snp.vcf.gz

tabix -p vcf "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.autosome_biallelic_het_snp.vcf.gz
