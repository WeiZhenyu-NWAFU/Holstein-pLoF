#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${PROJECT_ROOT:?Set PROJECT_ROOT; see SOFTWARE.md}"
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA; see SOFTWARE.md}"
## based on GATK v4.3.0.0

sample=$1

gatk --java-options "-Xmx300G -XX:ParallelGCThreads=4  -Djava.io.tmpdir=tmp/${sample}" ASEReadCounter \
	-R "${REFERENCE_FASTA}" \
       	-I "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.mkdup.bam \
	-V "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.autosome_biallelic_het_snp.vcf.gz \
	-O "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.autosome_biallelic_het_snp.ASE.table
