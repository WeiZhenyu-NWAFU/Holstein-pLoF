#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${PROJECT_ROOT:?Set PROJECT_ROOT; see SOFTWARE.md}"
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA; see SOFTWARE.md}"
##based GATK) v4.3.0.0

sample=$1


gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" GenotypeGVCFs \
	-R "${REFERENCE_FASTA}" \
	-V "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.g.vcf.gz \
	-O "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.vcf.gz \
	-stand-call-conf 20.0



