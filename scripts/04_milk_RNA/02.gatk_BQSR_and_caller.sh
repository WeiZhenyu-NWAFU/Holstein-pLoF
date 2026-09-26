#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${PROJECT_ROOT:?Set PROJECT_ROOT; see SOFTWARE.md}"
: "${REFERENCE_FASTA:?Set REFERENCE_FASTA; see SOFTWARE.md}"
##based GATK) v4.3.0.0
	
sample=$1


##Step 3. Split N trim and reassign mapping qualities.
gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" SplitNCigarReads \
	--spark-runner LOCAL \
	-R "${REFERENCE_FASTA}" \
	-I "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.mkdup.bam \
	-O "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.cigar.bam \
	--create-output-bam-index true

## Step 4. Base recalibration (BQSR).
gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" BaseRecalibrator \
	--spark-runner LOCAL \
	-R "${REFERENCE_FASTA}" \
	-I "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.cigar.bam \
	-O "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.bqsr.table \
	--known-sites "${PROJECT_ROOT}"/01.data/01.WGS/574_sample_var_for_RNA_BQSR/milk.574sample.clean.snp.recode.vcf.gz \
	--known-sites "${PROJECT_ROOT}"/01.data/01.WGS/574_sample_var_for_RNA_BQSR/milk.574sample.clean.indel.recode.vcf.gz

gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" ApplyBQSR \
	--spark-runner LOCAL \
	-R "${REFERENCE_FASTA}" \
	-I "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.cigar.bam \
	-O "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.bqsr.bam \
	--bqsr-recal-file "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.bqsr.table

## Step 5. Run the haplotypecaller.
gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" HaplotypeCaller \
	-R "${REFERENCE_FASTA}" \
       	-I "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.bqsr.bam \
	-O "${PROJECT_ROOT}"/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.g.vcf.gz \
	--dont-use-soft-clipped-bases \
	--emit-ref-confidence GVCF \
	-stand-call-conf 20
