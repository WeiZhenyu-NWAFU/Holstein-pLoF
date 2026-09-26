##based GATK) v4.3.0.0
	
sample=$1

source /storage/public/home/2020060185/.bashrc

##Step 3. Split N trim and reassign mapping qualities.
gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" SplitNCigarReads \
	--spark-runner LOCAL \
	-R /storage/public/home/2022060207/01.Dairycattle/00.Genome/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa \
	-I /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.mkdup.bam \
	-O /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.cigar.bam \
	--create-output-bam-index true

## Step 4. Base recalibration (BQSR).
gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" BaseRecalibrator \
	--spark-runner LOCAL \
	-R /storage/public/home/2022060207/01.Dairycattle/00.Genome/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa \
	-I /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.cigar.bam \
	-O /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.bqsr.table \
	--known-sites /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/01.WGS/574_sample_var_for_RNA_BQSR/milk.574sample.clean.snp.recode.vcf.gz \
	--known-sites /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/01.WGS/574_sample_var_for_RNA_BQSR/milk.574sample.clean.indel.recode.vcf.gz

gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" ApplyBQSR \
	--spark-runner LOCAL \
	-R /storage/public/home/2022060207/01.Dairycattle/00.Genome/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa \
	-I /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.cigar.bam \
	-O /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.bqsr.bam \
	--bqsr-recal-file /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.bqsr.table

## Step 5. Run the haplotypecaller.
gatk --java-options "-Xmx300G -XX:ParallelGCThreads=8 -Djava.io.tmpdir=tmp/${sample}" HaplotypeCaller \
	-R /storage/public/home/2022060207/01.Dairycattle/00.Genome/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa \
       	-I /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.bqsr.bam \
	-O /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.g.vcf.gz \
	--dont-use-soft-clipped-bases \
	--emit-ref-confidence GVCF \
	-stand-call-conf 20
