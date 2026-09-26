## based on GATK v4.3.0.0

sample=$1

source /storage/public/home/2020060185/.bashrc
gatk --java-options "-Xmx300G -XX:ParallelGCThreads=4  -Djava.io.tmpdir=tmp/${sample}" ASEReadCounter \
	-R /storage/public/home/2022060207/01.Dairycattle/00.Genome/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa \
       	-I /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.mkdup.bam \
	-V /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.autosome_biallelic_het_snp.vcf.gz \
	-O /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.autosome_biallelic_het_snp.ASE.table
