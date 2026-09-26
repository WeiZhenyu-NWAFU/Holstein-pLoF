## based GATK) v4.3.0.0

sample=$1

source /storage/public/home/2020060185/.bashrc

gatk --java-options "-Xmx300G -XX:ParallelGCThreads=1  -Djava.io.tmpdir=tmp/${sample}" SelectVariants \
        --variant  /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.vcf.gz \
        --output  /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.biallenic_snp.vcf.gz \
        --select-type-to-include  SNP \
        --exclude-non-variants true \
        --remove-unused-alternates true \
	--restrict-alleles-to BIALLELIC

/storage/public/apps/software/bcftools/bcftools-1.17/bin/bcftools view \
	--regions 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 \
	-i 'GT="RA"' \
	-o /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.autosome_biallelic_het_snp.vcf.gz \
	-Oz \
	/storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.biallenic_snp.vcf.gz

tabix -p vcf /storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/01.data/02.Transcriptome/533_sample/milk/${sample}/${sample}.autosome_biallelic_het_snp.vcf.gz
