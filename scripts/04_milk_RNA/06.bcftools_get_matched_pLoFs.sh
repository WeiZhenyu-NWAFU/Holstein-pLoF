sample=$1

/storage/public/apps/software/bcftools/bcftools-1.17/bin/bcftools isec \
	-c none \
	-p ${sample}/matched_pLoF \
	/storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/574_sample/02.annotation/04.venn_collect/v10_sLOF_a4/upper_bound.pLoF.vcf.gz \
	${sample}/${sample}.vcf.gz
