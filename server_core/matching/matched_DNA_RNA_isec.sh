sample=$1

/storage/public/apps/software/bcftools/bcftools-1.17/bin/bcftools isec \
	-c none \
	-p lung-${sample}/matched_pLoF \
	/storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/533_sample/pLoF_variant/533_indvduals_vcf/533_sample_private_variants_vcf/liuwei-chongcexu-${sample}/liuwei-chongcexu-${sample}.vcf.gz \
	lung-${sample}/lung-${sample}.vcf.gz
