/storage/public/apps/software/conda/Miniconda3-py310_23.3.1-0/envs/gatk-4.0.5.1_py35/bin/gatk SelectVariants \
	--variant  upper_bound.pLoF.vcf.gz \
	--output   upper_bound.pLoF.filter.vcf.gz \
	--max-indel-size 50 \
	--exclude-non-variants true \
        --remove-unused-alternates true

