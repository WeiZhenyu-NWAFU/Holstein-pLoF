sample=$1
gatk ASEReadCounter \
	--reference /storage/public/home/2022060207/01.Dairycattle/00.Genome/Bos_taurus.ARS-UCD1.2.dna_sm.toplevel.fa \
	--input $sample.bam \
	--variant sites.vcf.gz \
	--output $sample.ASE.table \
	--showHidden
