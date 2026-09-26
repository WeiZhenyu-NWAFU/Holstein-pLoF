for i in /storage/public/home/yuwang/liuanguo/master_thesis/lof_in_574_and_533_sample/574_sample/01.vcf_split/snp/*.gz;do vep --cache --species bos_taurus -i $i -o `basename $i .vcf.gz`.vep_anno.vcf;done

