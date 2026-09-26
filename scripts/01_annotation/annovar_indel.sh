file=$1
annotate_variation.pl \
 /storage/public/home/yuwang/liuanguo/master_thesis/lof_in_574_and_533_sample/574_sample/01.vcf_split/indel/avinput/$file \
 /storage/public/home/yuwang/biosoftware/annovar/cattledb/ \
 --buildver ARS-UCD1.2 \
 --outfile `basename $file .avinput`
