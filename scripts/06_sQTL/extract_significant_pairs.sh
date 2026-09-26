#!/usr/bin/env bash
set -euo pipefail
tis=${1:?Usage: bash extract_significant_pairs.sh milk OUTPUT_DIR THREADS}
output_dir=${2:?Set completed map_cis.sh output directory}
threads=${3:?Set thread count}
[[ -d "$output_dir" ]] || { echo 'Missing result directory' >&2; exit 1; }
[[ ! -e "${output_dir}/${tis}.pval_nominal_threshold" && ! -e "${output_dir}/${tis}.cis_qtl_pairs.sig.txt.gz" ]] || { echo 'Output already exists' >&2; exit 1; }
zcat ${output_dir}/${tis}.cis_qtl.txt.gz | awk 'NR==1||$NF<=0.05' | csvtk -t -j $threads cut -f 'pheno_id,pval_g1_threshold' > ${output_dir}/${tis}.pval_nominal_threshold
zcat ${output_dir}/${tis}.cis_qtl_pairs.chr1.txt.gz | csvtk join -t -j $threads -f 'pheno_id;pheno_id' - ${output_dir}/${tis}.pval_nominal_threshold | awk '$7<=$8' | gzip -c > ${output_dir}/${tis}.cis_qtl_pairs.sig.txt.gz
for chr in chr{2..29}
do
	echo $chr
	zcat ${output_dir}/${tis}.cis_qtl_pairs.${chr}.txt.gz | csvtk join -t -j $threads -f 'pheno_id;pheno_id' - ${output_dir}/${tis}.pval_nominal_threshold | awk '$7<=$8' | sed '1d' | gzip -c >> ${output_dir}/${tis}.cis_qtl_pairs.sig.txt.gz
done
