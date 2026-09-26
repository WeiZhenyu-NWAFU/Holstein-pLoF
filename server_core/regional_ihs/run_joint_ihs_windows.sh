#!/usr/bin/env bash
set -euo pipefail

analysis_dir=/storage/public/home/2020110005/weizhenyu/01Dairycattle/09.Evolution/iHS
source_file=/storage/public/home/2022060207/wangtao/04.1000+our/03.iHS/04.ihs_norm/02.merge_norm/header_all_merge_ChrAuto_snps.selscan.ihs.out.100bins.norm.txt
output_dir=${analysis_dir}/window_joint_norm_20260822
python_bin=/storage/public/apps/software/conda/Miniconda3-py310_23.3.1-0/bin/python3

mkdir -p "${output_dir}"
"${python_bin}" "${analysis_dir}/build_joint_ihs_windows.py" \
  --input "${source_file}" \
  --outdir "${output_dir}" \
  --windows 50000 100000 200000 \
  --extreme-threshold 2.0 \
  --min-snps 20 \
  --snp-bins 10
