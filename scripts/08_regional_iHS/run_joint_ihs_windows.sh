#!/usr/bin/env bash
set -euo pipefail

analysis_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source_file=${IHS_INPUT:?Set normalized iHS input file}
output_dir=${IHS_OUTDIR:?Set a new output directory}
python_bin=${PYTHON_BIN:-python3}

[[ ! -e "$output_dir" ]] || { echo "Output directory already exists" >&2; exit 1; }
mkdir -p "${output_dir}"
"${python_bin}" "${analysis_dir}/build_joint_ihs_windows.py" \
  --input "${source_file}" \
  --outdir "${output_dir}" \
  --windows 50000 100000 200000 \
  --extreme-threshold 2.0 \
  --min-snps 20 \
  --snp-bins 10
