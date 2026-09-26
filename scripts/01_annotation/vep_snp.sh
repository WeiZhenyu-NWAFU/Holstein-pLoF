#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${PROJECT_ROOT:?Set PROJECT_ROOT; see SOFTWARE.md}"
for i in "${PROJECT_ROOT}"/574_sample/01.vcf_split/snp/*.gz;do vep --cache --species bos_taurus -i $i -o `basename $i .vcf.gz`.vep_anno.vcf;done

