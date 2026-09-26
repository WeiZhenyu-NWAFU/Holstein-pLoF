#!/usr/bin/env bash
# Activate the required tools before running; see SOFTWARE.md.
: "${PROJECT_ROOT:?Set PROJECT_ROOT; see SOFTWARE.md}"
: "${ANNOVAR_DB:?Set ANNOVAR_DB; see SOFTWARE.md}"
file=$1
annotate_variation.pl \
 "${PROJECT_ROOT}"/574_sample/01.vcf_split/indel/avinput/$file \
 "${ANNOVAR_DB}"/ \
 --buildver ARS-UCD1.2 \
 --outfile `basename $file .avinput`
