#!/usr/bin/env bash
set -euo pipefail
# Inputs: aligned BAM with STAR/WASP tags, phased VCF, mappability blacklist,
# and gene-feature BED. Activate the compatible phASER/Python 2 environment first.
: "${PHASER_DIR:?Set phASER installation directory}"
: "${SAMPLE:?Set VCF sample ID}"
: "${BAM:?Set aligned RNA BAM}"
: "${VCF:?Set phased population VCF}"
: "${BLACKLIST:?Set haplotype-count blacklist BED}"
: "${FEATURES:?Set gene-feature BED}"
: "${OUTDIR:?Set a new sample output directory}"
[[ ! -e "$OUTDIR" ]] || { echo 'Output directory already exists' >&2; exit 1; }
mkdir -p "$OUTDIR"
# The source filter rejects vW codes 2–7; absence of a tag is not excluded.
# This filter consumes existing WASP tags; it does not run WASP remapping.
samtools view -h -q 255 "$BAM" | awk '$0 !~ /vW:i:[2-7]/' | samtools view -b > "$OUTDIR/filtered.bam"
samtools index "$OUTDIR/filtered.bam"
gatk AddOrReplaceReadGroups -I "$OUTDIR/filtered.bam" -O "$OUTDIR/readgroups.bam" \
  -RGID 4 -RGLB lib1 -RGPL illumina -RGPU run -RGSM 20 \
  -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate
gatk MarkDuplicates -I "$OUTDIR/readgroups.bam" -O "$OUTDIR/duplicates_marked.bam" \
  -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT --READ_NAME_REGEX null \
  -M "$OUTDIR/duplicate_metrics.txt"
mkdir "$OUTDIR/tmp"
python2 "$PHASER_DIR/phaser/phaser.py" \
  --vcf "$VCF" --bam "$OUTDIR/duplicates_marked.bam" \
  --paired_end 1 --mapq 255 --baseq 10 --sample "$SAMPLE" \
  --haplo_count_blacklist "$BLACKLIST" --threads 1 --pass_only 0 \
  --temp_dir "$OUTDIR/tmp" --o "$OUTDIR/$SAMPLE"
python2 "$PHASER_DIR/phaser_gene_ae/phaser_gene_ae.py" \
  --haplotypic_counts "$OUTDIR/$SAMPLE.haplotypic_counts.txt" \
  --features "$FEATURES" --o "$OUTDIR/${SAMPLE}_phaser.gene_ae.txt"
