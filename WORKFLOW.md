# Core commands and inputs

Paths below are relative to this repository. Shell scripts retain their analysis options and require local paths/environment configuration. They are command templates for the corresponding stage, not a scheduler-managed end-to-end pipeline.

## 1. Annotation and NMD

Run the VEP/ANNOVAR SNP and indel scripts in `scripts/01_annotation/` with the corresponding VCF or avinput files and bovine databases. VEP output format must be checked explicitly; a `.vcf` filename alone does not ensure VCF output. SnpEff annotations also contributed to collection, but its recovered invocation has an unresolved `-stats` argument and is not supplied as a runnable script.

Combine the tool-specific position lists after compatible coordinate/allele representation has been established:

```sh
cat annovar.all.pLoF.chr_site.position snpeff.all.pLoF.chr_site.position vep.all.pLoF.chr_site.position | sort -k1,1V -u > upper_bound.chr_site.list
```

`filter_indel_length.sh` uses GATK SelectVariants with `--max-indel-size 50`. This is an indel-length filter, not a filter for complex loci within a 50-bp neighborhood.

Sequence-based annotation uses the VEP NMD plugin:

```sh
vep --cache --species bos_taurus -i upper_bound.pLoF.vcf.gz --plugin NMD
```

Transcript-level predictions and observed ASE are distinct. This command does not perform final transcript reconciliation or ASE-based classification.

## 2. Genome–transcriptome matching

`scripts/02_genome_transcriptome_matching/matched_DNA_RNA_isec.sh` uses bcftools `isec -c none` for the same individual's DNA and RNA files. The supplied script is the lung example; corresponding inputs must be configured for other tissues. A shared variant record is not, by itself, proof of concordant genotypes.

## 3. ASE read counts

`scripts/03_ASE/gatk_ASEReadCounter.sh` consumes RNA BAMs, the reference and a site VCF. It returns allelic read counts, not the complete final ASE classifier. Coverage/eligibility filtering and final inference are separate steps.

## 4. Milk RNA processing

Run the scripts in `scripts/04_milk_RNA/` in numerical order, beginning with existing duplicate-marked RNA BAMs:

1. `02`: SplitNCigarReads, BaseRecalibrator, ApplyBQSR and HaplotypeCaller.
2. `03`: GenotypeGVCFs.
3. `04`: select autosomal biallelic heterozygous SNPs.
4. `05`: ASEReadCounter.
5. `06`: intersect RNA variants with the genomic pLoF catalogue.

These commands use RNA-derived genotypes. The DNA-derived known-sites used in BQSR are not independent per-cow DNA genotypes. The read-count command uses a duplicate-marked BAM, while HaplotypeCaller uses a recalibrated BAM; these input choices are preserved.

## 5. Homozygote-depletion screen

```sh
python scripts/05_homozygote_depletion/screen_homozygotes.py --input genotype_counts.tsv.gz --output homozygote_screen.tsv
```

Input columns: `variant`, `chrom`, `cohort`, `variant_class`, `n_refhom`, `n_het`, `n_althom`, `n_called`, `call_rate`. Chromosomes must be coded `1`–`29` for autosomal records. Each variant/cohort/class combination must occur once. For pLoF rows, ALT must identify the pLoF allele; do not substitute minor-allele counts. Supply all tested sites, not just preselected significant candidates.

The core calculation is unchanged from the study script: estimate ALT frequency from genotype counts, calculate HWE-expected ALT homozygotes, calculate the implemented binomial likelihood-ratio statistic and chi-square-tail P value, and apply BH within cohort and variant class. Depletion requires observed < expected and FDR < 0.05. Informative absence requires observed = 0, expected >= 3 and call rate >= 0.95; it is distinct from significant depletion. X-linked loci are excluded. This is a genotype-count screen, not a genotype-likelihood model. The CLI/input validation are a packaging refactor; candidate integration, enrichment and table/figure generation are deliberately omitted.

## 6. Regional iHS

```sh
python scripts/06_regional_iHS/build_joint_ihs_windows.py --input normalized_iHS.txt --outdir results/regional_iHS --windows 50000 100000 200000 --extreme-threshold 2.0 --min-snps 20 --snp-bins 10
```

Input header: `LocusID Position Freq_Derived iHH1 iHH0 iHS_Raw iHS_Norm Is_Sig` (tab-separated; see parser). Supply autosomes 1–29. The script uses non-overlapping 1-based windows, fraction with absolute normalized iHS >= 2, minimum 20 SNPs and 10 SNP-density strata. The 100-kb scale is primary; 50/200 kb provide sensitivity. Upstream phasing/selscan normalization and downstream pLoF-to-window integration are separate.

## 7. Relate / CLUES

Run `scripts/07_CLUES/02.sampleBranchLength.sh` followed by `03.clues_infer.sh` after supplying the required genealogy inputs and configuring the Relate/CLUES installations. The repository does not include complete genealogy construction, ancestral-state harmonization or neutral calibration. Inference output is not automatically a calibrated selection-support classification.

## Execution scope

Script reorganization does not constitute a new analysis. No research data were processed while preparing this repository. Configure environments, reference/database releases and resource settings before execution; the scripts have not been validated end to end on an independent installation.
