# Core commands and inputs

Paths below are relative to this repository. Shell scripts retain their analysis options; configure paths and environments as described in SOFTWARE.md. User-specific server paths and environment-activation commands have been replaced with configuration variables and tools on PATH. They are command templates for the corresponding stage, not a scheduler-managed end-to-end pipeline.

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

## 3. Allele-specific expression

### 3.1. GATK site-level counts

`scripts/03_ASE/GATK/gatk_ASEReadCounter.sh` consumes RNA BAMs, the reference and a site VCF. It returns allelic read counts, not the complete final ASE classifier. Coverage/eligibility filtering and final inference are separate steps.

### 3.2. phASER haplotype-based ASE

Configure `PHASER_DIR`, `SAMPLE`, `BAM`, `VCF`, `BLACKLIST`, `FEATURES` and a new `OUTDIR`, then run:

```sh
bash scripts/03_ASE/phASER/haplotype_counts.sh
```

The command retains MAPQ 255 alignments and removes alignments with WASP tags vW=2–7, adds read groups, marks duplicates, and runs phASER with paired-end mode, MAPQ 255, base quality 10, one thread, a haplotype-count blacklist and `--pass_only 0`. This last option disables a PASS-only restriction; it does not assert that every site passes quality control. `phaser_gene_ae.py` produces gene-level allelic counts. Python 2 and the compatible phASER environment are required. The source read-group constants (including RGSM 20) are preserved; phASER uses the explicit VCF `--sample` identifier.

For cis-variant aFC, configure `PHASER_DIR`, `GW_BED`, `VCF`, `PAIRS`, `SAMPLE_MAP` and `OUTPUT`, then run `scripts/03_ASE/phASER/cis_allelic_fold_change.sh`. The pair-file header is `gene_id, var_id, var_contig, var_pos, var_ref, var_alt` (tab-separated); the sample-map header is `vcf_sample, bed_sample`. The GW-phased BED must be prepared and indexed beforehand. Aggregation and pair-table assembly are intentionally omitted; preserve phase eligibility and allele alignment when preparing these inputs. aFC runs with 20 threads as in the supplied command.

## 4. Milk RNA processing

Run the scripts in `scripts/04_milk_RNA/` in numerical order, beginning with existing duplicate-marked RNA BAMs:

1. `02`: SplitNCigarReads, BaseRecalibrator, ApplyBQSR and HaplotypeCaller.
2. `03`: GenotypeGVCFs.
3. `04`: select autosomal biallelic heterozygous SNPs.
4. `05`: ASEReadCounter.
5. `06`: intersect RNA variants with the genomic pLoF catalogue.

These commands use RNA-derived genotypes. The DNA-derived known-sites used in BQSR are not independent per-cow DNA genotypes. The read-count command uses a duplicate-marked BAM, while HaplotypeCaller uses a recalibrated BAM; these input choices are preserved.

## 5–6. Milk cis-eQTL and cis-sQTL

The supplied mapping commands identify OmiGA **v1.1.3-beta.2+250831**. The eQTL genotype-conversion command uses PLINK **1.90b6.21**, `--maf 0.01`, `--chr-set 29` and `--keep-allele-order`. It is provided in `scripts/05_eQTL/prepare_genotypes.sh` with `VCF` and a new `OUTDIR` as environment variables. The sQTL source consumes prepared PLINK inputs; do not infer its upstream filters from the eQTL script alone.

Use separate eQTL/sQTL working directories. Required relative inputs are:

- `milk/genotypes/milk.{bed,bim,fam}`;
- `milk/covFile/milk.peer30.tsv`;
- eQTL: `milk/phenotypes/milk.expression.bed.gz`;
- sQTL: `milk/phenotypes/milk.leafcutter.sorted.bed.gz` and `milk/phenotypes/milk.leafcutter.phenotype_groups_omiga.txt`.

The main invocations, using the configured absolute repository path, are:

```sh
# From the eQTL work directory
bash "$REPO/scripts/05_eQTL/map_cis.sh" milk 30 20 results_eqtl
bash "$REPO/scripts/05_eQTL/extract_significant_pairs.sh" milk results_eqtl 4
# From the separate sQTL work directory
bash "$REPO/scripts/06_sQTL/map_cis.sh" milk 30 20 results_sqtl
bash "$REPO/scripts/06_sQTL/extract_significant_pairs.sh" milk results_sqtl 1
```

Both mapping scripts retain `--calcu-variant-threshold`, `--force-double-precision`, cis and cis_independent modes; sQTL additionally uses `--pheno-group`. The recorded invocations use 30 PEER factors and 20 mapping threads. The scripts do not explicitly set a cis-window, MAC threshold or phenotype normalization procedure; these must not be inferred from filenames or documented as explicit command options.

Significant-pair extraction preserves the source filters: the final cis-summary field <=0.05, the named `pval_g1_threshold` joined by `pheno_id`, and pair-table field 7 <= joined field 8 for chromosomes 1–29. These positional filters depend on the original OmiGA output schema; check headers before use with another version. No independent phenotype-normalization, LeafCutter or PEER-estimation pipeline is included.

These modules are concise, path-configurable adaptations of the supplied commands. Analysis flags are retained; automatic recursive cleanup, verbose status messages, sample-specific retry lists and summary/plotting code are omitted. New output directories are required to protect existing results. The sample-flow path inconsistency in the supplied BAM preprocessing commands is replaced by one explicit input/output chain; this packaging change has not been validated by rerunning the analyses.
## 7. Homozygote-depletion screen

```sh
python scripts/07_homozygote_depletion/screen_homozygotes.py --input genotype_counts.tsv.gz --output homozygote_screen.tsv
```

Input columns: `variant`, `chrom`, `cohort`, `variant_class`, `n_refhom`, `n_het`, `n_althom`, `n_called`, `call_rate`. Chromosomes must be coded `1`–`29` for autosomal records. Each variant/cohort/class combination must occur once. For pLoF rows, ALT must identify the pLoF allele; do not substitute minor-allele counts. Supply all tested sites, not just preselected significant candidates.

The core calculation is unchanged from the study script: estimate ALT frequency from genotype counts, calculate HWE-expected ALT homozygotes, calculate the implemented binomial likelihood-ratio statistic and chi-square-tail P value, and apply BH within cohort and variant class. Depletion requires observed < expected and FDR < 0.05. Informative absence requires observed = 0, expected >= 3 and call rate >= 0.95; it is distinct from significant depletion. X-linked loci are excluded. This is a genotype-count screen, not a genotype-likelihood model. The CLI/input validation are a packaging refactor; candidate integration, enrichment and table/figure generation are deliberately omitted.

## 8. Regional iHS

```sh
python scripts/08_regional_iHS/build_joint_ihs_windows.py --input normalized_iHS.txt --outdir results/regional_iHS --windows 50000 100000 200000 --extreme-threshold 2.0 --min-snps 20 --snp-bins 10
```

Input header: `LocusID Position Freq_Derived iHH1 iHH0 iHS_Raw iHS_Norm Is_Sig` (tab-separated; see parser). Supply autosomes 1–29. The script uses non-overlapping 1-based windows, fraction with absolute normalized iHS >= 2, minimum 20 SNPs and 10 SNP-density strata. The 100-kb scale is primary; 50/200 kb provide sensitivity. Upstream phasing/selscan normalization and downstream pLoF-to-window integration are separate.

## 9. Relate / CLUES

Run `scripts/09_CLUES/02.sampleBranchLength.sh` followed by `03.clues_infer.sh` after supplying the required genealogy inputs and configuring the Relate/CLUES installations. The repository does not include complete genealogy construction, ancestral-state harmonization or neutral calibration. Inference output is not automatically a calibrated selection-support classification.

## Execution scope

Script reorganization does not constitute a new analysis. No research data were processed while preparing this repository. Configure environments, reference/database releases and resource settings before execution; the scripts have not been validated end to end on an independent installation.

