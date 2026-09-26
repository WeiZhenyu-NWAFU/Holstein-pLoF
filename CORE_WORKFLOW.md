# Core workflow: commands, software and purpose

This is a compact guide to recovered project commands, not a new analysis or a one-command pipeline. Original scripts retain their source paths and options. Configure paths, environments and scheduler resources in a separate working copy before use; never run them over archived results or on a shared login node. They were read but not executed during curation.

## Workflow map

| Stage | Software / entry point | Work performed and required inputs |
|---|---|---|
| Genomic annotation | VEP and ANNOVAR: `server_core/annotation/`; SnpEff was also used in source commands | Annotate SNP/indel inputs against bovine annotations. Annotation databases and input VCF/avinput files are external. |
| Annotation collection | Source command excerpt below; GATK `filter_indel_length.sh` | Combine tool-specific position lists; apply maximum indel length of 50 bp. This is a length filter, not a 50-bp neighborhood/complex-locus filter. |
| DNA–RNA matching | bcftools `matching/matched_DNA_RNA_isec.sh` | Intersect per-individual DNA and RNA variant files with `isec -c none`; intersection alone does not establish genotype concordance. The lung script is the recovered example, not a newly reconstructed all-tissue pipeline. |
| Allelic read counts | GATK `ase/gatk_ASEReadCounter.sh` | Extract allele counts at supplied sites from RNA alignments. Site eligibility and post-count filtering must be supplied separately. |
| ASE tests | R `ase/repeated_95CI_with_p_value.R`, `ase/single_binom_test_2side.R` | Existing Repeated t-test/interval and Single binomial procedures. These are component scripts, not the complete final site-classification pipeline. Repeated records need not be distinct animals. |
| Sequence-based NMD annotation | VEP NMD plugin; exact command below | Annotate transcript consequences. This is distinct from observed ASE; positional classes and transcript ambiguity require downstream interpretation. |
| Milk RNA calls and ASE | GATK and bcftools, `server_core/milk/02` through `06` | Start from duplicate-marked RNA BAMs; split/recalibrate, call GVCFs, genotype, select autosomal biallelic heterozygous SNPs, count reads and intersect with genomic pLoFs. |
| Regional iHS | Python `regional_ihs/build_joint_ihs_windows.py`; shell runner | Consume already normalized iHS; compute extreme-SNP fractions and SNP-density-stratified empirical tails at 50/100/200 kb. The 100-kb scale is primary. No new phasing or selscan run is performed by this script. |
| Genealogy-based inference | Relate / CLUES: `server_core/clues/` | Sample branch lengths and perform locus-level inference from prepared genealogy inputs. Neutral calibration and upstream genealogy construction are not included. |
| Homozygote screen | Existing `analysis/section5_rebuild_20260921/audit_recompute.py` | Analyze supplied genotype counts and integrate existing evidence; this is not genotype-likelihood calling. Previously uploaded enrichment and one Fig.5 plot script remain available; no further summary or plotting scripts are added. |

## Important input distinctions

- The milk scripts infer genotypes from RNA. DNA-derived known-sites used for BQSR are not per-cow DNA genotype validation. The source ASE command uses a duplicate-marked BAM, whereas HaplotypeCaller uses the recalibrated BAM; these choices have not been silently changed.
- Original VEP output filenames ending in `.vcf` do not establish output format without the corresponding option. Check the actual output schema before downstream parsing.
- The NMD command below references a historical upstream annotation input. Its old record counts are not current manuscript counts. These commands do not redefine the final pLoF/mpLoF sets.
- The regional script uses non-overlapping 1-based windows, `|normalized iHS| >= 2`, at least 20 SNPs and 10 SNP-density strata, as specified in its supplied runner. Variant-to-window integration is a separate step.
- Source scripts contain site-specific paths and resource requests. Comments mentioning software versions are provenance, not verification of a runnable environment. No raw data or credentials are included.

## Scope intentionally left out

No additional plotting, intermediate summaries or complete command diaries are included. Upstream fastp/STAR/WGS-QC records contain incomplete or questionable command syntax and have not been promoted into a runnable pipeline. The recovered SnpEff command also needs its `-stats` argument checked. These are not silently repaired or represented as validated production commands.

Final phASER/eQTL/sQTL production scripts, full ASE postprocessing, CLUES calibration and complete genotype calling remain outside this compact subset. Their absence here does not mean those analyses were not performed. No thresholds, results or manuscript definitions were changed.

## Selected original command excerpts

The following are documentation excerpts, not an executable all-in-one script. Source paths are relative to `/storage/public/home/2020110005/liuanguo/master_thesis/lof_in_574_and_533_sample/`. Source hashes and line numbers identify the records inspected.

### 574_sample/02.annotation/01.vep/nmd_prediction/vep_cmd

Source SHA-256: `79f9957d6d97029a48e45f5b53af6b9fba517c171d36878d221b667229811388`

Line 2:

```sh
vep --cache --species bos_taurus -i upper_bound.pLoF.vcf.gz --plugin NMD
```

### 574_sample/02.annotation/04.venn_collect/v10_sLOF_a4/cmd

Source SHA-256: `a03469fdf8fc1dc57cb7974be813b9b0bcc4241942d58e3a3779a64de0d2fffb`

Line 10:

```sh
cat vep.indel.pLoF.chr_site.position snpeff.indel.pLoF.chr_site.position annovar.indel.pLoF.chr_site.position|sort -k1,1V -u > upper_bound.pLoF.indel.chr_site
```

Line 12:

```sh
cat vep.snp.pLoF.chr_site.position snpeff.snp.pLoF.chr_site.position annovar.snp.pLoF.chr_site.position |sort -k1,1V -u > upper_bound.pLoF.snp.chr_site
```

Line 15:

```sh
cat annovar.all.pLoF.chr_site.position snpeff.all.pLoF.chr_site.position vep.all.pLoF.chr_site.position |sort -k1,1V -u > upper_bound.chr_site.list
```
