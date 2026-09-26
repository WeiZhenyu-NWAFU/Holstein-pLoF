# Holstein pLoF: core processing workflow

Core commands and scripts for the study of putative loss-of-function variation in Holstein cattle. The repository follows the manuscript sequence: genomic annotation, transcriptomic evidence, and population-genetic analyses. Modules are grouped by purpose rather than constituting a single execution order. It is organized by analysis stage and documents the software, inputs and purpose of each step.

| Stage | Software | Entry point | Purpose |
|---|---|---|---|
| 1. Variant annotation | VEP, ANNOVAR, GATK; SnpEff in annotation collection | `scripts/01_annotation/` | Annotate SNPs and indels; apply the indel-length filter |
| 2. Genome–transcriptome matching | bcftools | `scripts/02_genome_transcriptome_matching/` | Intersect matched-individual DNA and RNA variant files |
| 3. ASE: GATK | GATK ASEReadCounter | `scripts/03_ASE/GATK/` | Obtain site-level reference and alternative read counts |
| 3. ASE: phASER | samtools, GATK, phASER | `scripts/03_ASE/phASER/` | Haplotype/gene allelic counts and cis-variant aFC |
| 4. Milk RNA processing | GATK, bcftools | `scripts/04_milk_RNA/` | Recalibrate RNA alignments, call/genotype variants, select heterozygous SNPs and count allelic reads |
| 5. cis-eQTL | PLINK, OmiGA, csvtk | `scripts/05_eQTL/` | Genotype conversion, cis mapping, independent signals and significant pairs |
| 6. cis-sQTL | OmiGA, csvtk | `scripts/06_sQTL/` | Map prepared LeafCutter phenotypes with phenotype groups |
| 7. Homozygote depletion | Python, NumPy, pandas | `scripts/07_homozygote_depletion/` | Screen autosomal genotype counts against HWE expectations |
| 8. Regional iHS | Python; normalized selscan input | `scripts/08_regional_iHS/` | Construct regional extreme-iHS fractions and empirical tail probabilities |
| 9. Genealogical selection inference | Relate, CLUES | `scripts/09_CLUES/` | Sample branch lengths and infer locus-level selection parameters |

Sequence-based NMD annotation and annotation-list collection commands are included in [WORKFLOW.md](WORKFLOW.md).

## Use

Read [WORKFLOW.md](WORKFLOW.md) for input requirements, software and command order. Configure reference resources, paths and compute resources before running the shell scripts in a separate working directory. Run computational jobs through the local HPC scheduler, not on a shared login node.

The repository contains the core technical workflow, not datasets, figure-generation scripts, enrichment analyses or manuscript-table assembly. Raw alignment/calling, phenotype normalization/PEER estimation, count-matrix assembly, full ASE classification and CLUES neutral calibration are not bundled. The QTL/phASER examples document the supplied milk workflows; other-tissue production commands are not inferred. It is not a one-command reproduction of every manuscript result.

The Python dependencies can be installed with `python -m pip install -r requirements.txt`; the regional-iHS script requires Python 3.9 or later. External bioinformatics tools and reference databases must be installed separately. Dependency versions are not pinned as a validated environment has not yet been packaged. There is no license grant or archived DOI release at present.
