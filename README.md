# Holstein pLoF: core processing workflow

Core commands and scripts for the study of putative loss-of-function variation in Holstein cattle. The repository is organized by analysis stage and documents the software, inputs and purpose of each step.

| Stage | Software | Entry point | Purpose |
|---|---|---|---|
| 1. Variant annotation | VEP, ANNOVAR, GATK; SnpEff in annotation collection | `scripts/01_annotation/` | Annotate SNPs and indels; apply the indel-length filter |
| 2. Genome–transcriptome matching | bcftools | `scripts/02_genome_transcriptome_matching/` | Intersect matched-individual DNA and RNA variant files |
| 3. Allelic read counting | GATK ASEReadCounter | `scripts/03_ASE/` | Obtain reference and alternative read counts at supplied sites |
| 4. Milk RNA processing | GATK, bcftools | `scripts/04_milk_RNA/` | Recalibrate RNA alignments, call/genotype variants, select heterozygous SNPs and count allelic reads |
| 5. Homozygote depletion | Python, NumPy, pandas | `scripts/05_homozygote_depletion/` | Screen autosomal genotype counts against HWE expectations |
| 6. Regional iHS | Python; normalized selscan input | `scripts/06_regional_iHS/` | Construct regional extreme-iHS fractions and empirical tail probabilities |
| 7. Genealogical selection inference | Relate, CLUES | `scripts/07_CLUES/` | Sample branch lengths and infer locus-level selection parameters |

Sequence-based NMD annotation and annotation-list collection commands are included in [WORKFLOW.md](WORKFLOW.md).

## Use

Read [WORKFLOW.md](WORKFLOW.md) for input requirements, software and command order. Configure reference resources, paths and compute resources before running the shell scripts in a separate working directory. Run computational jobs through the local HPC scheduler, not on a shared login node.

The repository contains the core technical workflow, not datasets, figure-generation scripts, enrichment analyses or manuscript-table assembly. Raw alignment/calling, final QTL/phASER workflows, full ASE classification and CLUES neutral calibration are not bundled. It is not a one-command reproduction of every manuscript result.

The Python dependencies can be installed with `python -m pip install -r requirements.txt`; the regional-iHS script requires Python 3.9 or later. External bioinformatics tools and reference databases must be installed separately. Dependency versions are not pinned as a validated environment has not yet been packaged. There is no license grant or archived DOI release at present.
