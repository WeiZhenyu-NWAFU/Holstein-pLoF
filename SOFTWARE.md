# Software and path configuration

These are core analysis commands for a Linux/HPC environment. Run jobs through your scheduler. Activate the appropriate tools before each module; scripts do not activate another user's environment. Use separate working directories and preserve input data. Some original commands write predictable output names, so do not rerun them over archived results.

## Dependencies

| Module | Dependency | Version information available |
|---|---|---|
| Annotation | VEP and bovine cache; ANNOVAR and bovine database; VEP NMD plugin | Exact releases not recorded in the supplied commands |
| Indel-length filtering | GATK | 4.0.5.1 identified in the original executable path |
| DNA–RNA matching / milk SNP selection | bcftools; tabix | bcftools 1.17 identified in the command path; tabix version not recorded |
| Milk RNA processing | GATK / Java | GATK 4.3.0.0 stated in source-script comments |
| GATK ASE | GATK / Java | General ASE example does not specify a release; milk example states 4.3.0.0 |
| phASER | Python 2, phASER, samtools, GATK | phASER commit and samtools release not recorded |
| QTL genotype preparation | PLINK | 1.90b6.21 identified in command path |
| eQTL / sQTL mapping | OmiGA | v1.1.3-beta.2+250831 identified in command path |
| Significant-pair extraction | csvtk, gzip, awk | Exact releases not recorded |
| Homozygote depletion | Python 3, NumPy, pandas | Python packages listed in requirements.txt; no validated version lock |
| Regional iHS | Python >=3.9 standard library; precomputed normalized iHS | No additional Python dependencies |
| Genealogical inference | Relate and CLUES | Relate v1.2.1 identified in command path; CLUES commit not recorded |

These are source-level version records, not a claim that one combined environment was tested. In particular, use separate environments for Python 2 phASER and Python 3 analysis. Install third-party tools under their own terms; their source code and databases are not redistributed here.

## Configuration variables

Export only the variables needed for the module you run:

```sh
export PROJECT_ROOT=/path/to/project
export REFERENCE_FASTA=/path/to/ARS-UCD1.2.fa
export ANNOVAR_DB=/path/to/annovar_database
export RELATE_DIR=/path/to/relate
export GENEALOGY_DIR=/path/to/prepared_holstein_genealogies
export CLUES_DIR=/path/to/clues
export IHS_INPUT=/path/to/normalized_iHS.txt
export IHS_OUTDIR=/path/to/new_iHS_output
```

`PROJECT_ROOT` represents the existing relative layout expected by annotation, matching and milk scripts, including `574_sample/`, `533_sample/` and `01.data/`. Configure or adapt this layout to your own inputs without changing allele definitions. Relative filenames and analysis flags are retained from the source commands.

The regional-iHS runner locates its Python script beside itself; `PYTHON_BIN` optionally overrides `python3`. Relate/CLUES scripts retain the `Chr${chrom}_Holstein_reinfer` genealogy naming and expect a prepared `timeBins.txt` and `clues_infer/` output layout. Other per-module variables and command arguments are documented in WORKFLOW.md.

## Checks and scope

Shell syntax, Python parsing and command-line help were checked during packaging. This does not validate biological results or demonstrate end-to-end reproduction. No research dataset was rerun for publication of this repository. Assembly-compatible reference files, indexes, phenotype/covariate matrices and phased genotype inputs must be supplied separately.
