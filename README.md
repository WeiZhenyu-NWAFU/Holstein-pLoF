# Holstein pLoF study: core analysis scripts

Initial code subset for *Matched multi-tissue transcriptomics and population genetics characterize putative loss-of-function variation in Holstein cattle*.

## Contents

- `analysis/section5_rebuild_20260921/audit_recompute.py`: existing autosomal genotype-count HWE screen, BH correction, additional-Holstein evaluation and integration of precomputed regional-iHS/CLUES results.
- `analysis/section5_rebuild_20260921/enrichment.py`: matched synonymous controls, GO/Reactome enrichment and matched sensitivity analyses.
- `figure5/scripts/plot_figure5.R`: Figure 5 plotting from an aggregate source table, without statistical fitting.

These are unchanged copies of existing project scripts. They have not been rerun as part of repository preparation. This initial release contains **code and documentation only**, no research data. It is not a complete WGS-to-manuscript workflow.

## Inputs and dependencies

See `INPUTS.md` for the original input layout. All listed analysis inputs and `figure5/source_tables/source_plot_data.tsv` must be supplied separately; none is included here. The figure input is an aggregate plotting table, not raw sequencing data.

The Python scripts require NumPy and pandas. The R script uses grid and Cairo-capable graphics, Arial and a Windows English UTF-8 locale setting. Exact production software versions remain to be documented; version numbers have not been guessed. Local font/locale availability may require configuration.

After providing the required inputs in a separate working copy, the original command order is:

```sh
python analysis/section5_rebuild_20260921/audit_recompute.py
python analysis/section5_rebuild_20260921/enrichment.py
Rscript figure5/scripts/plot_figure5.R
```

The Python scripts create/overwrite derived files under `analysis/section5_rebuild_20260921/tables/`; do not execute over archived results. For plotting, create `figure5/source_tables/` and supply `source_plot_data.tsv` first. Plot outputs are written under `figure5/`. The transformation from analysis output to the plotting input is not included in this subset. The R script reproduces the archived plot, not subsequent author-made typography edits.

## Scope and interpretation

The HWE script implements the existing genotype-count method, not a new genotype-likelihood analysis. It consumes existing iHS/CLUES results; it does not perform phasing, regional scans, genealogy inference or neutral calibration. Historical 14-site annotations are provenance fields, not a reinstated priority tier. Annotation, ASE/NMD, QTL and upstream selection workflows will be documented separately as available.

Software licensing and a permanent archived release remain subject to author approval. This repository should not yet be described as a complete reproducible pipeline or a DOI-archived release.
