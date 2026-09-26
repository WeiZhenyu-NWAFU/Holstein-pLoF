# Input dependencies of the existing analysis scripts

Paths below are relative to the repository's `analysis/` directory. They retain the original directory labels to document source dependencies, not to designate every historical file as a current published result. These inputs are **not bundled** in this initial code-only subset.

## audit_recompute.py

- `fig4_final/data/variant_level.charlier_lrt.tsv.gz`: genotype-count/statistical source, including cohorts and variant classes.
- `manuscript_supporting_information_final/phase1/variant_master.tsv`: annotation/allele master used by the existing script.
- `fig4_v3_final/data/ARS-USD1.2_ENSBTAG-gene_symbol`: source gene-symbol mapping (original filename retained).
- `fig4_v3_final/data/larger_holstein_replication_1148.raw.tsv`: allele-specific external cohort summary.
- `fig4_v2/qc/regional_iHS_100kb_pLoF_annotation.tsv`: precomputed regional-iHS results.
- `fig4_v2/qc/clues_all_evaluable_logLR_definition.tsv`: precomputed CLUES results.
- `fig4_v3_final/qc/stringent14_complete_candidate_annotation.tsv`: existing technical annotation of the historical 14-site set; not a current candidate-tier definition.

Important outputs include `01_depletion_candidate_evidence_matrix.tsv`, `01_autosomal_synonymous.tsv.gz`, `02_larger_Holstein_evaluation.tsv`, `04_regional_iHS_frozen.tsv` and `04_CLUES_frozen.tsv` under `section5_rebuild_20260921/tables/`.

## enrichment.py

Requires `section5_rebuild_20260921/tables/01_depletion_candidate_evidence_matrix.tsv` and `01_autosomal_synonymous.tsv.gz` from the preceding script, plus:

- `fig4_v3_final/data/ARS-USD1.2_ENSBTAG-gene_symbol`
- `fig4_v3_final/data/go-basic.obo`
- `fig4_v3_final/data/goa_cow.gaf.gz`
- `fig4_v3_final/data/Ensembl2Reactome.txt`

The exact resource releases and redistribution permissions must be documented before supplying these annotation resources. Do not substitute current downloads and claim they are the original versions. Main outputs are matched-synonymous effects and matches, GO/Reactome enrichment, matched sensitivities, and gene-set/background membership tables.

## plot_figure5.R

The input `figure5/source_tables/source_plot_data.tsv` is required but is not included in this code-only repository. Its columns and per-panel records are the archived plotting input; there is no statistical fitting in the plotting script. No individual-level sequencing data are included.
