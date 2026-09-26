# Selected HPC analysis scripts

Five small analysis scripts copied from the project's `lof_in_574_and_533_sample` directory. Source code is unchanged apart from filenames used to make the repository easier to navigate. `SOURCE_MANIFEST.tsv` records the original relative paths and SHA-256 checksums. Files were read only; no server analysis or script execution was performed during curation.

## ASE

- `ase/gatk_ASEReadCounter.sh`: counts allele-supporting reads from a BAM and a supplied `sites.vcf.gz`. The reference path must be configured. This command does not by itself establish DNA heterozygosity or implement the entire eligibility workflow.
- `ase/repeated_95CI_with_p_value.R`: accepts an input file of RAR values and an output path; computes the mean, t-based 95% CI and a one-sample test against 0.5. Example: `Rscript repeated_95CI_with_p_value.R input_RAR.txt output.tsv`. Repeated observations are not automatically independent animals. Constant or insufficient observations can make the original t-test undefined; the script does not implement a special-case rescue or the final directional classification.
- `ase/single_binom_test_2side.R`: reads `uniq.snp.matched_pLoF.stat`, using `ref_count` and `total_count` for two-sided binomial tests against 0.5. Its original absolute output path must be changed in a working copy before use. This is a per-record test, not a complete multiple-testing or final effect-classification procedure.

The source `reads_count_filter/cmd` records total-read filtering with `$6>=10`, followed by 632 sites, 498 Repeated and 134 Single. Its `uniq -d`/`uniq -u` operations require the original site-grouped input order. These are source-recorded counts, not newly calculated values. The command notebooks also contain intermediate/obsolete steps; they are not uploaded wholesale as executable pipelines.

The source `matched_pLoF_ref_ratio/cmd` calls the Repeated test and records two zero-variance exceptions handled manually. Consequently these scripts alone are not sufficient to recreate all final ASE/NMD classifications. The source notebook labels a mean column as `median`; the actual R code calculates a mean. Neither issue has been silently corrected or used to change published results here.

## Relate / CLUES

- `clues/02.sampleBranchLength.sh`: chromosome and position arguments; invokes the existing Relate branch-length sampling script with mutation rate `1.26e-8`, seed `1` and `200` samples.
- `clues/03.clues_infer.sh`: chromosome and position arguments; invokes CLUES inference with `--burnin 100`, `timeBins.txt`, sampled times and the specified coalescent input.

These two files document the recorded commands only. They do not contain initial genealogy estimation, phasing, allele polarization, neutral simulations, calibration or the final 211-site eligibility definition. Original software/reference paths and the existing genealogy inputs must be configured by the user. Use `bash script.sh ...` after configuring a separate working copy; do not run on archived project outputs.

## Not included yet

Sequence-based NMD positional prediction, the complete annotation pipeline, milk-specific classification, QTL mapping and calibrated neutral selection workflows. Incomplete annotation stubs and an NMD file-discovery utility were inspected but deliberately not presented as production analysis scripts.
