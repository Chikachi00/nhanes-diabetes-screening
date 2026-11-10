# Hypertension ascertainment and Figure 5

These post hoc additions accompany manuscript v11. They do not replace the main equation, its self-report hypertension input, or the development-derived R0/R1 thresholds. They do not add independent external validation.

Restore the public code (Data Sheet 2) and companion inputs (Data Sheet 3) to one directory. The original raw manifest downloader supplies the BPQ files. Six additional BPX files and their current archived SHA256 hashes are in `reproducibility/additional_raw_manifest_v11.csv`. The companion archive also supplies these six public raw files for offline reproduction.

From the restored root, with R packages `haven` and the original project dependencies:

```sh
python reproducibility/download_nhanes.py raw
python reproducibility/download_bpx_v11.py raw
Rscript analysis/reviewer_response_v11/hypertension_and_budget.R . raw v11_reproduced 1000
python analysis/reviewer_response_v11/figure5.py . v11_reproduced v11_figures
```

The R script reads `prepared.rds`, the original archived coefficients/thresholds, original evaluation predictions and the original Rao-Wu multiplier columns. It reuses their matched replicate order. A refitted main-equation prediction is reconstructed for every replicate; no model is reselected. The script computes 1,000 paired evaluation comparisons, common-sample and strict-BP-sample metrics, descriptive hypertension measures and a budget curve. It writes its analytic rules before numerical estimates. Runtime is several minutes on Windows.

Mean BP uses all valid available readings (1-4); a single available reading is retained. Medication skips are negative only after an explicit upstream no. An affirmative composite component is sufficient; a negative composite requires all components to be known negative. These are single-visit survey measurements, not confirmed clinical diagnoses. The activity mapping remains different from the original Bang publication.

The primary common evaluation sample is 11,624 participants with 348 HbA1c positives. The strict BP sample has 11,548 participants with 344 positives. Within the common sample the fixed-rule R0 sensitivity gain is 10.48 pp (normal 95% 2.57 to 18.39) against history-only Bang and 7.66 pp (-0.34 to 15.66) against composite Bang. The matched-budget gain against composite Bang has a normal interval of 0.04 to 15.28 and a percentile interval of -0.23 to 14.79. These differing bounds are both retained; the manuscript does not claim stable evidence of an advantage across interval methods.

Figure 5 uses the original main sample of 11,848, not the smaller hypertension sensitivity sample. Panel A shows post hoc budgets and pointwise percentile bands; Panel B uses percentile intervals for nested R1-R0 increments; Panel C uses paired normal intervals for changes in the high-minus-low access sensitivity gap. Earlier evaluation-informed selection is not covered by these intervals.

`results/` contains aggregate exports and replicate metrics only. The script additionally emits `participant_mapping_PRIVATE.csv` for local audit. That identifier-level file belongs in the companion data archive and must not be copied into the public code deposit.

All 1,000 replicates completed. An independent Python reconstruction read the raw SAS files, rebuilt Bang scores and checked point values plus five selected replicate columns (1,488 comparisons); maximum error was 6.83e-13. Recalculated summary intervals agreed within 9.67e-13. This is computational verification, not external validation.
