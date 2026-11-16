# v11 release

Adds post hoc hypertension-ascertainment sensitivity analyses and Figure 5. The main equation and thresholds are unchanged. All 1,000 original matched Rao-Wu draws and saved refitted coefficients were used. R 4.6.1 completed the analysis; Python independently read the raw SAS files and reconstructed 1,488 point and selected-replicate metrics, with maximum error below 7e-13. Summary intervals replayed within 1e-12.

The comparator closer to original-study hypertension ascertainment attenuated the fixed-rule R0 sensitivity gain. Main and supplementary prose now state this dependency. No external validation cohort was added. The aggregate results include both normal and percentile intervals. Participant mappings and raw BPX files remain in the companion archive, not this public code deposit.

Six source URLs and SHA256 values are in the additional v11 raw manifest. No software licence or DOI has been assigned. The original archived scientific outputs remain present without modification.

## Repository update, 16 November 2026

Renamed analysis folders to neutral names (`comparator_calibration`, `hypertension_budget`), added the ranking, integer-budget and cycle-stability analyses in `analysis/ranking_budget_stability/` (Supplements S26–S27; Figures S6–S7), and added `tools/render_production_figures.py`. Scientific results are unchanged.
