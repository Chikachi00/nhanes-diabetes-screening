# Reviewer-response analyses (manuscript v9)

Scripts added for the v9 revision. They reuse the frozen v8 inputs and do not change the equation, cut-offs or archived estimates. Run each from the repository root after restoring the companion participant-level inputs to `reproducibility/main_inputs/` (`prepared.rds`, `full_base.csv`) and, where noted, downloading the public XPT files to `./raw` with `python reproducibility/download_nhanes.py raw`.

| Script | Command | Output in `results/` | Manuscript |
|---|---|---|---|
| `comparator_decomp.R` | `Rscript analysis/reviewer_response_v9/comparator_decomp.R . out_s21 1000` (needs `raw`) | `S21_comparator_point.csv`, `S21_comparator_replicates.csv`, `S21_comparator_summary.csv`, `S21_same_info_coefficients.csv` | Section 3.4, Supplementary Table S21 |
| `calib_access.R` | `Rscript analysis/reviewer_response_v9/calib_access.R . out_cd` (needs `raw`) | `C1_subgroup_OE_Taylor.csv`, `C2_near_cutoff_calibration.csv`, `D_selfreport_by_access.csv`, `D2_items_among_positives.csv` | Sections 3.7, 4.4, Supplementary Tables S22, S23 |
| `htn_age.R` | `Rscript analysis/reviewer_response_v9/htn_age.R . out_d3` | `D3_htn_famhx_by_age_access.csv` | Section 4.4, Supplementary Table S23 |
| `r1_added.R` | `Rscript analysis/reviewer_response_v9/r1_added.R . out_f1` | `F1_R1_added_positivity.csv` | Section 3.7 (older adults referred by R1 but not R0) |
| `check_std.R` | `Rscript analysis/reviewer_response_v9/check_std.R .` | printed check | Supplementary Material S15 |

R0 results under FPG, either-test and any-of-three references (Section 3.8, Table 4B) are read from the existing `supplement/tables/S8_reference_performance.csv` (strategy `S0` / `S0_minus_Bang`).

For the v10 preparation, all five scripts were run again from this layout. The comparator was rerun for all 1,000 archived replicates: point estimates, coefficients, replicate metrics and summary intervals matched the supplied v9 exports within 1e-9. The other exported tables matched the supplied values within 1e-10, and the age-standardisation check reproduced the archived bootstrap standard error. These are computational reproduction checks, not independent external validation.
