# v10 code release

This update adds the five v9 reviewer-response R scripts and their ten aggregate result CSVs to the original core analysis deposit. The original main-equation coefficients, cut-offs, sample and result tables remain unchanged.

The comparator script was rerun for all 1,000 paired archived resampling replicates on Windows. Point estimates, coefficient exports, all replicate values and interval summaries agreed with the supplied v9 results within 1e-9. Subgroup calibration, near-cut-off positivity, self-report summaries, age/access summaries and positivity among older adults added by R1 were rerun and matched the supplied values. The age-standardisation check reproduced the archived bootstrap standard error. These are computational reproduction checks, not independent external validation.

Windows-compatible process selection and progress reporting are the only changes to the new comparator implementation. The `haven` dependency is now explicit. No software licence or DOI has been assigned in this update.
