# NHANES diabetes screening analysis

This repository accompanies the authors' v11 preparation for Frontiers in Public Health. The main equation and cut-offs remain frozen from v8; the reviewer-response analyses were added in v9. The v10 update clarifies reporting and enables the new comparator script to run on Windows. It contains the frozen main equation, model development and evaluation code, supplementary analyses, aggregate results, resampling metrics, and figure source data. The manuscript compares referral rules under operational constraints. Its score cutoffs are operational cutpoints, not validated absolute-risk thresholds.

## What runs directly

From the extracted repository root, with Python 3 and the packages in `requirements.txt`:

```sh
python -m pip install -r requirements.txt
python tools/check_release.py
python reproducibility/predict_locked.py examples/people_synthetic.csv predictions.csv
python tools/render_figures.py reproduced_figures
```

The prediction example contains synthetic people only. Required inputs are age, waist, height, male, family history, hypertension and inactivity. Waist and height use the same units; binary indicators must be observed 0/1. Unknown values must not be coded as zero. The output is research arithmetic, not a clinical diagnosis or a calibrated individual risk estimate.

`tools/check_release.py` verifies source hashes, Python syntax, frozen parameters, synthetic predictions and all normal/bootstrap-SE and percentile intervals against the archived 1,000-replicate metric files. It does not refit R models.

## Main analysis and required data

```sh
python reproducibility/download_nhanes.py raw
Rscript reproducibility/reproduce_main.R reproduced_main
# Optional complete refit of the 1,000 archived development/evaluation replicates:
Rscript reproducibility/reproduce_main.R reproduced_main --refit-all
```

The download command retrieves 63 public CDC NHANES XPT files and verifies the frozen hashes. Downloading those files alone does not build the frozen main-analysis cache. `reproduce_main.R` additionally needs `reproducibility/main_inputs/prepared.rds`, `main_inputs/full_base.csv`, `main_archive/evaluation_predictions.csv`, and `main_archive/development_oof.csv`. These participant-level inputs are supplied in the companion scientific reproduction archive, not this GitHub-ready package. Restore the same relative paths from that archive before running the main command. The matching development/evaluation CSV exports and public raw data are also in the separate dataset package. The original preparation/development implementation is preserved at `reproducibility/main_inputs/source_snapshot/refine_v5_equation.R`; its historical filesystem defaults need adaptation before standalone execution.

`reproduce_main.R` uses base R. Its default mode refits the development point model once and replays archived metric intervals. `--refit-all` performs the more expensive full replicate refits. Earlier manual candidate selection remains outside those intervals.

## Supplementary entry points

Run from the repository root; replace `.` in the fourth R argument below with a custom R package library only if needed.

```sh
Rscript analysis/rerun_s17.R raw reproduced_s17 supplement/code/external_validation_v5_WTPH2YR.R . supplement/tables
Rscript analysis/verify_s17_exports.R reproduced_s17
Rscript analysis/verify_s17_predictions.R raw reproduced_s17 analysis/predict_s17_derived.R
Rscript supplementary_500/joint_cutoff_v5_audit.R raw reproduced_joint
Rscript supplementary_500/access_gap_v5_audit.R raw reproduced_access
Rscript fext/fext_unknown_sensitivity.R . raw reproduced_fext supplement/code/denoise_features_v5.R supplement/tables
Rscript fext/validate_saved_fext.R . reproduced_fext supplement/tables raw
Rscript fext/replay_denoise_500.R . raw reproduced_denoise supplement/code/denoise_features_v5.R supplement/tables
Rscript analysis/local_calibration/local_calibration_v5.R . reproduced_local_calibration
python analysis/discrimination_v7/auc_v7.py --root . --evaluation PATH_TO_DATASET_PACKAGE/evaluation_2013_2018.csv
```

Local calibration requires the companion frozen main inputs. The AUC replay requires the participant-level evaluation CSV; it uses the included paired PSU multipliers and full-pipeline coefficient archives. Some fext replay/validation commands also require saved R model/replicate objects from the companion archive, restored to their original paths. `R_packages.txt` lists packages referenced by the archived R scripts; install them in R before rerunning supplementary analyses. No package version lock is claimed. During packaging, the main point model was refitted once in R 4.6.1 after restoring the companion inputs to a separate copy. The point metrics, predictions and all archived intervals passed verification. All 26 R sources parsed, and the portable five-figure renderer ran successfully. Full replicate refits and all supplementary R analyses were not rerun during packaging.

## File map

| Directory | Contents |
|---|---|
| `reproducibility/` | Download/hash checks, frozen main parameters, main point refit and interval replay |
| `analysis/` | S17 later-period implementation, saved-model checks, local calibration, S19 paired AUC replay |
| `fext/` | Family-history missing/unknown sensitivity and saved denoising replay |
| `supplementary_500/` | Audited 500-replicate access and joint-cutoff analyses |
| `supplement/code/` | Original supplementary model, missing-data, recalibration and evaluation implementations |
| `supplement/tables/` | Frozen aggregate tables for the manuscript supplements |
| `figures/source_data/` | Frozen aggregate inputs for five figures |
| `archive/` | Original figure and cohort-export scripts, with historical path assumptions |
| `tools/` | Portable figure renderer and release verification |
| `analysis/reviewer_response_v9/` | v9 reviewer-response scripts and their aggregate outputs (see its README) |

Legacy `S0/S1/S4` names in the code correspond to manuscript `R0/R1/R4`. Historical filenames ending in v5 or v7 identify code lineage; they do not mean that an older manuscript should replace v8. `SOURCE_PROVENANCE.json` records unchanged source files and deliberate exclusions. All legacy scientific scripts retain their original bytes. The portable figure wrapper changes input/output paths only; v8 submission images remain authoritative because font and raster settings can vary.

The original `auc_v7.py` writes historical provenance prose about retaining the old overall interval. In the v8 manuscript the overall equation AUC uses the full-pipeline interval. Consult the numeric `auc_summary.csv` and v8 text for the current reporting choice.

## Publication and reuse status

Repository: [Chikachi00/nhanes-diabetes-screening](https://github.com/Chikachi00/nhanes-diabetes-screening). This deposit contains core analysis code, frozen aggregate outputs and resampling metrics. Participant-level derived inputs and saved R objects are supplied separately in the companion reproduction archive. The project source code and its associated software documentation are licensed under the MIT License (see LICENSE). Cite the software using CITATION.cff. This licence does not relicense NHANES data, manuscript text, or third-party dependencies; no DOI is asserted here. Author declarations and institutional ethics documentation are not part of this code deposit.

`.gitattributes` preserves the archived file bytes so SHA-256 checks work after cloning on Windows.


## v10 execution note

The new analyses in `analysis/reviewer_response_v9/` additionally require the R package `haven` (now listed in `R_packages.txt`). `comparator_decomp.R` uses one process on Windows and four on platforms supporting forking; the fitting and resampling calculations are unchanged. Progress is printed every 50 replicates. The companion participant-level inputs remain necessary for these scripts. Internal numerical reproduction is distinct from independent external validation. See `RELEASE_NOTES_v10.md` for the recorded checks.


## v11 hypertension ascertainment and Figure 5

`analysis/reviewer_response_v11/` adds six-cycle public BPX ascertainment, common-sample and strict-BP-sample sensitivity comparisons, 1,000 paired replicate exports and the Figure 5 builder. The main self-report-based equation and R0/R1 thresholds remain unchanged. Composite hypertension ascertainment improves Bang performance; the fixed-rule R0 sensitivity gain falls to 7.66 pp (normal 95% -0.34 to 15.66) in the common sample. Both normal and percentile equal-budget intervals are supplied because their lower bounds differ around zero. See the subdirectory README for reproduction commands and source-file provenance. Independent numerical checks are computational verification, not external validation.
