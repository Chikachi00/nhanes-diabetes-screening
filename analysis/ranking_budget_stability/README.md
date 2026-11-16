# v11 merged ranking, budget and cycle checks

Added locally on 11 November 2026. Included in the GitHub repository and in Data Sheet 2. The original six-factor equation, self-report hypertension input and R0/R1 development cut-offs are unchanged.

These are post hoc additions following the earlier v11 findings. They reuse 1,000 paired archived development parameter rows and evaluation Rao–Wu multipliers; no new model is fitted. They do not adjust earlier model/strategy selection or constitute independent external validation.

Restore the merged Data Sheet 2 and the existing Data Sheet 3 into one root, preserving their paths. The latter already contains `analysis/hypertension_budget/participant_mapping.csv`. Prepare a local reference directory and run:

```text
python analysis/ranking_budget_stability/prepare_reference.py . v11_reference
Rscript analysis/ranking_budget_stability/supplementary_check.R . v11_reference reproduced_supplement 1000
python analysis/ranking_budget_stability/verify_independently.py . v11_reference reproduced_supplement independent_verification.json
python analysis/ranking_budget_stability/make_figures.py reproduced_supplement reproduced_figures
```

The reference directory contains participant-level mapping and must not be uploaded as public code. The provided `results/` contains only aggregates. Figure A corresponds to manuscript Figure S6, and Figure B to Figure S7. Generated PNG, PDF and SVG are the figure sources; manuscript TIFFs are LZW-compressed versions of the 450 dpi PNGs.

The main common sample is 11,624 people with 348 HbA1c positives; the strict BP sample is 11,548/344. Equation minus composite Bang AUC is 0.052318 (normal interval 0.031845 to 0.072791; percentile 0.031740 to 0.072053). The AUC differences retained positive directions across the three cycles and each omitted-cycle domain. Operational gains vary more.

All integer Bang cut-offs are retained. Matching includes zero referral, preserves full probability ties and favours lower referral at equal distance. At the composite ≥5 budget, sensitivity gain is 7.659471 pp; the percentile and simultaneous intervals include zero. At ≥6, the budget is approximately 19%, the sensitivity gain is 17.252002 pp (simultaneous 3.488117 to 31.015888), and extra positives are 3.258423 per 1,000 (simultaneous 0.621582 to 5.895264). This does not justify replacing the main cut-off or claiming gains throughout the budget range.

The 96-contrast simultaneous family covers both samples, both Bang hypertension definitions, all 12 integer cut-offs, and sensitivity/positives per 1,000. Cycle and equation-grid intervals remain pointwise. The files report normal and percentile intervals and valid replicate counts. Twelve PPV columns at no referral are undefined, not zero; all other 1,371 columns are complete in every replicate.

Independent Python reconstruction checked 8,226 values and 72 expected missing values at the point estimate and five fixed replicate indices. Largest numerical errors were below 1e-9. Original v11 fixed-rule values were also reproduced. These are computational checks, not clinical validation.
