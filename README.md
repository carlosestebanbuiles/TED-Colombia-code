# Thyroid eye disease in Colombian referral centers — analysis code

Analysis code for:

> Camargo González J, Aristizábal-Henao N, Pradierna-Espinoza S, et al. *Clinical Presentation and Recorded Follow-up Findings of Thyroid Eye Disease in Colombian Referral Centers: A Multicenter Retrospective Registry Study.* Therapeutic Advances in Ophthalmology (under review, manuscript OED-26-0109).

The code derives every variable from the registry spreadsheet and reproduces all numbers, tables and figures in the article and its supplementary material.

## Data

The patient-level registry data are **not** included: the participating institutions have not authorized public release. A pseudonymized version of the data is available from the corresponding author (Carlos E. Builes-Montaño, esteban.builes@udea.edu.co) on reasonable request, subject to institutional approval.

To make the code executable without the registry data, `data/synthetic_BD.xlsx` is a **synthetic** file with the same structure as the registry (two header rows, the same 155 columns, the same value formats and similar missingness). Each column was resampled independently, identifiers were replaced and free-text fields were replaced by generic phrases, so no synthetic row corresponds to a real patient and **results obtained from the synthetic file are meaningless**. Its only purpose is to let reviewers and readers run, inspect and modify the full pipeline.

## Running

Requirements: R ≥ 4.3 with `tidyverse`, `readxl`, `writexl`, `lme4`, `geepack`, `patchwork`, `scales` (and optionally `logistf`, see below).

```bash
export LC_ALL=C.UTF-8
export TED_BD=data/synthetic_BD.xlsx     # path to any spreadsheet built from the registry template

Rscript 01_derive.R              # import and harmonize; derive the analytic data set
Rscript 02_verify.R              # with TED_VERIFY=1 on the registry data: checks 26 published numbers
Rscript 03_revision_analyses.R   # all analyses -> outputs/revision_results.xlsx, outputs/03_log.txt
Rscript 04_figures.R             # Figures 1–4 and S1 (PNG 300 dpi, TIFF 600 dpi)
Rscript 05_tables_1_2.R          # Tables 1 and 2
```

### Running on other data

The code runs on any spreadsheet built with the same registry template (EUGOGO Database Template, Colombia: two header rows and the variable names of the glossary), whatever its number of patients or centres. Three environment variables adapt the study-specific settings:

| Variable | Default | Meaning |
|---|---|---|
| `TED_BD` | `BD.xlsx` | Path to the registry spreadsheet |
| `TED_BASELINE_CENTRES` | `Centre A` | Comma-separated centres in which the registry entry date is the date of the baseline visit, so that the baseline-to-last-visit interval can be computed |
| `TED_EXCLUDE_CENTRES` | `Centre B,Centre D` | Centres excluded jointly in the sensitivity analysis of incomplete documentation |
| `TED_CENTRE_MAP` | `config/centre_map.csv` | Optional file mapping the centre names in the data to pseudonyms (columns `centre_code,centre`) |
| `TED_VERIFY` | `0` | Set to `1` only with the registry data, to check the published numbers |

Example with another data set:

```bash
export TED_BD=my_registry.xlsx TED_BASELINE_CENTRES="Hospital X" TED_EXCLUDE_CENTRES="Hospital Z"
Rscript 01_derive.R && Rscript 03_revision_analyses.R && Rscript 04_figures.R && Rscript 05_tables_1_2.R
```

The code was tested on the registry data (147 patients, five centres) and on a synthetic data set with 90 patients and three differently named centres. Sections of the analysis that need a minimum number of events (for example the regression models) may return wide or unstable estimates in small data sets.

## Files

| File | Content |
|---|---|
| `00_helpers.R` | Text cleaning, yes/no harmonization, decimal commas and censored values, mixed date formats, formatting helpers |
| `00_firth.R` | Firth-penalized logistic regression with profile penalized-likelihood confidence intervals and likelihood-ratio tests; wraps `logistf::logistf()` when installed, otherwise uses a self-contained implementation that gives identical results |
| `01_derive.R` | Variable derivation: activity (recorded CAS or clinician's field), item-derived CAS, operational EUGOGO severity, treatment documentation windows, follow-up outcomes |
| `02_verify.R` | Reproduction check against the published numbers (registry data only) |
| `03_revision_analyses.R` | Sections A–H: activity measurement, severity denominators, follow-up and attrition, DON, treatment windows, center heterogeneity (fixed effects, random intercept, GEE, leave-one-center-out), paired effect sizes, protocol deviations, sex-specific exophthalmometry thresholds |
| `04_figures.R`, `05_tables_1_2.R` | Figures and Tables 1–2 |
| `config/centre_map.csv` | **Not included.** Maps institution names in the registry to the pseudonyms Centre A–E used in the article. Without it, centre names are used as they appear in the data (the synthetic file already uses pseudonyms) |

## Notes

- Quantiles: Tables 1 and 2 use type-2 quantiles (as produced by `gtsummary` in the original analysis); all other summaries use R's default (type 7). The last digit of a few Table 1 values can differ from the article by half-unit rounding.
- GEE models are fitted but not reported because, with five clusters, the sandwich variance estimator is degenerate.

## License

MIT (code). The synthetic data file is released under CC0.
