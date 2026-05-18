# OutcomeWAS Pipeline

OutcomeWAS estimates `analyte -> outcome` associations within randomized trial
datasets while retaining `TREATMENT_GROUP` as an adjustment covariate. It is
the companion to TreatmentWAS, which handled `treatment -> analyte`.

The current implementation supports:

- continuous outcomes via `OUTCOME`
- time-to-event outcomes via `OUTCOME_TIME` and `OUTCOME_STATUS`
- analyte `change` and analyte `level` models
- pooled, male-only, and female-only analyses
- Proteomics, Metabolomics, and DNAm inputs
- multiple follow-ups, handled as separate models for each nonzero `FU`
- DNAm-specific filtered BH correction
- parallel execution and optional checkpointing

The main technical reference is [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md).

## Public API

Load the package with:

```r
source("main.R")
```

Main functions:

- `FAST_outcome_WAS()`
- `FAST_outcome_WAS_reports()`

## Quick Example

```r
source("main.R")

results <- FAST_outcome_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("agebl", "mbmi"),
  n_cores = 4
)

reports <- FAST_outcome_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("agebl", "mbmi")
)
```

## Required Outcome Schema

`pheno` must contain exactly one of:

- continuous outcome: `OUTCOME`
- time-to-event outcome: `OUTCOME_TIME` and `OUTCOME_STATUS`

Outcome values must be constant within `SUBJECT_ID`.

## Modeling Summary

For each nonzero follow-up `FU = k`, OutcomeWAS fits a separate set of models
using subjects who have both baseline and that follow-up. It does not use mixed
models.

Continuous outcome:

```r
OUTCOME ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

Time-to-event outcome:

```r
survival::Surv(OUTCOME_TIME, OUTCOME_STATUS) ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

`analysis_change` uses `analyte = FU_value - baseline_value`.

`analysis_level` uses `analyte = FU_value`.

## Return Shape

`FAST_outcome_WAS()` returns:

```r
list(
  analysis_change = list(all = ..., male = ..., female = ...),
  analysis_level  = list(all = ..., male = ..., female = ...)
)
```

Each stratum contains:

- `coefficients`: all model terms, with `FU`
- `outcome_effects`: analyte-only summary table, with `FU`

For TTE runs:

- `EFFECT_SIZE` is log hazard ratio
- `HAZARD_RATIO = exp(EFFECT_SIZE)`

For DNAm runs:

- both output tables include `BH_P_VALUE_FILTERED`

## Reports

`FAST_outcome_WAS_reports()` returns:

```r
list(
  pheno_summary = ...,
  variable_summaries = ...,
  outcome_reports = ...
)
```

`outcome_reports$analysis_sample_summary` is one row per nonzero `FU`, so the
report object matches the per-follow-up analysis design.

## Testing

Run the regression test suite with:

```r
Rscript test_comprehensive.R
```

Current coverage includes:

- single-FU continuous and TTE
- multi-FU continuous and TTE
- non-DNAm and DNAm paths

## Notes

- `FU` must be encoded as consecutive integers starting at `0`.
- Additional covariates are taken from the FU analysis row.
- DNAm inputs should be M-values.
- Proteomics and Metabolomics inputs should be log2-transformed before running.
