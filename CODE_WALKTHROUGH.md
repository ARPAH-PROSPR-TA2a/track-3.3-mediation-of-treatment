# OutcomeWAS Code Walkthrough

OutcomeWAS estimates `analyte -> outcome` associations inside a trial, while
retaining `TREATMENT_GROUP` as an adjustment covariate. The current
implementation supports:

- continuous outcomes via `OUTCOME`
- time-to-event outcomes via `OUTCOME_TIME` and `OUTCOME_STATUS`
- analyte `change` models and analyte `level` models
- pooled, male-only, and female-only analyses
- Proteomics, Metabolomics, and DNAm inputs
- multiple follow-up visits, handled as separate models per follow-up

Unlike TreatmentWAS, OutcomeWAS does **not** use mixed models. If a
study has `FU = 0, 1, 2`, the pipeline fits one set of models comparing
baseline to FU1 and a second set comparing baseline to FU2.

## File Map

- `main.R`: public API
- `validation_helpers.R`: input validation and phenotype/omics harmonization
- `analysis_helpers.R`: per-analyte model fitting, FU looping, BH correction,
  DNAm filtered correction, checkpointing
- `reporting_helpers.R`: QC summaries and outcome summaries
- `test_comprehensive.R`: regression test suite covering single-FU and
  multi-FU, continuous and TTE, non-DNAm and DNAm

## Public API

The package exposes two entry points.

### `FAST_outcome_WAS()`

This is the analysis function. It:

1. validates `omics_type`
2. validates and subsets `pheno`
3. validates and harmonizes `omics`
4. loads DNAm probe manifests when `omics_type == "DNAm"`
5. runs stratified `change` analysis
6. runs stratified `level` analysis
7. returns both result trees

### `FAST_outcome_WAS_reports()`

This is the reporting function. It runs the same validation and harmonization
steps, then produces QC summaries instead of model fits.

## Phenotype Contract

`pheno` is one row per sample. Required columns are:

- `SAMPLE_ID`
- `SUBJECT_ID`
- `FU`
- `TREATMENT_GROUP`
- `FEMALE`
- one outcome schema:
  - continuous: `OUTCOME`
  - TTE: `OUTCOME_TIME`, `OUTCOME_STATUS`

Additional covariates are optional, but if they are named in
`additional_covariates` they must exist in `pheno`.

### Follow-up Encoding

`FU` must satisfy all of the following:

- integer-valued after coercion
- non-negative
- contain baseline `0`
- contain at least one nonzero follow-up, and specifically `1`
- be consecutive integers from `0` through `max(FU)`

That means raw visit encodings like `0, 3, 6, 12` must be recoded to
`0, 1, 2, 3` before running the package.

The requirement that `FU == 1` exists is intentional. The first follow-up is
treated as the minimum valid longitudinal structure for OutcomeWAS.

## Outcome-Type Detection

Outcome type is inferred automatically in `.detect_outcome_type()`:

- `OUTCOME` only: continuous
- `OUTCOME_TIME` + `OUTCOME_STATUS`: TTE
- both schemas present: error
- incomplete TTE schema: error
- neither schema present: error

There is no explicit `outcome_type` argument in the public API.

## Validation Flow

`FAST_outcome_WAS()` and `FAST_outcome_WAS_reports()` both call the same
validation stack.

### 1. `omics_type` validation

Accepted values are:

- `DNAm`
- `Proteomics`
- `Metabolomics`

This also prints a reminder message about expected input scale.

### 2. Phenotype validation

`.validate_pheno()` performs the following checks and transformations.

#### Core column checks

- confirms `pheno` is a `data.frame` or matrix
- checks all required columns exist
- checks `additional_covariates` is `NULL` or character

#### Binary field checks

- `FEMALE` must contain only `0/1`
- `TREATMENT_GROUP` must contain only `0/1`
- both treatment arms must be present somewhere in the dataset

Both fields are converted to factors if they are not already factors.

#### Sample uniqueness and replicate handling

- `SAMPLE_ID` must be unique globally
- duplicated `SUBJECT_ID/FU` pairs are not fatal
- if duplicated `SUBJECT_ID/FU` pairs are found, only the first occurrence is
  kept and a warning is emitted

This means technical replicates at the same visit are silently reduced to the
first row after warning. That behavior is inherited from the scaffold and is
worth keeping in mind.

#### Outcome validation

For continuous outcomes:

- `OUTCOME` must be numeric

For TTE outcomes:

- `OUTCOME_TIME` must be numeric
- `OUTCOME_TIME >= 0`
- `OUTCOME_STATUS` must be binary `0/1`
- `OUTCOME_STATUS` is coerced to integer

#### Missing covariates and missing outcomes

For each additional covariate:

- rows with missing values are dropped at the sample level
- a message is emitted with the number of dropped samples

For outcome fields:

- rows with missing or incomplete outcome data are dropped
- the message reports the number of dropped samples

The intent is to be strict and visible rather than allowing partner sites to
run partially missing outcome data without noticing.

#### Subject-level constancy

The following columns must be constant within `SUBJECT_ID`:

- `TREATMENT_GROUP`
- `FEMALE`
- `OUTCOME`, or `OUTCOME_TIME` and `OUTCOME_STATUS`

If not, validation stops.

This reflects the OutcomeWAS assumption that the outcome is subject-level and
shared across the person’s analyte rows.

#### Baseline/follow-up completeness

Subjects are retained if and only if they have:

- at least one baseline row (`FU == 0`)
- at least one nonzero follow-up row

Subjects missing either side are dropped with a message.

This is a dataset-level screen. Later, during analysis, each FU-specific model
uses the subset of subjects who have both baseline and that specific FU.

#### Sex-specific subsets

The validator returns:

- `all`
- `male`
- `female`

If the dataset contains only one sex, male/female subsets are returned as
`NULL` and a warning is emitted.

### 3. Omics validation

`.validate_omics()` expects:

- an `ANALYTE_NAME` column
- all other columns numeric
- column names matching `pheno$SAMPLE_ID`

It then:

- intersects omics columns with phenotype sample IDs
- reports counts of shared / omics-only / pheno-only samples
- warns if any analytes contain missing values
- warns if any analytes have near-zero variance
- returns aligned omics for `all`, `male`, and `female`

No analytes are dropped for NA or low variance at validation time.

### 4. DNAm-specific validation

If `omics_type == "DNAm"`:

- `Data/FAST_epicv1_epicv2_probe_list.rds` is loaded as the full probe list
- `Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds` is loaded as the
  filtered probe list
- coverage against both lists is checked and reported
- the omics table is subset to probes in the full probe list

The filtered probe list is used later only for the additional filtered BH
correction column.

## Analysis Design

### Overall structure

The analysis is run twice:

- `response_type = "change"`
- `response_type = "level"`

For each response type, `.run_stratified_analysis()` runs:

- `all`
- `male`
- `female`

For each stratum, `.perform_analysis()` loops over every nonzero follow-up
level and fits a separate set of analyte-wise models.

### FU-specific subject inclusion

For a given FU level `k`, the model includes a subject if that subject has:

- one baseline row
- one row at FU `k`

The subject does **not** need to have all intermediate follow-ups. For
example, a subject with FU0 and FU2 but no FU1 is eligible for the FU2 model.

### Analysis dataset construction

For a given stratum and FU:

1. `pheno_baseline_all` is defined as all baseline rows in that stratum
2. `pheno_analysis` is defined as rows for the target FU
3. subjects present in both are intersected
4. baseline omics values are matched by `SUBJECT_ID`
5. a model frame is built with:
   - outcome columns from the FU row
   - `analyte_baseline`
   - `analyte`
   - `TREATMENT_GROUP`
   - `FEMALE`
   - additional covariates

The package uses the FU-row covariates, not baseline covariates. The intended
use is that these are effectively constant subject-level variables, but the
implementation takes them from the analysis row.

### `change` versus `level`

For each analyte:

- `change`: `analyte = FU_value - baseline_value`
- `level`: `analyte = FU_value`

In both cases, `analyte_baseline` is included as a covariate.

This is the direct analog of the TreatmentWAS baseline-adjusted formulation.

## Model Formulas

The model builder always starts with:

- analyte term of interest: `analyte`
- baseline adjustment: `analyte_baseline`
- trial adjustment: `TREATMENT_GROUP`
- sex adjustment in pooled analyses: `FEMALE`
- optional user-provided covariates

Before fitting, `.drop_uninformative_covariates()` removes any adjustment term
with only one observed value in the current stratum/FU. This matters most for:

- `FEMALE` in sex-stratified analyses
- covariates that are constant within a sex stratum or FU subset

### Continuous outcome

The fitted model is:

```r
OUTCOME ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

using `lm()`.

### Time-to-event outcome

The fitted model is:

```r
survival::Surv(OUTCOME_TIME, OUTCOME_STATUS) ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

using `survival::coxph(..., ties = "efron")`.

If a stratum/FU has zero events, that FU-specific analysis returns `NULL` for
that stratum with a warning.

At the analyte level, any fit failure also returns `NULL` for that analyte and
emits a warning. Common reasons would be non-estimable coefficients or Cox
instability.

## Parallelism and Checkpointing

Within each FU/stratum, analytes are split into batches. For each batch:

- analyte-wise work is parallelized with `furrr::future_map()`
- `n_cores > 1` uses a `future::multisession` plan
- `n_cores == 1` uses `future::sequential`

Checkpointing is optional. If `checkpoint_dir` is provided, batch files are
written under:

```text
{checkpoint_dir}/{response_type}/{stratum}/FU{k}/batch_{b}.rds
```

This makes checkpoint reuse FU-specific.

## Output Objects

`FAST_outcome_WAS()` returns:

```r
list(
  analysis_change = list(all = ..., male = ..., female = ...),
  analysis_level  = list(all = ..., male = ..., female = ...)
)
```

Each non-`NULL` stratum contains:

- `coefficients`
- `outcome_effects`

### `coefficients`

This is one row per analyte, per FU, per model term. Core columns are:

- `ANALYTE_NAME`
- `FU`
- `COEFFICIENT`
- `N_OBS`
- `EFFECT_SIZE`
- `SE`
- `P_VALUE`
- `BH_P_VALUE`

TTE runs also include:

- `N_EVENTS`
- `HAZARD_RATIO`

DNAm runs also include:

- `BH_P_VALUE_FILTERED`

`EFFECT_SIZE` is:

- the linear-model coefficient for continuous outcomes
- the log-hazard ratio for TTE outcomes

### `outcome_effects`

This is one row per analyte per FU, keeping only the analyte term from the full
coefficient table. Core columns are:

- `ANALYTE_NAME`
- `FU`
- `EFFECT_SIZE`
- `SE`
- `P_VALUE`
- `BH_P_VALUE`

TTE runs also include:

- `HAZARD_RATIO`

DNAm runs also include:

- `BH_P_VALUE_FILTERED`

The point of `outcome_effects` is to provide the compact analyte-level summary
table that most downstream users will inspect first.

## Multiple Testing Correction

The grouping for BH correction is part of the statistical definition and is
important.

### `coefficients`

`BH_P_VALUE` is computed separately within each:

- `FU`
- `COEFFICIENT`

That means the analyte term at FU1 is corrected separately from the analyte
term at FU2, and separately from the baseline or treatment adjustment terms.

### `outcome_effects`

`BH_P_VALUE` is computed separately within each:

- `FU`

Since `outcome_effects` contains only the analyte term, no coefficient-name
stratification is needed there.

### DNAm filtered correction

For DNAm only, `BH_P_VALUE_FILTERED` is added:

- on `coefficients`, separately within `FU x COEFFICIENT`
- on `outcome_effects`, separately within `FU`

but only for analytes in the filtered probe list.

## Reporting Objects

`FAST_outcome_WAS_reports()` returns:

```r
list(
  pheno_summary = ...,
  variable_summaries = ...,
  outcome_reports = ...
)
```

### `pheno_summary`

This is a simple count table by `FU` and `FEMALE`, with:

- `N_SUBJECTS`
- `N_CONTROL`
- `N_TREATMENT`
- `N_SAMPLES`

### `variable_summaries`

This is returned separately for `all`, `male`, and `female`.

Within each stratum, summaries are keyed by `FU` and treatment arm:

- `omics_FU0_Tx0`
- `omics_FU2_Tx1`
- `covariates_FU1_Tx0`
- etc.

Each omics entry is one row per analyte with non-missing count and summary
statistics. Each covariate entry is one row per named covariate with a compact
typed summary.

### `outcome_reports`

This contains:

- `outcome_type`
- `analysis_sample_summary`
- `outcome_summary`

#### `analysis_sample_summary`

This is one row per nonzero FU with:

- `FU`
- `N_SUBJECTS`
- `N_BASELINE_SAMPLES`
- `N_FOLLOWUP_SAMPLES`

This table describes the eligible sample size for each FU-specific model set.

#### `outcome_summary`

Because the outcome is subject-level, this summary is built from a deduplicated
one-row-per-subject frame, not from all longitudinal rows. Without that step,
subjects with more follow-up rows would be over-counted.

For continuous outcomes, the table reports summary statistics for:

- `all`
- `control`
- `treatment`

For TTE outcomes, it reports:

- subject counts
- event counts
- censoring counts
- event rate
- summary statistics of follow-up time

## Testing

`test_comprehensive.R` currently covers:

- single-FU continuous + non-DNAm
- single-FU TTE + non-DNAm
- multi-FU continuous + non-DNAm
- multi-FU TTE + non-DNAm
- multi-FU continuous + DNAm
- multi-FU TTE + DNAm

The tests verify:

- top-level result structure
- presence and correctness of `FU` columns
- sex-stratified result trees
- TTE hazard-ratio columns
- DNAm filtered BH columns
- report structure
- exact consistency between `outcome_effects` and the analyte rows of
  `coefficients`

## Current Limitations and Open Points

- Outcome type is inferred only from the canonical column names. There is no
  aliasing layer.
- The pipeline assumes the outcome is constant within subject across visits.
- Covariates are taken from the FU row rather than being explicitly frozen at
  baseline.
- Duplicate `SUBJECT_ID/FU` rows are handled by keeping the first row after a
  warning, rather than requiring the user to resolve them upstream.
- There is no explicit minimum-event threshold for Cox models beyond the
  stratum/FU-level zero-event guard.
- Subjects contribute independently to each FU-specific model if they have the
  required baseline/FU pair.

Those are the main implementation decisions to review if we want to tighten the
spec further.
