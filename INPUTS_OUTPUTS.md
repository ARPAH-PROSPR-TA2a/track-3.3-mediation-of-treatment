# FAST OutcomeWAS: Inputs and Outputs

This document is the compact schema reference. For the implementation details
and control flow, see [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md).

## Inputs

### `pheno`

A data frame with one row per sample.

Required columns for all runs:

| Column | Type | Constraints |
|:---|:---:|:---|
| `SAMPLE_ID` | character | Must be unique |
| `SUBJECT_ID` | character | Repeated across visits for the same person |
| `FU` | factor or integer-coded numeric | Must be consecutive integers starting at `0`; must include `0` and `1` |
| `TREATMENT_GROUP` | factor or integer-coded numeric | Binary `0/1`; both arms must be present |
| `FEMALE` | factor or integer-coded numeric | Binary `0/1` |

Outcome columns must follow exactly one schema:

| Outcome Type | Required Columns | Constraints |
|:---|:---|:---|
| Continuous | `OUTCOME` | Numeric; constant within `SUBJECT_ID` |
| Time-to-event | `OUTCOME_TIME`, `OUTCOME_STATUS` | `OUTCOME_TIME` numeric and `>= 0`; `OUTCOME_STATUS` binary `0/1`; both constant within `SUBJECT_ID` |

Additional covariates named in `additional_covariates` must:

- exist in `pheno`
- be numeric, factor, or logical

Validation behavior that affects row retention:

- duplicated `SUBJECT_ID/FU` rows are reduced to the first occurrence with a
  warning
- samples missing additional covariates are dropped
- samples with missing or incomplete outcome data are dropped
- any subject missing baseline or missing all follow-up rows is dropped

Important follow-up rule:

- multi-follow-up data is supported
- OutcomeWAS fits separate models for each nonzero `FU`
- subjects are eligible for FU `k` if they have both baseline and FU `k`

### `omics`

A data frame with one row per analyte and one column per sample.

| Column | Type | Notes |
|:---|:---:|:---|
| `ANALYTE_NAME` | character | Must be present |
| Sample columns | numeric | Names must match `pheno$SAMPLE_ID` values |

The omics table is intersected with phenotype sample IDs. Extra omics columns
and extra phenotype samples are allowed but excluded from analysis.

### `omics_type`

Must be one of:

- `"Proteomics"`
- `"Metabolomics"`
- `"DNAm"`

DNAm runs trigger:

- probe-manifest coverage checks
- subsetting to the full probe list
- filtered-probe BH correction

### `additional_covariates`

Optional character vector of phenotype column names to include in every model.

These covariates are taken from the FU analysis row, not explicitly frozen at
baseline.

### `n_cores`

`FAST_outcome_WAS()` only. Number of workers used for analyte-level parallel
model fitting.

### `checkpoint_dir`

`FAST_outcome_WAS()` only. Optional directory for checkpoint files.

Checkpoint batches are written under:

```text
{checkpoint_dir}/{response_type}/{stratum}/FU{k}/batch_{b}.rds
```

### `checkpoint_batch_size`

`FAST_outcome_WAS()` only. Number of analytes per checkpoint batch.

## Model Construction

OutcomeWAS builds one model per:

- response type: `change` or `level`
- stratum: `all`, `male`, `female`
- nonzero follow-up `FU`
- analyte

### Subject inclusion for FU `k`

A subject is included in the FU `k` analysis if the subject has:

- one baseline sample at `FU == 0`
- one sample at `FU == k`

The subject does not need every intermediate FU.

### Analysis variables

For each eligible subject and analyte:

- `analyte_baseline`: analyte at baseline
- `analyte_level`: analyte at FU `k`
- `analyte_change`: `analyte_level - analyte_baseline`
- outcome variables from the phenotype table
- `TREATMENT_GROUP`, `FEMALE`, and additional covariates from the FU `k` row

### `analysis_change`

Continuous:

```r
OUTCOME ~ analyte_change + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

TTE:

```r
survival::Surv(OUTCOME_TIME, OUTCOME_STATUS) ~ analyte_change + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

### `analysis_level`

Continuous:

```r
OUTCOME ~ analyte_level + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

TTE:

```r
survival::Surv(OUTCOME_TIME, OUTCOME_STATUS) ~ analyte_level + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

### Notes on covariates

- sex-stratified analyses usually drop `FEMALE` because it becomes constant
- any other adjustment covariate with one observed value in the current
  stratum/FU subset is dropped automatically

## Outputs

### `FAST_outcome_WAS()`

Returns:

```r
list(
  analysis_change = list(all = ..., male = ..., female = ...),
  analysis_level  = list(all = ..., male = ..., female = ...)
)
```

Each stratum contains:

- `$coefficients`
- `$outcome_effects`

If a sex-specific subset is unavailable, that stratum is `NULL`.

If a TTE stratum/FU has zero events, that FU-specific analysis is skipped for
that stratum and may lead to `NULL` results there.

#### `$coefficients`

One row per analyte × follow-up × model term.

Common columns:

| Column | Meaning |
|:---|:---|
| `ANALYTE_NAME` | Analyte identifier |
| `FU` | Follow-up level for this model |
| `COEFFICIENT` | Model term name |
| `N_OBS` | Number of observations used in the fit |
| `EFFECT_SIZE` | Coefficient estimate |
| `SE` | Standard error |
| `P_VALUE` | Raw p-value |
| `BH_P_VALUE` | BH-adjusted p-value within `FU x COEFFICIENT` |

TTE-only additional columns:

| Column | Meaning |
|:---|:---|
| `N_EVENTS` | Number of events in the fit |
| `HAZARD_RATIO` | `exp(EFFECT_SIZE)` |

DNAm-only additional column:

| Column | Meaning |
|:---|:---|
| `BH_P_VALUE_FILTERED` | BH correction within the filtered probe set, grouped by `FU x COEFFICIENT` |

Interpretation of `EFFECT_SIZE`:

- continuous outcomes: linear-model coefficient
- TTE outcomes: log hazard ratio

#### `$outcome_effects`

One row per analyte × follow-up, keeping only the analyte term itself.

Common columns:

| Column | Meaning |
|:---|:---|
| `ANALYTE_NAME` | Analyte identifier |
| `FU` | Follow-up level for this model |
| `EFFECT_SIZE` | Analyte coefficient estimate |
| `SE` | Standard error |
| `P_VALUE` | Raw p-value |
| `BH_P_VALUE` | BH-adjusted p-value within `FU` |

TTE-only additional column:

| Column | Meaning |
|:---|:---|
| `HAZARD_RATIO` | `exp(EFFECT_SIZE)` |

DNAm-only additional column:

| Column | Meaning |
|:---|:---|
| `BH_P_VALUE_FILTERED` | BH correction within the filtered probe set, grouped by `FU` |

`outcome_effects` is a strict subset summary of the `COEFFICIENT == "analyte"`
rows from `$coefficients`.

### `FAST_outcome_WAS_reports()`

Returns:

```r
list(
  pheno_summary = ...,
  variable_summaries = ...,
  outcome_reports = ...
)
```

#### `pheno_summary`

Counts by `FU` and `FEMALE`, including:

- `N_SUBJECTS`
- `N_CONTROL`
- `N_TREATMENT`
- `N_SAMPLES`

#### `variable_summaries`

Sex-stratified summaries with entries keyed by `FU x TREATMENT_GROUP`, for
example:

- `omics_FU0_Tx0`
- `omics_FU2_Tx1`
- `covariates_FU1_Tx0`

Each omics entry is one row per analyte with:

- `N_NONMISSING`
- `MEAN`
- `MEDIAN`
- `SD`
- `MIN`
- `MAX`

Each covariate entry is one row per covariate with:

- `COVARIATE_NAME`
- `TYPE`
- `N_NA`
- `SUMMARY`

#### `outcome_reports`

Contains:

- `outcome_type`
- `analysis_sample_summary`
- `outcome_summary`

##### `analysis_sample_summary`

One row per nonzero FU:

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up level |
| `N_SUBJECTS` | Subjects with both baseline and this FU |
| `N_BASELINE_SAMPLES` | Baseline samples contributing to this FU analysis |
| `N_FOLLOWUP_SAMPLES` | Follow-up samples contributing to this FU analysis |

##### `outcome_summary`

This summary is computed on a deduplicated one-row-per-subject dataset so the
subject-level outcome is not double-counted across repeated FU rows.

Continuous outcome columns:

- `GROUP`
- `N_SUBJECTS`
- `MEAN`
- `MEDIAN`
- `SD`
- `MIN`
- `MAX`

TTE outcome columns:

- `GROUP`
- `N_SUBJECTS`
- `N_EVENTS`
- `N_CENSORED`
- `EVENT_RATE`
- `TIME_MEAN`
- `TIME_MEDIAN`
- `TIME_MIN`
- `TIME_MAX`

## Validation and Failure Behavior

- invalid required columns or invalid field types stop the run
- non-consecutive `FU` encoding stops the run
- non-constant within-subject outcome, sex, or treatment stops the run
- analytes with fit failures are returned as `NULL` internally and omitted from
  the final bound tables, with warnings
- a TTE stratum/FU with zero events returns `NULL` for that FU-specific
  analysis with a warning

## Testing

The main regression test is `test_comprehensive.R`. It currently exercises:

- single-FU continuous + non-DNAm
- single-FU TTE + non-DNAm
- multi-FU continuous + non-DNAm
- multi-FU TTE + non-DNAm
- multi-FU continuous + DNAm
- multi-FU TTE + DNAm
