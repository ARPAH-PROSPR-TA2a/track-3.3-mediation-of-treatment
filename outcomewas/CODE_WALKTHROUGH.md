# OutcomeWAS Pipeline: Code Walkthrough

<style>
table {
  border-collapse: collapse;
  width: 100%;
  margin: 1em 0 1.5em 0;
}

th,
td {
  border: 1px solid #d0d7de;
  padding: 0.45em 0.65em;
  vertical-align: top;
}

th {
  background: #f6f8fa;
  font-weight: 600;
}

pre {
  margin: 1em 0 1.5em 0;
}

.note {
  border-left: 4px solid #6f42c1;
  background: #f6f3fb;
  padding: 0.75em 1em;
  margin: 1em 0 1.5em 0;
}
</style>

This walkthrough documents the current behavior of the OutcomeWAS pipeline in
`outcomewas/main.R` and its pipeline-local helper files. The pipeline estimates `analyte -> outcome`
associations inside randomized trial datasets while retaining
`TREATMENT_GROUP` as an adjustment covariate.

OutcomeWAS exposes two public functions:

- `FAST_outcome_WAS()`: runs the statistical analyses with parallelization and
  optional checkpointing.
- `FAST_outcome_WAS_reports()`: generates QC and data summary reports without
  fitting models.

<div class="note">
Unlike TreatmentWAS, OutcomeWAS does not use mixed models. Multi-follow-up data
are handled as separate baseline-to-follow-up analyses for each nonzero `FU`.
</div>

## Table of Contents

1. [File Structure](#file-structure)
2. [Main Functions](#main-functions)
3. [Accepted Inputs](#accepted-inputs)
4. [Validation Flow](#validation-flow)
5. [High-Level Pipeline Flow](#high-level-pipeline-flow)
6. [Analysis Design](#analysis-design)
7. [Model Formulas](#model-formulas)
8. [Parallelization and Checkpointing](#parallelization-and-checkpointing)
9. [Multiple Testing Correction](#multiple-testing-correction)
10. [DNAm Probe Sets](#dnam-probe-sets)
11. [Reporting Pipeline](#reporting-pipeline)
12. [Results Output](#results-output)

---

## File Structure

```text
outcomewas/
  main.R                     Public API: FAST_outcome_WAS(), FAST_outcome_WAS_reports()
  R/
    validation_helpers.R     Input validation and phenotype/omics harmonization
    analysis_helpers.R       Model fitting, FU looping, BH correction, checkpointing
    reporting_helpers.R      QC summaries and outcome reports
    plotting_helpers.R       QQ and volcano plots from outcome_effects
  Data/                      OutcomeWAS DNAm probe lists
  Examples/                  Pipeline-specific example inputs, outputs, and scripts
tests/outcomewas/
  test_comprehensive.R       Main regression suite
  test_parallel_checkpoint.R Parallelization and checkpointing tests
```

Function locations:

| File | Key functions |
|:---|:---|
| `outcomewas/main.R` | `FAST_outcome_WAS()`, `FAST_outcome_WAS_reports()` |
| `outcomewas/R/validation_helpers.R` | `.detect_outcome_type()`, `.validate_omics_type()`, `.validate_pheno()`, `.validate_omics()`, `.validate_dnam_probe_coverage()`, `.subset_omics_list()` |
| `outcomewas/R/analysis_helpers.R` | `.perform_continuous_analysis()`, `.perform_tte_analysis()`, `.perform_analysis()`, `.run_stratified_analysis()`, `.apply_multiple_testing_correction()`, `.add_filtered_bh_correction()` |
| `outcomewas/R/reporting_helpers.R` | `.generate_reports()`, `.create_pheno_data_report()`, `.create_omics_data_report()`, `.create_addx_covariate_report()`, `.create_analysis_sample_summary()`, `.create_continuous_outcome_report()`, `.create_tte_outcome_report()` |
| `outcomewas/R/plotting_helpers.R` | `plot_qq()`, `plot_volcano()`, `generate_all_plots()` |

---

## Main Functions

### `FAST_outcome_WAS()`

File: `outcomewas/main.R`

```r
FAST_outcome_WAS <- function(pheno,
                             omics,
                             omics_type = "Proteomics",
                             additional_covariates = NULL,
                             n_cores = NULL,
                             checkpoint_dir = NULL,
                             checkpoint_batch_size = 2000L)
```

What it does:

1. Resolves `n_cores`; defaults to `max(1, parallel::detectCores() - 1)`.
2. Sets a `future` plan for analyte-level parallelization and restores the
   previous plan on exit.
3. Validates `omics_type`.
4. Validates `pheno`, detects outcome type, and creates sex-specific subsets.
5. Validates and aligns `omics` to the retained phenotype samples.
6. For DNAm, loads probe lists, checks coverage, and subsets to the full probe
   list.
7. Runs stratified `change` analysis.
8. Runs stratified `level` analysis.
9. Returns both result trees.

Return shape:

```r
list(
  analysis_change = list(all = ..., male = ..., female = ...),
  analysis_level  = list(all = ..., male = ..., female = ...)
)
```

### `FAST_outcome_WAS_reports()`

File: `outcomewas/main.R`

```r
FAST_outcome_WAS_reports <- function(pheno,
                                     omics,
                                     omics_type = "Proteomics",
                                     additional_covariates = NULL)
```

What it does:

1. Runs the same validation and harmonization stack as `FAST_outcome_WAS()`.
2. Applies the same DNAm probe-list handling when `omics_type == "DNAm"`.
3. Calls `.generate_reports()`.
4. Returns QC summaries and outcome summaries.

Reports and analysis are intentionally separate. Reports are cheap and
descriptive; model fitting can be long-running, parallelized, and checkpointed.

---

## Accepted Inputs

The full user-facing input/output contract is in `INPUTS_OUTPUTS.md`. This
section records the implementation contract.

### Phenotype Data

`pheno` must be a `data.frame` or matrix with one row per sample.

Required columns:

| Column | Requirement |
|:---|:---|
| `SAMPLE_ID` | Globally unique sample ID |
| `SUBJECT_ID` | Subject ID repeated across visits |
| `FU` | Integer-valued visit index; consecutive from `0`; must include `0` and `1` |
| `TREATMENT_GROUP` | Binary `0/1`; both treatment arms must be present |
| `FEMALE` | Binary `0/1` |

Exactly one outcome schema must be present:

| Outcome type | Required columns |
|:---|:---|
| Continuous | `OUTCOME` |
| Time-to-event | `OUTCOME_TIME`, `OUTCOME_STATUS`; `OUTCOME_TIME` is assumed to be measured in years |

Additional covariates are optional. If supplied, each named column must exist
and be numeric, factor, or logical.

### Omics Data

`omics` must be a data frame with:

| Column | Requirement |
|:---|:---|
| `ANALYTE_NAME` | One feature identifier per row |
| Sample columns | Numeric columns named by `pheno$SAMPLE_ID` |

Extra omics samples and extra phenotype samples are allowed. The validator
intersects them and reports sample counts.

### Omics Type

Accepted values:

- `Proteomics`
- `Metabolomics`
- `DNAm`

DNAm triggers probe-list loading, probe coverage checks, full-probe-list
subsetting, and filtered-probe BH correction.

---

## Validation Flow

`FAST_outcome_WAS()` and `FAST_outcome_WAS_reports()` call the same validation
stack.

### Outcome Type Detection

Function: `.detect_outcome_type()`

| Columns present | Result |
|:---|:---|
| `OUTCOME` only | continuous |
| `OUTCOME_TIME` and `OUTCOME_STATUS` only | time-to-event |
| both schemas | error |
| incomplete time-to-event schema | error |
| no outcome schema | error |

There is no public `outcome_type` argument. The schema determines the model
family.

### Omics Type Validation

Function: `.validate_omics_type()`

- Enforces the accepted `omics_type` values.
- Prints preprocessing reminders, such as DNAm M-values and log2-transformed
  Proteomics/Metabolomics.
- Accepts `"other"` for non-omics analytes without printing a preprocessing
  reminder.

### Phenotype Validation

Function: `.validate_pheno()`

Core checks:

- Confirms `pheno` is a `data.frame` or matrix.
- Checks all required columns.
- Checks `additional_covariates` is `NULL` or character.
- Converts `FU`, `TREATMENT_GROUP`, and `FEMALE` to factors when needed.

Follow-up checks:

- `FU` must be integer-valued after coercion.
- `FU` must be non-negative.
- `FU == 0` must exist.
- `FU == 1` must exist.
- FU levels must be consecutive integers from `0` through `max(FU)`.

Binary field checks:

- `FEMALE` must contain only `0/1`.
- `TREATMENT_GROUP` must contain only `0/1`.
- Both treatment arms must be present.

Duplicate handling:

- `SAMPLE_ID` must be unique globally.
- Duplicate `SUBJECT_ID/FU` pairs are not fatal.
- If duplicates exist, only the first row is kept and a warning is emitted.

Outcome checks:

| Outcome type | Checks |
|:---|:---|
| Continuous | `OUTCOME` must be numeric |
| Time-to-event | `OUTCOME_TIME` must be numeric and non-negative; validation warns that it is assumed to be measured in years; `OUTCOME_STATUS` must be binary `0/1` and is coerced to integer |

Missing-data behavior:

- Samples missing any requested additional covariate are dropped.
- Samples missing required outcome data are dropped.
- Dropped sample counts are reported by messages.

Subject-level constancy:

- `TREATMENT_GROUP` must be constant within `SUBJECT_ID`.
- `FEMALE` must be constant within `SUBJECT_ID`.
- Outcome values must be constant within `SUBJECT_ID`.

Subject retention:

- A subject must have at least one baseline row (`FU == 0`).
- A subject must have at least one nonzero follow-up row.
- This is a dataset-level screen; FU-specific inclusion happens later.

Sex strata:

- The validator returns `all`, `male`, and `female`.
- If the dataset is single-sex, both sex-specific subsets are set to `NULL`.

### Omics Validation

Function: `.validate_omics()`

- Requires `ANALYTE_NAME`.
- Requires all measurement columns to be numeric.
- Intersects omics sample columns with retained phenotype sample IDs.
- Reports shared, omics-only, and pheno-only sample counts.
- Warns on analytes with missing values.
- Warns on analytes with near-zero variance.
- Does not drop analytes for missingness or low variance at validation time.
- Returns aligned omics tables for `all`, `male`, and `female`.

### DNAm Validation

Functions: `.validate_dnam_probe_coverage()`, `.subset_omics_list()`

For DNAm:

1. `outcomewas/Data/FAST_epicv1_epicv2_probe_list.rds` is loaded as the full probe list.
2. `outcomewas/Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds` is loaded as the
   filtered probe list.
3. Coverage against both lists is checked and reported.
4. Omics tables are subset to probes in the full probe list.
5. The filtered probe list is retained for `BH_P_VALUE_FILTERED`.

---

## High-Level Pipeline Flow

```text
FAST_outcome_WAS()
|
|-- resolve n_cores
|-- set future::plan()
|-- .validate_omics_type()
|-- .validate_pheno()       -> pheno_list
|-- .validate_omics()       -> omics_list
|
|-- [DNAm] load probe lists, validate coverage, subset omics_list
|
|-- .run_stratified_analysis(response_type = "change")
|     |-- all
|     |-- male
|     `-- female
|
`-- .run_stratified_analysis(response_type = "level")
      |-- all
      |-- male
      `-- female

FAST_outcome_WAS_reports()
|
|-- same validation and DNAm setup
`-- .generate_reports()
```

---

## Analysis Design

OutcomeWAS runs a separate model set for every nonzero follow-up. If a dataset
contains `FU = 0, 1, 2`, there is one model set for baseline-to-FU1 and another
for baseline-to-FU2.

The analysis dimensions are:

| Dimension | Values |
|:---|:---|
| response type | `change`, `level` |
| stratum | `all`, `male`, `female` |
| follow-up | every nonzero `FU` |
| analyte | every retained row of `omics` |

### FU-Specific Subject Inclusion

For a given follow-up `FU = k`, a subject is included only if they have:

- one retained baseline row at `FU == 0`
- one retained follow-up row at `FU == k`

The subject does not need all intermediate visits. A subject with FU0 and FU2
but no FU1 is eligible for the FU2 model.

### Analysis Dataset Construction

For each stratum and follow-up:

1. `pheno_baseline_all` is all baseline rows in the stratum.
2. `pheno_analysis` is all rows for the target follow-up.
3. Subjects present in both are intersected.
4. Baseline omics values are matched by `SUBJECT_ID`.
5. A model frame is built with the outcome, analyte values, treatment, sex, and
   requested covariates.

The implementation takes `TREATMENT_GROUP`, `FEMALE`, and additional covariates
from the follow-up analysis row. The intended use is that these are effectively
subject-level covariates, but they are not explicitly frozen at baseline.

### Change Versus Level

For each analyte and FU-specific subject:

| Analysis | `analyte` term |
|:---|:---|
| `analysis_change` | `FU_value - baseline_value` |
| `analysis_level` | `FU_value` |

Both include `analyte_baseline = baseline_value` as an adjustment covariate.

---

## Model Formulas

The model builder starts with:

- analyte term of interest: `analyte`
- baseline analyte adjustment: `analyte_baseline`
- trial adjustment: `TREATMENT_GROUP`
- sex adjustment: `FEMALE`
- optional user-provided covariates

Before fitting, `.drop_uninformative_covariates()` removes adjustment terms with
only one observed value in the current stratum/FU model frame. This usually
drops `FEMALE` in sex-stratified analyses.

### Continuous Outcome

Function: `.perform_continuous_analysis()`

Model:

```r
OUTCOME ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

Fit method:

```r
lm(...)
```

Output:

- `coefficients`: all rows from `summary(fit)$coefficients`
- `outcome_effects`: the `analyte` coefficient only

### Time-to-Event Outcome

Function: `.perform_tte_analysis()`

Model:

```r
survival::Surv(OUTCOME_TIME, OUTCOME_STATUS) ~
  analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

Fit method:

```r
survival::coxph(..., ties = "efron")
```

Output:

- `coefficients`: all rows from `summary(fit)$coefficients`
- `outcome_effects`: the `analyte` coefficient only
- `HAZARD_RATIO`: `exp(EFFECT_SIZE)`
- `N_EVENTS`: number of events in the fitted model

Failure behavior:

- If a stratum/FU has zero events, the entire FU-specific model set is skipped
  for that stratum.
- If an individual analyte fit fails, that analyte returns `NULL`, emits a
  warning, and is omitted from the final bound table.

---

## Parallelization and Checkpointing

Parallelization is configured in `FAST_outcome_WAS()`:

| `n_cores` | future plan |
|:---|:---|
| `1` | `future::sequential` |
| `> 1` | `future::multisession` |

Within each response type, stratum, and follow-up:

1. Analytes are split into batches of `checkpoint_batch_size`.
2. Each batch is processed with `furrr::future_map()`.
3. If checkpointing is enabled, each completed batch is saved as an `.rds`.
4. After all batches finish, non-`NULL` analyte results are row-bound.

Checkpoint path:

```text
{checkpoint_dir}/{response_type}/{stratum}/FU{k}/batch_{b}.rds
```

Example:

```text
checkpoints/change/all/FU1/batch_1.rds
checkpoints/change/all/FU2/batch_1.rds
checkpoints/level/female/FU1/batch_3.rds
```

Checkpoint reuse is FU-specific. A rerun with the same checkpoint directory
loads completed batch files and skips recomputation.

Important assumption:

- Checkpoints are tied to analyte ordering and `checkpoint_batch_size`.
- Do not change the omics data or batch size between a run and resume.

---

## Multiple Testing Correction

Function: `.apply_multiple_testing_correction()`

Benjamini-Hochberg correction is applied after all batches for a stratum and
response type are assembled.

| Output table | BH grouping |
|:---|:---|
| `coefficients` | within each `FU x COEFFICIENT` |
| `outcome_effects` | within each `FU` |

This means, for example, that analyte effects at FU1 are corrected separately
from analyte effects at FU2, and separately from baseline or treatment
adjustment terms in the full coefficient table.

---

## DNAm Probe Sets

DNAm-specific behavior is split across `outcomewas/main.R` and
`outcomewas/R/analysis_helpers.R`.

In `outcomewas/main.R`:

1. Load full and filtered probe lists from `outcomewas/Data/`.
2. Validate overlap between incoming omics probes and both reference lists.
3. Subset the omics tables to the full probe list.

In `outcomewas/R/analysis_helpers.R`:

1. Fit models on the full retained probe set.
2. Add `BH_P_VALUE_FILTERED` for analytes in the filtered probe list.
3. Leave `BH_P_VALUE_FILTERED` as `NA` outside the filtered set.

Filtered BH grouping:

| Output table | Filtered BH grouping |
|:---|:---|
| `coefficients` | within each `FU x COEFFICIENT` |
| `outcome_effects` | within each `FU` |

This allows full-set and filtered-set significance thresholds to be compared
without refitting models.

---

## Reporting Pipeline

Function: `.generate_reports()`

`FAST_outcome_WAS_reports()` returns:

```r
list(
  pheno_summary = ...,
  variable_summaries = ...,
  outcome_reports = ...
)
```

Reports are descriptive and independent of response type.

### `pheno_summary`

Function: `.create_pheno_data_report()`

One row per `FU x FEMALE` cell:

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up level |
| `FEMALE` | Sex indicator |
| `N_SUBJECTS` | Unique subjects in the cell |
| `N_CONTROL` | Unique control subjects in the cell |
| `N_TREATMENT` | Unique treatment subjects in the cell |
| `N_SAMPLES` | Raw sample-row count |

### `variable_summaries`

Built separately for `all`, `male`, and `female`.

Within each stratum, summaries are keyed by observed `FU x TREATMENT_GROUP`
cells:

```r
reports$variable_summaries$all$omics_FU0_Tx0
reports$variable_summaries$all$omics_FU1_Tx1
reports$variable_summaries$all$covariates_FU0_Tx0
```

Omics summaries contain:

- `ANALYTE_NAME`
- `N_NONMISSING`
- `MEAN`
- `MEDIAN`
- `SD`
- `MIN`
- `MAX`

Covariate summaries contain:

- `COVARIATE_NAME`
- `TYPE`
- `N_NA`
- `SUMMARY`

### `outcome_reports`

Contains:

- `outcome_type`
- `analysis_sample_summary`
- `outcome_summary`

`analysis_sample_summary` has one row per nonzero FU:

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up level |
| `N_SUBJECTS` | Subjects with both baseline and this FU |
| `N_BASELINE_SAMPLES` | Baseline samples contributing to this FU-specific model set |
| `N_FOLLOWUP_SAMPLES` | Follow-up samples contributing to this FU-specific model set |

`outcome_summary` is built from a deduplicated one-row-per-subject frame so
repeated follow-up rows do not double-count subject-level outcomes.

For continuous outcomes, groups are:

- `all`
- `control`
- `treatment`

For time-to-event outcomes, the summary includes subject counts, event counts,
censoring counts, event rate, and follow-up-time summaries.

---

## Results Output

### Analysis Results

`FAST_outcome_WAS()` returns:

```r
results$analysis_change$all
results$analysis_change$male
results$analysis_change$female

results$analysis_level$all
results$analysis_level$male
results$analysis_level$female
```

Each non-`NULL` stratum contains:

```r
list(
  coefficients = data.frame(...),
  outcome_effects = data.frame(...)
)
```

`coefficients` columns:

- `ANALYTE_NAME`
- `FU`
- `COEFFICIENT`
- `N_OBS`
- `EFFECT_SIZE`
- `SE`
- `P_VALUE`
- `BH_P_VALUE`
- `N_EVENTS` for time-to-event outcomes
- `HAZARD_RATIO` for time-to-event outcomes
- `BH_P_VALUE_FILTERED` for DNAm

`outcome_effects` columns:

- `ANALYTE_NAME`
- `FU`
- `EFFECT_SIZE`
- `SE`
- `P_VALUE`
- `BH_P_VALUE`
- `HAZARD_RATIO` for time-to-event outcomes
- `BH_P_VALUE_FILTERED` for DNAm

`outcome_effects` should match the `COEFFICIENT == "analyte"` rows from
`coefficients`, with only the compact analyte-effect columns retained.

### Report Results

`FAST_outcome_WAS_reports()` returns:

```r
reports$pheno_summary
reports$variable_summaries
reports$outcome_reports
```

The report object is intended to describe input data and outcome availability;
it does not contain model results.
