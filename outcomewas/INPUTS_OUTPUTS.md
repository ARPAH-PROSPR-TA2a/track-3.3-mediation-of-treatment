# FAST OutcomeWAS: Inputs and Outputs

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
</style>

This document describes exactly what the OutcomeWAS pipeline expects as input
and what it returns as output. For implementation details, see
[CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md).

---

## Inputs

### `pheno` - Phenotype Data

`pheno` is a data frame with one row per sample. A subject can appear multiple
times, once per follow-up visit.

Required columns for all runs:

| Column | Type | Requirement |
|:---|:---|:---|
| `SAMPLE_ID` | character | Unique sample identifier; no duplicates |
| `SUBJECT_ID` | character | Subject identifier shared across visits |
| `FU` | factor or integer-coded numeric | Visit index; baseline is `0`; values must be consecutive integers and include `0` and `1` |
| `TREATMENT_GROUP` | factor or integer-coded numeric | Binary treatment assignment; `0` = control, `1` = treatment; both arms must be present |
| `FEMALE` | factor or integer-coded numeric | Binary sex indicator; `0` = male, `1` = female |

Outcome columns must follow exactly one of these two schemas:

| Outcome type | Required columns | Requirement |
|:---|:---|:---|
| Continuous | `OUTCOME` | Numeric; constant within each `SUBJECT_ID` |
| Time-to-event | `OUTCOME_TIME`, `OUTCOME_STATUS` | `OUTCOME_TIME` is numeric, non-negative, and measured in years; `OUTCOME_STATUS` is binary `0/1`; both are constant within each `SUBJECT_ID` |

Optional additional covariates:

- Named through the `additional_covariates` argument.
- Must exist as columns in `pheno`.
- Must be numeric, factor, or logical.
- Missing values are allowed, but samples missing any requested covariate are
  dropped before analysis.

Example continuous-outcome phenotype data:

```text
SAMPLE_ID  SUBJECT_ID  FU  TREATMENT_GROUP  FEMALE  OUTCOME  agebl  mbmi
s001       subj_01     0   1                1       12.4     55     28.1
s002       subj_01     1   1                1       12.4     55     28.3
s003       subj_02     0   0                0       10.1     62     31.0
s004       subj_02     1   0                0       10.1     62     30.7
s005       subj_03     0   1                0       15.0     48     26.2
s006       subj_03     2   1                0       15.0     48     26.4
```

Example time-to-event outcome columns, with time measured in years:

```text
SUBJECT_ID  OUTCOME_TIME  OUTCOME_STATUS
subj_01     4.2           1
subj_02     5.0           0
subj_03     2.8           1
```

Important row-retention behavior:

- Duplicate `SUBJECT_ID/FU` rows are reduced to the first row, with a warning.
- Samples with missing outcome data are dropped.
- Samples with missing requested covariates are dropped.
- Subjects must have baseline and at least one follow-up sample after filtering.
- For a specific follow-up `FU = k`, a subject is analyzed only if they have
  both baseline and `FU = k`.

---

### `omics` - Omics Data

`omics` is a data frame with one row per analyte and one column per sample.

| Column | Type | Requirement |
|:---|:---|:---|
| `ANALYTE_NAME` | character | Feature identifier, such as protein, metabolite, or CpG probe |
| Sample ID columns | numeric | One column per sample, named exactly as in `pheno$SAMPLE_ID` |

Example omics data:

```text
ANALYTE_NAME  s001    s002    s003    s004    s005    s006
protein_A     1.204   1.318   0.987   1.052   1.401   1.289
protein_B     0.523   0.489   0.601   0.578   0.512   0.531
protein_C     2.017   2.103   1.889   1.944   2.201   2.088
```

Sample matching:

- Samples present in `omics` but not in `pheno` are excluded.
- Samples present in `pheno` but not in `omics` are excluded.
- Analytes with missing values or near-zero variance generate warnings but are
  not dropped during validation.

---

### Other Input Arguments

- `omics_type`: Character string indicating the data type. Must be one of
  `"Proteomics"`, `"Metabolomics"`, `"DNAm"`, or `"other"`. Use `"other"`
  for non-omics analytes; it does not emit a preprocessing reminder. DNAm runs
  also check probe-list coverage, subset to the full probe list, and add a
  filtered-probe BH correction column.
- `additional_covariates`: Optional character vector naming extra phenotype
  columns to include as model covariates. These covariates are taken from the
  follow-up row used for each analysis.
- `n_cores`: `FAST_outcome_WAS()` only. Number of workers for parallel model
  fitting. Defaults to `max(1, parallel::detectCores() - 1)`, with
  `future::availableCores()` as the fallback when detection is unavailable.
  Set to `1` for serial runs.
- `checkpoint_dir`: `FAST_outcome_WAS()` only. Optional directory for resumable
  batch checkpoints. If set, checkpoint files are written under:

```text
checkpoint_dir/
  change/all/FU1/batch_1.rds
  change/male/FU1/batch_1.rds
  level/all/FU2/batch_3.rds
  ...
```

- `checkpoint_batch_size`: `FAST_outcome_WAS()` only. Number of analytes per
  checkpoint batch. Default: `2000`.

---

## Models Used

OutcomeWAS runs two analyte definitions for every eligible nonzero follow-up:

- `analysis_change`: analyte value at follow-up minus analyte value at baseline.
- `analysis_level`: analyte value at follow-up.

Both analyses adjust for the baseline analyte value.

For a continuous outcome, the model is:

```r
OUTCOME ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

For a time-to-event outcome, the model is:

```r
survival::Surv(OUTCOME_TIME, OUTCOME_STATUS) ~
  analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

Models are fit separately for each:

- analyte
- nonzero follow-up level
- analysis type: `change` and `level`
- stratum: `all`, `male`, and `female` when available

OutcomeWAS does not use mixed models. If a study has `FU = 0, 1, 2`, the
pipeline fits one set of models for baseline-to-FU1 and another set for
baseline-to-FU2.

In sex-stratified analyses, `FEMALE` is usually dropped automatically because it
is constant within the stratum. Other covariates with only one observed value in
a specific model subset are also dropped.

---

## Outputs

### `FAST_outcome_WAS()`

Returns:

```r
list(
  analysis_change = list(all = ..., male = ..., female = ...),
  analysis_level  = list(all = ..., male = ..., female = ...)
)
```

Each non-`NULL` stratum contains two tables:

- `coefficients`: all model terms.
- `outcome_effects`: analyte term only; this is the primary result table for
  most downstream review.

Sex-specific strata are `NULL` when the subset is unavailable. For
time-to-event analyses, a follow-up/stratum with zero events is skipped with a
warning.

---

### `$coefficients`

One row per analyte, follow-up, and model term.

Example:

```text
ANALYTE_NAME  FU  COEFFICIENT        N_OBS  EFFECT_SIZE  SE     P_VALUE  BH_P_VALUE
protein_A     1   analyte            178    0.312        0.067  0.000    0.002
protein_A     1   analyte_baseline   178    0.794        0.041  0.000    0.000
protein_A     1   TREATMENT_GROUP1   178    0.105        0.061  0.086    0.210
protein_B     1   analyte            178   -0.028        0.061  0.647    0.712
```

Common columns:

| Column | Meaning |
|:---|:---|
| `ANALYTE_NAME` | Analyte identifier |
| `FU` | Follow-up level for this model |
| `COEFFICIENT` | Model term |
| `N_OBS` | Number of observations used in the fit |
| `EFFECT_SIZE` | Coefficient estimate |
| `SE` | Standard error |
| `P_VALUE` | Raw p-value |
| `BH_P_VALUE` | Benjamini-Hochberg adjusted p-value within each `FU x COEFFICIENT` group |

Additional columns:

| Column | Present when | Meaning |
|:---|:---|:---|
| `N_EVENTS` | Time-to-event outcome | Number of events in the fitted model |
| `HAZARD_RATIO` | Time-to-event outcome | `exp(EFFECT_SIZE)` |
| `BH_P_VALUE_FILTERED` | DNAm | BH correction within the filtered probe set; `NA` for probes outside that set |

`EFFECT_SIZE` is interpreted as:

- Continuous outcome: linear-model coefficient.
- Time-to-event outcome: log hazard ratio.

---

### `$outcome_effects`

One row per analyte and follow-up. This table keeps only the analyte term from
`coefficients`.

Example:

```text
ANALYTE_NAME  FU  EFFECT_SIZE  SE     P_VALUE  BH_P_VALUE
protein_A     1   0.312        0.067  0.000    0.002
protein_B     1  -0.028        0.061  0.647    0.712
protein_C     2   0.198        0.077  0.011    0.044
```

Common columns:

| Column | Meaning |
|:---|:---|
| `ANALYTE_NAME` | Analyte identifier |
| `FU` | Follow-up level for this model |
| `EFFECT_SIZE` | Analyte coefficient estimate |
| `SE` | Standard error |
| `P_VALUE` | Raw p-value |
| `BH_P_VALUE` | Benjamini-Hochberg adjusted p-value within each `FU` |

Additional columns:

| Column | Present when | Meaning |
|:---|:---|:---|
| `HAZARD_RATIO` | Time-to-event outcome | `exp(EFFECT_SIZE)` |
| `BH_P_VALUE_FILTERED` | DNAm | BH correction within the filtered probe set; `NA` for probes outside that set |

---

### Optional annotated proteomics tables

Protein annotation is a separate post-processing step for `outcome_effects`.
It does not change the object returned by `FAST_outcome_WAS()`, any
`coefficients` table, or the object returned by `FAST_outcome_WAS_reports()`.

```r
get_annotated_outcome_effects(
  results,
  protein_annotation,
  analysis_type = "change",
  group = "all",
  fu = NULL,
  fu_labels = NULL
)

write_annotated_outcome_effect_tables(
  results,
  protein_annotation,
  output_rds,
  fu_labels = NULL
)
```

`results` must be an OutcomeWAS result object. `analysis_type` is `"change"`
or `"level"`; `group` is `"all"`, `"male"`, or `"female"`. `fu` may be one
or more follow-up values or `NULL` for every available follow-up. Follow-ups are
not restricted to `1:4`.

`protein_annotation` is a data frame requiring these columns:

| Column | Role |
|:---|:---|
| `AptName` | Unique translation key matched to `outcome_effects$ANALYTE_NAME` |
| `SomaId` | SomaScan identifier |
| `TargetFullName` | Full target name |
| `Target` | Short target name |
| `UniProt` | UniProt identifier supplied by the translation table |
| `EntrezGeneID` | Entrez gene identifier |
| `EntrezGeneSymbol` | Entrez gene symbol |
| `Uniprot_Unique` | Uniqueness-qualified UniProt value |
| `Symbol_Unique` | Uniqueness-qualified gene-symbol value |

Missing required columns are errors. `AptName` must be non-missing, non-empty,
and unique; violations are errors because they would make the annotation join
ambiguous. The join preserves the OutcomeWAS row order and row count. Analytes
without an `AptName` match remain in the output with `NA` annotation fields and
produce a warning. If no selected analyte matches, the function stops.
Translation rows unused by the selected result table are allowed.

`fu_labels`, when supplied, is a named character vector whose names are the
character representations of `FU`, for example:

```r
c("1" = "3mo", "2" = "6mo", "3" = "12mo", "4" = "24mo")
```

It must cover every FU returned by the call. It supplies the `FU_LABEL` field
and the follow-up portion of names in the writer's returned list. This explicit
mapping supports arbitrary follow-up values without imposing CALERIE's
four-visit convention. When `fu_labels = NULL`, `FU_LABEL` is omitted and writer
names use `FU` plus the integer value, such as `change_all_FU1`.

For continuous outcomes with `fu_labels` supplied, each annotated table has
exactly these 15 columns, in the same order as Track 1.1.1 annotated
treatment-effect tables:

```text
ANALYTE_NAME
SomaId
TargetFullName
Target
UniProt
EntrezGeneID
EntrezGeneSymbol
Uniprot_Unique
Symbol_Unique
FU
FU_LABEL
EFFECT_SIZE
SE
P_VALUE
BH_P_VALUE
```

Without `fu_labels`, the continuous table omits `FU_LABEL` and has 14 columns.
For time-to-event outcomes, the same contract also contains `HAZARD_RATIO`
immediately after `EFFECT_SIZE`. The writer returns the named list invisibly
after saving it to `output_rds`; the output directory must already exist.

---

### `FAST_outcome_WAS_reports()`

Returns:

```r
list(
  pheno_summary = ...,
  variable_summaries = ...,
  outcome_reports = ...
)
```

### `pheno_summary`

Counts by follow-up and sex.

```text
FU  FEMALE  N_SUBJECTS  N_CONTROL  N_TREATMENT  N_SAMPLES
0   0       85          42         43           85
0   1       93          46         47           93
1   0       83          41         42           83
1   1       91          45         46           91
```

Columns:

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up level |
| `FEMALE` | Sex indicator |
| `N_SUBJECTS` | Number of unique subjects |
| `N_CONTROL` | Number of unique control subjects |
| `N_TREATMENT` | Number of unique treatment subjects |
| `N_SAMPLES` | Number of sample rows |

### `variable_summaries`

Sex-stratified summaries of omics values and requested covariates, grouped by
follow-up and treatment arm.

Example entries:

```text
reports$variable_summaries$all$omics_FU0_Tx0
reports$variable_summaries$all$omics_FU1_Tx1
reports$variable_summaries$all$covariates_FU0_Tx0
reports$variable_summaries$female$covariates_FU2_Tx1
```

Each omics entry contains:

- `ANALYTE_NAME`
- `N_NONMISSING`
- `MEAN`
- `MEDIAN`
- `SD`
- `MIN`
- `MAX`

Each covariate entry contains:

- `COVARIATE_NAME`
- `TYPE`
- `N_NA`
- `SUMMARY`

### `outcome_reports`

Contains:

- `outcome_type`: `"continuous"` or `"tte"`.
- `analysis_sample_summary`: sample size available for each FU-specific model
  set.
- `outcome_summary`: subject-level outcome summary.

`analysis_sample_summary` has one row per nonzero follow-up:

| Column | Meaning |
|:---|:---|
| `FU` | Follow-up level |
| `N_SUBJECTS` | Subjects with both baseline and this follow-up |
| `N_BASELINE_SAMPLES` | Baseline samples contributing to this follow-up analysis |
| `N_FOLLOWUP_SAMPLES` | Follow-up samples contributing to this follow-up analysis |

`outcome_summary` is computed after reducing the phenotype data to one row per
subject, so repeated follow-up rows do not double-count the same outcome.

Continuous outcome summary columns:

- `GROUP`
- `N_SUBJECTS`
- `MEAN`
- `MEDIAN`
- `SD`
- `MIN`
- `MAX`

Time-to-event outcome summary columns:

- `GROUP`
- `N_SUBJECTS`
- `N_EVENTS`
- `N_CENSORED`
- `EVENT_RATE`
- `TIME_MEAN`
- `TIME_MEDIAN`
- `TIME_MIN`
- `TIME_MAX`

---

## Failure Behavior

- Invalid required columns or invalid field types stop the run.
- Non-consecutive `FU` encoding stops the run.
- Non-constant within-subject outcome, sex, or treatment stops the run.
- Time-to-event validation warns that `OUTCOME_TIME` is assumed to be measured
  in years.
- Analyte-level model failures are omitted from the final result tables, with
  warnings.
- Time-to-event follow-up/stratum combinations with zero events are skipped,
  with warnings.
