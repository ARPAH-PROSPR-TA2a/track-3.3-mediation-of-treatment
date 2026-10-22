# OutcomeWAS Pipeline

This pipeline, written in R, provides an implementation of omics-wide
association studies for subject-level outcomes in randomized trial datasets.
OutcomeWAS estimates `analyte -> outcome` associations while retaining
`TREATMENT_GROUP` as an adjustment covariate. It is the companion to
TreatmentWAS, which handled `treatment -> analyte`.

The pipeline supports continuous outcomes, time-to-event outcomes, Proteomics,
Metabolomics, and DNA Methylation inputs. It fits separate baseline-adjusted
models for each nonzero follow-up level and automatically stratifies results by
sex when data permits. Binary calls remain unchanged. Multi-arm studies use the
joint model with categorical treatment adjustment and one common
analyte-to-outcome slope across arms. This is the same model used for binary
studies.

The pipeline exposes `FAST_outcome_WAS()` for the statistical analysis and
`FAST_outcome_WAS_reports()` for QC and data summary reports. Optional
proteomics post-processing returns an annotated copy of the statistical results
with protein metadata on both `coefficients` and `outcome_effects`; reports are
unchanged.

## Installation

Clone the repository:

```bash
git clone https://github.com/ARPAH-PROSPR-TA2a/track-3.3-mediation-of-treatment
```

In your R script, set your working directory to the repository directory and
source the main module:

```r
setwd("path/to/Track3.3")
source(file.path("outcomewas", "main.R"), chdir = TRUE)
```

Run OutcomeWAS from the repository root. The pipeline automatically sources
its implementation files from `outcomewas/R/`.

Core analysis dependencies include `future`, `furrr`, and `survival`. Plotting
helpers additionally require `qqman` and `ggplot2`.

## FAST_outcome_WAS()

The main function for conducting omics-wide outcome association analyses.

### Quick Example

```r
# Load your phenotype and omics data
pheno <- read.csv("my_phenotype.csv")
omics <- read.csv("my_omics_data.csv")

# Run the association analysis
results <- FAST_outcome_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("agebl", "mbmi"),
  n_cores = 8,
  checkpoint_dir = "checkpoints",
  verbose = TRUE
)

# View analysis results
head(results$analysis_change$all$coefficients)
head(results$analysis_change$all$outcome_effects)

head(results$analysis_level$all$coefficients)
head(results$analysis_level$all$outcome_effects)
```

### Multiple treatment arms

`TREATMENT_GROUP` remains the sole assignment column. Use `0` for control and
positive integer codes for active arms; codes need not be consecutive and are
preserved as categories. Supply a label for every observed code when there is
more than one active arm:

```r
treatment_labels <- c(`0` = "Control", `1` = "Diet", `2` = "Exercise")
results <- FAST_outcome_WAS(
  pheno, omics, omics_type = "other", n_cores = 1,
  treatment_labels = treatment_labels
)
reports <- FAST_outcome_WAS_reports(
  pheno, omics, omics_type = "other",
  treatment_labels = treatment_labels
)
```

The result has one common beta per `ANALYTE_NAME × FU` within each analysis and
stratum. The treatment-adjustment coefficients describe arm differences; they
are not the beta path. Different interventions can have different effects on an
analyte (their alpha paths), while the analyte's modeled association with the
outcome (beta) is shared. Mediation therefore computes `alpha_1 * beta`,
`alpha_2 * beta`, and so on.

There is no treatment-by-analyte interaction option. Outcome-effect rows have
no `ARM`, `CONTRAST`, `TREATMENT`, or `REFERENCE` fields. See the synthetic
[three-arm example](Examples/run_multiarm_example.R).

### Parameters

- **`pheno`** (`data.frame`): Phenotype data with required columns. See
  [Data Format Requirements](#data-format-requirements).
- **`omics`** (`data.frame`): Omics data with one row per analyte and one
  column per sample. See [Data Format Requirements](#data-format-requirements).
- **`omics_type`** (`character`): Type of omics data being analyzed. Options:
  `"Proteomics"`, `"Metabolomics"`, `"DNAm"`, or `"other"` for non-omics
  analytes. Default: `"Proteomics"`.
- **`additional_covariates`** (`character` vector, optional): Phenotype columns
  to include as additional model covariates. These columns must be numeric,
  factor, or logical.
- **`treatment_labels`** (named character vector, optional): Appended to both
  public function signatures. Names must match every observed treatment code
  exactly, including `0`; labels must be unique, nonmissing, and nonblank.
  Required for multiple active arms. Binary defaults are `Control` and
  `Treatment`.
- **`n_cores`** (`integer`, optional): Number of workers used for analyte-level
  parallel model fitting. Defaults to `max(1, parallel::detectCores() - 1)`,
  with `future::availableCores()` as the fallback when detection is unavailable.
  Set to `1` to run serially.
- **`checkpoint_dir`** (`character`, optional): Directory for per-batch
  checkpoints. If `NULL`, checkpointing is disabled. If provided, completed
  batches are reused when the same run is resumed.
- **`checkpoint_batch_size`** (`integer`): Number of analytes per checkpoint
  batch. Default: `2000`. Only relevant when `checkpoint_dir` is set.
- **`verbose`** (`logical`): If `TRUE`, emit `[3.3]` progress for validation,
  response/stratum/FU work, checkpoint reuse, worker PIDs, failures, timing, and
  completion. Default: `FALSE`. Set `options(track33.progress_log = "path")`
  to append timestamped progress to a file as well as the console.

### Return Value

A list with the two analysis trees and a `metadata` element:

- **`$analysis_change`**: Change-score analysis, where the analyte term is
  follow-up minus baseline.
- **`$analysis_level`**: Level analysis, where the analyte term is the follow-up
  value.

Each contains results stratified by sex:

- **`$all`**: Results from the full dataset.
- **`$male`**: Results from the male subset, if available.
- **`$female`**: Results from the female subset, if available.

Each non-`NULL` stratum contains:

- **`$coefficients`**: Data frame of all model coefficients for each analyte and
  follow-up.
- **`$outcome_effects`**: One common effect per `ANALYTE_NAME × FU`, including
  an unavailable row when the analyte slope cannot be estimated. `STATUS` is
  `ok`, `not_estimable`, or `fit_error`; unavailable rows retain their identity
  with `NA` inference. A zero-event survival cell remains visible.

`metadata` records the treatment code-to-label mapping, treatment contrasts,
outcome type, additional covariates, common analyte slope, and BH family.
`analyte_slope = "common_across_treatment_arms"` describes the fitted model.
The treatment contrast mapping describes the study; it does not duplicate beta
rows by arm.

For continuous outcomes, `EFFECT_SIZE` is the linear-model coefficient. For
time-to-event outcomes, `EFFECT_SIZE` is the log hazard ratio and
`HAZARD_RATIO = exp(EFFECT_SIZE)`.

**DNAm only**: `BH_P_VALUE_FILTERED` is added to `coefficients` and
`outcome_effects` for the pre-specified filtered probe set. `outcome_effects`
also contains logical `IN_FILTERED_SET`, calculated independently of inference.
A filtered probe remains a member even if its effect is unavailable. BH uses
available p-values within each analysis, stratum, and FU. Coefficient-table
correction additionally separates model terms within FU.

## Optional Proteomics Annotation

`get_annotated_results()` returns a complete annotated copy of an OutcomeWAS
result object. It preserves the nested `analysis_change`/`analysis_level` and
`all`/`male`/`female` topology, including `NULL` strata and any other fields,
while annotating both `coefficients` and `outcome_effects` in each non-`NULL`
stratum. All follow-ups remain together in their original tables.

```r
protein_annotation <- read.csv(
  "CALERIE_cleaned_protein_translation_table.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
fu_labels <- c("1" = "3mo", "2" = "6mo", "3" = "12mo", "4" = "24mo")

annotated_results <- get_annotated_results(
  results = results,
  protein_annotation = protein_annotation,
  fu_labels = fu_labels
)

saveRDS(annotated_results, "Proteomics_3.3OWAS_INF_results_annotated.rds")
```

The translation table is joined from each result table's `ANALYTE_NAME` to its
unique, non-missing `AptName` column. It must contain the annotation fields used
by Track 1.1.1: `AptName`, `SomaId`, `TargetFullName`, `Target`, `UniProt`,
`EntrezGeneID`, `EntrezGeneSymbol`, `Uniprot_Unique`, and `Symbol_Unique`.
Missing required columns or duplicate, missing, or empty `AptName` values are
errors. Source column names must be unique and cannot reuse an output annotation
field or `FU_LABEL`. Partial coverage retains unmatched result rows with `NA`
annotations and produces a warning; zero coverage within any FU is an error.
Extra translation rows are allowed. Annotation preserves row count and order.

For a continuous outcome with `fu_labels` supplied, annotated
`outcome_effects` retain the following fields:

```text
ANALYTE_NAME, SomaId, TargetFullName, Target, UniProt, EntrezGeneID,
EntrezGeneSymbol, Uniprot_Unique, Symbol_Unique, FU, FU_LABEL,
EFFECT_SIZE, SE, P_VALUE, STATUS, BH_P_VALUE
```

Time-to-event tables additionally preserve `HAZARD_RATIO` in source-column
order. Follow-up values are not restricted to the CALERIE `1:4`
convention. `fu_labels` is an optional named character vector whose unique FU
names map to nonblank labels; label values may repeat. With `fu_labels = NULL`,
`FU_LABEL` is omitted.

In `coefficients`, the same eight annotation columns are inserted immediately
after `ANALYTE_NAME`, followed by `FU`, optional `FU_LABEL`, and all remaining
source columns in their original order. Thus the standard continuous order is:

```text
ANALYTE_NAME, SomaId, TargetFullName, Target, UniProt, EntrezGeneID,
EntrezGeneSymbol, Uniprot_Unique, Symbol_Unique, FU, FU_LABEL,
COEFFICIENT, N_OBS, EFFECT_SIZE, SE, P_VALUE, BH_P_VALUE
```

The function does not mutate `results` and performs no file I/O. Reports from
`FAST_outcome_WAS_reports()` are separate and untouched. The CALERIE runners
save the original results and annotated results separately as
`Proteomics_3.3OWAS_INF_results_annotated.rds` and
`Proteomics_3.3OWAS_MetS_results_annotated.rds`.

## FAST_outcome_WAS_reports()

Generates QC and data summary reports. Takes the same `pheno`, `omics`,
`omics_type`, `additional_covariates`, and `treatment_labels` arguments as
`FAST_outcome_WAS()` and runs the same input validation, but does not fit
models.

### Quick Example

```r
reports <- FAST_outcome_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("agebl", "mbmi"),
  verbose = TRUE
)
```

`verbose = TRUE` reports validation and summary-stage progress without changing
the report calculations or return value.

Reports run serially. Omics summaries use bounded numeric-matrix chunks and
preallocated output vectors to avoid repeated data-frame row extraction and
column copying; statistics, sample selection, and output ordering are unchanged.

### Return Value

A list containing:

- **`$pheno_summary`**: Study-level counts by `FU` and `FEMALE`.
- **`$variable_summaries`**: Sex-stratified omics and covariate summaries by
  `FU x TREATMENT_GROUP` cell.
- **`$outcome_reports`**: Outcome type, FU-specific analysis sample sizes, and
  subject-level outcome summaries.
- **`$arm_counts`**: Labeled subject/sample counts for each treatment arm,
  follow-up, and sex.
- **`$metadata`**: The same treatment and model contract as the analysis output.

`N_TREATMENT` pools all active arms. Variable-summary names preserve the codes
(e.g., `omics_FU1_Tx2`); the label mapping identifies them in the reports.

## Checkpointing

For large analyses, passing `checkpoint_dir` enables resumable execution:
analytes are processed in batches, and each completed batch is saved to disk.
If the run is interrupted, re-running the same call skips completed batches and
continues from where it left off.

```r
results <- FAST_outcome_WAS(
  pheno,
  omics,
  omics_type = "DNAm",
  n_cores = 16,
  checkpoint_dir = "checkpoints",
  checkpoint_batch_size = 2000
)
```

Checkpoints are organized by response type, stratum, and follow-up:

```text
checkpoints/
  change/all/FU1/batch_1.rds
  change/all/FU2/batch_1.rds
  change/male/FU1/batch_1.rds
  level/all/FU1/batch_1.rds
  ...
```

Checkpoint files are tied to the inputs, outcome, covariates, visit coding,
analyte ordering, and `checkpoint_batch_size`. Cached batches are loaded by
file existence without checking those settings. Use a new checkpoint directory
when any of them changes, including treatment labels, switching
from M-values to beta values, or updating across this result-schema change. Save/resume
assumes unchanged data and settings; it does not fingerprint the run.

## Data Format Requirements

### Phenotype Data

Phenotype data must be a data frame with one row per sample.

Required columns for all runs:

| Column | Type | Description |
|:---|:---|:---|
| `SAMPLE_ID` | character | Unique sample identifier |
| `SUBJECT_ID` | character | Subject identifier shared across visits |
| `FU` | factor or integer-coded numeric | Visit index; baseline is `0` |
| `TREATMENT_GROUP` | factor or integer-coded numeric | Categorical assignment: `0` = control, positive integers = active arms |
| `FEMALE` | factor or integer-coded numeric | Binary sex indicator, `0/1` |

Outcome columns must follow exactly one schema:

| Outcome Type | Required Columns |
|:---|:---|
| Continuous | `OUTCOME` |
| Time-to-event | `OUTCOME_TIME`, `OUTCOME_STATUS` |

For time-to-event outcomes, `OUTCOME_TIME` must be measured in years. The
validator warns about this assumption but cannot verify the unit from the data
alone.

Data requirements:

- `FU` must be consecutive integers starting at `0` and must include `1`.
- Every retained subject must have baseline and at least one follow-up sample.
- `TREATMENT_GROUP`, `FEMALE`, and outcome values must be constant within
  `SUBJECT_ID`.
- `TREATMENT_GROUP` must contain control and at least one active arm.
- For multiple active arms, `treatment_labels` must name every observed code.
- `SAMPLE_ID` values must be unique.
- Duplicate `SUBJECT_ID/FU` rows are reduced to the first row with a warning.
- Additional covariates must be numeric, factor, or logical. Samples with
  missing additional covariates are dropped.

Example continuous-outcome structure:

```text
SAMPLE_ID   SUBJECT_ID  FU  TREATMENT_GROUP  FEMALE  OUTCOME  agebl  mbmi
sample_001  subj_001    0   1                1       12.4     55     28.1
sample_002  subj_001    1   1                1       12.4     55     28.3
sample_003  subj_002    0   0                0       10.1     62     31.0
sample_004  subj_002    1   0                0       10.1     62     30.7
```

### Omics Data

DNAm beta values and M-values are used on the supplied scale; the pipeline
does not transform or standardize them. The CALERIE beta runners use the same
raw beta matrix and DNAm covariate files as the Track 1.1.1 beta runner. See
[CALERIE DNAm beta runs](INPUTS_OUTPUTS.md#calerie-dnam-beta-runs) for paths,
outputs, and launch commands.

For CALERIE, the INF/MetS outcome is fixed at 24 months and is modeled from
both 12-month (`FU = 1`) and 24-month (`FU = 2`) methylation, with baseline
adjustment. Archive the previous 24-month-only beta output directory before
rerunning: its `FU1` checkpoints contain 24-month fits and cannot be reused.

Omics data must be a data frame with:

| Column | Type | Description |
|:---|:---|:---|
| `ANALYTE_NAME` | character | Feature identifier |
| Sample ID columns | numeric | One column per sample, named as in `pheno$SAMPLE_ID` |

Data requirements:

- `ANALYTE_NAME` is required.
- All measurement columns must be numeric.
- Omics sample columns are intersected with phenotype `SAMPLE_ID` values.
- Extra omics samples and extra phenotype samples are allowed but excluded.
- Analytes with missing values or near-zero variance generate warnings but are
  not dropped during validation.

Example structure:

```text
ANALYTE_NAME  sample_001  sample_002  sample_003  sample_004
cg00000029    0.602       0.515       0.684       0.598
cg00000103    0.456       0.468       0.412       0.401
cg00000109    0.721       0.735       0.691       0.702
```

## Analysis Methods

OutcomeWAS fits separate models for each nonzero follow-up level. For a given
`FU = k`, a subject is included if they have both baseline and FU `k`; they do
not need every intermediate follow-up. OutcomeWAS does not use mixed models.

For continuous outcomes:

```r
OUTCOME ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

For time-to-event outcomes:

```r
survival::Surv(OUTCOME_TIME, OUTCOME_STATUS) ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + [additional covariates]
```

`analysis_change` uses `analyte = FU_value - baseline_value`.

`analysis_level` uses `analyte = FU_value`.

Both analyses adjust for `analyte_baseline`. `TREATMENT_GROUP` is categorical,
and all arms enter one joint model. The analyte slope is common across arms
within each analyte/FU/analysis/stratum fit. No treatment-by-analyte or
baseline-analyte interactions are added. The formulas and inference match the
ordinary binary workflow.

For continuous outcomes, beta is a linear slope. For survival outcomes, it is a
log-hazard slope under proportional hazards. These adjusted associations alone
do not establish causal mediation. Temporal ordering, compatible populations
and scales, and mediator–outcome confounding assumptions require a
study-specific justification. OutcomeWAS does not construct a post-mediator
risk window or perform landmark analysis.

Beta uses each fit's usable observations; filtering can reduce which arms are
represented. A non-estimable beta remains visible with unavailable inference.
Downstream mediation must consider the analysis sample and common-slope
assumption.

## Plotting Results

The pipeline includes plotting helpers for QQ and volcano plots. Source the
pipeline-specific plotting module and call `generate_all_plots()` on the full
results object or a stratum-level result object.

```r
source(file.path("outcomewas", "R", "plotting_helpers.R"))

generate_all_plots(results, analysis = "analysis_change")
generate_all_plots(results, figures_dir = "my_figures", analysis = "analysis_level")
```

For DNAm results, both full and filtered probe set plots use explicit probe
membership. Plots display available common beta estimates. A cell with no
estimable results receives a clearly labeled empty plot.

## Examples and Tests

Run examples and tests from the repository root:

```bash
Rscript outcomewas/Examples/run_proteomics_example.R
Rscript outcomewas/Examples/run_metabolomics_example.R
Rscript outcomewas/Examples/run_DNAm_example.R
Rscript outcomewas/Examples/run_multiarm_example.R /tmp/outcomewas_multiarm
Rscript tests/outcomewas/test_multiarm.R
Rscript tests/outcomewas/test_multiarm_reports.R
Rscript tests/run_tests.R
```

## Optional Practice Data

Large local practice datasets may be kept under `outcomewas/PracticeData/`,
with an optional archive at `outcomewas/PracticeData.zip`. These paths are
ignored by Git and are not required by the test suite; tests use the tracked
fixtures under `outcomewas/Examples/ExampleData/` and generated synthetic
multi-arm fixtures. The multi-arm suites cover pooled LM/Cox references,
unavailable results, labeled reports, filtered membership, and batch resume.

## More Information

- [INPUTS_OUTPUTS.md](INPUTS_OUTPUTS.md): compact input and output schema
  reference.
- [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md): implementation details, control
  flow, and modeling decisions.
