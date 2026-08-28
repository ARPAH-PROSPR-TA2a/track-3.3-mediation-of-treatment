# OutcomeWAS Pipeline

This pipeline, written in R, provides an implementation of omics-wide
association studies for subject-level outcomes in randomized trial datasets.
OutcomeWAS estimates `analyte -> outcome` associations while retaining
`TREATMENT_GROUP` as an adjustment covariate. It is the companion to
TreatmentWAS, which handled `treatment -> analyte`.

The pipeline supports continuous outcomes, time-to-event outcomes, Proteomics,
Metabolomics, and DNA Methylation inputs. It fits separate baseline-adjusted
models for each nonzero follow-up level and automatically stratifies results by
sex when data permits.

The pipeline exposes `FAST_outcome_WAS()` for the statistical analysis and
`FAST_outcome_WAS_reports()` for QC and data summary reports. Optional
proteomics post-processing functions add protein names to compact
`outcome_effects` tables without changing either core result object.

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
  checkpoint_dir = "checkpoints"
)

# View analysis results
head(results$analysis_change$all$coefficients)
head(results$analysis_change$all$outcome_effects)

head(results$analysis_level$all$coefficients)
head(results$analysis_level$all$outcome_effects)
```

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
- **`n_cores`** (`integer`, optional): Number of workers used for analyte-level
  parallel model fitting. Defaults to `max(1, parallel::detectCores() - 1)`,
  with `future::availableCores()` as the fallback when detection is unavailable.
  Set to `1` to run serially.
- **`checkpoint_dir`** (`character`, optional): Directory for per-batch
  checkpoints. If `NULL`, checkpointing is disabled. If provided, completed
  batches are reused when the same run is resumed.
- **`checkpoint_batch_size`** (`integer`): Number of analytes per checkpoint
  batch. Default: `2000`. Only relevant when `checkpoint_dir` is set.

### Return Value

A list with two top-level elements:

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
- **`$outcome_effects`**: Compact analyte-effect table, keeping only the
  analyte term from `coefficients`.

For continuous outcomes, `EFFECT_SIZE` is the linear-model coefficient. For
time-to-event outcomes, `EFFECT_SIZE` is the log hazard ratio and
`HAZARD_RATIO = exp(EFFECT_SIZE)`.

**DNAm only**: `BH_P_VALUE_FILTERED` is added to `coefficients` and
`outcome_effects` for the pre-specified filtered probe set.

## Optional Proteomics Annotation

`get_annotated_outcome_effects()` traverses every available analysis, stratum,
and follow-up in an OutcomeWAS result object and returns a named list of
annotated `outcome_effects` tables. Saving that list is an explicit caller
decision.

```r
protein_annotation <- read.csv(
  "CALERIE_cleaned_protein_translation_table.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
fu_labels <- c("1" = "3mo", "2" = "6mo", "3" = "12mo", "4" = "24mo")

annotated_tables <- get_annotated_outcome_effects(
  results = results,
  protein_annotation = protein_annotation,
  fu_labels = fu_labels
)

saveRDS(annotated_tables, "annotated_outcome_effect_tables.rds")
```

The translation table is joined from `outcome_effects$ANALYTE_NAME` to its
unique, non-missing `AptName` column. It must contain the annotation fields used
by Track 1.1.1: `AptName`, `SomaId`, `TargetFullName`, `Target`, `UniProt`,
`EntrezGeneID`, `EntrezGeneSymbol`, `Uniprot_Unique`, and `Symbol_Unique`.
Missing required columns or duplicate, missing, or empty `AptName` values are
errors. Source column names must be unique and cannot reuse an output annotation
field or `FU_LABEL`. Partial coverage retains unmatched result analytes with
`NA` annotations and produces a warning; zero coverage is an error. Extra
translation rows are allowed. Annotation preserves the original result row
count and order.

For a continuous outcome with `fu_labels` supplied, each annotated table has the
same 15-column contract as the Track 1.1.1 annotated treatment-effect tables:

```text
ANALYTE_NAME, SomaId, TargetFullName, Target, UniProt, EntrezGeneID,
EntrezGeneSymbol, Uniprot_Unique, Symbol_Unique, FU, FU_LABEL,
EFFECT_SIZE, SE, P_VALUE, BH_P_VALUE
```

Time-to-event tables additionally contain `HAZARD_RATIO` immediately after
`EFFECT_SIZE`. Follow-up values are not restricted to the CALERIE `1:4`
convention; `fu_labels` is an explicit named character vector used to populate
`FU_LABEL` and output table names. With `fu_labels = NULL`, `FU_LABEL` is omitted
and list names use suffixes such as `FU1`.

This function annotates only compact `outcome_effects`. It does not modify the
input `results` object, its `coefficients` tables, or any object returned by
`FAST_outcome_WAS_reports()`.

## FAST_outcome_WAS_reports()

Generates QC and data summary reports. Takes the same `pheno`, `omics`,
`omics_type`, and `additional_covariates` arguments as `FAST_outcome_WAS()` and
runs the same input validation, but does not fit models.

### Quick Example

```r
reports <- FAST_outcome_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "Proteomics",
  additional_covariates = c("agebl", "mbmi")
)
```

### Return Value

A list with three top-level elements:

- **`$pheno_summary`**: Study-level counts by `FU` and `FEMALE`.
- **`$variable_summaries`**: Sex-stratified omics and covariate summaries by
  `FU x TREATMENT_GROUP` cell.
- **`$outcome_reports`**: Outcome type, FU-specific analysis sample sizes, and
  subject-level outcome summaries.

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

Checkpoint files are tied to a specific analyte ordering and
`checkpoint_batch_size`. Do not change the `omics` data or batch size between a
run and its resume.

## Data Format Requirements

### Phenotype Data

Phenotype data must be a data frame with one row per sample.

Required columns for all runs:

| Column | Type | Description |
|:---|:---|:---|
| `SAMPLE_ID` | character | Unique sample identifier |
| `SUBJECT_ID` | character | Subject identifier shared across visits |
| `FU` | factor or integer-coded numeric | Visit index; baseline is `0` |
| `TREATMENT_GROUP` | factor or integer-coded numeric | Binary treatment assignment, `0/1` |
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
- `TREATMENT_GROUP` must contain both arms.
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

In both cases, `analyte_baseline` is included as an adjustment covariate.

## Plotting Results

The pipeline includes plotting helpers for QQ and volcano plots. Source the
pipeline-specific plotting module and call `generate_all_plots()` on the full
results object or a stratum-level result object.

```r
source(file.path("outcomewas", "R", "plotting_helpers.R"))

generate_all_plots(results, analysis = "analysis_change")
generate_all_plots(results, figures_dir = "my_figures", analysis = "analysis_level")
```

For DNAm results, both full and filtered probe set plots are generated when
filtered probe results are available.

## Examples and Tests

Run examples and tests from the repository root:

```bash
Rscript outcomewas/Examples/run_proteomics_example.R
Rscript outcomewas/Examples/run_metabolomics_example.R
Rscript outcomewas/Examples/run_DNAm_example.R
Rscript tests/run_tests.R
```

## Optional Practice Data

Large local practice datasets may be kept under `outcomewas/PracticeData/`,
with an optional archive at `outcomewas/PracticeData.zip`. These paths are
ignored by Git and are not required by the test suite; tests use the tracked
fixtures under `outcomewas/Examples/ExampleData/`.

## More Information

- [INPUTS_OUTPUTS.md](INPUTS_OUTPUTS.md): compact input and output schema
  reference.
- [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md): implementation details, control
  flow, and modeling decisions.
