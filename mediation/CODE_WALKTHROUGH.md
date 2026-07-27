# FAST Mediation Pipeline: Code Walkthrough

This walkthrough documents the product-of-coefficients mediation pipeline in
`mediation/`. The pipeline is much smaller than Track 1.1.1 or OutcomeWAS
because it fits no models. It combines two sets of model results that have
already been produced.

The central distinction is:

- `FAST_mediation()` validates and routes analysis/stratum pairs.
- `.compute_mediation()` performs the complete mediation calculation for one
  analysis/stratum pair.

## Table of Contents

1. [Pipeline Mental Model](#pipeline-mental-model)
2. [File Structure](#file-structure)
3. [Orchestration](#orchestration)
4. [Validation](#validation)
5. [Mediation Computation](#mediation-computation)
6. [Returned Result Object](#returned-result-object)

---

## Pipeline Mental Model

For a given analyte and follow-up:

- Track 1.1.1 supplies `alpha`, the treatment-to-analyte coefficient.
- OutcomeWAS supplies `beta`, the analyte-to-outcome coefficient conditional on
  treatment and the other OutcomeWAS covariates.
- Mediation calculates the indirect effect `alpha * beta`.

The complete flow is:

```text
Track 1.1.1 result object                 OutcomeWAS result object
          |                                        |
          `------------- FAST_mediation() ----------'
                              |
                    validate and route pairs
                              |
              call .compute_mediation() for each shared
                    analysis type and stratum
                              |
                    match alpha/beta rows
                              |
             calculate indirect effects + Sobel inference
```

> **Important:** This pipeline uses fitted summary results only. It does not
> read subject-level phenotype or omics data, and it does not refit either
> upstream model.

---

## File Structure

```text
mediation/
  main.R
  R/
    validation_helpers.R
    analysis_helpers.R
```

| File | Responsibility |
|:---|:---|
| `mediation/main.R` | Public entry point, input-level validation, analysis/stratum traversal, and result metadata |
| `mediation/R/validation_helpers.R` | Definition and validation of the upstream result contract |
| `mediation/R/analysis_helpers.R` | Row alignment, indirect-effect calculation, Sobel inference, and BH correction |

`main.R` sources both helper files and exposes one public function:
`FAST_mediation()`.

---

## Orchestration

File: `mediation/main.R`

```r
FAST_mediation <- function(treatment_results,
                           outcome_results,
                           outcome_type = c("continuous", "tte"))
```

The arguments are:

| Argument | Expected value |
|:---|:---|
| `treatment_results` | Result object returned by Track 1.1.1 `FAST_omics_WAS()` |
| `outcome_results` | Result object returned by OutcomeWAS `FAST_outcome_WAS()` |
| `outcome_type` | `"continuous"` or `"tte"`; determines outcome-table validation and TTE output |

The function itself does little statistical work. Its job is to establish the
correct pair of upstream entries for each call to `.compute_mediation()`.

In execution order it:

1. Resolves `outcome_type` with `match.arg()`.
2. Validates the complete treatment result object.
3. Validates the complete outcome result object against `outcome_type`.
4. Creates output branches for `analysis_change` and `analysis_level`.
5. Within each analysis, finds the union of supported stratum names supplied by
   the two inputs.
6. Calls `.compute_mediation()` when the same stratum is non-`NULL` in both
   inputs.
7. Warns and returns `NULL` for a stratum that is non-`NULL` in only one input.
8. Attaches `outcome_type`, `method`, and `inference` metadata.

Supported analysis types and strata are defined once in
`validation_helpers.R`:

```r
.supported_analyses <- c("analysis_change", "analysis_level")
.supported_strata <- c("all", "male", "female")
```

The important result of this orchestration is that `.compute_mediation()`
always receives one treatment entry and one outcome entry for the same analysis
type and stratum.

---

## Validation

File: `mediation/R/validation_helpers.R`

Validation runs before any analysis/stratum pair is calculated. It establishes
the structural and numeric assumptions used by `.compute_mediation()`.

### Result-Object Structure

`.validate_results_object()` requires both `analysis_change` and
`analysis_level`. Each analysis must name at least one member of the supported
stratum set, but may omit the others. Unknown stratum names fail validation.

Within a non-`NULL` stratum:

- Track 1.1.1 must provide `treatment_effects`.
- OutcomeWAS must provide `outcome_effects`.

The broader upstream `coefficients` tables are not used by mediation.

### Effect-Table Structure

Both compact effect tables require:

| Column | Validation |
|:---|:---|
| `ANALYTE_NAME` | Non-missing and non-empty |
| `FU` | Positive integer |
| `EFFECT_SIZE` | Finite numeric value |
| `SE` | Finite numeric value greater than zero |

`ANALYTE_NAME x FU` must be unique within each table. Duplicate keys are
rejected rather than resolved implicitly during the join.

For time-to-event outcomes, `outcome_effects` must contain `HAZARD_RATIO`. For
continuous outcomes it must not contain that column. This makes an incorrect
`outcome_type` argument fail early.

After validation, the core calculation can assume that every path estimate has
one key, one finite coefficient, and one positive standard error.

---

## Mediation Computation

File: `mediation/R/analysis_helpers.R`

```r
.compute_mediation <- function(treatment_entry,
                               outcome_entry,
                               analysis,
                               stratum,
                               outcome_type)
```

This function performs the full calculation for one analysis type and one
stratum.

### Prepare and Align the Inputs

- Select `treatment_effects` and `outcome_effects`; the broader upstream
  coefficient tables are never used.
- Rename their estimates and SEs to `ALPHA_EFFECT`, `ALPHA_SE`, `BETA_EFFECT`,
  and `BETA_SE`.
- Match rows by `ANALYTE_NAME x FU` using an inner join.
- Record matched and unmatched row counts in `join_diagnostics`.
- Stop for zero overlap; warn and continue for partial overlap.
- If `BH_P_VALUE_FILTERED` is present, require it in both inputs and require the
  same matched rows to identify the filtered probe set. Only its non-`NA`
  pattern is used; the upstream BH values are not reused.

### Calculate the Indirect Effects

For each matched row, let:

- `a = ALPHA_EFFECT`;
- `b = BETA_EFFECT`;
- `SE(a) = ALPHA_SE`;
- `SE(b) = BETA_SE`.

The product-of-coefficients estimate is:

```text
INDIRECT_EFFECT = a * b
```

The first-order Sobel standard error is:

```text
INDIRECT_SE = sqrt(b^2 * SE(a)^2 + a^2 * SE(b)^2)
```

The function then calculates:

```text
Z_VALUE = INDIRECT_EFFECT / INDIRECT_SE
P_VALUE = 2 * pnorm(-abs(Z_VALUE))
```

Confidence intervals are fixed at 95%:

```text
CI_LOWER = INDIRECT_EFFECT - qnorm(0.975) * INDIRECT_SE
CI_UPPER = INDIRECT_EFFECT + qnorm(0.975) * INDIRECT_SE
```

The Sobel variance does not include covariance between `a` and `b`; that
covariance is not available in the two upstream summary-result objects.

If both estimated paths are exactly zero, both the numerator and denominator of
the z statistic are zero. The implementation defines that singular case as
`Z_VALUE = 0` and `P_VALUE = 1`.

### Correct and Return the Results

- Apply BH correction to the mediation p-values within each `FU`. Because this
  function handles one analysis/stratum pair, the complete grouping is
  `analysis type x stratum x FU`.
- For DNAm, separately correct within filtered probes and each `FU`; other
  probes receive `NA` in `BH_P_VALUE_FILTERED`.
- For TTE outcomes, add
  `INDIRECT_HAZARD_RATIO = exp(INDIRECT_EFFECT)`.
- Order rows by `ANALYTE_NAME` and `FU` and return:

```r
list(
  mediation_effects = mediation,
  join_diagnostics = diagnostics
)
```

---

## Returned Result Object

For inputs containing all three strata, the result is accessed as:

```r
results$analysis_change$all
results$analysis_change$male
results$analysis_change$female

results$analysis_level$all
results$analysis_level$male
results$analysis_level$female
```

Only strata named by at least one input appear. A stratum unavailable in one of
the two inputs is present as `NULL` after its warning.

Each calculated stratum contains `mediation_effects` and `join_diagnostics`.
`mediation_effects` retains the four component values used in the calculation:

- `ANALYTE_NAME`, `FU`
- `ALPHA_EFFECT`, `ALPHA_SE`
- `BETA_EFFECT`, `BETA_SE`

The remaining columns are computed by this pipeline:

- `INDIRECT_EFFECT`, `INDIRECT_SE`
- `Z_VALUE`, `P_VALUE`, `CI_LOWER`, `CI_UPPER`, `BH_P_VALUE`
- `BH_P_VALUE_FILTERED` when filtered DNAm inputs are supplied
- `INDIRECT_HAZARD_RATIO` for time-to-event outcomes

The complete object has:

```r
attr(results, "outcome_type")  # "continuous" or "tte"
attr(results, "method")        # "product_of_coefficients"
attr(results, "inference")     # "sobel"
```
