# Mediation: Inputs and Outputs

## Purpose

`FAST_mediation()` performs large-scale, single-mediator
product-of-coefficients analysis from two existing result objects:

- Track 1.1.1 estimates the treatment-to-analyte path, `alpha`.
- OutcomeWAS estimates the analyte-to-outcome path conditional on treatment,
  `beta`.

The function does not fit models or read subject-level phenotype or omics data.
It joins summary results, calculates indirect effects, performs Sobel inference,
and applies mediation-specific multiple-testing correction.

## Public Function

```r
FAST_mediation(
  treatment_results,
  outcome_results,
  outcome_type = c("continuous", "tte")
)
```

### Arguments

- `treatment_results`: complete `FAST_omics_WAS()` result object from Track
  1.1.1.
- `outcome_results`: complete `FAST_outcome_WAS()` result object.
- `outcome_type`: `"continuous"` for linear outcome models or `"tte"` for
  Cox proportional-hazards outcome models.

## Required Upstream Structure

Both inputs must contain `analysis_change` and `analysis_level`. Within each
analysis, they may contain any non-empty subset of the supported strata:

```text
analysis_change/
  all/
  male/
  female/
analysis_level/
  all/
  male/
  female/
```

For example, a single-sex dataset may provide only `all`. A stratum may also be
present with value `NULL`. Matching non-`NULL` strata are analyzed. A stratum
available in only one input produces a warning and a `NULL` output entry; the
remaining matching strata continue. Within each available stratum:

- Track 1.1.1 must provide `$treatment_effects`.
- OutcomeWAS must provide `$outcome_effects`.

Each effect table requires:

| Column | Meaning |
|:---|:---|
| `ANALYTE_NAME` | Mediator identifier |
| `FU` | Positive integer follow-up level |
| `EFFECT_SIZE` | Path coefficient estimate |
| `SE` | Strictly positive standard error |

For time-to-event outcomes, OutcomeWAS `$outcome_effects` must also contain
`HAZARD_RATIO`, confirming that `EFFECT_SIZE` is on the log-hazard scale.

## Alignment Preconditions

Results are joined within the same analysis type and sex stratum using
`ANALYTE_NAME` and `FU`. The caller is responsible for ensuring that the two
upstream analyses use compatible:

- treatment coding;
- analyte identifiers, transformations, and units;
- follow-up definitions;
- study population and inclusion rules;
- baseline adjustment and additional covariates.

The summary outputs cannot prove those provenance conditions. Duplicate keys
or zero overlap stop the analysis. Partial overlap is allowed but produces a
warning and explicit join diagnostics.

## Statistical Calculation

For each matched analyte and follow-up:

```text
indirect effect = alpha * beta
```

The conventional first-order Sobel standard error is:

```text
sqrt(beta^2 * SE_alpha^2 + alpha^2 * SE_beta^2)
```

The pipeline calculates a standard-normal z statistic, two-sided p-value, and
symmetric 95% normal-theory confidence interval. If both estimated paths are
exactly zero, the first-order standard error is zero; the implementation defines
`Z_VALUE = 0` and `P_VALUE = 1`.

No component-path significance filter is applied before mediation testing.
Mediation p-values receive Benjamini-Hochberg correction separately within each
analysis type, sex stratum, and `FU`.

If both upstream tables contain `BH_P_VALUE_FILTERED`, their non-`NA` patterns
must identify the same filtered DNAm probes. The mediation output then includes
`BH_P_VALUE_FILTERED`, corrected within filtered probes separately by `FU`.

## Output Structure

The returned object contains the union of the stratum names supplied by the two
upstream result objects. For the complete three-stratum case:

```text
analysis_change/
  all/
    mediation_effects
    join_diagnostics
  male/
  female/
analysis_level/
  ...
```

Each `mediation_effects` table contains:

| Column | Meaning |
|:---|:---|
| `ANALYTE_NAME`, `FU` | Join key |
| `ALPHA_EFFECT`, `ALPHA_SE` | Treatment-to-analyte estimate and SE used in the calculation |
| `BETA_EFFECT`, `BETA_SE` | Analyte-to-outcome estimate and SE used in the calculation |
| `INDIRECT_EFFECT` | `ALPHA_EFFECT * BETA_EFFECT` |
| `INDIRECT_SE` | Sobel delta-method standard error |
| `Z_VALUE`, `P_VALUE` | Sobel test statistic and raw p-value |
| `CI_LOWER`, `CI_UPPER` | 95% normal-theory confidence interval |
| `BH_P_VALUE` | Mediation BH correction within `FU` |
| `BH_P_VALUE_FILTERED` | Filtered-probe mediation BH correction, when applicable |
| `INDIRECT_HAZARD_RATIO` | `exp(INDIRECT_EFFECT)` for time-to-event outcomes |

`join_diagnostics` reports treatment rows, outcome rows, matched rows, and rows
found only in either input.

The returned object has `outcome_type`, `method = "product_of_coefficients"`,
and `inference = "sobel"` attributes.
