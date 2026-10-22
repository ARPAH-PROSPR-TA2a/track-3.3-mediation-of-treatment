# FAST Mediation Pipeline: Code Walkthrough

Mediation combines fitted summary results; it does not refit either upstream
model or read subject-level data. `FAST_mediation()` validates and routes
analysis/stratum pairs. `.compute_mediation()` performs the calculation for one
such pair.

## Files and Entry Point

| File | Responsibility |
|:---|:---|
| `mediation/main.R` | Public API, normalization/validation, traversal, provenance |
| `mediation/R/validation_helpers.R` | Upstream table and model contracts |
| `mediation/R/analysis_helpers.R` | Row alignment, indirect effects, Sobel inference, BH |

```r
FAST_mediation(
  treatment_results,
  outcome_results,
  outcome_type = c("continuous", "tte")
)
```

Both inputs contain `analysis_change` and `analysis_level`. The supported
strata are `all`, `male`, and `female`; each analysis must name at least one,
but may omit the rest or hold `NULL` entries. The function traverses the union
of supplied supported strata, calculates matching non-`NULL` pairs, and warns
and returns `NULL` when a stratum is available in only one input.

## Validate and Normalize Inputs

The input tables are `$treatment_effects` from Track 1.1.1 and
`$outcome_effects` from OutcomeWAS. Coefficient tables are not used.

| Table | Unique identity within analysis/stratum |
|:---|:---|
| Alpha | `ANALYTE_NAME × FU × CONTRAST` |
| Beta | `ANALYTE_NAME × FU` |

Codes and labels remain separate. Alpha `CONTRAST` identifies the active arm
against control (`k_vs_0`), and `TREATMENT`/`REFERENCE` supply readable labels.
Beta is the common analyte slope across categorical treatment arms. Outcome
tables with an `ARM` column are incompatible and rejected.

All tables require a nonempty analyte identifier, positive integer FU, numeric
`EFFECT_SIZE` and `SE`, and a recognized availability state. An available row
requires finite estimates and a strictly positive SE. Explicit `not_estimable`
or `fit_error` rows require `NA` effect and SE. Invalid available values
and duplicate keys stop the run.

Time-to-event beta tables must contain `HAZARD_RATIO`; continuous tables must
not. This identifies the expected effect scale and catches an inconsistent
`outcome_type` argument.

Complete legacy binary inputs normalize missing contrast fields to `1_vs_0`
and missing statuses to `ok`. They retain the existing requirement for finite
effects and positive SEs. This fallback is not evidence that old output used a
multi-arm model. Metadata is optional for compatibility. Updated outputs
record `analyte_slope = "common_across_treatment_arms"`, and supplied treatment
mappings validate codes and labels. Previous interaction-model outputs are
rejected. Retain upstream provenance because table shape cannot establish
which participants or treatment arms contributed to an old fit.

DNAm membership uses `IN_FILTERED_SET` independently of inference. Matched
inputs must agree on membership. Legacy nonmissing-`BH_P_VALUE_FILTERED`
membership is usable only when every row has available inference; ambiguous
inputs require regeneration. Neither upstream p-values nor upstream adjusted
p-values replace the new mediation tests.

## Align Alpha With the Intended Beta

For each alpha row, preserve `ANALYTE_NAME`, `FU`, `CONTRAST`, `TREATMENT`,
`REFERENCE`, its estimate/SE, and its status.

Match `ANALYTE_NAME × FU`; several alpha contrasts use the same beta row
without losing their identities. The treatment coefficients in the OutcomeWAS
model are adjustment terms, not distinct beta paths.

The output is based on alpha rows. A genuinely missing beta warns but retains
the alpha row with `BETA_STATUS = "not_estimable"`; beta and indirect inference
are `NA`. Explicit unavailable upstream rows are retained in the same way,
while available component estimates remain visible. Zero overlap stops the
pair rather than returning a wholly unmatched calculation.

`join_diagnostics` distinguishes alpha contrast counts from distinct beta-row
counts. `matched_rows` counts matched alpha rows; `matched_outcome_rows` counts
the distinct beta rows they use.

## Calculate the Product and Its Inference

For `k_vs_0`, let `a` be that treatment effect on the mediator and `b` the
common analyte-to-outcome coefficient:

```text
INDIRECT_EFFECT = a * b
INDIRECT_SE = sqrt(b^2 * SE(a)^2 + a^2 * SE(b)^2)
Z_VALUE = INDIRECT_EFFECT / INDIRECT_SE
P_VALUE = 2 * pnorm(-abs(Z_VALUE))
CI_LOWER = INDIRECT_EFFECT - qnorm(0.975) * INDIRECT_SE
CI_UPPER = INDIRECT_EFFECT + qnorm(0.975) * INDIRECT_SE
```

The Sobel method and fixed 95% normal intervals are unchanged. For two exactly
zero available estimates, the implementation defines `Z_VALUE = 0` and
`P_VALUE = 1`. No significance screening of alpha or beta precedes these tests.

Calculations run only for available paths. Preserve `ALPHA_STATUS` and
`BETA_STATUS`; overall `STATUS` is `fit_error` if either path has that state,
otherwise `not_estimable` if either is unavailable, otherwise `ok`. Derived
inference is `NA` for unavailable rows.

### Estimand and Covariance

For continuous outcomes, `alpha_k * beta` gives the same pure and total natural
indirect effect under the common-slope model. Interventions have separate alpha
effects; beta is shared. The outcome model is unchanged from the binary case,
with categorical treatment adjustment and no treatment-by-analyte or
baseline-analyte interactions.

This interpretation requires causal identification assumptions beyond treatment
randomization: mediator/outcome timing, compatible populations and scales,
correct models, and no unmeasured or treatment-induced mediator–outcome
confounding.

The existing Sobel expression omits covariance between alpha and beta.
Indirect effects for different arms are correlated through their common beta,
joint alpha model, and shared controls. Comparing them requires
cross-contrast covariance not exported by these tables.

For Cox outcomes, `INDIRECT_HAZARD_RATIO = exp(INDIRECT_EFFECT)` is retained for
compatibility. It is an approximation, not a general causal indirect hazard
ratio; the product interpretation needs additional conditions such as rare
events, proportional hazards, and a justified post-mediator risk window.
OutcomeWAS does not construct that window or implement landmarking.

## Apply Multiple-Testing Correction

Within each analysis/stratum pair:

- Apply BH to available mediation p-values within `FU × CONTRAST`.
- For DNAm, repeat within `IN_FILTERED_SET` members in each `FU × CONTRAST`.
  Nonmembers and unavailable inference receive `NA` filtered p-values.
- Preserve all alpha identities and readable labels in the result.

These are separate contrast-specific families, not a combined FDR guarantee
across interventions.

## Returned Object

Each calculated stratum returns:

```r
list(
  mediation_effects = mediation,
  join_diagnostics = diagnostics
)
```

`mediation_effects` includes:

- `ANALYTE_NAME`, `FU`, `CONTRAST`, `TREATMENT`, `REFERENCE`;
- `ALPHA_EFFECT`, `ALPHA_SE`, `ALPHA_STATUS`;
- `BETA_EFFECT`, `BETA_SE`, `BETA_STATUS`, overall `STATUS`;
- `INDIRECT_EFFECT`, `INDIRECT_SE`, `Z_VALUE`, `P_VALUE`, `CI_LOWER`, `CI_UPPER`,
  `BH_P_VALUE`;
- `IN_FILTERED_SET` and `BH_P_VALUE_FILTERED` for DNAm;
- `INDIRECT_HAZARD_RATIO` for time-to-event outcomes.

The complete object retains attributes `outcome_type`,
`method = "product_of_coefficients"`, and `inference = "sobel"`, with
`attr(results, "metadata")` describing the common analyte slope, treatment mapping,
estimand, effect scale, BH family, and omitted covariance. See
[INPUTS_OUTPUTS.md](INPUTS_OUTPUTS.md) for the public contract and
[Examples/README.md](Examples/README.md) for runnable examples.
