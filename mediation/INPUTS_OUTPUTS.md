# Mediation: Inputs and Outputs

## Public Function

`FAST_mediation()` combines existing Track 1.1.1 treatment-to-analyte estimates
(`alpha`) with OutcomeWAS analyte-to-outcome estimates (`beta`). It does not
read subject-level data or fit models.

```r
FAST_mediation(
  treatment_results,
  outcome_results,
  outcome_type = c("continuous", "tte")
)
```

- `treatment_results`: result object from Track 1.1.1 `FAST_omics_WAS()`.
- `outcome_results`: result object from OutcomeWAS `FAST_outcome_WAS()`.
- `outcome_type`: `"continuous"` for linear models or `"tte"` for Cox models.

OutcomeWAS fits the existing common-beta model across all categorical
treatment arms. Mediation uses its analyte coefficient for every treatment
contrast at the same analyte and follow-up; it does not fit or select an
outcome model.

## Upstream Structure

Both inputs must contain `analysis_change` and `analysis_level`. Each analysis
may contain a nonempty subset of `all`, `male`, and `female`; strata may be
`NULL`. Matching non-`NULL` strata are analyzed. A stratum available in only one
input produces a warning and a `NULL` output entry.

Within an available stratum, use `$treatment_effects` from Track 1.1.1 and
`$outcome_effects` from OutcomeWAS. The broader coefficient tables are unused.

Common effect-table fields:

| Column | Contract |
|:---|:---|
| `ANALYTE_NAME` | Nonmissing, nonempty mediator identifier |
| `FU` | Positive integer follow-up code |
| `EFFECT_SIZE` | Numeric path coefficient; finite for available inference |
| `SE` | Numeric standard error; finite and strictly positive for available inference |
| `STATUS` | `ok`, `not_estimable`, or `fit_error`; unavailable rows require `NA` effect and SE |

Additional identity fields:

| Table | Unique key | Label fields |
|:---|:---|:---|
| Treatment effects | `ANALYTE_NAME × FU × CONTRAST` | `TREATMENT`, `REFERENCE` |
| Outcome effects | `ANALYTE_NAME × FU` | No per-row arm label |

`CONTRAST` retains the upstream active-versus-control identifier, e.g.,
`1_vs_0`. Outcome beta is the common analyte slope across arms, not a
treatment contrast. Outcome tables containing an `ARM` column are incompatible
and rejected.

Time-to-event outcome effects must also contain `HAZARD_RATIO`, identifying
`EFFECT_SIZE` as a log-hazard slope. Continuous outcome tables must not contain
that field. Unavailable inference must be declared by status with both `EFFECT_SIZE` and
`SE` set to `NA`; malformed available estimates and duplicate keys are errors.

Complete legacy binary tables without contrast/status fields remain supported:
alpha identity normalizes to `1_vs_0` and missing statuses to `ok`. This
compatibility does not make incomplete legacy inference valid or infer a
multi-arm model from an old binary output.

Metadata is optional for compatibility. Updated OutcomeWAS output describes
`analyte_slope = "common_across_treatment_arms"`. Supplied treatment-code/label
mappings must agree with row labels and requested arms. Previous interaction
model outputs are incompatible and must be regenerated with the common-beta
model. Retain upstream metadata: table shape alone cannot prove which study
arms contributed to an old beta fit.

### Follow-up labels

Labels are read from each upstream result's `metadata$fu_labels` and
from `FU_LABEL` columns in its effect tables. There is no separate mediation
`fu_labels` argument. Numeric FU codes remain the join keys; labels never
remap visits.

For matching numeric FU codes:

| Label availability | Behavior |
|:---|:---|
| Neither input supplies labels | Continue without label fields or a label warning |
| Only one input supplies labels | Inherit those labels and warn once per call that label agreement could not be checked |
| Both supply matching labels | Continue with those labels and no label warning |
| Both supply conflicting labels | Error; resolve the disagreement before mediation |

Input validation and partial joins follow these rules:

- Metadata maps use unique non-negative integer code names and nonmissing,
  nonblank character labels. Baseline `0` is allowed but not required because
  mediation effect tables contain only follow-ups. An explicitly supplied map
  must cover that source's observed effect-table visits.
- When present, `FU_LABEL` is a nonmissing, nonblank character column. Each
  code must have one consistent label across that source's analyses and strata;
  table labels must agree with metadata when both are supplied.
- The two inputs must agree exactly on shared codes, including `0` when both
  metadata maps contain it. `"12mo"` and `"Month 12"` are a conflict;
  neither source takes precedence and no textual equivalence is inferred.
- Either source can supply labels. Their known mappings are combined without
  requiring identical follow-up sets. Partial joins keep their existing rows
  and statuses; a retained legacy visit with no known label gets `NA`.

Only presentation fields change. Products, Sobel standard errors, confidence
intervals, hazard-ratio approximations, BH values, and join diagnostics are
unchanged.

## Alignment and Availability

The caller must ensure compatible treatment coding, analyte identifiers/scales,
follow-up definitions, study populations, inclusion rules, baseline adjustment,
and additional covariates. Summary tables cannot prove those conditions.

For each alpha row, mediation preserves its contrast and labels, then matches
the beta with the same `ANALYTE_NAME × FU`. Several alpha contrasts can use the
same beta without losing their identities.

Within matching strata, alpha rows are retained even if their beta is absent.
Missing beta rows warn and receive `BETA_STATUS = "not_estimable"`, with `NA`
beta and indirect inference. Zero overlap is an error. An explicit unavailable
alpha or beta row likewise remains visible; available component estimates are
retained. Overall `STATUS` uses this precedence: `fit_error`, then
`not_estimable`, then `ok`. A failed beta affects every matching contrast.

For DNAm, both inputs must agree on the filtered probe set. Use logical
`IN_FILTERED_SET`, which describes membership independently of available
inference. An included probe with a failed fit remains a member. The legacy
nonmissing-`BH_P_VALUE_FILTERED` convention is accepted only when all rows
have available inference; regenerate ambiguous inputs with explicit membership.
Upstream adjusted p-values are not reused for mediation inference.

## Statistical Contract

For contrast `k_vs_0`:

```text
indirect effect = alpha_k * beta
```

The alpha effect can differ by intervention; beta is common across treatment
arms within each analyte/FU/analysis/stratum. For continuous outcomes, the pure
and total natural indirect effects coincide under this linear model without
treatment-by-mediator interactions.

A causal interpretation additionally requires appropriate temporal ordering,
compatible populations and scales, correct model specification, and mediation
identification assumptions, including no unmeasured mediator–outcome
confounding or treatment-induced mediator–outcome confounding. Treatment
randomization alone does not establish these assumptions.

The existing first-order Sobel standard error remains:

```text
sqrt(beta^2 * SE_alpha_k^2 + alpha_k^2 * SE_beta^2)
```

The pipeline returns a standard-normal z statistic, two-sided p-value, and symmetric 95% interval.
When both available path estimates are exactly zero, it defines `Z_VALUE = 0`
and `P_VALUE = 1`. No component-path significance filter is applied.

Sobel variance omits alpha–beta covariance. Separate indirect effects are
correlated through their shared beta, joint alpha model, and shared controls.
Comparing interventions requires covariance absent from these summary tables.
Significance for one arm alone does not establish a difference between arms.

For survival outcomes, beta is a Cox log-hazard slope, and the output retains
`INDIRECT_HAZARD_RATIO = exp(INDIRECT_EFFECT)` for compatibility. This is a
product-method approximation, not a general causal indirect hazard ratio.
Causal interpretation requires additional conditions, including a rare-event
approximation for this method, proportional hazards, and an appropriate
post-mediator risk window. OutcomeWAS does not construct that window or perform
landmarking; the caller must address timing.

BH correction uses available mediation p-values within each **analysis ×
stratum × FU × CONTRAST**. Filtered BH uses members of `IN_FILTERED_SET` within
the same families. There is no combined FDR guarantee across contrasts.

## Output

Each matching analysis/stratum contains:

```r
list(
  mediation_effects = data.frame(...),
  join_diagnostics = list(...)
)
```

The result preserves the union of supplied supported strata; unavailable
whole-stratum pairs remain `NULL` after a warning.

| Field | Meaning |
|:---|:---|
| `ANALYTE_NAME`, `FU`, `CONTRAST` | Unique indirect-effect identity |
| `FU_LABEL` | When any upstream visit labels exist: display label beside `FU`; `NA` for unknown legacy visits |
| `TREATMENT`, `REFERENCE` | Readable active and control labels from alpha |
| `ALPHA_EFFECT`, `ALPHA_SE`, `ALPHA_STATUS` | Treatment-to-analyte path |
| `BETA_EFFECT`, `BETA_SE`, `BETA_STATUS` | Common analyte-to-outcome path |
| `STATUS` | Overall availability |
| `INDIRECT_EFFECT`, `INDIRECT_SE` | Product estimate and Sobel SE |
| `Z_VALUE`, `P_VALUE` | Sobel statistic and raw p-value |
| `CI_LOWER`, `CI_UPPER` | Fixed 95% normal-theory interval |
| `BH_P_VALUE` | BH within FU and contrast, separately by analysis and stratum |
| `IN_FILTERED_SET`, `BH_P_VALUE_FILTERED` | DNAm membership and filtered BH, when applicable |
| `INDIRECT_HAZARD_RATIO` | Exponentiated product for time-to-event outcomes |

Unavailable rows preserve identity and any available component estimates;
derived inference is `NA`.

`join_diagnostics` contains:

| Field | Count |
|:---|:---|
| `treatment_rows` | Alpha contrast rows requested, including unavailable alpha |
| `outcome_rows` | All beta rows supplied |
| `matched_rows` | Alpha rows with a matching beta key, regardless of status |
| `treatment_only_rows` | Alpha rows whose beta key is absent; retained as placeholders |
| `outcome_only_rows` | Beta rows without a matching alpha analyte/FU |
| `matched_outcome_rows` | Distinct beta rows matched; can be fewer than `matched_rows` because contrasts share beta |

The returned object retains attributes `outcome_type`,
`method = "product_of_coefficients"`, and `inference = "sobel"`.
`attr(results, "metadata")` records:

- `schema_version`, `analyte_slope = "common_across_treatment_arms"`,
  `treatment_labels`, and `treatments`;
- `estimand`: `pure_and_total_natural_indirect_effect_under_shared_beta` for
  continuous outcomes or `cox_product_of_coefficients_approximation`;
- `effect_scale`: `outcome` or `log_hazard`;
- `bh_family`: `BH within analysis, stratum, FU and contrast`;
- `alpha_beta_covariance = "not_included"` and
  `cross_contrast_covariance = "not_available"`;
- `survival_interpretation`, the approximation caveat for survival output.

When any follow-up labels are supplied, this metadata attribute also includes
`fu_labels`, the resolved map of known codes (including baseline if provided).
