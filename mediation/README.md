# Mediation

This product-of-coefficients pipeline combines treatment-to-analyte estimates
from FAST Track 1.1.1 with analyte-to-outcome estimates from OutcomeWAS. It
returns a separate indirect effect for each analyte, follow-up, and treatment
contrast. Inference retains the conventional Sobel delta-method standard error
and normal approximation.

Ordinary binary calls remain unchanged:

```r
source(file.path("mediation", "main.R"))

mediation_results <- FAST_mediation(
  treatment_results = treatment_results,
  outcome_results = outcome_results,
  outcome_type = "continuous"
)
```

Both inputs are result objects from their respective pipelines, with matching
analysis types and compatible study populations, analyte scales, and follow-up
definitions. Sex-specific strata are optional; an object containing only `all`
is valid. Complete legacy binary tables remain supported.

For multiple arms, Track 1.1.1 supplies one alpha per `CONTRAST`, such as
`1_vs_0` and `2_vs_0`. OutcomeWAS supplies one common beta per analyte and
follow-up, adjusted for all categorical treatment arms in the same model.
Mediation pairs each contrast with that beta: `alpha_1 * beta`,
`alpha_2 * beta`, and so on. The treatment coefficients in the outcome model
are adjustment terms, not additional beta paths.

Different interventions can change an analyte by different amounts. The model
assumes that a unit change in that analyte has the same outcome association
across arms. This extends the ordinary binary workflow without changing its
outcome model or adding treatment interactions.

For continuous outcomes, pure and total natural indirect effects coincide
under this common-slope model. A causal interpretation still requires timing
and mediation identification assumptions. For survival outcomes, the
exponentiated product is retained as an approximation; it is not generally a
causal indirect hazard ratio.

Outputs preserve `CONTRAST`, readable treatment/reference labels, and
availability statuses. Unavailable paths remain visible with `NA` indirect
inference. BH correction is within analysis, stratum, FU, and contrast; filtered
DNAm correction uses explicit `IN_FILTERED_SET` membership.

Follow-up labels are inherited from upstream `metadata$fu_labels` and/or
effect-table `FU_LABEL` columns. Labeled output adds `FU_LABEL` beside the
numeric `FU` and stores the resolved map in `attr(results, "metadata")$fu_labels`.
No additional argument is required. If neither input supplies labels, mediation
continues without a label warning or added label fields. If only one input
supplies labels, mediation inherits them and warns once that agreement could
not be checked. When both supply labels, shared codes must agree exactly,
including baseline metadata; conflicts are errors. Numeric FU codes remain the
join keys, and partial legacy inputs retain unknown labels as `NA`.
See [the label rules](INPUTS_OUTPUTS.md#follow-up-labels).

The pipeline uses summary results only and does not source or refit either
upstream pipeline. Sobel inference ignores alpha–beta covariance, and the
summary tables do not supply covariance for comparing interventions.

Run `Rscript tests/mediation/test_comprehensive.R` for the primary scenarios,
including visible follow-up labels, and `Rscript tests/mediation/test_multiarm.R`
for contrast and inference checks.
`Rscript tests/mediation/test_upstream_multiarm.R` additionally fits the synthetic
data through the local Track 1.1.1 checkout, then OutcomeWAS and mediation; set
`TRACK111_DIR` if that checkout is not the sibling directory. All suites are
registered in `Rscript tests/run_tests.R`.

The statistical interpretation follows the no-treatment-interaction case in
[Valeri and VanderWeele's mediation framework](https://pubmed.ncbi.nlm.nih.gov/23379553/). The additional
limitations of Cox products are discussed in [VanderWeele's survival mediation
paper](https://pubmed.ncbi.nlm.nih.gov/21642779/).

Documentation:

- [INPUTS_OUTPUTS.md](INPUTS_OUTPUTS.md): input/output contract and interpretation.
- [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md): validation, joins, and calculations.
- [Examples/README.md](Examples/README.md): paired binary and synthetic multi-arm
  examples.
- [CALERIE DNAm beta / INF runner](CALERIE/DNAm/INF_24mo/README.md): production
  inputs, checks, resource guidance, and outputs.
- [CALERIE DNAm beta / MetS runner](CALERIE/DNAm/MetS_24mo/README.md): the same
  workflow for the 24-month MetS outcome.
