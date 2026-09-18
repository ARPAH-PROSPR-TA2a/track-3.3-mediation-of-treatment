# Mediation

This product-of-coefficients mediation pipeline combines the
treatment-to-analyte estimates from FAST
Track 1.1.1 with the analyte-to-outcome estimates from OutcomeWAS. The indirect
effect for each analyte and follow-up is the product of coefficients
`alpha * beta`; inference uses the conventional Sobel delta-method standard
error and normal approximation.

```r
source(file.path("mediation", "main.R"))

mediation_results <- FAST_mediation(
  treatment_results = treatment_results,
  outcome_results = outcome_results,
  outcome_type = "continuous"
)
```

Both inputs must be complete results objects from their respective pipelines.
The function aligns results within each analysis type and sex stratum using
`ANALYTE_NAME` and `FU`. Sex-specific strata are optional; a dataset containing
only `all` is valid. See `INPUTS_OUTPUTS.md` for the complete contract,
statistical formulas, and interpretation boundaries.

The pipeline is self-contained and does not source either upstream pipeline.
The example input-generation test separately verifies live interoperability
with Track 1.1.1 and OutcomeWAS when a Track 1.1.1 checkout is available.

Documentation:

- [INPUTS_OUTPUTS.md](INPUTS_OUTPUTS.md): complete user-facing data contract.
- [CODE_WALKTHROUGH.md](CODE_WALKTHROUGH.md): implementation flow, calculations,
  validation, and output structure.
- [Examples/README.md](Examples/README.md): reproducible paired example.
- [CALERIE DNAm beta / INF runner](CALERIE/DNAm/INF_24mo/README.md): production
  input paths, preflight checks, resource guidance, launch command and outputs.
