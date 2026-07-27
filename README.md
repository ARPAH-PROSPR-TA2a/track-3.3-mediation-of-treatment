# FAST Track 3.3: Mediation of Treatment Effects

This repository contains three independent analysis pipelines:

- [OutcomeWAS](outcomewas/README.md): implemented omics-wide association
  analyses for subject-level outcomes.
- [Mediation](mediation/README.md): product-of-
  coefficients mediation using Track 1.1.1 treatment effects and OutcomeWAS
  outcome effects, with Sobel inference.
- [DACT Mediation](dact_mediation/README.md): reserved for the DACT mediation
  pipeline; implementation has not started.

The pipelines are separate user surfaces and do not source or run one another
implicitly. Mediation explicitly accepts result objects produced
by Track 1.1.1 and OutcomeWAS.

## Current Entry Point

Source the entry point for the pipeline you intend to run. For OutcomeWAS:

```r
source(file.path("outcomewas", "main.R"))

results <- FAST_outcome_WAS(...)
reports <- FAST_outcome_WAS_reports(...)
```

For mediation:

```r
source(file.path("mediation", "main.R"))

mediation_results <- FAST_mediation(
  treatment_results,
  outcome_results,
  outcome_type = "continuous"
)
```

The repository root does not expose a default analysis entrypoint. Source the
specific pipeline you intend to run.

## Tests

Run the complete implemented test suite from the repository root:

```bash
Rscript tests/run_tests.R
```

Pipeline-specific tests live under `tests/<pipeline>/`.
