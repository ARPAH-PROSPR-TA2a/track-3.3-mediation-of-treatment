# FAST Track 3.3: Mediation of Treatment Effects

This repository contains three independent analysis pipelines:

- [OutcomeWAS](outcomewas/README.md): implemented omics-wide association
  analyses for subject-level outcomes.
- [Traditional Mediation](traditional_mediation/README.md): reserved for the
  traditional mediation pipeline; implementation has not started.
- [DACT Mediation](dact_mediation/README.md): reserved for the DACT mediation
  pipeline; implementation has not started.

The pipelines are separate user surfaces. They are not stages of a single
workflow, and one pipeline's output is not an implicit input to another.

## Current Entry Point

OutcomeWAS is the currently implemented pipeline. Run from the repository root:

```r
source(file.path("outcomewas", "main.R"))

results <- FAST_outcome_WAS(...)
reports <- FAST_outcome_WAS_reports(...)
```

The repository root does not expose a default analysis entrypoint. Source the
specific pipeline you intend to run.

## Tests

Run the complete implemented test suite from the repository root:

```bash
Rscript tests/run_tests.R
```

Pipeline-specific tests live under `tests/<pipeline>/`.
