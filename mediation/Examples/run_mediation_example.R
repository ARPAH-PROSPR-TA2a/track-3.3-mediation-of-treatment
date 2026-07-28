source(file.path("mediation", "main.R"))

input_dir <- file.path("mediation", "Examples", "ExampleInputs")
result_dir <- file.path("mediation", "Examples", "ExampleResults")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

treatment_results <- readRDS(
  file.path(input_dir, "proteomics_treatment_results.rds")
)
outcome_results <- readRDS(
  file.path(input_dir, "proteomics_outcome_results.rds")
)

mediation_results <- FAST_mediation(
  treatment_results = treatment_results,
  outcome_results = outcome_results,
  outcome_type = "continuous"
)

saveRDS(
  mediation_results,
  file.path(result_dir, "proteomics_mediation_results.rds")
)

effects <- mediation_results$analysis_change$all$mediation_effects
cat("Mediation proteomics example complete\n")
cat("  Matched analyte-FU rows: ", nrow(effects), "\n", sep = "")
cat("  Raw P < 0.05:            ", sum(effects$P_VALUE < 0.05), "\n", sep = "")
cat("  BH P < 0.05:             ", sum(effects$BH_P_VALUE < 0.05), "\n", sep = "")
