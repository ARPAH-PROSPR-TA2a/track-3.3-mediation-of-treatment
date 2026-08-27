repo <- "~/FAST/GitHub/track-3.3"
setwd(repo)

source(file.path(repo, "mediation", "main.R"))

input_dir <- "~/FAST/Outputs/"
result_dir <- "~/FAST/Outputs/3.3/Proteomics_3.3Med_MetS_INF"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

treatment_results <- readRDS(
  file.path(input_dir, "1.1.1", "Proteomics_1.1.1_Output", "results_full_Proteomics_annotated.rds")
)

outcome_files <- list(
  MetS_24mo = file.path(input_dir, "3.3", "Proteomics_3.3OWAS_MetS", "Proteomics_3.3OWAS_MetS_results.rds"),
  INF_24mo = file.path(input_dir, "3.3", "Proteomics_3.3OWAS_INF", "Proteomics_3.3OWAS_INF_results.rds")
)

mediation_results_list <- lapply(names(outcome_files), function(outcome_name) {
  outcome_results <- readRDS(outcome_files[[outcome_name]])
  
  mediation_results <- FAST_mediation(
    treatment_results = treatment_results,
    outcome_results = outcome_results,
    outcome_type = "continuous"
  )
  
  saveRDS(
    mediation_results,
    file.path(result_dir, paste0("mediation_results_", outcome_name, ".rds"))
  )
  
  effects <- mediation_results$analysis_change$all$mediation_effects
  
  cat("\n", outcome_name, "\n", sep = "")
  cat("  Matched analyte-FU rows: ", nrow(effects), "\n", sep = "")
  cat("  Raw P < 0.05:            ", sum(effects$P_VALUE < 0.05), "\n", sep = "")
  cat("  BH P < 0.05:             ", sum(effects$BH_P_VALUE < 0.05), "\n", sep = "")
  
  mediation_results
})

names(mediation_results_list) <- names(outcome_files)