source(file.path("outcomewas", "main.R"))
source(file.path("outcomewas", "R", "plotting_helpers.R"))

pheno <- readRDS("outcomewas/Examples/ExampleData/pheno_example.rds")
omics <- readRDS("outcomewas/Examples/ExampleData/metabolomics_log2.rds")
dir.create("outcomewas/Examples/ExampleResults", showWarnings = FALSE)
figures_dir <- "outcomewas/Examples/ExampleFigures/Metabolomics"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

additional_covariates <- c("agebl", "agevis", "mbmi")

results <- FAST_outcome_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = "Metabolomics",
  additional_covariates = additional_covariates,
  n_cores = 1
)

reports <- FAST_outcome_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = "Metabolomics",
  additional_covariates = additional_covariates
)

saveRDS(results, "outcomewas/Examples/ExampleResults/metabolomics_results.rds")
saveRDS(reports, "outcomewas/Examples/ExampleResults/metabolomics_reports.rds")

generate_all_plots(results$analysis_change, figures_dir = file.path(figures_dir, "analysis_change"))
generate_all_plots(results$analysis_level, figures_dir = file.path(figures_dir, "analysis_level"))

cat("Metabolomics OutcomeWAS example complete\n")
cat("  Pheno samples: ", nrow(pheno), "\n", sep = "")
cat("  Metabolomic analytes: ", nrow(omics), "\n", sep = "")
cat("  Change outcome effects: ", nrow(results$analysis_change$all$outcome_effects), "\n", sep = "")
cat("  Level outcome effects: ", nrow(results$analysis_level$all$outcome_effects), "\n", sep = "")
cat("  Report summary rows: ", nrow(reports$pheno_summary), "\n", sep = "")
cat("  Figures: ", figures_dir, "\n", sep = "")
