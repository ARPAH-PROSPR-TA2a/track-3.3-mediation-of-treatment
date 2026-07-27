args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(getwd(), mustWork = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else {
  file.path(repo_root, "mediation", "Examples", "ExampleInputs")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

default_track111 <- file.path(dirname(repo_root), "Track1.1.1")
track111_dir <- Sys.getenv("TRACK111_DIR", unset = default_track111)
if (!dir.exists(track111_dir)) {
  stop(
    "Track 1.1.1 checkout not found at '", track111_dir, "'. ",
    "Set TRACK111_DIR to its repository path."
  )
}
track111_dir <- normalizePath(track111_dir, mustWork = TRUE)

data_dir <- file.path(repo_root, "mediation", "Examples", "ExampleData")
pheno <- readRDS(file.path(data_dir, "pheno_example.rds"))
omics <- readRDS(file.path(data_dir, "proteomics_log2.rds"))
additional_covariates <- c("agebl", "agevis", "mbmi")

old_wd <- setwd(track111_dir)
source("main.R")
treatment_results <- suppressMessages(
  FAST_omics_WAS(
    pheno = pheno,
    omics = omics,
    omics_type = "Proteomics",
    additional_covariates = additional_covariates,
    n_cores = 1
  )
)
setwd(repo_root)
saveRDS(
  treatment_results,
  file.path(output_dir, "proteomics_treatment_results.rds")
)

source(file.path("outcomewas", "main.R"))
outcome_results <- suppressMessages(
  FAST_outcome_WAS(
    pheno = pheno,
    omics = omics,
    omics_type = "Proteomics",
    additional_covariates = additional_covariates,
    n_cores = 1
  )
)
saveRDS(
  outcome_results,
  file.path(output_dir, "proteomics_outcome_results.rds")
)

cat("Mediation upstream inputs generated\n")
cat("  Track 1.1.1: ", track111_dir, "\n", sep = "")
cat("  Output:      ", output_dir, "\n", sep = "")
