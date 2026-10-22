# Run from the repository root:
#   Rscript outcomewas/Examples/run_multiarm_example.R [output_directory]
# All data are simulated. Without an output directory, files use R's tempdir().
source(file.path("outcomewas", "main.R"))

set.seed(3303)
n_subjects <- 96L
n_analytes <- 4L
subject_id <- sprintf("synthetic_%03d", seq_len(n_subjects))
arm <- rep(c(0L, 1L, 3L), each = n_subjects / 3L)
female <- rep(c(0L, 1L), length.out = n_subjects)
age <- rnorm(n_subjects, mean = 55, sd = 8)
baseline <- matrix(rnorm(n_analytes * n_subjects), nrow = n_analytes)
change <- matrix(rnorm(n_analytes * n_subjects), nrow = n_analytes)
change[1L, ] <- change[1L, ] + ifelse(arm == 1L, 0.5, ifelse(arm == 3L, 0.9, 0))
followup <- baseline + change

# Interventions shift the marker differently; its outcome slope is common to all arms.
# Each subject's outcome is repeated across their two sample rows.
outcome <- 0.8 * change[1L, ] + 0.3 * baseline[1L, ] +
  0.4 * (arm == 1L) - 0.2 * (arm == 3L) + 0.02 * age +
  rnorm(n_subjects)
pheno <- data.frame(
  SAMPLE_ID = c(paste0(subject_id, "_FU0"), paste0(subject_id, "_FU1")),
  SUBJECT_ID = rep(subject_id, times = 2L),
  FU = factor(rep(c(0L, 1L), each = n_subjects)),
  TREATMENT_GROUP = rep(arm, times = 2L),
  FEMALE = factor(rep(female, times = 2L)),
  OUTCOME = rep(outcome, times = 2L),
  age = rep(age, times = 2L),
  stringsAsFactors = FALSE
)
omics <- data.frame(
  ANALYTE_NAME = paste0("synthetic_marker_", seq_len(n_analytes)),
  cbind(baseline, followup), check.names = FALSE
)
names(omics)[-1L] <- pheno$SAMPLE_ID

# Nonconsecutive active-arm codes are preserved as categories.
treatment_labels <- c(`0` = "Control", `1` = "Diet", `3` = "Exercise")
results <- FAST_outcome_WAS(
  pheno, omics, omics_type = "other", additional_covariates = "age",
  n_cores = 1, treatment_labels = treatment_labels
)
reports <- FAST_outcome_WAS_reports(
  pheno, omics, omics_type = "other", additional_covariates = "age",
  treatment_labels = treatment_labels
)

# One common analyte slope per analyte/FU, adjusted for every treatment arm.
print(results$analysis_change$all$outcome_effects)
print(results$metadata)
print(reports$arm_counts)
print(reports$outcome_reports$outcome_summary)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1L]] else file.path(tempdir(), "outcomewas_multiarm")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, "synthetic_multiarm_results.rds"))
saveRDS(reports, file.path(output_dir, "synthetic_multiarm_reports.rds"))
cat("Synthetic three-arm OutcomeWAS example complete\n")
cat("  Outputs: ", normalizePath(output_dir), "\n", sep = "")
if (!length(args)) {
  cat("  Supply output_directory to retain files after this R session ends.\n")
}
