# From the repository root:
# Rscript mediation/Examples/run_multiarm_example.R [output_directory]
# All data are simulated. Fits a joint alpha model, common-slope OutcomeWAS,
# and contrast-aware mediation without research data or a sibling checkout.
# An explicit directory retains artifacts after this R session.
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

outcome_results <- results
source(file.path("mediation", "main.R"))

# Demonstrate the Track 1.1.1 single-follow-up alpha model, jointly over all arms.
# In a real analysis, read the saved Track 1.1.1 treatment_results instead.
fit_alpha <- function(response, stratum) {
  keep <- if (stratum == "all") rep(TRUE, n_subjects) else {
    female == if (stratum == "male") 0L else 1L
  }
  rows <- lapply(seq_len(n_analytes), function(j) {
    md <- data.frame(
      M = if (response == "change") change[j, keep] else followup[j, keep],
      M_baseline = baseline[j, keep], treatment = factor(arm[keep]),
      female = female[keep], age = age[keep]
    )
    formula <- if (stratum == "all") {
      M ~ treatment + M_baseline + female + age
    } else {
      M ~ treatment + M_baseline + age
    }
    coefficients <- coef(summary(lm(formula, md)))
    active <- c("1", "3")
    terms <- paste0("treatment", active)
    data.frame(
      ANALYTE_NAME = omics$ANALYTE_NAME[j], FU = 1L,
      CONTRAST = paste0(active, "_vs_0"),
      TREATMENT = unname(treatment_labels[active]),
      REFERENCE = unname(treatment_labels["0"]),
      EFFECT_SIZE = coefficients[terms, "Estimate"],
      SE = coefficients[terms, "Std. Error"],
      P_VALUE = coefficients[terms, "Pr(>|t|)"], STATUS = "ok",
      row.names = NULL
    )
  })
  effects <- do.call(rbind, rows)
  effects$BH_P_VALUE <- NA_real_
  for (contrast in unique(effects$CONTRAST)) {
    index <- effects$CONTRAST == contrast
    effects$BH_P_VALUE[index] <- p.adjust(effects$P_VALUE[index], "BH")
  }
  list(treatment_effects = effects)
}
treatment_results <- setNames(lapply(c("change", "level"), function(response) {
  setNames(lapply(c("all", "male", "female"), function(stratum) {
    fit_alpha(response, stratum)
  }), c("all", "male", "female"))
}), c("analysis_change", "analysis_level"))
treatment_results$metadata <- list(treatment_labels = treatment_labels)

mediation_results <- FAST_mediation(treatment_results, outcome_results)
effects <- mediation_results$analysis_change$all$mediation_effects
stopifnot(nrow(effects) == n_analytes * 2L,
          setequal(effects$CONTRAST, c("1_vs_0", "3_vs_0")),
          all(effects$STATUS == "ok"))
# Each intervention has its own alpha, paired with the same analyte/FU beta.
stopifnot(all(vapply(split(effects$BETA_EFFECT, effects$ANALYTE_NAME),
                     function(x) length(unique(x)) == 1L, logical(1))))
print(effects[c("ANALYTE_NAME", "CONTRAST", "TREATMENT", "REFERENCE",
                "ALPHA_EFFECT", "BETA_EFFECT", "INDIRECT_EFFECT", "STATUS")])
saveRDS(treatment_results, file.path(output_dir, "synthetic_multiarm_treatment_results.rds"))
saveRDS(mediation_results, file.path(output_dir, "synthetic_multiarm_mediation_results.rds"))
saveRDS(list(
  synthetic = TRUE, seed = 3303L, treatment_labels = treatment_labels,
  outcome_model = outcome_results$metadata,
  mediation = attr(mediation_results, "metadata"), session_info = sessionInfo()
), file.path(output_dir, "synthetic_multiarm_provenance.rds"))
cat("Synthetic joint multi-arm mediation complete; outputs: ", output_dir, "\n", sep = "")
