source("main.R")


run_checks <- function(checks) {
  all_pass <- TRUE
  for (name in names(checks)) {
    status <- if (checks[[name]]) "PASS" else "FAIL"
    cat(status, " ", name, "\n", sep = "")
    if (!checks[[name]]) all_pass <- FALSE
  }
  all_pass
}


# =============================================================================
# OMICS TYPE VALIDATION
# =============================================================================

other_omics_output <- capture.output(
  other_omics_result <- .validate_omics_type("other"),
  type = "message"
)
invalid_omics_errors <- tryCatch(
  {
    .validate_omics_type("invalid")
    FALSE
  },
  error = function(e) TRUE
)

omics_type_validation_pass <- run_checks(list(
  "Other omics type passes validation" = is.null(other_omics_result),
  "Other omics type is silent" = length(other_omics_output) == 0L,
  "Invalid omics type errors" = invalid_omics_errors
))
cat("\n")


print_results_summary <- function(results, label = "") {
  if (label != "") cat(label, "\n", sep = "")
  if (!is.null(results$coefficients)) {
    coefs <- results$coefficients
    cat("  Coefficients rows: ", nrow(coefs), "\n")
    cat("  Analytes:          ", length(unique(coefs$ANALYTE_NAME)), "\n")
    cat("  FU levels:         ", paste(sort(unique(coefs$FU)), collapse = ", "), "\n")
    cat("  Columns:           ", paste(colnames(coefs), collapse = ", "), "\n")
  }
  if (!is.null(results$outcome_effects)) {
    oe <- results$outcome_effects
    cat("  Outcome effects:   ", nrow(oe), " rows\n")
    cat("  FU levels:         ", paste(sort(unique(oe$FU)), collapse = ", "), "\n")
  }
}


validate_outcome_effects <- function(results, expect_hazard_ratio = FALSE) {
  coefs <- results$coefficients
  oe <- results$outcome_effects

  analyte_coefs <- coefs[coefs$COEFFICIENT == "analyte", ]
  merged <- merge(
    oe,
    analyte_coefs,
    by = c("ANALYTE_NAME", "FU"),
    suffixes = c("_oe", "_coef")
  )

  checks <- list(
    "All outcome effects map to analyte coefficients" = nrow(merged) == nrow(oe),
    "Effect size matches analyte coefficient" = nrow(merged) == nrow(oe) &&
      all(abs(merged$EFFECT_SIZE_oe - merged$EFFECT_SIZE_coef) < 1e-10),
    "SE matches analyte coefficient" = nrow(merged) == nrow(oe) &&
      all(abs(merged$SE_oe - merged$SE_coef) < 1e-10),
    "P-value matches analyte coefficient" = nrow(merged) == nrow(oe) &&
      all(abs(merged$P_VALUE_oe - merged$P_VALUE_coef) < 1e-10)
  )

  if (expect_hazard_ratio) {
    checks[["Hazard ratio matches analyte coefficient"]] <- nrow(merged) == nrow(oe) &&
      all(abs(merged$HAZARD_RATIO_oe - merged$HAZARD_RATIO_coef) < 1e-10)
  }

  checks
}


check_top_level_structure <- function(results) {
  list(
    "Has analysis_change and analysis_level" = all(c("analysis_change", "analysis_level") %in% names(results)),
    "analysis_change has all/male/female" = all(c("all", "male", "female") %in% names(results$analysis_change)),
    "analysis_level has all/male/female" = all(c("all", "male", "female") %in% names(results$analysis_level))
  )
}


check_analysis_structure <- function(analysis, expected_fu_levels,
                                     expect_hazard_ratio = FALSE,
                                     expect_filtered_bh = FALSE) {
  coefs <- analysis$all$coefficients
  oe <- analysis$all$outcome_effects

  checks <- list(
    "Coefficients exist" = !is.null(coefs),
    "Coefficients have rows" = !is.null(coefs) && nrow(coefs) > 0,
    "Coefficients have FU" = !is.null(coefs) && "FU" %in% colnames(coefs),
    "Coefficients FU levels correct" = !is.null(coefs) &&
      identical(sort(unique(coefs$FU)), expected_fu_levels),
    "Coefficients have N_OBS" = !is.null(coefs) && "N_OBS" %in% colnames(coefs),
    "Coefficients have BH correction" = !is.null(coefs) && "BH_P_VALUE" %in% colnames(coefs),
    "Outcome effects exist" = !is.null(oe),
    "Outcome effects have rows" = !is.null(oe) && nrow(oe) > 0,
    "Outcome effects FU levels correct" = !is.null(oe) &&
      identical(sort(unique(oe$FU)), expected_fu_levels),
    "Outcome effects have BH correction" = !is.null(oe) && "BH_P_VALUE" %in% colnames(oe),
    "Sex stratification works" = !is.null(analysis$male$coefficients) && !is.null(analysis$female$coefficients)
  )

  if (expect_hazard_ratio) {
    checks[["Coefficients have HAZARD_RATIO"]] <- "HAZARD_RATIO" %in% colnames(coefs)
    checks[["Outcome effects have HAZARD_RATIO"]] <- "HAZARD_RATIO" %in% colnames(oe)
  }

  if (expect_filtered_bh) {
    checks[["Coefficients have BH_P_VALUE_FILTERED"]] <- "BH_P_VALUE_FILTERED" %in% colnames(coefs)
    checks[["Outcome effects have BH_P_VALUE_FILTERED"]] <- "BH_P_VALUE_FILTERED" %in% colnames(oe)
  }

  checks
}


check_reports_structure <- function(reports, outcome_type, expected_fu_levels,
                                    has_covariates = TRUE) {
  vs <- reports$variable_summaries
  or <- reports$outcome_reports

  all_keys <- names(vs$all)
  omics_keys <- all_keys[startsWith(all_keys, "omics_")]
  covar_keys <- all_keys[startsWith(all_keys, "covariates_")]

  expected_baseline_keys <- c("omics_FU0_Tx0", "omics_FU0_Tx1")
  expected_followup_keys <- as.vector(outer(
    paste0("omics_FU", expected_fu_levels),
    c("_Tx0", "_Tx1"),
    paste0
  ))

  checks <- list(
    "Has pheno_summary" = !is.null(reports$pheno_summary),
    "Pheno summary has rows" = !is.null(reports$pheno_summary) && nrow(reports$pheno_summary) > 0,
    "Has variable_summaries" = !is.null(vs),
    "variable_summaries has all stratum" = !is.null(vs$all),
    "Baseline omics entries present" = all(expected_baseline_keys %in% omics_keys),
    "Follow-up omics entries present" = all(expected_followup_keys %in% omics_keys),
    "Omics entries have rows" = length(omics_keys) > 0 && all(sapply(omics_keys, function(k) nrow(vs$all[[k]]) > 0)),
    "Covariate entries present" = !has_covariates || length(covar_keys) > 0,
    "Has outcome_reports" = !is.null(or),
    "Outcome type matches" = !is.null(or) && identical(or$outcome_type, outcome_type),
    "Has analysis_sample_summary" = !is.null(or) && !is.null(or$analysis_sample_summary),
    "analysis_sample_summary FU levels correct" = !is.null(or) &&
      identical(sort(or$analysis_sample_summary$FU), expected_fu_levels),
    "Has outcome_summary" = !is.null(or) && !is.null(or$outcome_summary) && nrow(or$outcome_summary) > 0
  )

  checks
}


create_subject_level_outcomes <- function(pheno_df, omics_df, preferred_fu = 1L) {
  fu_num <- as.integer(as.character(pheno_df$FU))
  pheno_baseline <- pheno_df[fu_num == 0, ]
  pheno_followup <- pheno_df[fu_num != 0, ]

  followup_priority <- ifelse(
    as.integer(as.character(pheno_followup$FU)) == preferred_fu,
    -1L,
    as.integer(as.character(pheno_followup$FU))
  )
  pheno_followup <- pheno_followup[order(pheno_followup$SUBJECT_ID, followup_priority), ]
  pheno_followup <- pheno_followup[!duplicated(pheno_followup$SUBJECT_ID), ]

  subjects <- intersect(pheno_baseline$SUBJECT_ID, pheno_followup$SUBJECT_ID)
  pheno_baseline <- pheno_baseline[match(subjects, pheno_baseline$SUBJECT_ID), ]
  pheno_followup <- pheno_followup[match(subjects, pheno_followup$SUBJECT_ID), ]

  signal_baseline <- as.numeric(omics_df[1, pheno_baseline$SAMPLE_ID])
  signal_followup <- as.numeric(omics_df[1, pheno_followup$SAMPLE_ID])
  signal_change <- signal_followup - signal_baseline

  treatment <- as.numeric(as.character(pheno_followup$TREATMENT_GROUP))
  bmi_scaled <- as.numeric(scale(pheno_followup$mbmi))
  signal_scaled <- as.numeric(scale(signal_change))

  set.seed(20260518)
  continuous_outcome <- 0.9 * signal_scaled + 0.4 * treatment + 0.25 * bmi_scaled + rnorm(length(signal_scaled), sd = 0.7)

  lp <- 0.7 * signal_scaled + 0.35 * treatment + 0.2 * bmi_scaled
  event_rate <- 0.08 * exp(lp)
  event_time <- rexp(length(lp), rate = pmax(event_rate, 1e-4))
  censor_time <- rexp(length(lp), rate = 0.06)
  outcome_status <- as.integer(event_time <= censor_time)
  outcome_time <- pmin(event_time, censor_time)

  data.frame(
    SUBJECT_ID = pheno_followup$SUBJECT_ID,
    OUTCOME = continuous_outcome,
    OUTCOME_TIME = outcome_time,
    OUTCOME_STATUS = outcome_status,
    stringsAsFactors = FALSE
  )
}


attach_subject_level_outcome <- function(pheno_df, subject_outcomes, type = c("continuous", "tte")) {
  type <- match.arg(type)
  out <- merge(pheno_df, subject_outcomes, by = "SUBJECT_ID", all.x = TRUE, sort = FALSE)
  out <- out[match(pheno_df$SAMPLE_ID, out$SAMPLE_ID), ]

  if (type == "continuous") {
    out$OUTCOME_TIME <- NULL
    out$OUTCOME_STATUS <- NULL
  } else {
    out$OUTCOME <- NULL
  }

  out
}


restrict_to_complete_subjects <- function(pheno_df) {
  fu_num <- as.integer(as.character(pheno_df$FU))
  subjects_with_baseline <- unique(pheno_df$SUBJECT_ID[fu_num == 0])
  subjects_with_followup <- unique(pheno_df$SUBJECT_ID[fu_num > 0])
  complete_subjects <- intersect(subjects_with_baseline, subjects_with_followup)
  pheno_df[pheno_df$SUBJECT_ID %in% complete_subjects, ]
}


cat("OutcomeWAS Comprehensive Test Suite\n")
cat("===================================\n\n")

pheno_raw <- readRDS("PracticeData/pheno_example.rds")
omics_raw <- readRDS("PracticeData/synth_small_betas.rds")

if (any(sapply(pheno_raw, function(x) inherits(x, "haven_labelled")))) {
  for (col in names(pheno_raw)) {
    if (inherits(pheno_raw[[col]], "haven_labelled")) {
      raw_values <- as.vector(pheno_raw[[col]])
      numeric_attempt <- suppressWarnings(as.numeric(raw_values))
      if (!all(is.na(numeric_attempt))) {
        pheno_raw[[col]] <- numeric_attempt
      } else {
        pheno_raw[[col]] <- as.character(raw_values)
      }
    }
  }
}

pheno_raw$SUBJECT_ID <- as.character(pheno_raw$SUBJECT_ID)
pheno_raw$FU <- factor(as.integer(pheno_raw$FU))
pheno_raw$FEMALE <- factor(as.integer(pheno_raw$FEMALE))
pheno_raw$TREATMENT_GROUP <- factor(1 - pheno_raw$CONTROL_STATUS)
pheno_raw$CONTROL_STATUS <- NULL

required_base_cols <- c("SAMPLE_ID", "FU", "SUBJECT_ID", "FEMALE", "TREATMENT_GROUP",
                        "agebl", "agevis", "ethnic", "race3", "mbmi")
pheno <- pheno_raw[complete.cases(pheno_raw[, required_base_cols]), ]
pheno <- pheno[!duplicated(paste(pheno$SUBJECT_ID, pheno$FU)), ]

analyte_names <- rownames(omics_raw)
sample_names <- colnames(omics_raw)
omics_full <- as.data.frame(omics_raw)
colnames(omics_full) <- sample_names
omics_full <- cbind(ANALYTE_NAME = analyte_names, omics_full, stringsAsFactors = FALSE)

shared_samples <- intersect(pheno$SAMPLE_ID, sample_names)
pheno <- pheno[pheno$SAMPLE_ID %in% shared_samples, ]
pheno <- restrict_to_complete_subjects(pheno)

pheno_single_fu <- pheno[pheno$FU %in% c(0, 1), ]
pheno_single_fu <- restrict_to_complete_subjects(pheno_single_fu)
pheno_multi_fu <- pheno

omics_full <- omics_full[, c("ANALYTE_NAME", shared_samples), drop = FALSE]

full_probes <- readRDS("Data/FAST_epicv1_epicv2_probe_list.rds")
filtered_probes <- readRDS("Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")

filtered_available <- omics_full[omics_full$ANALYTE_NAME %in% filtered_probes, ]
full_only_available <- omics_full[omics_full$ANALYTE_NAME %in% full_probes &
                                    !omics_full$ANALYTE_NAME %in% filtered_probes, ]

omics_dnam <- rbind(
  filtered_available[1:min(25, nrow(filtered_available)), ],
  full_only_available[1:min(75, nrow(full_only_available)), ]
)
omics_non_dnam <- omics_dnam[1:min(40, nrow(omics_dnam)), ]

subject_outcomes <- create_subject_level_outcomes(pheno_multi_fu, omics_non_dnam, preferred_fu = 1L)
pheno_single_continuous <- attach_subject_level_outcome(pheno_single_fu, subject_outcomes, type = "continuous")
pheno_single_tte <- attach_subject_level_outcome(pheno_single_fu, subject_outcomes, type = "tte")
pheno_multi_continuous <- attach_subject_level_outcome(pheno_multi_fu, subject_outcomes, type = "continuous")
pheno_multi_tte <- attach_subject_level_outcome(pheno_multi_fu, subject_outcomes, type = "tte")

single_fu_levels <- 1L
multi_fu_levels <- sort(unique(as.integer(as.character(pheno_multi_fu$FU))))
multi_fu_levels <- multi_fu_levels[multi_fu_levels != 0]

additional_covariates <- c("agebl", "agevis", "ethnic", "race3", "mbmi")

cat("Data prepared\n")
cat("  Single-FU subjects: ", length(unique(pheno_single_fu$SUBJECT_ID)), "\n")
cat("  Multi-FU subjects:  ", length(unique(pheno_multi_fu$SUBJECT_ID)), "\n")
cat("  Single-FU samples:  ", nrow(pheno_single_fu), "\n")
cat("  Multi-FU samples:   ", nrow(pheno_multi_fu), "\n")
cat("  Multi-FU levels:    ", paste(multi_fu_levels, collapse = ", "), "\n")
cat("  Non-DNAm analytes:  ", nrow(omics_non_dnam), "\n")
cat("  DNAm analytes:      ", nrow(omics_dnam), "\n")
cat("  TTE events:         ", sum(subject_outcomes$OUTCOME_STATUS), "\n\n")

# =============================================================================
# TEST 1: Single FU continuous + non-DNAm
# =============================================================================

cat("TEST 1: Single FU continuous + non-DNAm\n")
cat("=======================================\n\n")

results_test1 <- FAST_outcome_WAS(
  pheno = pheno_single_continuous,
  omics = omics_non_dnam,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1
)
reports_test1 <- FAST_outcome_WAS_reports(
  pheno = pheno_single_continuous,
  omics = omics_non_dnam,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

cat("Top-Level Structure Checks\n")
test1_top_pass <- run_checks(check_top_level_structure(results_test1))
cat("\n")

print_results_summary(results_test1$analysis_change$all, "Change Results")
print_results_summary(results_test1$analysis_level$all, "Level Results")
cat("\n")

cat("Change Structural Checks\n")
test1_change_struct_pass <- run_checks(check_analysis_structure(results_test1$analysis_change, expected_fu_levels = single_fu_levels))
cat("\n")

cat("Level Structural Checks\n")
test1_level_struct_pass <- run_checks(check_analysis_structure(results_test1$analysis_level, expected_fu_levels = single_fu_levels))
cat("\n")

cat("Reports Structural Checks\n")
test1_reports_pass <- run_checks(check_reports_structure(reports_test1, outcome_type = "continuous", expected_fu_levels = single_fu_levels))
cat("\n")

cat("Change Outcome Effect Validation\n")
test1_change_oe_pass <- run_checks(validate_outcome_effects(results_test1$analysis_change$all))
cat("\n")

cat("Level Outcome Effect Validation\n")
test1_level_oe_pass <- run_checks(validate_outcome_effects(results_test1$analysis_level$all))
cat("\n")

test1_pass <- test1_top_pass && test1_change_struct_pass && test1_level_struct_pass &&
  test1_reports_pass && test1_change_oe_pass && test1_level_oe_pass

# =============================================================================
# TEST 2: Single FU TTE + non-DNAm
# =============================================================================

cat("TEST 2: Single FU TTE + non-DNAm\n")
cat("================================\n\n")

results_test2 <- FAST_outcome_WAS(
  pheno = pheno_single_tte,
  omics = omics_non_dnam,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1
)
reports_test2 <- FAST_outcome_WAS_reports(
  pheno = pheno_single_tte,
  omics = omics_non_dnam,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

cat("Top-Level Structure Checks\n")
test2_top_pass <- run_checks(check_top_level_structure(results_test2))
cat("\n")

print_results_summary(results_test2$analysis_change$all, "Change Results")
print_results_summary(results_test2$analysis_level$all, "Level Results")
cat("\n")

cat("Change Structural Checks\n")
test2_change_struct_pass <- run_checks(check_analysis_structure(results_test2$analysis_change, expected_fu_levels = single_fu_levels, expect_hazard_ratio = TRUE))
cat("\n")

cat("Level Structural Checks\n")
test2_level_struct_pass <- run_checks(check_analysis_structure(results_test2$analysis_level, expected_fu_levels = single_fu_levels, expect_hazard_ratio = TRUE))
cat("\n")

cat("Reports Structural Checks\n")
test2_reports_pass <- run_checks(check_reports_structure(reports_test2, outcome_type = "tte", expected_fu_levels = single_fu_levels))
cat("\n")

cat("Change Outcome Effect Validation\n")
test2_change_oe_pass <- run_checks(validate_outcome_effects(results_test2$analysis_change$all, expect_hazard_ratio = TRUE))
cat("\n")

cat("Level Outcome Effect Validation\n")
test2_level_oe_pass <- run_checks(validate_outcome_effects(results_test2$analysis_level$all, expect_hazard_ratio = TRUE))
cat("\n")

test2_pass <- test2_top_pass && test2_change_struct_pass && test2_level_struct_pass &&
  test2_reports_pass && test2_change_oe_pass && test2_level_oe_pass

# =============================================================================
# TEST 3: Multi-FU continuous + non-DNAm
# =============================================================================

cat("TEST 3: Multi-FU continuous + non-DNAm\n")
cat("======================================\n\n")

results_test3 <- FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_non_dnam,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1
)
reports_test3 <- FAST_outcome_WAS_reports(
  pheno = pheno_multi_continuous,
  omics = omics_non_dnam,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

cat("Top-Level Structure Checks\n")
test3_top_pass <- run_checks(check_top_level_structure(results_test3))
cat("\n")

print_results_summary(results_test3$analysis_change$all, "Change Results")
print_results_summary(results_test3$analysis_level$all, "Level Results")
cat("\n")

cat("Change Structural Checks\n")
test3_change_struct_pass <- run_checks(check_analysis_structure(results_test3$analysis_change, expected_fu_levels = multi_fu_levels))
cat("\n")

cat("Level Structural Checks\n")
test3_level_struct_pass <- run_checks(check_analysis_structure(results_test3$analysis_level, expected_fu_levels = multi_fu_levels))
cat("\n")

cat("Reports Structural Checks\n")
test3_reports_pass <- run_checks(check_reports_structure(reports_test3, outcome_type = "continuous", expected_fu_levels = multi_fu_levels))
cat("\n")

cat("Change Outcome Effect Validation\n")
test3_change_oe_pass <- run_checks(validate_outcome_effects(results_test3$analysis_change$all))
cat("\n")

cat("Level Outcome Effect Validation\n")
test3_level_oe_pass <- run_checks(validate_outcome_effects(results_test3$analysis_level$all))
cat("\n")

test3_pass <- test3_top_pass && test3_change_struct_pass && test3_level_struct_pass &&
  test3_reports_pass && test3_change_oe_pass && test3_level_oe_pass

# =============================================================================
# TEST 4: Multi-FU TTE + non-DNAm
# =============================================================================

cat("TEST 4: Multi-FU TTE + non-DNAm\n")
cat("===============================\n\n")

results_test4 <- FAST_outcome_WAS(
  pheno = pheno_multi_tte,
  omics = omics_non_dnam,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1
)
reports_test4 <- FAST_outcome_WAS_reports(
  pheno = pheno_multi_tte,
  omics = omics_non_dnam,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates
)

cat("Top-Level Structure Checks\n")
test4_top_pass <- run_checks(check_top_level_structure(results_test4))
cat("\n")

print_results_summary(results_test4$analysis_change$all, "Change Results")
print_results_summary(results_test4$analysis_level$all, "Level Results")
cat("\n")

cat("Change Structural Checks\n")
test4_change_struct_pass <- run_checks(check_analysis_structure(results_test4$analysis_change, expected_fu_levels = multi_fu_levels, expect_hazard_ratio = TRUE))
cat("\n")

cat("Level Structural Checks\n")
test4_level_struct_pass <- run_checks(check_analysis_structure(results_test4$analysis_level, expected_fu_levels = multi_fu_levels, expect_hazard_ratio = TRUE))
cat("\n")

cat("Reports Structural Checks\n")
test4_reports_pass <- run_checks(check_reports_structure(reports_test4, outcome_type = "tte", expected_fu_levels = multi_fu_levels))
cat("\n")

cat("Change Outcome Effect Validation\n")
test4_change_oe_pass <- run_checks(validate_outcome_effects(results_test4$analysis_change$all, expect_hazard_ratio = TRUE))
cat("\n")

cat("Level Outcome Effect Validation\n")
test4_level_oe_pass <- run_checks(validate_outcome_effects(results_test4$analysis_level$all, expect_hazard_ratio = TRUE))
cat("\n")

test4_pass <- test4_top_pass && test4_change_struct_pass && test4_level_struct_pass &&
  test4_reports_pass && test4_change_oe_pass && test4_level_oe_pass

# =============================================================================
# TEST 5: Multi-FU continuous + DNAm
# =============================================================================

cat("TEST 5: Multi-FU continuous + DNAm\n")
cat("==================================\n\n")

results_test5 <- FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_dnam,
  omics_type = "DNAm",
  additional_covariates = additional_covariates,
  n_cores = 1
)
reports_test5 <- FAST_outcome_WAS_reports(
  pheno = pheno_multi_continuous,
  omics = omics_dnam,
  omics_type = "DNAm",
  additional_covariates = additional_covariates
)

cat("Top-Level Structure Checks\n")
test5_top_pass <- run_checks(check_top_level_structure(results_test5))
cat("\n")

cat("Change Structural Checks\n")
test5_change_struct_pass <- run_checks(check_analysis_structure(results_test5$analysis_change, expected_fu_levels = multi_fu_levels, expect_filtered_bh = TRUE))
cat("\n")

cat("Level Structural Checks\n")
test5_level_struct_pass <- run_checks(check_analysis_structure(results_test5$analysis_level, expected_fu_levels = multi_fu_levels, expect_filtered_bh = TRUE))
cat("\n")

cat("Reports Structural Checks\n")
test5_reports_pass <- run_checks(check_reports_structure(reports_test5, outcome_type = "continuous", expected_fu_levels = multi_fu_levels))
cat("\n")

cat("Change Outcome Effect Validation\n")
test5_change_oe_pass <- run_checks(validate_outcome_effects(results_test5$analysis_change$all))
cat("\n")

cat("Level Outcome Effect Validation\n")
test5_level_oe_pass <- run_checks(validate_outcome_effects(results_test5$analysis_level$all))
cat("\n")

test5_pass <- test5_top_pass && test5_change_struct_pass && test5_level_struct_pass &&
  test5_reports_pass && test5_change_oe_pass && test5_level_oe_pass

# =============================================================================
# TEST 6: Multi-FU TTE + DNAm
# =============================================================================

cat("TEST 6: Multi-FU TTE + DNAm\n")
cat("===========================\n\n")

results_test6 <- FAST_outcome_WAS(
  pheno = pheno_multi_tte,
  omics = omics_dnam,
  omics_type = "DNAm",
  additional_covariates = additional_covariates,
  n_cores = 1
)
reports_test6 <- FAST_outcome_WAS_reports(
  pheno = pheno_multi_tte,
  omics = omics_dnam,
  omics_type = "DNAm",
  additional_covariates = additional_covariates
)

cat("Top-Level Structure Checks\n")
test6_top_pass <- run_checks(check_top_level_structure(results_test6))
cat("\n")

cat("Change Structural Checks\n")
test6_change_struct_pass <- run_checks(check_analysis_structure(results_test6$analysis_change, expected_fu_levels = multi_fu_levels, expect_hazard_ratio = TRUE, expect_filtered_bh = TRUE))
cat("\n")

cat("Level Structural Checks\n")
test6_level_struct_pass <- run_checks(check_analysis_structure(results_test6$analysis_level, expected_fu_levels = multi_fu_levels, expect_hazard_ratio = TRUE, expect_filtered_bh = TRUE))
cat("\n")

cat("Reports Structural Checks\n")
test6_reports_pass <- run_checks(check_reports_structure(reports_test6, outcome_type = "tte", expected_fu_levels = multi_fu_levels))
cat("\n")

cat("Change Outcome Effect Validation\n")
test6_change_oe_pass <- run_checks(validate_outcome_effects(results_test6$analysis_change$all, expect_hazard_ratio = TRUE))
cat("\n")

cat("Level Outcome Effect Validation\n")
test6_level_oe_pass <- run_checks(validate_outcome_effects(results_test6$analysis_level$all, expect_hazard_ratio = TRUE))
cat("\n")

test6_pass <- test6_top_pass && test6_change_struct_pass && test6_level_struct_pass &&
  test6_reports_pass && test6_change_oe_pass && test6_level_oe_pass

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n")
cat("FINAL SUMMARY\n")
cat("=============\n")
cat("Omics type validation:                      ", if (omics_type_validation_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Test 1 (Single-FU continuous + non-DNAm): ", if (test1_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Test 2 (Single-FU TTE + non-DNAm):        ", if (test2_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Test 3 (Multi-FU continuous + non-DNAm):  ", if (test3_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Test 4 (Multi-FU TTE + non-DNAm):         ", if (test4_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Test 5 (Multi-FU continuous + DNAm):      ", if (test5_pass) "PASS" else "FAIL", "\n", sep = "")
cat("Test 6 (Multi-FU TTE + DNAm):             ", if (test6_pass) "PASS" else "FAIL", "\n", sep = "")

all_tests_pass <- all(omics_type_validation_pass, test1_pass, test2_pass,
                      test3_pass, test4_pass, test5_pass, test6_pass)
cat("\nOverall: ", if (all_tests_pass) "ALL TESTS PASSED" else "SOME TESTS FAILED", "\n", sep = "")

if (!all_tests_pass) {
  stop("One or more OutcomeWAS comprehensive tests failed.")
}
