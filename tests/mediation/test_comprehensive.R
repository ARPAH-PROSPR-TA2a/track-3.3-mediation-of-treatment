source(file.path("mediation", "main.R"))


expect_true <- function(value, label) {
  if (!isTRUE(value)) stop("FAIL: ", label)
  cat("PASS ", label, "\n", sep = "")
}


expect_equal <- function(actual, expected, label, tolerance = 1e-10) {
  equal <- isTRUE(all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE))
  if (!equal) {
    stop("FAIL: ", label, "\nExpected: ", paste(expected, collapse = ", "),
         "\nActual: ", paste(actual, collapse = ", "))
  }
  cat("PASS ", label, "\n", sep = "")
}


expect_error <- function(expression, pattern, label) {
  message <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  expect_true(!is.null(message) && grepl(pattern, message), label)
}


expect_warning <- function(expression, pattern, label) {
  messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(warning) {
      messages <<- c(messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl(pattern, messages)), label)
  value
}


make_effects <- function(outcome_type = "continuous", filtered = FALSE) {
  analyte <- c("A", "B", "C", "D")
  fu <- c(1L, 1L, 2L, 2L)
  alpha <- data.frame(
    ANALYTE_NAME = analyte,
    FU = fu,
    EFFECT_SIZE = c(0.5, -0.2, 0.3, 0),
    SE = c(0.1, 0.1, 0.15, 0.2),
    stringsAsFactors = FALSE
  )
  beta <- data.frame(
    ANALYTE_NAME = analyte,
    FU = fu,
    EFFECT_SIZE = c(0.4, 0.1, -0.25, 0),
    SE = c(0.2, 0.05, 0.1, 0.3),
    stringsAsFactors = FALSE
  )
  if (outcome_type == "tte") beta$HAZARD_RATIO <- exp(beta$EFFECT_SIZE)

  if (filtered) {
    alpha$BH_P_VALUE_FILTERED <- c(0.001, NA, 0.03, NA)
    beta$BH_P_VALUE_FILTERED <- c(0.04, NA, 0.02, NA)
  }

  list(
    treatment = list(treatment_effects = alpha),
    outcome = list(outcome_effects = beta)
  )
}


make_results <- function(entry, sex_strata = TRUE) {
  branch <- list(
    all = entry,
    male = if (sex_strata) entry else NULL,
    female = if (sex_strata) entry else NULL
  )
  list(analysis_change = branch, analysis_level = branch)
}


make_all_only_results <- function(entry) {
  branch <- list(all = entry)
  list(analysis_change = branch, analysis_level = branch)
}


cat("Mediation comprehensive tests\n\n")

fixture <- make_effects()
treatment_results <- make_results(fixture$treatment)
outcome_results <- make_results(fixture$outcome)
results <- FAST_mediation(treatment_results, outcome_results, "continuous")
effects <- results$analysis_change$all$mediation_effects
row_a <- effects[effects$ANALYTE_NAME == "A" & effects$FU == 1L, ]

expected_effect <- 0.5 * 0.4
expected_se <- sqrt(0.4^2 * 0.1^2 + 0.5^2 * 0.2^2)
expected_z <- expected_effect / expected_se
expected_p <- 2 * pnorm(-abs(expected_z))

expect_equal(row_a$INDIRECT_EFFECT, expected_effect, "Product-of-coefficients estimate")
expect_equal(row_a$INDIRECT_SE, expected_se, "Sobel standard error")
expect_equal(row_a$Z_VALUE, expected_z, "Sobel z statistic")
expect_equal(row_a$P_VALUE, expected_p, "Sobel p-value")
expect_equal(row_a$CI_LOWER, expected_effect - qnorm(0.975) * expected_se,
             "Normal-theory lower confidence limit")
expect_equal(row_a$CI_UPPER, expected_effect + qnorm(0.975) * expected_se,
             "Normal-theory upper confidence limit")
expect_true(!any(grepl("^DIRECT_", names(effects))),
            "Output is limited to indirect effects")
expect_true(
  !any(c(
    "ALPHA_P_VALUE", "ALPHA_BH_P_VALUE",
    "BETA_P_VALUE", "BETA_BH_P_VALUE"
  ) %in% names(effects)),
  "Unused upstream p-values are not copied"
)
expect_true(identical(attr(results, "method"), "product_of_coefficients"),
            "Method metadata")
expect_true(identical(attr(results, "inference"), "sobel"), "Inference metadata")
expect_true(results$analysis_change$all$join_diagnostics$matched_rows == 4L,
            "Complete join diagnostics")

all_only_results <- FAST_mediation(
  make_all_only_results(fixture$treatment),
  make_all_only_results(fixture$outcome),
  "continuous"
)
expect_true(
  identical(names(all_only_results$analysis_change), "all"),
  "All-only inputs do not require sex strata"
)

fu1 <- effects$FU == 1L
expect_equal(
  effects$BH_P_VALUE[fu1],
  p.adjust(effects$P_VALUE[fu1], method = "BH"),
  "BH correction is applied within FU"
)

row_d <- effects[effects$ANALYTE_NAME == "D", ]
expect_equal(row_d$INDIRECT_EFFECT, 0, "Zero-path indirect effect")
expect_equal(row_d$Z_VALUE, 0, "Singular zero-path z statistic")
expect_equal(row_d$P_VALUE, 1, "Singular zero-path p-value")

filtered_fixture <- make_effects(filtered = TRUE)
filtered_results <- FAST_mediation(
  make_results(filtered_fixture$treatment, sex_strata = FALSE),
  make_results(filtered_fixture$outcome, sex_strata = FALSE),
  "continuous"
)
filtered_effects <- filtered_results$analysis_change$all$mediation_effects
expect_true("BH_P_VALUE_FILTERED" %in% names(filtered_effects), "Filtered DNAm BH column")
expect_true(all(!is.na(filtered_effects$BH_P_VALUE_FILTERED[filtered_effects$ANALYTE_NAME %in% c("A", "C")])),
            "Filtered probes receive mediation BH values")
expect_true(all(is.na(filtered_effects$BH_P_VALUE_FILTERED[filtered_effects$ANALYTE_NAME %in% c("B", "D")])),
            "Non-filtered probes retain NA")

tte_fixture <- make_effects(outcome_type = "tte")
tte_results <- FAST_mediation(
  make_results(tte_fixture$treatment, sex_strata = FALSE),
  make_results(tte_fixture$outcome, sex_strata = FALSE),
  "tte"
)
tte_effects <- tte_results$analysis_change$all$mediation_effects
expect_equal(tte_effects$INDIRECT_HAZARD_RATIO, exp(tte_effects$INDIRECT_EFFECT),
             "Time-to-event mediated path hazard ratio")

partial_outcome <- make_results(fixture$outcome, sex_strata = FALSE)
partial_outcome$analysis_change$all$outcome_effects <-
  partial_outcome$analysis_change$all$outcome_effects[-1, ]
partial <- expect_warning(
  FAST_mediation(
    make_results(fixture$treatment, sex_strata = FALSE),
    partial_outcome,
    "continuous"
  ),
  "Partial mediation join",
  "Partial joins emit diagnostics"
)
expect_true(partial$analysis_change$all$join_diagnostics$treatment_only_rows == 1L,
            "Partial join counts treatment-only rows")

duplicate_treatment <- make_results(fixture$treatment, sex_strata = FALSE)
duplicate_treatment$analysis_change$all$treatment_effects <- rbind(
  duplicate_treatment$analysis_change$all$treatment_effects,
  duplicate_treatment$analysis_change$all$treatment_effects[1, ]
)
expect_error(
  FAST_mediation(duplicate_treatment, make_results(fixture$outcome, FALSE), "continuous"),
  "duplicate ANALYTE_NAME x FU",
  "Duplicate keys fail validation"
)

invalid_se <- make_results(fixture$treatment, sex_strata = FALSE)
invalid_se$analysis_change$all$treatment_effects$SE[1] <- 0
expect_error(
  FAST_mediation(invalid_se, make_results(fixture$outcome, FALSE), "continuous"),
  "SE must be strictly positive",
  "Nonpositive standard errors fail validation"
)

filtered_mismatch <- make_effects(filtered = TRUE)
filtered_mismatch$outcome$outcome_effects$BH_P_VALUE_FILTERED[1] <- NA
expect_error(
  FAST_mediation(
    make_results(filtered_mismatch$treatment, FALSE),
    make_results(filtered_mismatch$outcome, FALSE),
    "continuous"
  ),
  "Filtered-probe membership disagrees",
  "Filtered membership disagreement fails"
)

stratum_mismatch_treatment <- make_results(fixture$treatment, sex_strata = TRUE)
stratum_mismatch_outcome <- make_results(fixture$outcome, sex_strata = TRUE)
stratum_mismatch_outcome$analysis_change["male"] <- list(NULL)
stratum_result <- expect_warning(
  FAST_mediation(
    stratum_mismatch_treatment,
    stratum_mismatch_outcome,
    "continuous"
  ),
  "available in only one input",
  "Mismatched strata warn and skip"
)
expect_true(is.null(stratum_result$analysis_change$male), "Mismatched stratum returns NULL")

unsupported_stratum <- make_all_only_results(fixture$treatment)
unsupported_stratum$analysis_change$femlae <- fixture$treatment
expect_error(
  FAST_mediation(
    unsupported_stratum,
    make_all_only_results(fixture$outcome),
    "continuous"
  ),
  "unsupported strata",
  "Unknown stratum names fail validation"
)

cat("\nAll mediation comprehensive tests passed.\n")
