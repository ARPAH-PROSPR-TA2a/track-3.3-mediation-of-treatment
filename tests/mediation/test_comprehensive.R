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


expect_warnings <- function(expression, pattern = NULL, label, count = NULL) {
  messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(warning) {
      messages <<- c(messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  matched <- if (is.null(pattern)) !length(messages) else any(grepl(pattern, messages))
  expect_true(matched && (is.null(count) || length(messages) == count), label)
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


make_results <- function(entry, sex_strata = TRUE, fu_labels = NULL) {
  if (!is.null(fu_labels)) {
    for (table in intersect(c("treatment_effects", "outcome_effects"), names(entry))) {
      entry[[table]]$FU_LABEL <- unname(fu_labels[as.character(entry[[table]]$FU)])
    }
  }
  branch <- list(
    all = entry,
    male = if (sex_strata) entry else NULL,
    female = if (sex_strata) entry else NULL
  )
  results <- list(analysis_change = branch, analysis_level = branch)
  if (!is.null(fu_labels)) results$metadata <- list(fu_labels = fu_labels)
  results
}


make_all_only_results <- function(entry) {
  branch <- list(all = entry)
  list(analysis_change = branch, analysis_level = branch)
}


cat("Mediation comprehensive tests\n\n")

fixture <- make_effects()
fu_labels <- c(`0` = "Baseline", `1` = "12mo", `2` = "24mo")
treatment_results <- make_results(fixture$treatment, fu_labels = fu_labels)
outcome_results <- make_results(fixture$outcome, fu_labels = fu_labels)
results <- expect_warnings(
  FAST_mediation(treatment_results, outcome_results, "continuous"),
  label = "Matching follow-up labels proceed without warnings"
)
effects <- results$analysis_change$all$mediation_effects
print(unique(effects[c("FU", "FU_LABEL")]))
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
  make_results(filtered_fixture$treatment, sex_strata = FALSE, fu_labels = fu_labels),
  make_results(filtered_fixture$outcome, sex_strata = FALSE, fu_labels = fu_labels),
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
  make_results(tte_fixture$treatment, sex_strata = FALSE, fu_labels = fu_labels),
  make_results(tte_fixture$outcome, sex_strata = FALSE, fu_labels = fu_labels),
  "tte"
)
tte_effects <- tte_results$analysis_change$all$mediation_effects
expect_equal(tte_effects$INDIRECT_HAZARD_RATIO, exp(tte_effects$INDIRECT_EFFECT),
             "Time-to-event mediated path hazard ratio")

partial_outcome <- make_results(fixture$outcome, sex_strata = FALSE)
partial_outcome$analysis_change$all$outcome_effects <-
  partial_outcome$analysis_change$all$outcome_effects[-1, ]
partial <- expect_warnings(
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
stratum_result <- expect_warnings(
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

# Follow-up label rules preserve numeric matching and inference.
strip_fu_labels <- function(result) {
  attr(result, "metadata")$fu_labels <- NULL
  for (analysis in c("analysis_change", "analysis_level")) {
    for (stratum in names(result[[analysis]])) {
      if (is.null(result[[analysis]][[stratum]])) next
      result[[analysis]][[stratum]]$mediation_effects$FU_LABEL <- NULL
    }
  }
  result
}
cases <- list(
  continuous = list(fixture = fixture, labeled = results, strata = TRUE, type = "continuous"),
  filtered = list(fixture = filtered_fixture, labeled = filtered_results, strata = FALSE, type = "continuous"),
  tte = list(fixture = tte_fixture, labeled = tte_results, strata = FALSE, type = "tte")
)
for (name in names(cases)) {
  case <- cases[[name]]
  legacy <- expect_warnings(
    FAST_mediation(make_results(case$fixture$treatment, case$strata),
                   make_results(case$fixture$outcome, case$strata), case$type),
    label = paste(name, "unlabeled inputs proceed without warnings")
  )
  expect_true(identical(strip_fu_labels(case$labeled), legacy),
              paste(name, "labels preserve every numeric column, status, BH value, and join diagnostic"))
  expect_true(is.null(attr(legacy, "metadata")$fu_labels) &&
                !"FU_LABEL" %in% names(legacy$analysis_change$all$mediation_effects),
              paste(name, "legacy inputs retain their output structure"))
  expect_true(identical(attr(case$labeled, "metadata")$fu_labels, fu_labels),
              paste(name, "metadata preserves the full visit map including baseline"))
  for (analysis in c("analysis_change", "analysis_level")) {
    for (stratum in names(case$labeled[[analysis]])) {
      table <- case$labeled[[analysis]][[stratum]]$mediation_effects
      if (is.null(table)) next
      expect_true(identical(names(table)[1:3], c("ANALYTE_NAME", "FU", "FU_LABEL")) &&
                    identical(table$FU_LABEL, unname(fu_labels[as.character(table$FU)])),
                  paste(name, analysis, stratum, "retains visit labels beside numeric FU"))
    }
  }
}

# A single labeled source is enough. Row-only sources need no baseline metadata.
alpha_only <- expect_warnings(
  FAST_mediation(treatment_results, make_results(fixture$outcome)),
  "only by alpha.*Label agreement could not be checked",
  "Alpha-only labels emit one warning across all analyses and strata", count = 1L
)
beta_only <- expect_warnings(
  FAST_mediation(make_results(fixture$treatment), outcome_results),
  "only by beta.*Label agreement could not be checked",
  "Beta-only labels emit one warning across all analyses and strata", count = 1L
)
expect_true(identical(alpha_only, results) && identical(beta_only, results),
            "Either labeled source supplies the full display mapping")
rows_only <- treatment_results
rows_only$metadata <- NULL
row_result <- expect_warnings(
  FAST_mediation(rows_only, make_results(fixture$outcome)),
  "only by alpha.*Label agreement could not be checked",
  "Row-only labels emit one inheritance warning", count = 1L
)
expect_true(identical(attr(row_result, "metadata")$fu_labels, fu_labels[-1]) &&
              identical(row_result$analysis_change$all$mediation_effects, effects),
            "FU_LABEL columns supply a mapping without baseline metadata")
metadata_only <- outcome_results
for (analysis in c("analysis_change", "analysis_level")) {
  for (stratum in c("all", "male", "female")) {
    metadata_only[[analysis]][[stratum]]$outcome_effects$FU_LABEL <- NULL
  }
}
metadata_result <- expect_warnings(
  FAST_mediation(make_results(fixture$treatment), metadata_only),
  "only by beta.*Label agreement could not be checked",
  "Metadata-only labels emit one inheritance warning", count = 1L
)
expect_true(identical(metadata_result, results),
            "Metadata supplies labels without repeated table columns")

conflicting <- outcome_results
conflicting$metadata$fu_labels["1"] <- "24mo"
expect_error(FAST_mediation(treatment_results, conflicting), "FU_LABEL disagrees with metadata",
              "A source's table labels must agree with its metadata")
conflicting <- make_results(fixture$outcome, fu_labels = c(`0` = "Baseline", `1` = "Month 12", `2` = "24mo"))
expect_error(FAST_mediation(treatment_results, conflicting), "Follow-up labels disagree",
              "Conflicting follow-up labels stop mediation even when numeric codes match")
conflicting <- outcome_results
conflicting$metadata$fu_labels["0"] <- "Pre-intervention"
expect_error(FAST_mediation(treatment_results, conflicting), "Follow-up labels disagree.*0",
              "Shared baseline metadata labels must agree")
conflicting <- rows_only
conflicting$analysis_level$male$treatment_effects$FU_LABEL[1] <- "Different"
expect_error(FAST_mediation(conflicting, outcome_results), "consistently across analyses and strata",
              "Rows must map each code consistently across the source")
incomplete <- outcome_results
incomplete$metadata$fu_labels <- fu_labels[1:2]
expect_error(FAST_mediation(treatment_results, incomplete), "cover every observed FU",
              "An explicit source map must cover its observed effect-table visits")
invalid <- outcome_results
invalid$metadata$fu_labels <- setNames(unname(fu_labels), c("0", "1", "1"))
expect_error(FAST_mediation(treatment_results, invalid), "unique non-negative integer FU codes",
              "Duplicate visit-map keys fail validation")

# Labels must not remove an otherwise valid legacy partial-join result.
fu1_beta <- fixture$outcome
fu1_beta$outcome_effects <- fu1_beta$outcome_effects[fu1_beta$outcome_effects$FU == 1L, ]
partial_labeled <- expect_warnings(FAST_mediation(
  make_results(fixture$treatment, FALSE), make_results(fu1_beta, FALSE, fu_labels[1:2])
), "Partial mediation join", "Partially known labels preserve partial joins")
partial_legacy <- suppressWarnings(FAST_mediation(
  make_results(fixture$treatment, FALSE), make_results(fu1_beta, FALSE)
))
partial_table <- partial_labeled$analysis_change$all$mediation_effects
expect_true(identical(strip_fu_labels(partial_labeled), partial_legacy) &&
              all(partial_table$FU_LABEL[partial_table$FU == 1L] == "12mo") &&
              all(is.na(partial_table$FU_LABEL[partial_table$FU == 2L])),
            "Unknown legacy visit labels remain NA without changing status or inference")

cat("\nAll mediation comprehensive tests passed.\n")
