source(file.path("mediation", "main.R"))

check <- function(value, label) {
  if (!isTRUE(value)) stop("FAIL: ", label, call. = FALSE)
  cat("PASS ", label, "\n", sep = "")
}
equal <- function(actual, expected) {
  isTRUE(all.equal(actual, expected, tolerance = 1e-12, check.attributes = FALSE))
}
reject <- function(expression, pattern, label) {
  message <- tryCatch({ force(expression); NULL }, error = function(error) conditionMessage(error))
  check(!is.null(message) && grepl(pattern, message), label)
}
warnings_from <- function(expression) {
  messages <- character()
  result <- withCallingHandlers(expression, warning = function(warning) {
    messages <<- c(messages, conditionMessage(warning))
    invokeRestart("muffleWarning")
  })
  list(result = result, warnings = messages)
}
wrap <- function(table, role, metadata = NULL) {
  entry <- setNames(list(table), if (role == "alpha") "treatment_effects" else "outcome_effects")
  result <- list(analysis_change = list(all = entry), analysis_level = list(all = entry))
  if (!is.null(metadata)) result$metadata <- metadata
  result
}
extract <- function(result) result$analysis_change$all$mediation_effects
key <- function(table) paste(table$ANALYTE_NAME, table$FU, table$CONTRAST, sep = "/")
labels <- c(`0` = "Usual care", `2` = "Diet", `7` = "Exercise")
base <- expand.grid(ANALYTE_NAME = LETTERS[1:4], FU = 1:2, stringsAsFactors = FALSE)
shared <- transform(base, EFFECT_SIZE = rep(c(0.5, -0.2, 0.1, 0), 2L),
                    SE = rep(c(0.1, 0.05, 0.12, 0.2), 2L), STATUS = "ok")
alpha <- do.call(rbind, lapply(c("2", "7"), function(arm) {
  transform(base, CONTRAST = paste0(arm, "_vs_0"), TREATMENT = unname(labels[arm]),
            REFERENCE = unname(labels["0"]), EFFECT_SIZE = if (arm == "2") {
              c(0.8, 0.3, -0.4, 0, 0.6, 0.1, -0.25, 0.2)
            } else c(-0.2, 0.8, 0.1, 0, 0.15, -0.8, 0.6, 0.1),
            SE = 0.1, STATUS = "ok")
}))
shared_metadata <- list(analyte_slope = "common_across_treatment_arms", treatment_labels = labels)
alpha_result <- wrap(alpha, "alpha", list(treatment_labels = labels))
shared_result <- wrap(shared[rev(seq_len(nrow(shared))), ], "beta", shared_metadata)

cat("Synthetic multi-arm mediation tests\n\n")
shared_mediation <- FAST_mediation(alpha_result, shared_result)
shared_effects <- extract(shared_mediation)
expected_beta <- shared$EFFECT_SIZE[match(paste(shared_effects$ANALYTE_NAME, shared_effects$FU),
                                        paste(shared$ANALYTE_NAME, shared$FU))]
check(nrow(shared_effects) == nrow(alpha) && !anyDuplicated(key(shared_effects)) &&
        equal(shared_effects$BETA_EFFECT, expected_beta) &&
        equal(shared_effects$INDIRECT_EFFECT, shared_effects$ALPHA_EFFECT * expected_beta),
      "Each alpha contrast joins the same shared beta by exact analyte/FU")
diagnostics <- shared_mediation$analysis_change$all$join_diagnostics
check(diagnostics$matched_rows == 16L && diagnostics$matched_outcome_rows == 8L &&
        diagnostics$treatment_only_rows == 0L && diagnostics$outcome_only_rows == 0L,
      "Shared-beta diagnostics distinguish contrast rows from unique matched beta rows")
effects <- shared_effects
check(equal(effects$BETA_EFFECT[effects$CONTRAST == "2_vs_0"],
            effects$BETA_EFFECT[effects$CONTRAST == "7_vs_0"]) &&
        !equal(effects$INDIRECT_EFFECT[effects$CONTRAST == "2_vs_0"],
               effects$INDIRECT_EFFECT[effects$CONTRAST == "7_vs_0"]),
      "Different alpha contrasts produce separate indirect effects using one common beta")
check(all(effects$TREATMENT[effects$CONTRAST == "2_vs_0"] == "Diet") &&
        all(effects$TREATMENT[effects$CONTRAST == "7_vs_0"] == "Exercise") &&
        all(effects$REFERENCE == "Usual care"), "Contrast labels survive alignment")
expected_se <- sqrt(effects$BETA_EFFECT^2 * effects$ALPHA_SE^2 + effects$ALPHA_EFFECT^2 * effects$BETA_SE^2)
check(equal(effects$INDIRECT_EFFECT, effects$ALPHA_EFFECT * effects$BETA_EFFECT) &&
        equal(effects$INDIRECT_SE, expected_se) &&
        equal(effects$CI_LOWER, effects$INDIRECT_EFFECT - qnorm(0.975) * expected_se) &&
        equal(effects$CI_UPPER, effects$INDIRECT_EFFECT + qnorm(0.975) * expected_se),
      "Common-beta mediation retains the exact Sobel formula and 95% normal confidence intervals")
for (fu in 1:2) {
  for (contrast in c("2_vs_0", "7_vs_0")) {
    rows <- effects$FU == fu & effects$CONTRAST == contrast
    check(equal(effects$BH_P_VALUE[rows], p.adjust(effects$P_VALUE[rows], "BH")),
          paste("BH uses a separate family for", contrast, "FU", fu))
  }
}
check(identical(attr(shared_mediation, "metadata")$estimand, "pure_and_total_natural_indirect_effect_under_shared_beta") &&
        identical(attr(shared_mediation, "metadata")$analyte_slope, "common_across_treatment_arms") &&
        is.null(attr(shared_mediation, "metadata")$beta_model) &&
        identical(attr(shared_mediation, "inference"), "sobel"),
      "Metadata describes the common slope and estimand without exposing a model option")
no_metadata <- FAST_mediation(wrap(alpha, "alpha"), wrap(shared, "beta"))
check(equal(extract(no_metadata), effects), "Common-beta tables remain accepted without metadata")
old_metadata <- FAST_mediation(alpha_result, wrap(shared, "beta",
                               list(beta_model = "shared_across_treatment_arms", treatment_labels = labels)))
check(equal(extract(old_metadata), effects), "Prior shared-beta metadata remains compatible")
legacy_shared <- FAST_mediation(wrap(alpha, "alpha"), wrap(shared[c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE")], "beta"))
check(equal(extract(legacy_shared), shared_effects), "Legacy shared beta tables accept multiple labeled alpha contrasts")

legacy_alpha <- alpha[alpha$CONTRAST == "2_vs_0", c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE")]
legacy_beta <- shared[c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE")]
legacy <- extract(FAST_mediation(wrap(legacy_alpha, "alpha"), wrap(legacy_beta, "beta")))
check(all(legacy$CONTRAST == "1_vs_0") && all(legacy$TREATMENT == "Treatment") &&
        all(legacy$REFERENCE == "Control") && all(legacy$STATUS == "ok") &&
        equal(legacy$INDIRECT_EFFECT, legacy$ALPHA_EFFECT * legacy$BETA_EFFECT),
      "Unlabeled binary inputs retain numerical behavior and receive familiar contrast labels")
binary <- extract(FAST_mediation(wrap(legacy_alpha, "alpha"), wrap(legacy_beta, "beta",
                 list(treatment_labels = c(`0` = "Usual care", `1` = "Intervention")))))
check(equal(binary$INDIRECT_EFFECT, legacy$INDIRECT_EFFECT),
      "Unlabeled legacy alpha does not falsely conflict with explicit beta metadata labels")

unavailable_alpha <- alpha
unavailable_alpha$STATUS[1L] <- "not_estimable"
unavailable_alpha[1L, c("EFFECT_SIZE", "SE")] <- NA_real_
unavailable_beta <- shared
failed <- unavailable_beta$ANALYTE_NAME == "B" & unavailable_beta$FU == 2L
unavailable_beta$STATUS[failed] <- "fit_error"
unavailable_beta[failed, c("EFFECT_SIZE", "SE")] <- NA_real_
missing <- unavailable_beta$ANALYTE_NAME == "C" & unavailable_beta$FU == 1L
unavailable_beta <- unavailable_beta[!missing, ]
partial <- warnings_from(FAST_mediation(wrap(unavailable_alpha, "alpha"), wrap(unavailable_beta, "beta")))
unavailable <- extract(partial$result)
missing_rows <- unavailable[unavailable$ANALYTE_NAME == "C" & unavailable$FU == 1L, ]
check(length(partial$warnings) == 2L && all(grepl("Partial mediation join", partial$warnings)) &&
        nrow(unavailable) == nrow(alpha) && nrow(missing_rows) == 2L &&
        all(missing_rows$BETA_STATUS == "not_estimable") &&
        all(missing_rows$STATUS == "not_estimable") && all(is.na(missing_rows$BETA_EFFECT)),
      "Absent common beta rows warn and preserve every requested alpha contrast")
check(unavailable$ALPHA_STATUS[key(unavailable) == "A/1/2_vs_0"] == "not_estimable" &&
        all(unavailable$BETA_STATUS[unavailable$ANALYTE_NAME == "B" & unavailable$FU == 2L] == "fit_error") &&
        all(unavailable$STATUS[unavailable$ANALYTE_NAME == "B" & unavailable$FU == 2L] == "fit_error") &&
        all(is.na(as.matrix(unavailable[unavailable$STATUS != "ok", c("INDIRECT_EFFECT", "INDIRECT_SE", "P_VALUE", "BH_P_VALUE")]))),
      "Unavailable path statuses propagate without losing contrast identities or producing inference")
all_unavailable <- alpha
all_unavailable$STATUS <- "not_estimable"
all_unavailable[c("EFFECT_SIZE", "SE")] <- NA_real_
empty_inference <- extract(FAST_mediation(wrap(all_unavailable, "alpha"), shared_result))
check(nrow(empty_inference) == nrow(alpha) && all(is.na(empty_inference$BH_P_VALUE)),
      "Entirely unavailable inference retains rows and empty BH families")

filtered_alpha <- unavailable_alpha
filtered_beta <- shared
filtered_alpha$IN_FILTERED_SET <- filtered_alpha$ANALYTE_NAME %in% c("A", "C")
filtered_beta$IN_FILTERED_SET <- filtered_beta$ANALYTE_NAME %in% c("A", "C")
filtered_alpha$BH_P_VALUE_FILTERED <- NA_real_
filtered_beta$BH_P_VALUE_FILTERED <- NA_real_
filtered <- extract(FAST_mediation(wrap(filtered_alpha, "alpha"), wrap(filtered_beta, "beta")))
check(all(filtered$IN_FILTERED_SET[filtered$ANALYTE_NAME == "A"]) &&
        all(is.na(filtered$BH_P_VALUE_FILTERED[!filtered$IN_FILTERED_SET])) &&
        is.na(filtered$BH_P_VALUE_FILTERED[key(filtered) == "A/1/2_vs_0"]),
      "Explicit membership is authoritative even when upstream filtered inference is unavailable")
for (group in unique(paste(filtered$FU, filtered$CONTRAST))) {
  rows <- paste(filtered$FU, filtered$CONTRAST) == group & filtered$IN_FILTERED_SET
  check(equal(filtered$BH_P_VALUE_FILTERED[rows], p.adjust(filtered$P_VALUE[rows], "BH")),
        paste("Filtered BH uses membership and available p-values in family", group))
}
old_filtered_beta <- shared
old_filtered_beta$BH_P_VALUE_FILTERED <- ifelse(shared$ANALYTE_NAME %in% c("A", "C"), 0.5, NA_real_)
mixed <- extract(FAST_mediation(wrap(filtered_alpha, "alpha"), wrap(old_filtered_beta, "beta")))
check(all(mixed$IN_FILTERED_SET == (mixed$ANALYTE_NAME %in% c("A", "C"))),
      "New explicit membership interoperates with complete legacy filtered-BH membership")
old_filtered_alpha <- alpha
old_filtered_alpha$BH_P_VALUE_FILTERED <- ifelse(alpha$ANALYTE_NAME %in% c("A", "C"), 0.7, NA_real_)
new_filtered_beta <- shared
new_filtered_beta$IN_FILTERED_SET <- shared$ANALYTE_NAME %in% c("A", "C")
mixed_reverse <- extract(FAST_mediation(wrap(old_filtered_alpha, "alpha"), wrap(new_filtered_beta, "beta")))
check(identical(mixed_reverse$IN_FILTERED_SET, mixed$IN_FILTERED_SET),
      "Complete legacy alpha membership interoperates with new beta membership")
bad <- filtered_alpha
bad$IN_FILTERED_SET <- NULL
reject(FAST_mediation(wrap(bad, "alpha"), wrap(filtered_beta, "beta")), "unavailable inference without IN_FILTERED_SET",
       "Ambiguous legacy membership with unavailable inference is rejected")
bad <- filtered_beta
bad$IN_FILTERED_SET[bad$ANALYTE_NAME == "A"] <- FALSE
reject(FAST_mediation(wrap(filtered_alpha, "alpha"), wrap(bad, "beta")), "membership disagrees", "Mismatched explicit membership is rejected")

tte_beta <- shared
tte_beta$HAZARD_RATIO <- exp(tte_beta$EFFECT_SIZE)
tte_result <- FAST_mediation(alpha_result, wrap(tte_beta, "beta", shared_metadata), "tte")
check(equal(extract(tte_result)$INDIRECT_HAZARD_RATIO, exp(effects$INDIRECT_EFFECT)) &&
        identical(attr(tte_result, "metadata")$estimand, "cox_product_of_coefficients_approximation") &&
        grepl("not a general causal", attr(tte_result, "metadata")$survival_interpretation),
      "Cox retains exp(product) with explicit approximation metadata")

reject(FAST_mediation(wrap(rbind(alpha, alpha[1L, ]), "alpha"), shared_result), "duplicate ANALYTE_NAME x FU x CONTRAST", "Duplicate alpha contrast keys fail")
reject(FAST_mediation(alpha_result, wrap(rbind(shared, shared[1L, ]), "beta")),
       "duplicate ANALYTE_NAME x FU", "Duplicate beta analyte/FU keys fail")
arm_beta <- transform(shared, ARM = "2", TREATMENT = "Diet")
reject(FAST_mediation(alpha_result, wrap(arm_beta, "beta")), "one common analyte slope",
       "Prior ARM beta tables cannot be silently treated as common slopes")
reject(FAST_mediation(alpha_result, wrap(arm_beta, "beta", shared_metadata)), "one common analyte slope",
       "Common-slope metadata cannot disguise ARM beta tables")
reject(FAST_mediation(alpha_result, wrap(shared, "beta", list(beta_model = "arm_specific"))),
       "arm_specific.*not supported", "Prior arm-specific metadata is rejected even without ARM columns")
reject(FAST_mediation(alpha_result, wrap(shared, "beta", list(analyte_slope = "different_by_arm"))),
       "analyte_slope must be", "Conflicting slope metadata is rejected")
bad <- shared_result
bad$analysis_level$all$outcome_effects <- arm_beta
reject(FAST_mediation(alpha_result, bad), "one common analyte slope", "An ARM table in any analysis is rejected")
for (contrast in c("0_vs_0", "02_vs_0", "2_vs_7", "Diet")) {
  bad <- alpha
  bad$CONTRAST[1L] <- contrast
  reject(FAST_mediation(wrap(bad, "alpha"), shared_result), "CONTRAST", paste("Invalid contrast fails:", contrast))
}
bad <- alpha
bad$TREATMENT[1L] <- "Wrong intervention"
reject(FAST_mediation(wrap(bad, "alpha"), shared_result), "consistently", "One arm cannot have conflicting labels")
bad <- shared_result
bad$metadata$treatment_labels["2"] <- "Wrong intervention"
reject(FAST_mediation(wrap(alpha, "alpha"), bad), "labels disagree", "Conflicting alpha/beta metadata labels fail")
bad <- alpha
bad$REFERENCE <- "Different control"
reject(FAST_mediation(wrap(bad, "alpha"), shared_result), "labels disagree", "Conflicting control labels fail")
bad <- shared_result
bad$metadata$treatment_labels <- c(`0` = "Usual care", `1` = "Other")
reject(FAST_mediation(alpha_result, bad), "requested alpha arm", "Wrong shared-beta treatment mapping fails")
bad <- legacy_alpha
bad$SE[1L] <- NA_real_
reject(FAST_mediation(wrap(bad, "alpha"), wrap(legacy_beta, "beta")), "finite numeric", "Legacy unavailable inference needs explicit STATUS")
bad <- alpha
bad$EFFECT_SIZE[1L] <- NA_real_
reject(FAST_mediation(wrap(bad, "alpha"), shared_result), "finite numeric", "Rows claiming ok must have finite estimates")
bad <- alpha
bad$STATUS[1L] <- "fit_error"
reject(FAST_mediation(wrap(bad, "alpha"), shared_result), "must be NA", "Unavailable status cannot retain contradictory finite path values")
bad <- alpha
bad$STATUS[1L] <- "unknown"
reject(FAST_mediation(wrap(bad, "alpha"), shared_result), "STATUS must be", "Unknown status values fail")
bad <- alpha
bad$FU[1L] <- 1.5
reject(FAST_mediation(wrap(bad, "alpha"), shared_result), "positive integer", "Fractional follow-up keys cannot be silently truncated")
bad <- shared
bad$ANALYTE_NAME <- paste0("other_", bad$ANALYTE_NAME)
reject(FAST_mediation(alpha_result, wrap(bad, "beta")), "No matching", "Complete key disagreement retains the existing no-overlap error")
check(equal(alpha_result, wrap(alpha, "alpha", list(treatment_labels = labels))) &&
        equal(shared_result, wrap(shared[rev(seq_len(nrow(shared))), ], "beta", shared_metadata)),
      "Input objects are unchanged by normalization and alignment")
cat("\nAll synthetic multi-arm mediation tests passed.\n")
