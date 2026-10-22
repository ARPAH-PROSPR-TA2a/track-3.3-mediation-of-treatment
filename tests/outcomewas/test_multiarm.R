source(file.path("outcomewas", "main.R"))
source(file.path("tests", "outcomewas", "synthetic_multiarm.R"))

check <- function(value, label) {
  if (!isTRUE(value)) stop("FAIL: ", label, call. = FALSE)
  cat("PASS ", label, "\n", sep = "")
}
equal <- function(actual, expected, tolerance = 1e-9) {
  isTRUE(all.equal(actual, expected, tolerance = tolerance, check.attributes = FALSE))
}
reject <- function(expression, pattern, label) {
  message <- tryCatch({ .quiet_multiarm(expression); NULL },
                      error = function(error) conditionMessage(error))
  check(!is.null(message) && grepl(pattern, message), label)
}
run_fixture <- function(fixture, labels = fixture$treatment_labels, ...) {
  .quiet_multiarm(FAST_outcome_WAS(
    fixture$pheno, fixture$omics, omics_type = "other",
    additional_covariates = "AGE", treatment_labels = labels, ...
  ))
}
analyses <- c("analysis_change", "analysis_level")
strata <- c("all", "male", "female")

check_references <- function(result, fixture, outcome_type, label) {
  for (response in c("change", "level")) {
    for (stratum in strata) {
      effects <- result[[paste0("analysis_", response)]][[stratum]]$outcome_effects
      check(nrow(effects) == 2L * nrow(fixture$omics) &&
              !anyDuplicated(effects[c("ANALYTE_NAME", "FU")]) &&
              !any(c("CONTRAST", "ARM") %in% names(effects)) && all(effects$STATUS == "ok"),
            paste(label, response, stratum, "has one shared beta per analyte/FU"))
      for (i in seq_len(nrow(effects))) {
        expected <- .manual_multiarm_outcome(
          fixture, response, stratum, effects$FU[i], effects$ANALYTE_NAME[i], outcome_type
        )
        if (!equal(as.numeric(effects[i, names(expected)]), unname(expected))) {
          stop("Independent model mismatch: ", label, "/", response, "/", stratum,
               "/", effects$ANALYTE_NAME[i], "/FU", effects$FU[i])
        }
      }
      for (fu in 1:2) {
        rows <- effects$FU == fu
        check(equal(effects$BH_P_VALUE[rows], p.adjust(effects$P_VALUE[rows], "BH")),
              paste(label, response, stratum, paste0("FU", fu), "BH family is unchanged"))
      }
      if (outcome_type == "tte") {
        check(equal(effects$HAZARD_RATIO, exp(effects$EFFECT_SIZE)),
              paste(label, response, stratum, "exports log-hazard beta and hazard ratio"))
      }
    }
  }
  cat("PASS ", label, " all estimates, SEs and p-values match independent pooled models\n", sep = "")
}

cat("OutcomeWAS synthetic multi-arm tests\n\n")
continuous <- .make_multiarm_outcome_fixture()
tte <- .make_multiarm_outcome_fixture("tte")
continuous_result <- run_fixture(continuous, n_cores = 1L)
tte_result <- run_fixture(tte, n_cores = 1L)
check_references(continuous_result, continuous, "continuous", "Multi-arm linear")
check_references(tte_result, tte, "tte", "Multi-arm Cox")

metadata <- continuous_result$metadata
check(identical(metadata$treatment_labels, continuous$treatment_labels) &&
        identical(metadata$treatments$ARM, c("2", "7")) &&
        identical(metadata$treatments$CONTRAST, c("2_vs_0", "7_vs_0")) &&
        identical(metadata$treatments$REFERENCE, rep("Usual care", 2L)) &&
        identical(metadata$outcome_type, "continuous") &&
        identical(metadata$additional_covariates, "AGE") &&
        identical(metadata$analyte_slope, "common_across_treatment_arms") &&
        !"beta_model" %in% names(metadata) &&
        identical(metadata$p_adjust, "BH within analysis, stratum and FU") &&
        identical(metadata$schema_version, "2.0"),
      "Metadata preserves nonconsecutive categorical codes, labels and shared-beta specification")
check(identical(tte_result$metadata$outcome_type, "tte"), "Cox metadata identifies outcome type")
check(!"beta_model" %in% names(formals(FAST_outcome_WAS)) &&
        !"beta_model" %in% names(formals(FAST_outcome_WAS_reports)),
      "Analysis and reports expose no slope-model choice")
check(identical(metadata$model_formula,
                "OUTCOME ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + AGE") &&
        identical(tte_result$metadata$model_formula,
                  "survival::Surv(OUTCOME_TIME, OUTCOME_STATUS) ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE + AGE"),
      "Metadata describes the fixed common-slope linear and Cox formulas")
for (outcome_type in c("continuous", "tte")) {
  result <- if (outcome_type == "continuous") continuous_result else tte_result
  for (analysis in analyses) {
    for (stratum in strata) {
      terms <- unique(result[[analysis]][[stratum]]$coefficients$COEFFICIENT)
      expected <- c("analyte", "analyte_baseline", "TREATMENT_GROUP2", "TREATMENT_GROUP7", "AGE")
      if (outcome_type == "continuous") expected <- c("(Intercept)", expected)
      if (stratum == "all") expected <- c(expected, "FEMALE1")
      check(setequal(terms, expected),
            paste(outcome_type, analysis, stratum,
                  "expands categorical treatment main effects without analyte interactions"))
    }
  }
}

binary <- .make_multiarm_outcome_fixture(arms = c("0", "1"), n_analytes = 3L)
binary_result <- run_fixture(binary, labels = NULL, n_cores = 1L)
check_references(binary_result, binary, "continuous", "Ordinary binary linear")
check(identical(binary_result$metadata$treatment_labels, c(`0` = "Control", `1` = "Treatment")),
      "Binary labels remain optional with familiar defaults")
binary_tte <- .make_multiarm_outcome_fixture("tte", arms = c("0", "1"), n_analytes = 3L)
binary_tte_result <- run_fixture(binary_tte, labels = NULL, n_cores = 1L)
check_references(binary_tte_result, binary_tte, "tte", "Ordinary binary Cox")
single_active <- .make_multiarm_outcome_fixture(arms = c("0", "7"), n_analytes = 1L)
single_result <- run_fixture(single_active, labels = NULL, n_cores = 1L)
check(identical(single_result$metadata$treatments$CONTRAST, "7_vs_0"),
      "A sole nonstandard positive treatment code remains categorical")

reject(run_fixture(continuous, labels = NULL, n_cores = 1L), "treatment_labels", "Multiple active arms require labels")
bad_labels <- list(
  missing_control = continuous$treatment_labels[-1L],
  extra = c(continuous$treatment_labels, `9` = "Unused"),
  duplicate_label = c(`0` = "Same", `2` = "Same", `7` = "Other"),
  blank = c(`0` = "Control", `2` = " ", `7` = "Other"),
  missing = c(`0` = "Control", `2` = NA_character_, `7` = "Other"),
  unnamed = unname(continuous$treatment_labels)
)
for (name in names(bad_labels)) {
  reject(run_fixture(continuous, labels = bad_labels[[name]], n_cores = 1L),
         "treatment_labels", paste("Rejects invalid treatment labels:", name))
}
for (code in c("-1", "1.5", "01", "A", NA_character_)) {
  bad <- continuous
  bad$pheno$TREATMENT_GROUP <- as.character(bad$pheno$TREATMENT_GROUP)
  bad$pheno$TREATMENT_GROUP[1L] <- code
  reject(run_fixture(bad, n_cores = 1L), "TREATMENT_GROUP",
         paste("Rejects invalid assignment code", code))
}
bad <- continuous
bad$pheno$TREATMENT_GROUP <- as.character(bad$pheno$TREATMENT_GROUP)
bad$pheno$TREATMENT_GROUP[1L] <- "2"
bad$pheno$AGE[1L] <- NA_real_
reject(run_fixture(bad, n_cores = 1L), "constant within SUBJECT_ID",
       "Assignment constancy is checked before covariate filtering")
for (keep in list(c("2", "7"), "0")) {
  bad <- continuous
  bad$pheno <- bad$pheno[bad$pheno$TREATMENT_GROUP %in% keep, ]
  reject(run_fixture(bad, n_cores = 1L), "control|active treatment",
         paste("Rejects missing study-wide control or active arm:", paste(keep, collapse = ",")))
}

# An outcome slope can remain estimable without every arm, including without
# control, because beta is a common slope rather than a treatment contrast.
for (outcome_type in c("continuous", "tte")) {
  missing_arms <- if (outcome_type == "continuous") continuous else tte
  pheno <- missing_arms$pheno
  missing_arms$omics[1L, pheno$SAMPLE_ID[pheno$TREATMENT_GROUP == "7"]] <- NA_real_
  missing_arms$omics[2L, pheno$SAMPLE_ID[pheno$TREATMENT_GROUP == "0"]] <- NA_real_
  missing_arms$omics[3L, pheno$SAMPLE_ID[pheno$TREATMENT_GROUP != "2"]] <- NA_real_
  missing_result <- run_fixture(missing_arms, n_cores = 1L)
  check_references(missing_result, missing_arms, outcome_type,
                   paste("Missing-arm", outcome_type))
}

filtered_arm <- continuous
filtered_arm$pheno$AGE[filtered_arm$pheno$TREATMENT_GROUP == "7"] <- NA_real_
filtered_result <- run_fixture(filtered_arm, n_cores = 1L)
check(identical(filtered_result$metadata$treatment_labels, continuous$treatment_labels) &&
        all(filtered_result$analysis_change$all$outcome_effects$STATUS == "ok"),
      "An arm removed by phenotype filtering retains its requested treatment metadata")

for (outcome_type in c("continuous", "tte")) {
  nuisance_alias <- .make_multiarm_outcome_fixture(outcome_type, n_analytes = 2L)
  nuisance_alias$pheno$AGE_COPY <- nuisance_alias$pheno$AGE
  nuisance_result <- .quiet_multiarm(FAST_outcome_WAS(
    nuisance_alias$pheno, nuisance_alias$omics, omics_type = "other",
    additional_covariates = c("AGE", "AGE_COPY"), n_cores = 1L,
    treatment_labels = nuisance_alias$treatment_labels
  ))
  for (analysis in analyses) {
    check(all(nuisance_result[[analysis]]$all$outcome_effects$STATUS == "ok"),
          paste(outcome_type, analysis, "a redundant nuisance term does not invalidate an identifiable beta"))
  }
  analyte_alias <- nuisance_alias
  baseline <- analyte_alias$pheno[analyte_alias$pheno$FU == "0", ]
  for (fu in 1:2) {
    followup <- analyte_alias$pheno[analyte_alias$pheno$FU == fu, ]
    matched_baseline <- baseline[match(followup$SUBJECT_ID, baseline$SUBJECT_ID), ]
    analyte_alias$omics[1L, followup$SAMPLE_ID] <-
      as.numeric(analyte_alias$omics[1L, matched_baseline$SAMPLE_ID]) + followup$AGE
  }
  aliased_result <- run_fixture(analyte_alias, n_cores = 1L)
  for (analysis in analyses) {
    effects <- aliased_result[[analysis]]$all$outcome_effects
    target <- effects$ANALYTE_NAME == analyte_alias$omics$ANALYTE_NAME[1L]
    check(all(effects$STATUS[target] == "not_estimable") &&
            all(is.na(effects$EFFECT_SIZE[target])) && all(effects$STATUS[!target] == "ok"),
          paste(outcome_type, analysis, "analyte aliased with later nuisance terms is not estimable"))
  }
}

unavailable <- continuous
unavailable$omics[2L, -1L] <- 1
unavailable$omics[3L, -1L] <- NA_real_
unavailable_result <- run_fixture(unavailable, n_cores = 1L)
for (analysis in analyses) {
  for (stratum in strata) {
    effects <- unavailable_result[[analysis]][[stratum]]$outcome_effects
    constant <- effects$ANALYTE_NAME == unavailable$omics$ANALYTE_NAME[2L]
    failed <- effects$ANALYTE_NAME == unavailable$omics$ANALYTE_NAME[3L]
    check(nrow(effects) == 2L * nrow(unavailable$omics) &&
            all(effects$STATUS[constant] == "not_estimable") &&
            all(effects$STATUS[failed] == "fit_error") &&
            all(is.na(as.matrix(effects[constant | failed, c("EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")]))),
          paste(analysis, stratum, "retains distinct unavailable and failed rows with NA inference"))
    for (fu in 1:2) {
      rows <- effects$FU == fu
      check(equal(effects$BH_P_VALUE[rows], p.adjust(effects$P_VALUE[rows], "BH")),
            paste(analysis, stratum, paste0("FU", fu), "BH excludes unavailable p-values"))
    }
  }
}
all_failed <- continuous
all_failed$omics[, -1L] <- NA_real_
# One observed follow-up value has no baseline partner. Validation must let it
# reach the fit-failure placeholder rather than evaluating an undefined variance.
one_sample <- all_failed$pheno$SAMPLE_ID[all_failed$pheno$FU == "1"][1L]
all_failed$omics[1L, one_sample] <- 0.25
all_failed_result <- run_fixture(all_failed, n_cores = 1L)
for (analysis in analyses) {
  for (stratum in strata) {
    entry <- all_failed_result[[analysis]][[stratum]]
    check(is.data.frame(entry$coefficients) && nrow(entry$coefficients) == 0L &&
            all(c("ANALYTE_NAME", "FU", "COEFFICIENT", "N_OBS", "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE") %in% names(entry$coefficients)) &&
            is.numeric(entry$coefficients$EFFECT_SIZE) &&
            nrow(entry$outcome_effects) == 2L * nrow(all_failed$omics) &&
            all(entry$outcome_effects$STATUS == "fit_error"),
          paste(analysis, stratum, "all-fit-failure tables retain keys and an empty typed coefficient table"))
  }
}
check(all(all_failed_result$analysis_change$all$outcome_effects$STATUS[
  all_failed_result$analysis_change$all$outcome_effects$ANALYTE_NAME == all_failed$omics$ANALYTE_NAME[1L]
] == "fit_error"), "A single observed omics value yields placeholders without aborting validation")
no_events <- tte
no_events$pheno$OUTCOME_STATUS <- 0L
no_event_result <- run_fixture(no_events, n_cores = 1L)
for (analysis in analyses) {
  for (stratum in strata) {
    entry <- no_event_result[[analysis]][[stratum]]
    check(is.data.frame(entry$coefficients) && nrow(entry$coefficients) == 0L &&
            nrow(entry$outcome_effects) == 2L * nrow(no_events$omics) &&
            all(entry$outcome_effects$STATUS == "not_estimable") &&
            all(is.na(entry$outcome_effects$HAZARD_RATIO)),
          paste(analysis, stratum, "zero-event Cox models retain unavailable beta rows"))
  }
}
per_probe_no_events <- tte
per_probe_no_events$omics[1L, tte$pheno$SAMPLE_ID[tte$pheno$OUTCOME_STATUS == 1L]] <- NA_real_
per_probe_result <- run_fixture(per_probe_no_events, n_cores = 1L)
for (analysis in analyses) {
  effects <- per_probe_result[[analysis]]$all$outcome_effects
  check(all(effects$STATUS[effects$ANALYTE_NAME == tte$omics$ANALYTE_NAME[1L]] == "not_estimable") &&
          all(effects$STATUS[effects$ANALYTE_NAME != tte$omics$ANALYTE_NAME[1L]] == "ok"),
        paste(analysis, "detects no events after analyte-specific missingness"))
}
tte_unavailable <- tte
tte_unavailable$omics[2L, -1L] <- 1
tte_unavailable$omics[3L, -1L] <- NA_real_
tte_unavailable_result <- run_fixture(tte_unavailable, n_cores = 1L)
for (analysis in analyses) {
  for (stratum in strata) {
    effects <- tte_unavailable_result[[analysis]][[stratum]]$outcome_effects
    constant <- effects$ANALYTE_NAME == tte$omics$ANALYTE_NAME[2L]
    failed <- effects$ANALYTE_NAME == tte$omics$ANALYTE_NAME[3L]
    check(nrow(effects) == 2L * nrow(tte$omics) &&
            all(effects$STATUS[constant] == "not_estimable") &&
            all(effects$STATUS[failed] == "fit_error") &&
            all(is.na(as.matrix(effects[constant | failed, c("EFFECT_SIZE", "SE", "P_VALUE", "HAZARD_RATIO", "BH_P_VALUE")]))),
          paste(analysis, stratum, "Cox failures and unavailable coefficients retain NA inference"))
  }
}

# The public probe lists supply identifiers only; all measurements are simulated.
full_probes <- readRDS("outcomewas/Data/FAST_epicv1_epicv2_probe_list.rds")
filtered_probes <- readRDS("outcomewas/Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
dnam <- unavailable
dnam$omics$ANALYTE_NAME <- c(head(intersect(full_probes, filtered_probes), 3L),
                            head(setdiff(full_probes, filtered_probes), 3L))
run_dnam <- function(n_cores = 1L, checkpoint_dir = NULL) .quiet_multiarm(FAST_outcome_WAS(
  dnam$pheno, dnam$omics, omics_type = "DNAm", additional_covariates = "AGE",
  n_cores = n_cores, checkpoint_dir = checkpoint_dir, checkpoint_batch_size = 2L,
  treatment_labels = dnam$treatment_labels
))
dnam_result <- run_dnam()
for (analysis in analyses) {
  for (stratum in strata) {
    effects <- dnam_result[[analysis]][[stratum]]$outcome_effects
    check(identical(effects$IN_FILTERED_SET, effects$ANALYTE_NAME %in% filtered_probes) &&
            all(effects$IN_FILTERED_SET[effects$STATUS != "ok"]) &&
            all(is.na(effects$BH_P_VALUE_FILTERED[effects$STATUS != "ok"])),
          paste(analysis, stratum, "DNAm membership survives unavailable inference"))
    for (fu in 1:2) {
      rows <- effects$FU == fu & effects$IN_FILTERED_SET
      check(equal(effects$BH_P_VALUE_FILTERED[rows], p.adjust(effects$P_VALUE[rows], "BH")),
            paste(analysis, stratum, paste0("FU", fu), "filtered BH uses explicit membership"))
    }
  }
}
check(!"IN_FILTERED_SET" %in% names(continuous_result$analysis_change$all$outcome_effects),
      "Non-DNAm results do not gain probe-membership columns")

local({
  checkpoint <- tempfile("outcomewas-multiarm-checkpoint-")
  on.exit(unlink(checkpoint, recursive = TRUE), add = TRUE)
  checkpointed <- run_dnam(n_cores = 2L, checkpoint_dir = checkpoint)
  check(equal(checkpointed, dnam_result), "Serial and two-worker DNAm results agree, including failures and metadata")
  files <- sort(list.files(checkpoint, "^batch_[0-9]+[.]rds$", recursive = TRUE, full.names = TRUE))
  check(length(files) == 36L, "Existing change/level x stratum x FU batch layout is retained")
  old_time <- as.POSIXct("2001-01-01", tz = "UTC")
  invisible(vapply(files, function(path) Sys.setFileTime(path, old_time), logical(1)))
  mtimes <- file.info(files)$mtime
  resumed <- run_dnam(checkpoint_dir = checkpoint)
  check(equal(resumed, dnam_result) && identical(file.info(files)$mtime, mtimes),
        "Full resume restores statuses, membership and metadata without rewriting batches")
  removed <- files[c(1L, length(files))]
  kept <- setdiff(files, removed)
  kept_times <- file.info(kept)$mtime
  unlink(removed)
  partial <- run_dnam(n_cores = 2L, checkpoint_dir = checkpoint)
  check(equal(partial, dnam_result) && all(file.exists(removed)) &&
          identical(file.info(kept)$mtime, kept_times),
        "Partial resume recreates missing batches and preserves all cached batches")
  check(length(list.files(checkpoint, "[.]tmp$", recursive = TRUE)) == 0L,
        "Completed resume leaves no temporary checkpoint files")
})
tte_parallel <- run_fixture(tte, n_cores = 2L)
check(equal(tte_parallel, tte_result), "Serial and two-worker multi-arm Cox results agree")
local({
  checkpoint <- tempfile("outcomewas-multiarm-cox-checkpoint-")
  on.exit(unlink(checkpoint, recursive = TRUE), add = TRUE)
  fresh <- run_fixture(tte_unavailable, n_cores = 2L, checkpoint_dir = checkpoint,
                       checkpoint_batch_size = 2L)
  files <- list.files(checkpoint, "^batch_[0-9]+[.]rds$", recursive = TRUE, full.names = TRUE)
  invisible(vapply(files, function(path) Sys.setFileTime(path, as.POSIXct("2001-01-01", tz = "UTC")), logical(1)))
  mtimes <- file.info(files)$mtime
  resumed <- run_fixture(tte_unavailable, n_cores = 1L, checkpoint_dir = checkpoint,
                         checkpoint_batch_size = 2L)
  check(length(files) == 36L && equal(fresh, tte_unavailable_result) &&
          equal(resumed, tte_unavailable_result) && identical(file.info(files)$mtime, mtimes),
        "Cox resume restores unavailable/failed rows and metadata without recomputing batches")
})

annotation <- data.frame(AptName = continuous$omics$ANALYTE_NAME, stringsAsFactors = FALSE)
for (column in .protein_annotation_columns) annotation[[column]] <- paste(column, annotation$AptName)
annotated <- get_annotated_results(unavailable_result, annotation, c(`1` = "12 months", `2` = "24 months"))
check(identical(annotated$metadata, unavailable_result$metadata) &&
        identical(annotated$analysis_change$all$outcome_effects$STATUS,
                  unavailable_result$analysis_change$all$outcome_effects$STATUS) &&
        nrow(annotated$analysis_change$all$outcome_effects) == nrow(unavailable_result$analysis_change$all$outcome_effects),
      "Protein annotation preserves treatment metadata, statuses and unavailable rows")

cat("\nAll OutcomeWAS synthetic multi-arm tests passed.\n")
