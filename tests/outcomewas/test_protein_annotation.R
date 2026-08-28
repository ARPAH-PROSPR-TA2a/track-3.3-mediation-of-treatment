source(file.path("outcomewas", "main.R"))

expect_true <- function(value, label) {
  if (!isTRUE(value)) stop("FAIL: ", label)
  cat("PASS ", label, "\n", sep = "")
}

expect_equal <- function(actual, expected, label) {
  if (!isTRUE(all.equal(actual, expected, check.attributes = TRUE))) {
    stop("FAIL: ", label, "\nExpected: ", paste(expected, collapse = ", "),
         "\nActual: ", paste(actual, collapse = ", "))
  }
  cat("PASS ", label, "\n", sep = "")
}

expect_error <- function(expression, pattern, label) {
  message <- tryCatch({ force(expression); NULL },
                      error = function(error) conditionMessage(error))
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

expect_rscript_success <- function(lines, label) {
  script <- tempfile("outcomewas-source-smoke-", fileext = ".R")
  writeLines(lines, script)
  output <- system2(file.path(R.home("bin"), "Rscript"), script,
                    stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  unlink(script)
  expect_true(
    identical(as.integer(status), 0L),
    paste0(label, if (length(output)) paste0(": ", paste(output, collapse = " | ")) else "")
  )
}

make_effects <- function(tte = FALSE, extra = FALSE) {
  value <- data.frame(
    ANALYTE_NAME = c("B", "A", "B", "A"),
    FU = c(1L, 1L, 5L, 5L),
    EFFECT_SIZE = c(0.4, -0.2, 0.8, -0.6),
    SE = c(0.10, 0.20, 0.15, 0.25),
    P_VALUE = c(0.01, 0.20, 0.03, 0.40),
    BH_P_VALUE = c(0.02, 0.20, 0.06, 0.40),
    stringsAsFactors = FALSE
  )
  if (tte) {
    value$HAZARD_RATIO <- exp(value$EFFECT_SIZE)
    value <- value[c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "HAZARD_RATIO",
                     "SE", "P_VALUE", "BH_P_VALUE")]
  }
  if (extra) value$FUTURE_OUTCOME_FIELD <- c("w", "x", "y", "z")
  value
}

make_coefficients <- function(tte = FALSE, extra = FALSE) {
  value <- data.frame(
    ANALYTE_NAME = rep(c("B", "A", "B", "A"), each = 2L),
    FU = rep(c(1L, 1L, 5L, 5L), each = 2L),
    COEFFICIENT = rep(c("analyte", "treatment"), 4L),
    N_OBS = rep(c(42L, 42L, 38L, 38L), each = 2L),
    EFFECT_SIZE = seq(0.1, 0.8, by = 0.1),
    SE = rep(c(0.1, 0.2), 4L),
    P_VALUE = seq(0.01, 0.08, by = 0.01),
    BH_P_VALUE = seq(0.02, 0.16, by = 0.02),
    stringsAsFactors = FALSE
  )
  if (tte) {
    value$N_EVENTS <- rep(c(12L, 12L, 9L, 9L), each = 2L)
    value$HAZARD_RATIO <- exp(value$EFFECT_SIZE)
    value <- value[c("ANALYTE_NAME", "FU", "COEFFICIENT", "N_OBS", "N_EVENTS",
                     "EFFECT_SIZE", "HAZARD_RATIO", "SE", "P_VALUE", "BH_P_VALUE")]
  }
  if (extra) value$FUTURE_COEFFICIENT_FIELD <- seq_len(nrow(value))
  value
}

make_entry <- function(coefficients, effects, note) {
  list(coefficients = coefficients, outcome_effects = effects,
       diagnostics = list(note = note, converged = TRUE))
}

make_results <- function(coefficients = make_coefficients(), effects = make_effects()) {
  list(
    analysis_change = list(
      all = make_entry(coefficients, effects, "change-all"),
      male = make_entry(coefficients, effects, "change-male"),
      female = NULL,
      report = data.frame(metric = "change-report", value = 1L)
    ),
    analysis_level = list(
      all = make_entry(coefficients, effects, "level-all"),
      male = make_entry(coefficients, effects, "level-male"),
      female = NULL,
      report = data.frame(metric = "level-report", value = 2L)
    ),
    reports = list(input_summary = data.frame(n = 84L)),
    run_metadata = list(version = "synthetic")
  )
}

make_annotation <- function(keys = c("A", "B", "UNUSED")) {
  data.frame(
    AptName = keys,
    SomaId = paste0("Soma_", keys),
    TargetFullName = paste("Full target", keys),
    Target = paste0("Target_", keys),
    UniProt = paste0("UniProt_", keys),
    EntrezGeneID = seq_along(keys) + 100L,
    EntrezGeneSymbol = paste0("Gene_", keys),
    Uniprot_Unique = paste0("UniqueUniProt_", keys),
    Symbol_Unique = paste0("UniqueSymbol_", keys),
    stringsAsFactors = FALSE
  )
}

annotation_columns <- c("SomaId", "TargetFullName", "Target", "UniProt",
                        "EntrezGeneID", "EntrezGeneSymbol", "Uniprot_Unique",
                        "Symbol_Unique")

cat("OutcomeWAS protein annotation tests\n\n")

expect_true(exists("get_annotated_results"), "New nested-results API is exported")
expect_true(!exists("get_annotated_outcome_effects"), "Old flat API is absent")
expect_true(!exists("write_annotated_outcome_effect_tables"), "Old writer API is absent")

source_assertion <- paste0(
  'stopifnot(exists("FAST_outcome_WAS"), exists("get_annotated_results"), ',
  '!exists("get_annotated_outcome_effects"), ',
  '!exists("write_annotated_outcome_effect_tables"))'
)
expect_rscript_success(
  c('source(file.path("outcomewas", "main.R"), chdir = TRUE)', source_assertion),
  "main.R sources the new API from a relative path with chdir = TRUE"
)
main_path_literal <- paste(capture.output(dput(normalizePath(file.path("outcomewas", "main.R")))),
                           collapse = "")
expect_rscript_success(
  c(paste0("source(", main_path_literal, ", chdir = TRUE)"), source_assertion),
  "main.R sources the new API from an absolute path with chdir = TRUE"
)

results <- make_results()
annotation <- make_annotation()
annotation$IGNORED_TRANSLATION_FIELD <- "not exported"
results_before <- unserialize(serialize(results, NULL))
annotation_before <- unserialize(serialize(annotation, NULL))
fu_labels <- c("1" = "3mo", "5" = "24mo")
annotated <- get_annotated_results(results, annotation, fu_labels)

track111_outcome_columns <- c(
  "ANALYTE_NAME", annotation_columns, "FU", "FU_LABEL",
  "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE"
)
expected_coefficient_columns <- c(
  "ANALYTE_NAME", annotation_columns, "FU", "FU_LABEL", "COEFFICIENT",
  "N_OBS", "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE"
)

expect_true(identical(names(annotated), names(results)), "Top-level topology is unchanged")
for (analysis_name in c("analysis_change", "analysis_level")) {
  expect_true(identical(names(annotated[[analysis_name]]), names(results[[analysis_name]])),
              paste(analysis_name, "topology is unchanged"))
  for (group in c("all", "male")) {
    output <- annotated[[analysis_name]][[group]]
    input <- results[[analysis_name]][[group]]
    expect_true(identical(names(output), names(input)),
                paste(analysis_name, group, "entry topology is unchanged"))
    expect_true(identical(names(output$outcome_effects), track111_outcome_columns),
                paste(analysis_name, group, "outcomes match the exact 15-column contract"))
    expect_true(identical(names(output$coefficients), expected_coefficient_columns),
                paste(analysis_name, group, "coefficients match the annotated schema"))
    expect_true(identical(output$diagnostics, input$diagnostics),
                paste(analysis_name, group, "other entry fields are unchanged"))
  }
  expect_true(is.null(annotated[[analysis_name]]$female),
              paste(analysis_name, "NULL strata stay NULL"))
  expect_true(identical(annotated[[analysis_name]]$report, results[[analysis_name]]$report),
              paste(analysis_name, "reports are not annotated"))
}
expect_true(identical(annotated$reports, results$reports), "Top-level reports are unchanged")
expect_true(identical(annotated$run_metadata, results$run_metadata), "Metadata is unchanged")
expect_true(identical(results, results_before), "Input result object is not mutated")
expect_true(identical(annotation, annotation_before), "Annotation input is not mutated")

source_effects <- results$analysis_change$all$outcome_effects
annotated_effects <- annotated$analysis_change$all$outcome_effects
expect_true(identical(annotated_effects$ANALYTE_NAME, source_effects$ANALYTE_NAME),
            "Outcome-effect row order is unchanged")
expect_true(nrow(annotated_effects) == nrow(source_effects),
            "Outcome-effect row count is unchanged")
expect_true(identical(annotated_effects$FU, source_effects$FU),
            "All outcome-effect FUs remain together and ordered")
expect_equal(annotated_effects$Target,
             c("Target_B", "Target_A", "Target_B", "Target_A"),
             "Outcome annotations follow exact AptName matches")
expect_equal(annotated_effects$FU_LABEL, c("3mo", "3mo", "24mo", "24mo"),
             "Outcome FU labels map by arbitrary FU value")
for (column in c("EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")) {
  expect_true(identical(annotated_effects[[column]], source_effects[[column]]),
              paste("Outcome source column is unchanged:", column))
}
expect_true(!"IGNORED_TRANSLATION_FIELD" %in% names(annotated_effects),
            "Non-contract annotation columns do not leak into results")

source_coefficients <- results$analysis_change$all$coefficients
annotated_coefficients <- annotated$analysis_change$all$coefficients
expect_true(identical(annotated_coefficients$ANALYTE_NAME,
                      source_coefficients$ANALYTE_NAME),
            "Coefficient row order is unchanged")
expect_true(nrow(annotated_coefficients) == nrow(source_coefficients),
            "Coefficient row count is unchanged")
expect_true(identical(annotated_coefficients$FU, source_coefficients$FU),
            "All coefficient FUs remain together and ordered")
expect_true(identical(annotated_coefficients$COEFFICIENT,
                      source_coefficients$COEFFICIENT),
            "Coefficient terms are unchanged")
expect_equal(annotated_coefficients$Target,
             rep(c("Target_B", "Target_A", "Target_B", "Target_A"), each = 2L),
             "Coefficient annotations repeat across model terms")
expect_equal(annotated_coefficients$FU_LABEL,
             rep(c("3mo", "3mo", "24mo", "24mo"), each = 2L),
             "Coefficient FU labels repeat across model terms")

without_labels <- get_annotated_results(results, annotation)
for (analysis_name in c("analysis_change", "analysis_level")) {
  for (group in c("all", "male")) {
    output <- without_labels[[analysis_name]][[group]]
    expect_true(!"FU_LABEL" %in% names(output$coefficients) &&
                  !"FU_LABEL" %in% names(output$outcome_effects),
                paste(analysis_name, group, "omits FU_LABEL by default"))
  }
}
duplicate_label_values <- get_annotated_results(
  results, annotation, fu_labels = c("1" = "visit", "5" = "visit")
)
expect_true(
  all(duplicate_label_values$analysis_change$all$outcome_effects$FU_LABEL == "visit"),
  "Duplicate FU label values are allowed as row metadata"
)

tte_results <- make_results(
  coefficients = make_coefficients(tte = TRUE, extra = TRUE),
  effects = make_effects(tte = TRUE, extra = TRUE)
)
tte <- get_annotated_results(tte_results, annotation, fu_labels)
expected_tte_outcome_columns <- c(
  track111_outcome_columns[1:12], "HAZARD_RATIO",
  track111_outcome_columns[13:15], "FUTURE_OUTCOME_FIELD"
)
expected_tte_coefficient_columns <- c(
  "ANALYTE_NAME", annotation_columns, "FU", "FU_LABEL", "COEFFICIENT",
  "N_OBS", "N_EVENTS", "EFFECT_SIZE", "HAZARD_RATIO", "SE", "P_VALUE",
  "BH_P_VALUE", "FUTURE_COEFFICIENT_FIELD"
)
expect_true(
  identical(names(tte$analysis_change$all$outcome_effects),
            expected_tte_outcome_columns),
  "TTE outcomes preserve hazard ratios and future fields"
)
expect_true(
  identical(names(tte$analysis_change$all$coefficients),
            expected_tte_coefficient_columns),
  "TTE coefficients preserve event counts, hazard ratios, and future fields"
)
expect_true(
  identical(tte$analysis_change$all$outcome_effects$HAZARD_RATIO,
            tte_results$analysis_change$all$outcome_effects$HAZARD_RATIO) &&
    identical(tte$analysis_change$all$coefficients$HAZARD_RATIO,
              tte_results$analysis_change$all$coefficients$HAZARD_RATIO),
  "TTE hazard-ratio values are unchanged in both tables"
)

partial <- expect_warning(
  get_annotated_results(results, make_annotation("A")),
  "unique (coefficients|outcome_effects) analytes were not found",
  "Partial annotation coverage warns"
)
partial_outcomes <- partial$analysis_change$all$outcome_effects
partial_coefficients <- partial$analysis_change$all$coefficients
expect_true(
  all(is.na(partial_outcomes$SomaId[partial_outcomes$ANALYTE_NAME == "B"])) &&
    all(is.na(partial_coefficients$SomaId[partial_coefficients$ANALYTE_NAME == "B"])),
  "Partially unmatched rows remain in both tables with NA annotations"
)
expect_error(
  get_annotated_results(results, make_annotation(c("X", "Y"))),
  "No unique coefficients analytes matched",
  "Zero annotation coverage errors"
)

per_fu_coefficients <- make_coefficients()
per_fu_coefficient_rows <- per_fu_coefficients$FU == 5L
per_fu_coefficients$ANALYTE_NAME[per_fu_coefficient_rows] <- paste0(
  "UNMATCHED_",
  per_fu_coefficients$ANALYTE_NAME[per_fu_coefficient_rows]
)
expect_error(
  get_annotated_results(make_results(coefficients = per_fu_coefficients),
                        make_annotation(c("A", "B")), fu_labels),
  "No unique coefficients analytes matched.*FU: 5",
  "Coefficient annotation rejects zero coverage within any FU"
)
per_fu_effects <- make_effects()
per_fu_effect_rows <- per_fu_effects$FU == 5L
per_fu_effects$ANALYTE_NAME[per_fu_effect_rows] <- paste0(
  "UNMATCHED_",
  per_fu_effects$ANALYTE_NAME[per_fu_effect_rows]
)
expect_error(
  get_annotated_results(make_results(effects = per_fu_effects),
                        make_annotation(c("A", "B")), fu_labels),
  "No unique outcome_effects analytes matched.*FU: 5",
  "Outcome annotation rejects zero coverage within any FU"
)

duplicate_annotation <- rbind(annotation, annotation[1L, ])
expect_error(get_annotated_results(results, duplicate_annotation),
             "AptName must be unique", "Duplicate annotation keys error")
blank_annotation <- annotation
blank_annotation$AptName[1L] <- " "
expect_error(get_annotated_results(results, blank_annotation),
             "missing or blank", "Blank annotation keys error")
missing_annotation_column <- annotation
missing_annotation_column$UniProt <- NULL
expect_error(get_annotated_results(results, missing_annotation_column),
             "missing required columns: UniProt", "Missing annotation columns error")
duplicate_annotation_columns <- annotation
names(duplicate_annotation_columns)[2L] <- "AptName"
expect_error(get_annotated_results(results, duplicate_annotation_columns),
             "column names must be unique", "Duplicate annotation column names error")
expect_error(get_annotated_results(results, "not a data frame"),
             "protein_annotation must be a data.frame",
             "Non-data-frame annotation errors")
expect_error(get_annotated_results(results, make_annotation(c("a", "b"))),
             "No unique coefficients analytes matched",
             "AptName matching is exact and case-sensitive")

reserved_coefficients <- results
reserved_coefficients$analysis_change$all$coefficients$SomaId <- "collision"
expect_error(get_annotated_results(reserved_coefficients, annotation),
             "coefficients contains reserved annotation columns: SomaId",
             "Reserved annotation columns in coefficients error")
reserved_outcomes <- results
reserved_outcomes$analysis_change$all$outcome_effects$FU_LABEL <- "collision"
expect_error(get_annotated_results(reserved_outcomes, annotation),
             "outcome_effects contains reserved annotation columns: FU_LABEL",
             "Reserved annotation columns in outcomes error")

duplicate_coefficient_columns <- results
names(duplicate_coefficient_columns$analysis_change$all$coefficients)[4L] <- "FU"
expect_error(get_annotated_results(duplicate_coefficient_columns, annotation),
             "coefficients column names must be unique",
             "Duplicate coefficient column names error")
duplicate_outcome_columns <- results
names(duplicate_outcome_columns$analysis_change$all$outcome_effects)[3L] <- "FU"
expect_error(get_annotated_results(duplicate_outcome_columns, annotation),
             "outcome_effects column names must be unique",
             "Duplicate outcome column names error")

duplicate_coefficient_keys <- results
duplicate_coefficient_keys$analysis_change$all$coefficients <- rbind(
  duplicate_coefficient_keys$analysis_change$all$coefficients,
  duplicate_coefficient_keys$analysis_change$all$coefficients[1L, ]
)
expect_error(get_annotated_results(duplicate_coefficient_keys, annotation),
             "unique ANALYTE_NAME x FU x COEFFICIENT",
             "Duplicate coefficient keys error")
duplicate_outcome_keys <- results
duplicate_outcome_keys$analysis_change$all$outcome_effects <- rbind(
  duplicate_outcome_keys$analysis_change$all$outcome_effects,
  duplicate_outcome_keys$analysis_change$all$outcome_effects[1L, ]
)
expect_error(get_annotated_results(duplicate_outcome_keys, annotation),
             "unique ANALYTE_NAME x FU rows", "Duplicate outcome keys error")

blank_result_ids <- results
blank_result_ids$analysis_change$all$coefficients$ANALYTE_NAME[1L] <- ""
expect_error(get_annotated_results(blank_result_ids, annotation),
             "must not contain missing or blank", "Blank result identifiers error")
invalid_fu_results <- results
invalid_fu_results$analysis_change$all$coefficients$FU[1L] <- 1.5
expect_error(get_annotated_results(invalid_fu_results, annotation),
             "positive integer", "Noninteger observed FU errors")
missing_coefficient_column <- results
missing_coefficient_column$analysis_change$all$coefficients$COEFFICIENT <- NULL
expect_error(get_annotated_results(missing_coefficient_column, annotation),
             "coefficients is missing required columns: COEFFICIENT",
             "Missing coefficient key column errors")
missing_coefficient_stat <- results
missing_coefficient_stat$analysis_change$all$coefficients$N_OBS <- NULL
expect_error(get_annotated_results(missing_coefficient_stat, annotation),
             "coefficients is missing required columns: N_OBS",
             "Missing coefficient contract column errors")
missing_outcome_column <- results
missing_outcome_column$analysis_change$all$outcome_effects$BH_P_VALUE <- NULL
expect_error(get_annotated_results(missing_outcome_column, annotation),
             "outcome_effects is missing required columns: BH_P_VALUE",
             "Missing outcome contract column errors")

normalized_duplicate_results <- results
normalized_duplicate_row <- normalized_duplicate_results$analysis_change$all$outcome_effects[1L, ]
normalized_duplicate_row$FU <- "01"
normalized_duplicate_results$analysis_change$all$outcome_effects$FU <- as.character(
  normalized_duplicate_results$analysis_change$all$outcome_effects$FU
)
normalized_duplicate_results$analysis_change$all$outcome_effects <- rbind(
  normalized_duplicate_results$analysis_change$all$outcome_effects,
  normalized_duplicate_row
)
expect_error(get_annotated_results(normalized_duplicate_results, annotation),
             "unique ANALYTE_NAME x FU rows",
             "Duplicate normalized FU keys error")

expect_error(get_annotated_results("not a result", annotation),
             "results must be a FAST_outcome_WAS result list",
             "Non-list result object errors")
expect_error(get_annotated_results(list(), annotation), "missing the required list",
             "Missing analysis structure errors")
missing_group <- results
missing_group$analysis_change <- missing_group$analysis_change[
  names(missing_group$analysis_change) != "male"
]
expect_error(get_annotated_results(missing_group, annotation), "missing group: male",
             "Missing group structure errors")
nonlist_group <- results
nonlist_group$analysis_change$all <- 1L
expect_error(get_annotated_results(nonlist_group, annotation), "must be a list or NULL",
             "Malformed non-NULL group errors")
missing_table <- results
missing_table$analysis_change$all <- missing_table$analysis_change$all[
  names(missing_table$analysis_change$all) != "outcome_effects"
]
expect_error(get_annotated_results(missing_table, annotation),
             "is missing outcome_effects", "Missing nested result table errors")

expect_error(get_annotated_results(results, annotation, c("1" = "3mo")),
             "missing labels for FU: 5", "Incomplete FU labels error")
expect_error(get_annotated_results(results, annotation, c("3mo", "24mo")),
             "nonblank name for every label", "Unnamed FU labels error")
duplicate_fu_names <- structure(c("first", "second"), names = c("01", "1"))
expect_error(get_annotated_results(results, annotation, duplicate_fu_names),
             "unique FU names", "Duplicate normalized FU names error")
expect_error(get_annotated_results(results, annotation, c("1" = "3mo", "5" = " ")),
             "only nonblank labels", "Blank FU label values error")

output_rds <- tempfile("annotated-results-", fileext = ".rds")
saveRDS(annotated, output_rds)
expect_true(identical(readRDS(output_rds), annotated),
            "Caller can save and reload the complete annotated result object")
unlink(output_rds)

cat("\nAll OutcomeWAS protein annotation tests passed.\n")
