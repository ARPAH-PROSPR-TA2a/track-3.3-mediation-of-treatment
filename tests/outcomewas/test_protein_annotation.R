source(file.path("outcomewas", "main.R"))


expect_true <- function(value, label) {
  if (!isTRUE(value)) stop("FAIL: ", label)
  cat("PASS ", label, "\n", sep = "")
}


expect_equal <- function(actual, expected, label) {
  if (!isTRUE(all.equal(actual, expected, check.attributes = TRUE))) {
    stop(
      "FAIL: ", label,
      "\nExpected: ", paste(expected, collapse = ", "),
      "\nActual: ", paste(actual, collapse = ", ")
    )
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


expect_rscript_success <- function(lines, label) {
  script <- tempfile("outcomewas-source-smoke-", fileext = ".R")
  writeLines(lines, script)
  output <- system2(
    file.path(R.home("bin"), "Rscript"),
    script,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  unlink(script)
  expect_true(
    identical(as.integer(status), 0L),
    paste0(label, if (length(output) > 0L) paste0(": ", paste(output, collapse = " | ")) else "")
  )
}


make_effects <- function(tte = FALSE, extra = FALSE) {
  effects <- data.frame(
    ANALYTE_NAME = c("B", "A", "B", "A"),
    FU = c(1L, 1L, 5L, 5L),
    EFFECT_SIZE = c(0.4, -0.2, 0.8, -0.6),
    SE = c(0.10, 0.20, 0.15, 0.25),
    P_VALUE = c(0.01, 0.20, 0.03, 0.40),
    BH_P_VALUE = c(0.02, 0.20, 0.06, 0.40),
    stringsAsFactors = FALSE
  )

  if (tte) {
    effects$HAZARD_RATIO <- exp(effects$EFFECT_SIZE)
    effects <- effects[
      c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "HAZARD_RATIO", "SE", "P_VALUE", "BH_P_VALUE")
    ]
  }
  if (extra) effects$FUTURE_FIELD <- c("w", "x", "y", "z")
  effects
}


make_results <- function(effects = make_effects()) {
  entry <- list(
    coefficients = data.frame(SENTINEL = 1L),
    outcome_effects = effects
  )
  list(
    analysis_change = list(all = entry, male = entry, female = NULL),
    analysis_level = list(all = entry, male = entry, female = NULL)
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


cat("OutcomeWAS protein annotation tests\n\n")

expect_rscript_success(
  c(
    'source(file.path("outcomewas", "main.R"), chdir = TRUE)',
    'stopifnot(exists("FAST_outcome_WAS"), exists("get_annotated_outcome_effects"))'
  ),
  "main.R sources from a relative path with chdir = TRUE"
)
main_path_literal <- paste(
  capture.output(dput(normalizePath(file.path("outcomewas", "main.R")))),
  collapse = ""
)
expect_rscript_success(
  c(
    paste0("source(", main_path_literal, ", chdir = TRUE)"),
    'stopifnot(exists("FAST_outcome_WAS"), exists("get_annotated_outcome_effects"))'
  ),
  "main.R sources from an absolute path with chdir = TRUE"
)

results <- make_results()
annotation <- make_annotation()
results_before <- unserialize(serialize(results, NULL))
annotation_before <- unserialize(serialize(annotation, NULL))
fu_labels <- c("1" = "3mo", "5" = "24mo")

annotated <- get_annotated_outcome_effects(
  results = results,
  protein_annotation = annotation,
  analysis_type = "change",
  group = "all",
  fu_labels = fu_labels
)

track111_columns <- c(
  "ANALYTE_NAME",
  "SomaId",
  "TargetFullName",
  "Target",
  "UniProt",
  "EntrezGeneID",
  "EntrezGeneSymbol",
  "Uniprot_Unique",
  "Symbol_Unique",
  "FU",
  "FU_LABEL",
  "EFFECT_SIZE",
  "SE",
  "P_VALUE",
  "BH_P_VALUE"
)

expect_true(
  identical(names(annotated), track111_columns),
  "Continuous annotated columns exactly match Track 1.1.1"
)
expect_true(
  identical(annotated$ANALYTE_NAME, results$analysis_change$all$outcome_effects$ANALYTE_NAME),
  "Annotation preserves result row order"
)
expect_true(nrow(annotated) == 4L, "Annotation preserves result row count")
expect_equal(annotated$Target, c("Target_B", "Target_A", "Target_B", "Target_A"),
             "Annotation values follow exact AptName matches")
expect_equal(annotated$FU_LABEL, c("3mo", "3mo", "24mo", "24mo"),
             "Custom FU labels map by FU")
for (column in c("FU", "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")) {
  expect_true(
    identical(annotated[[column]], results$analysis_change$all$outcome_effects[[column]]),
    paste("Source column is unchanged:", column)
  )
}
expect_true(identical(results, results_before), "Input result object is not mutated")
expect_true(identical(annotation, annotation_before), "Annotation input is not mutated")

without_labels <- get_annotated_outcome_effects(results, annotation)
expect_true(!"FU_LABEL" %in% names(without_labels), "FU_LABEL is omitted by default")

fu5 <- get_annotated_outcome_effects(results, annotation, fu = 5L, fu_labels = fu_labels)
expect_true(nrow(fu5) == 2L && all(fu5$FU == 5L), "FU filtering accepts values beyond four")
expect_true(
  identical(fu5$ANALYTE_NAME, c("B", "A")),
  "FU filtering preserves relative row order"
)

tte_results <- make_results(make_effects(tte = TRUE, extra = TRUE))
tte <- get_annotated_outcome_effects(tte_results, annotation, fu_labels = fu_labels)
expected_tte_columns <- append(track111_columns, "HAZARD_RATIO", after = 12L)
expected_tte_columns <- c(expected_tte_columns, "FUTURE_FIELD")
expect_true(
  identical(names(tte), expected_tte_columns),
  "TTE and future columns are preserved dynamically"
)
expect_true(
  identical(tte$HAZARD_RATIO, tte_results$analysis_change$all$outcome_effects$HAZARD_RATIO),
  "Hazard ratios are unchanged"
)

partial_annotation <- make_annotation("A")
partial <- expect_warning(
  get_annotated_outcome_effects(results, partial_annotation),
  "unique outcome analytes were not found",
  "Partial annotation coverage warns"
)
expect_true(
  all(is.na(partial$SomaId[partial$ANALYTE_NAME == "B"])),
  "Unmatched result rows are retained with NA annotation"
)

expect_error(
  get_annotated_outcome_effects(results, make_annotation(c("X", "Y"))),
  "No unique outcome analytes matched",
  "Zero annotation coverage errors"
)

duplicate_annotation <- rbind(annotation, annotation[1L, ])
expect_error(
  get_annotated_outcome_effects(results, duplicate_annotation),
  "AptName must be unique",
  "Duplicate annotation keys error"
)

blank_annotation <- annotation
blank_annotation$AptName[1L] <- " "
expect_error(
  get_annotated_outcome_effects(results, blank_annotation),
  "missing or blank",
  "Blank annotation keys error"
)

missing_annotation_column <- annotation
missing_annotation_column$UniProt <- NULL
expect_error(
  get_annotated_outcome_effects(results, missing_annotation_column),
  "missing required columns: UniProt",
  "Missing annotation columns error"
)

duplicate_results <- results
duplicate_results$analysis_change$all$outcome_effects <- rbind(
  duplicate_results$analysis_change$all$outcome_effects,
  duplicate_results$analysis_change$all$outcome_effects[1L, ]
)
expect_error(
  get_annotated_outcome_effects(duplicate_results, annotation),
  "unique ANALYTE_NAME x FU",
  "Duplicate OutcomeWAS keys error"
)

blank_results <- results
blank_results$analysis_change$all$outcome_effects$ANALYTE_NAME[1L] <- ""
expect_error(
  get_annotated_outcome_effects(blank_results, annotation),
  "must not contain missing or blank",
  "Blank result identifiers error"
)

expect_error(
  get_annotated_outcome_effects(list(), annotation),
  "missing the required list",
  "Malformed result structure errors"
)
expect_error(
  get_annotated_outcome_effects(results, annotation, fu = 2L),
  "Requested fu values are not present",
  "Absent requested FU errors"
)
expect_error(
  get_annotated_outcome_effects(results, annotation, fu = 1.5),
  "positive integer",
  "Noninteger requested FU errors"
)
expect_error(
  get_annotated_outcome_effects(results, annotation, fu_labels = c("1" = "3mo")),
  "missing labels for FU: 5",
  "Incomplete FU labels error"
)
expect_true(
  is.null(get_annotated_outcome_effects(results, annotation, group = "female")),
  "Unavailable strata return NULL"
)

output_rds <- tempfile("annotated-outcome-effects-", fileext = ".rds")
tables <- write_annotated_outcome_effect_tables(
  results,
  annotation,
  output_rds,
  fu_labels = fu_labels
)
expected_table_names <- c(
  "change_all_3mo", "change_all_24mo",
  "change_male_3mo", "change_male_24mo",
  "level_all_3mo", "level_all_24mo",
  "level_male_3mo", "level_male_24mo"
)
expect_true(identical(names(tables), expected_table_names), "Writer uses labeled table names")
expect_true(identical(readRDS(output_rds), tables), "Writer saves the returned table list")
expect_true(
  all(vapply(tables, function(table) identical(names(table), track111_columns), logical(1))),
  "Every continuous exported table uses the Track 1.1.1 schema"
)
expect_true(
  !any(grepl("female", names(tables))),
  "Writer skips NULL strata"
)

default_output_rds <- tempfile("annotated-outcome-effects-default-", fileext = ".rds")
default_tables <- write_annotated_outcome_effect_tables(
  results,
  annotation,
  default_output_rds
)
expect_true(
  all(grepl("_(FU1|FU5)$", names(default_tables))),
  "Writer uses FU<n> names when labels are absent"
)
expect_true(
  all(vapply(default_tables, function(table) !"FU_LABEL" %in% names(table), logical(1))),
  "Default exported tables omit FU_LABEL"
)

per_fu_results <- make_results(data.frame(
  ANALYTE_NAME = c("A", "UNMATCHED"),
  FU = c(1L, 5L),
  EFFECT_SIZE = c(0.1, 0.2),
  SE = c(0.01, 0.02),
  P_VALUE = c(0.1, 0.2),
  BH_P_VALUE = c(0.1, 0.2),
  stringsAsFactors = FALSE
))
per_fu_output_rds <- tempfile("annotated-outcome-effects-per-fu-", fileext = ".rds")
expect_error(
  write_annotated_outcome_effect_tables(
    per_fu_results,
    make_annotation("A"),
    per_fu_output_rds,
    fu_labels = fu_labels
  ),
  "No unique outcome analytes matched",
  "Writer rejects zero annotation coverage within an exported FU table"
)
expect_true(
  !file.exists(per_fu_output_rds),
  "Writer does not save a partial table list after per-FU coverage failure"
)

unlink(c(output_rds, default_output_rds, per_fu_output_rds))
cat("\nAll OutcomeWAS protein annotation tests passed.\n")
