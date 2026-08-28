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
    paste0(
      'stopifnot(exists("FAST_outcome_WAS"), ',
      'exists("get_annotated_outcome_effects"), ',
      '!exists("write_annotated_outcome_effect_tables"))'
    )
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
    paste0(
      'stopifnot(exists("FAST_outcome_WAS"), ',
      'exists("get_annotated_outcome_effects"), ',
      '!exists("write_annotated_outcome_effect_tables"))'
    )
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

expected_table_names <- c(
  "change_all_3mo", "change_all_24mo",
  "change_male_3mo", "change_male_24mo",
  "level_all_3mo", "level_all_24mo",
  "level_male_3mo", "level_male_24mo"
)

expect_true(is.list(annotated), "Annotation returns a list of tables")
expect_true(
  identical(names(annotated), expected_table_names),
  "Annotation traverses analyses, non-NULL strata, and observed FU values"
)
expect_true(
  all(vapply(
    annotated,
    function(table) identical(names(table), track111_columns),
    logical(1)
  )),
  "Every continuous table exactly matches the Track 1.1.1 schema"
)
expect_true(
  !any(grepl("female", names(annotated))),
  "Annotation skips NULL strata"
)

change_all_3mo <- annotated$change_all_3mo
source_3mo <- results$analysis_change$all$outcome_effects
source_3mo <- source_3mo[source_3mo$FU == 1L, , drop = FALSE]
expect_true(
  identical(change_all_3mo$ANALYTE_NAME, source_3mo$ANALYTE_NAME),
  "Each FU table preserves source row order"
)
expect_true(
  nrow(change_all_3mo) == nrow(source_3mo),
  "Each FU table preserves source row count"
)
expect_equal(
  change_all_3mo$Target,
  c("Target_B", "Target_A"),
  "Annotation values follow exact AptName matches"
)
expect_equal(
  change_all_3mo$FU_LABEL,
  c("3mo", "3mo"),
  "Custom FU labels map by FU"
)
for (column in c("FU", "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")) {
  expect_true(
    identical(change_all_3mo[[column]], source_3mo[[column]]),
    paste("Source column is unchanged:", column)
  )
}
expect_true(identical(results, results_before), "Input result object is not mutated")
expect_true(identical(annotation, annotation_before), "Annotation input is not mutated")

without_labels <- get_annotated_outcome_effects(results, annotation)
expect_true(
  all(grepl("_(FU1|FU5)$", names(without_labels))),
  "Default table names support arbitrary observed FU values"
)
expect_true(
  all(vapply(
    without_labels,
    function(table) !"FU_LABEL" %in% names(table),
    logical(1)
  )),
  "FU_LABEL is omitted from every table by default"
)

tte_results <- make_results(make_effects(tte = TRUE, extra = TRUE))
tte <- get_annotated_outcome_effects(tte_results, annotation, fu_labels = fu_labels)
expected_tte_columns <- append(track111_columns, "HAZARD_RATIO", after = 12L)
expected_tte_columns <- c(expected_tte_columns, "FUTURE_FIELD")
expect_true(
  all(vapply(
    tte,
    function(table) identical(names(table), expected_tte_columns),
    logical(1)
  )),
  "Every TTE table preserves hazard ratios and future columns dynamically"
)
expect_true(
  identical(
    tte$change_all_3mo$HAZARD_RATIO,
    tte_results$analysis_change$all$outcome_effects$HAZARD_RATIO[1:2]
  ),
  "Hazard ratios are unchanged"
)

partial_annotation <- make_annotation("A")
partial <- expect_warning(
  get_annotated_outcome_effects(results, partial_annotation),
  "unique outcome analytes were not found",
  "Partial annotation coverage warns"
)
expect_true(
  all(vapply(
    partial,
    function(table) all(is.na(table$SomaId[table$ANALYTE_NAME == "B"])),
    logical(1)
  )),
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

duplicate_annotation_columns <- annotation
names(duplicate_annotation_columns)[2L] <- "AptName"
expect_error(
  get_annotated_outcome_effects(results, duplicate_annotation_columns),
  "column names must be unique",
  "Duplicate annotation column names error"
)

reserved_column_results <- results
reserved_column_results$analysis_change$all$outcome_effects$UniProt <- "collision"
expect_error(
  get_annotated_outcome_effects(reserved_column_results, annotation),
  "reserved annotation columns: UniProt",
  "Reserved annotation columns in outcome effects error"
)

duplicate_result_columns <- results
names(duplicate_result_columns$analysis_change$all$outcome_effects)[3L] <- "FU"
expect_error(
  get_annotated_outcome_effects(duplicate_result_columns, annotation),
  "column names must be unique",
  "Duplicate outcome-effect column names error"
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

invalid_fu_results <- results
invalid_fu_results$analysis_change$all$outcome_effects$FU[1L] <- 1.5
expect_error(
  get_annotated_outcome_effects(invalid_fu_results, annotation),
  "positive integer",
  "Noninteger observed FU errors"
)

expect_error(
  get_annotated_outcome_effects(results, annotation, fu_labels = c("1" = "3mo")),
  "missing labels for FU: 5",
  "Incomplete FU labels error"
)

output_rds <- tempfile("annotated-outcome-effects-", fileext = ".rds")
saveRDS(annotated, output_rds)
expect_true(
  identical(readRDS(output_rds), annotated),
  "Caller can save and reload the returned table list"
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
expect_error(
  get_annotated_outcome_effects(
    per_fu_results,
    make_annotation("A"),
    fu_labels = fu_labels
  ),
  "No unique outcome analytes matched",
  "Annotation rejects zero coverage within any FU table"
)

unlink(output_rds)
cat("\nAll OutcomeWAS protein annotation tests passed.\n")
