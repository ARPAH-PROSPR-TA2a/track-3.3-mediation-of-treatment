source(file.path("outcomewas", "main.R"))

# Freeze the original helper so this test does not depend on Git at runtime.
reference_omics_report <- function(sample_ids, omics_df) {
  analyte_names <- omics_df$ANALYTE_NAME
  omics_numeric <- omics_df[, setdiff(names(omics_df), "ANALYTE_NAME"), drop = FALSE]
  omics_numeric <- omics_numeric[, colnames(omics_numeric) %in% sample_ids, drop = FALSE]
  report <- data.frame(
    ANALYTE_NAME = analyte_names, N_NONMISSING = NA_integer_, MEAN = NA_real_,
    MEDIAN = NA_real_, SD = NA_real_, MIN = NA_real_, MAX = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(analyte_names)) {
    analyte_values <- as.numeric(omics_numeric[i, ])
    report$N_NONMISSING[i] <- sum(!is.na(analyte_values))
    report$MEAN[i] <- mean(analyte_values, na.rm = TRUE)
    report$MEDIAN[i] <- median(analyte_values, na.rm = TRUE)
    report$SD[i] <- sd(analyte_values, na.rm = TRUE)
    report$MIN[i] <- min(analyte_values, na.rm = TRUE)
    report$MAX[i] <- max(analyte_values, na.rm = TRUE)
  }
  report
}

capture_report <- function(fun, sample_ids, omics) {
  warnings <- character()
  value <- withCallingHandlers(fun(sample_ids, omics), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

check_report <- function(label, sample_ids, omics) {
  expected <- capture_report(reference_omics_report, sample_ids, omics)
  actual <- capture_report(.create_omics_data_report, sample_ids, omics)
  # identical checks numeric types, NA versus NaN, attributes, and row/column order.
  stopifnot(identical(actual$value, expected$value),
            identical(actual$warnings, expected$warnings))
  cat("PASS ", label, ": exact report and ordered warning equivalence\n", sep = "")
  invisible(actual)
}

edge_values <- rbind(
  c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_),
  rep(NaN, 6L), rep(7, 6L), rep(Inf, 6L), rep(-Inf, 6L),
  c(Inf, -Inf, 3, NA, NaN, 5), c(NA, NA, 8, NaN, NA, NA),
  c(-1e16, 1, 1e16, 2, 3, 4), c(0.1, 0.7, NA, 0.3, NaN, 0.9)
)
colnames(edge_values) <- c("sample-C", "sample-A", "sample-F", "sample-B", "sample-E", "sample-D")
omics <- data.frame(ANALYTE_NAME = c("z", "a", "z", "inf", "negative-inf",
                                     "mixed", "singleton", "cancellation", "betas"),
                    edge_values, check.names = FALSE)
rownames(omics) <- paste0("input-row-", rev(seq_len(nrow(omics))))
check_report("double and exceptional values", rev(colnames(edge_values)), omics)
check_report("shuffled, repeated, and unknown sample IDs",
             c("sample-D", "unknown", "sample-C", "sample-D", NA, "ANALYTE_NAME"), omics)
check_report("singleton sample selection", "sample-F", omics)
empty <- check_report("no matching samples", c("unknown", NA), omics)
stopifnot(length(empty$warnings) == 2L * nrow(omics))
check_report("empty sample ID vector", character(), omics)
check_report("one analyte", colnames(edge_values), omics[3L, , drop = FALSE])
factor_omics <- omics
factor_omics$ANALYTE_NAME <- factor(factor_omics$ANALYTE_NAME,
                                   levels = rev(unique(factor_omics$ANALYTE_NAME)))
check_report("factor analyte identifiers", colnames(edge_values), factor_omics)

integer_omics <- data.frame(ANALYTE_NAME = c("third", "first", "second"),
                            C = c(1L, NA_integer_, 8L), A = c(3L, NA_integer_, 8L),
                            B = c(5L, NA_integer_, 8L))
check_report("integer sample columns", c("B", "C", "A"), integer_omics)
check_report("integer singleton selection", "A", integer_omics)
check_report("integer columns without matching samples", "unknown", integer_omics)
integer_omics$B <- as.double(integer_omics$B) + 0.25
check_report("mixed integer and double columns", c("A", "B", "C"), integer_omics)

set.seed(20260913)
chunk_values <- matrix(rnorm(2007L * 5L), nrow = 2007L)
chunk_values[2000L, ] <- NA_real_
chunk_values[2001L, ] <- NaN
chunk_values[2007L, ] <- 0
colnames(chunk_values) <- paste0("sample-", c(4L, 2L, 5L, 1L, 3L))
chunk_omics <- data.frame(ANALYTE_NAME = paste0("probe-", sample.int(2007L)),
                          chunk_values, check.names = FALSE)
chunk_report <- check_report("2007 analytes across the matrix chunk boundary",
                              colnames(chunk_values), chunk_omics)
stopifnot(length(chunk_report$warnings) == 4L)

# Exercise the public reports entry point without fitting any models.
subject <- rep(seq_len(8L), each = 3L)
pheno <- data.frame(
  SAMPLE_ID = paste0("sample-", seq_along(subject)), SUBJECT_ID = paste0("subject-", subject),
  FU = factor(rep(0:2, 8L)), FEMALE = factor((subject - 1L) %% 2L),
  TREATMENT_GROUP = factor(((subject - 1L) %/% 2L) %% 2L), OUTCOME = subject / 3
)
public_values <- matrix(rnorm(5L * nrow(pheno)), nrow = 5L,
                         dimnames = list(NULL, rev(pheno$SAMPLE_ID)))
public_omics <- data.frame(ANALYTE_NAME = paste0("probe-", c(3L, 1L, 5L, 2L, 4L)),
                           public_values, check.names = FALSE)
reports <- FAST_outcome_WAS_reports(pheno, public_omics, omics_type = "other")
for (stratum in c("all", "male", "female")) {
  subset <- pheno
  if (stratum != "all") subset <- subset[subset$FEMALE == as.integer(stratum == "female"), ]
  for (fu in 0:2) {
    for (tx in 0:1) {
      sample_ids <- subset$SAMPLE_ID[subset$FU == fu & subset$TREATMENT_GROUP == tx]
      expected <- reference_omics_report(sample_ids, public_omics)
      actual <- reports$variable_summaries[[stratum]][[paste0("omics_FU", fu, "_Tx", tx)]]
      stopifnot(identical(actual, expected))
    }
  }
}
cat("PASS public reports: exact summaries for every stratum, follow-up, and treatment cell\n")
