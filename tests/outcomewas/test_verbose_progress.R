source(file.path("outcomewas", "main.R"))


expect_true <- function(condition, label) {
  if (!isTRUE(condition)) {
    stop("FAIL ", label, call. = FALSE)
  }
  cat("PASS ", label, "\n", sep = "")
}


capture_conditions <- function(expression) {
  messages <- character()
  warnings <- character()

  value <- withCallingHandlers(
    force(expression),
    message = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleMessage")
    },
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, messages = messages, warnings = warnings)
}


has_progress_line <- function(messages) {
  any(startsWith(messages, "[3.3] "))
}


has_message <- function(messages, text) {
  any(grepl(text, messages, fixed = TRUE))
}


make_test_data <- function() {
  set.seed(3303)

  n_subjects <- 20L
  subject_ids <- sprintf("subject_%02d", seq_len(n_subjects))
  female <- rep(0:1, each = n_subjects / 2L)
  treatment <- rep(rep(0:1, each = n_subjects / 4L), 2L)
  outcome <- rnorm(n_subjects)

  pheno <- data.frame(
    SAMPLE_ID = as.vector(rbind(
      paste0(subject_ids, "_baseline"),
      paste0(subject_ids, "_followup")
    )),
    FU = factor(rep(0:1, times = n_subjects)),
    SUBJECT_ID = rep(subject_ids, each = 2L),
    FEMALE = factor(rep(female, each = 2L)),
    TREATMENT_GROUP = factor(rep(treatment, each = 2L)),
    OUTCOME = rep(outcome, each = 2L),
    stringsAsFactors = FALSE
  )

  omics_values <- rbind(
    matrix(rnorm(3L * nrow(pheno)), nrow = 3L),
    rep(1, nrow(pheno))
  )
  colnames(omics_values) <- pheno$SAMPLE_ID
  omics <- data.frame(
    ANALYTE_NAME = c("analyte_1", "analyte_2", "analyte_3", "constant_failure"),
    omics_values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  list(pheno = pheno, omics = omics)
}


run_verbose_progress_tests <- function() {
  original_progress_log <- getOption("track33.progress_log")
  work_dir <- tempfile("track33-verbose-progress-")
  dir.create(work_dir)
  on.exit(options(track33.progress_log = original_progress_log), add = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  test_data <- make_test_data()
  analysis_args <- list(
    pheno = test_data$pheno,
    omics = test_data$omics,
    omics_type = "Proteomics",
    n_cores = 2L,
    checkpoint_batch_size = 2L
  )

  default_run <- capture_conditions(do.call(
    FAST_outcome_WAS,
    c(analysis_args, list(checkpoint_dir = file.path(work_dir, "default-checkpoints")))
  ))
  quiet_run <- capture_conditions(do.call(
    FAST_outcome_WAS,
    c(
      analysis_args,
      list(
        checkpoint_dir = file.path(work_dir, "quiet-checkpoints"),
        verbose = FALSE
      )
    )
  ))

  expect_true(!has_progress_line(default_run$messages),
              "Default analysis emits no [3.3] progress lines")
  expect_true(!has_progress_line(quiet_run$messages),
              "verbose = FALSE analysis emits no [3.3] progress lines")
  expect_true(identical(default_run$value, quiet_run$value),
              "Default and verbose = FALSE analysis results are identical")

  progress_log <- file.path(work_dir, "progress.log")
  options(track33.progress_log = progress_log)
  verbose_checkpoint_dir <- file.path(work_dir, "verbose-checkpoints")

  verbose_run <- capture_conditions(do.call(
    FAST_outcome_WAS,
    c(
      analysis_args,
      list(checkpoint_dir = verbose_checkpoint_dir, verbose = TRUE)
    )
  ))

  expect_true(identical(default_run$value, verbose_run$value),
              "Verbose analysis preserves results")
  expect_true(
    has_message(verbose_run$messages, "[3.3] analysis validating inputs and preparing data (serial)") &&
      has_message(verbose_run$messages, "[3.3] omics quality checks: 0/4 analytes (serial)"),
    "Verbose analysis logs validation"
  )
  expect_true(
    has_message(verbose_run$messages, "[3.3] change/all starting") &&
      has_message(verbose_run$messages, "[3.3] change/all: 1 follow-up level to process") &&
      has_message(verbose_run$messages, "[3.3] change/all/FU1 starting") &&
      has_message(verbose_run$messages, "[3.3] level/female/FU1 starting"),
    "Verbose analysis logs response, stratum, and follow-up"
  )
  expect_true(
    has_message(verbose_run$messages, "2 batches (0 cached, 2 pending)") &&
      has_message(verbose_run$messages, "dispatching 2 analytes across up to 2 workers"),
    "Verbose analysis logs checkpoint status and dispatch"
  )

  pid_lines <- verbose_run$messages[
    grepl("worker PIDs observed: ", verbose_run$messages, fixed = TRUE)
  ]
  pid_payloads <- trimws(sub("^.*worker PIDs observed: ", "", pid_lines))
  valid_pid_payloads <- grepl("^[0-9]+(, [0-9]+)*$", pid_payloads)
  parsed_pid_counts <- vapply(
    pid_payloads[valid_pid_payloads],
    function(payload) length(unique(strsplit(payload, ", ", fixed = TRUE)[[1L]])),
    integer(1)
  )
  expect_true(
    length(pid_lines) > 0L && all(valid_pid_payloads) && any(parsed_pid_counts >= 2L),
    "Verbose analysis logs exact worker PIDs from both workers"
  )
  expect_true(
    has_message(verbose_run$messages, "1 newly attempted analyte failed") &&
      has_message(verbose_run$messages, "model batches complete; applying BH correction (serial)"),
    "Verbose analysis logs failures and BH correction"
  )
  expect_true(has_message(verbose_run$messages, "[3.3] analysis complete ("),
              "Verbose analysis logs completion")

  cached_run <- capture_conditions(do.call(
    FAST_outcome_WAS,
    c(
      analysis_args,
      list(checkpoint_dir = verbose_checkpoint_dir, verbose = TRUE)
    )
  ))
  expect_true(identical(verbose_run$value, cached_run$value),
              "Cached analysis preserves results")
  expect_true(
    has_message(cached_run$messages, "2 batches (2 cached, 0 pending)") &&
      has_message(cached_run$messages, "all batches cached; no workers launched"),
    "Second analysis reports fully cached execution"
  )
  expect_true(
    !has_message(cached_run$messages, "dispatching") &&
      !has_message(cached_run$messages, "worker PIDs observed"),
    "Fully cached analysis does not report worker dispatch"
  )

  report_args <- list(
    pheno = test_data$pheno,
    omics = test_data$omics,
    omics_type = "Proteomics"
  )
  default_reports <- capture_conditions(do.call(FAST_outcome_WAS_reports, report_args))
  quiet_reports <- capture_conditions(do.call(
    FAST_outcome_WAS_reports,
    c(report_args, list(verbose = FALSE))
  ))
  verbose_reports <- capture_conditions(do.call(
    FAST_outcome_WAS_reports,
    c(report_args, list(verbose = TRUE))
  ))

  expect_true(
    !has_progress_line(default_reports$messages) &&
      !has_progress_line(quiet_reports$messages),
    "Default and verbose = FALSE reports emit no [3.3] progress lines"
  )
  expect_true(
    identical(default_reports$value, quiet_reports$value) &&
      identical(default_reports$value, verbose_reports$value),
    "Verbose reporting preserves results"
  )
  expect_true(
    has_message(verbose_reports$messages, "[3.3] reports validating inputs (serial)") &&
      has_message(verbose_reports$messages, "[3.3] reports variable summaries/all starting (serial)") &&
      has_message(verbose_reports$messages, "[3.3] reports outcome summary (serial)") &&
      has_message(verbose_reports$messages, "[3.3] reports complete ("),
    "Verbose reporting logs major stages and completion"
  )

  expect_true(file.exists(progress_log), "Progress log file is written")
  progress_lines <- readLines(progress_log, warn = FALSE)
  timestamp_pattern <- paste0(
    "^\\[3\\.3\\] [0-9]{4}-[0-9]{2}-[0-9]{2} ",
    "[0-9]{2}:[0-9]{2}:[0-9]{2} \\| "
  )
  expect_true(
    length(progress_lines) > 0L && all(grepl(timestamp_pattern, progress_lines)),
    "Progress log contains timestamped [3.3] lines"
  )
  expect_true(
    has_message(progress_lines, "analysis complete (") &&
      has_message(progress_lines, "reports complete ("),
    "Progress log records analysis and report completion"
  )

  cat("\nAll OutcomeWAS verbose progress tests passed.\n")
}


run_verbose_progress_tests()
