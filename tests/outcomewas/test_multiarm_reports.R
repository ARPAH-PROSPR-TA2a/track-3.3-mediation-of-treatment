source(file.path("outcomewas", "main.R"))
source(file.path("outcomewas", "R", "plotting_helpers.R"))
source(file.path("tests", "outcomewas", "synthetic_multiarm.R"))

run_multiarm_report_tests <- function() {
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  fixture <- .make_multiarm_outcome_fixture(n_per_cell = 4L, n_analytes = 2L)
  reports <- quiet(FAST_outcome_WAS_reports(
    fixture$pheno, fixture$omics, omics_type = "other",
    treatment_labels = fixture$treatment_labels
  ))
  counts <- reports$arm_counts
  stopifnot(
    identical(names(counts), c("FU", "FEMALE", "ARM", "TREATMENT", "IS_CONTROL",
                              "N_SUBJECTS", "N_SAMPLES")),
    nrow(counts) == 18L,
    all(counts$N_SUBJECTS == 4L), all(counts$N_SAMPLES == 4L),
    identical(counts$TREATMENT, unname(fixture$treatment_labels[counts$ARM])),
    identical(counts$IS_CONTROL, counts$ARM == "0"),
    all(reports$pheno_summary$N_CONTROL == 4L),
    all(reports$pheno_summary$N_TREATMENT == 8L)
  )
  summary <- reports$outcome_reports$outcome_summary
  stopifnot(
    identical(summary$GROUP, c("all", "arm_0", "arm_2", "arm_7")),
    identical(summary$ARM, c(NA_character_, "0", "2", "7")),
    identical(summary$TREATMENT[-1L], unname(fixture$treatment_labels)),
    identical(summary$N_SUBJECTS, c(24L, 8L, 8L, 8L))
  )
  for (stratum in c("all", "male", "female")) {
    expected_keys <- as.vector(outer(0:2, c("0", "2", "7"),
                                    function(fu, arm) paste0("omics_FU", fu, "_Tx", arm)))
    stopifnot(setequal(names(reports$variable_summaries[[stratum]]), expected_keys))
  }
  cat("PASS multi-arm report counts, labels, and variable-summary keys\n")

  missing_cell <- fixture$pheno[!(fixture$pheno$FEMALE == 1 &
                                  fixture$pheno$FU == 2 &
                                  fixture$pheno$TREATMENT_GROUP == 7), ]
  missing_reports <- quiet(FAST_outcome_WAS_reports(
    missing_cell, fixture$omics, omics_type = "other",
    treatment_labels = fixture$treatment_labels
  ))
  absent <- with(missing_reports$arm_counts, FU == 2 & FEMALE == 1 & ARM == "7")
  stopifnot(sum(absent) == 1L,
            missing_reports$arm_counts$N_SUBJECTS[absent] == 0L,
            missing_reports$arm_counts$N_SAMPLES[absent] == 0L,
            !"omics_FU2_Tx7" %in% names(missing_reports$variable_summaries$female))
  cat("PASS absent arm cells retain labeled zero counts\n")

  binary <- .make_multiarm_outcome_fixture(arms = c("0", "1"),
                                          n_per_cell = 4L, n_analytes = 2L)
  binary_reports <- quiet(FAST_outcome_WAS_reports(binary$pheno, binary$omics,
                                                  omics_type = "other"))
  binary_summary <- binary_reports$outcome_reports$outcome_summary
  stopifnot(identical(binary_summary$GROUP, c("all", "control", "treatment")),
            identical(binary_summary$N_SUBJECTS, c(16L, 8L, 8L)),
            all(binary_reports$pheno_summary$N_TREATMENT == 4L),
            all(binary_reports$pheno_summary$N_CONTROL == 4L),
            identical(binary_summary$TREATMENT, c(NA_character_, "Control", "Treatment")))
  subject_data <- binary$pheno[binary$pheno$FU == 0, ]
  expected_means <- c(mean(subject_data$OUTCOME),
                      mean(subject_data$OUTCOME[subject_data$TREATMENT_GROUP == 0]),
                      mean(subject_data$OUTCOME[subject_data$TREATMENT_GROUP == 1]))
  stopifnot(isTRUE(all.equal(binary_summary$MEAN, expected_means)))
  cat("PASS ordinary binary report groups and values remain familiar\n")

  tte <- .make_multiarm_outcome_fixture(outcome_type = "tte",
                                       n_per_cell = 4L, n_analytes = 2L)
  tte_reports <- quiet(FAST_outcome_WAS_reports(
    tte$pheno, tte$omics, omics_type = "other", treatment_labels = tte$treatment_labels
  ))
  tte_summary <- tte_reports$outcome_reports$outcome_summary
  subjects <- tte$pheno[tte$pheno$FU == 0, ]
  stopifnot(identical(tte_summary$GROUP, summary$GROUP),
            identical(tte_summary$ARM, summary$ARM),
            tte_summary$N_EVENTS[1L] == sum(subjects$OUTCOME_STATUS),
            sum(tte_summary$N_EVENTS[-1L]) == tte_summary$N_EVENTS[1L],
            sum(tte_summary$N_CENSORED[-1L]) == tte_summary$N_CENSORED[1L])
  treatments <- .prepare_treatments(tte$pheno, tte$treatment_labels)$treatments
  empty_arm <- .create_tte_outcome_report(tte$pheno[tte$pheno$TREATMENT_GROUP != 7, ],
                                         treatments)
  stopifnot(empty_arm$N_SUBJECTS[empty_arm$ARM %in% "7"] == 0L,
            is.na(empty_arm$TIME_MEAN[empty_arm$ARM %in% "7"]))
  cat("PASS survival summaries include every labeled arm and empty groups\n")

  effects <- data.frame(
    ANALYTE_NAME = c("estimated", "failed_filtered", "excluded"), FU = 1L,
    EFFECT_SIZE = c(0.2, NA_real_, -0.3), SE = c(0.1, NA_real_, 0.1),
    P_VALUE = c(0.01, NA_real_, 0.02), BH_P_VALUE = c(0.02, NA_real_, 0.02),
    BH_P_VALUE_FILTERED = c(NA_real_, NA_real_, 0.02),
    IN_FILTERED_SET = c(FALSE, TRUE, FALSE),
    STATUS = c("ok", "fit_error", "ok"), stringsAsFactors = FALSE
  )
  unavailable <- effects[2L, ]
  unavailable$FU <- 2L
  unavailable$STATUS <- "not_estimable"
  results <- list(all = list(outcome_effects = effects),
                  male = list(outcome_effects = unavailable), female = NULL)
  stopifnot(identical(.outcome_filtered_membership(effects), c(FALSE, TRUE, FALSE)))
  filtered <- .outcome_plot_rows(results, "all", 1L, "filtered")
  stopifnot(identical(filtered$ANALYTE_NAME, "failed_filtered"),
            nrow(filtered) == 1L, !any(.outcome_plot_estimable(filtered)))
  failed_plot <- plot_volcano(results, probe_set = "filtered")
  stopifnot(identical(ggplot_build(failed_plot)$data[[1L]]$label, "No estimable results"))
  no_members <- results
  no_members$all$outcome_effects$IN_FILTERED_SET <- FALSE
  no_members_plot <- plot_volcano(no_members, probe_set = "filtered")
  stopifnot(identical(ggplot_build(no_members_plot)$data[[1L]]$label,
                      "No probes in filtered set"))
  legacy <- effects
  legacy$IN_FILTERED_SET <- NULL
  stopifnot(identical(.outcome_filtered_membership(legacy), c(FALSE, FALSE, TRUE)))
  legacy$BH_P_VALUE_FILTERED <- c(Inf, -0.1, 1.1)
  stopifnot(!any(.outcome_filtered_membership(legacy)))
  cat("PASS filtered membership remains independent of inference, with legacy fallback\n")

  work_dir <- tempfile("outcomewas-multiarm-plots-")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  original_qq <- plot_qq
  original_volcano <- plot_volcano
  on.exit(assign("plot_qq", original_qq, envir = .GlobalEnv), add = TRUE)
  on.exit(assign("plot_volcano", original_volcano, envir = .GlobalEnv), add = TRUE)
  calls <- character()
  assign("plot_qq", function(results, stratum = "all", fu = 1, probe_set = "full", ...) {
    calls <<- c(calls, paste("qq", stratum, fu, probe_set, sep = "/"))
    original_qq(results, stratum = stratum, fu = fu, probe_set = probe_set, ...)
  }, envir = .GlobalEnv)
  assign("plot_volcano", function(results, stratum = "all", fu = 1, probe_set = "full", ...) {
    calls <<- c(calls, paste("volcano", stratum, fu, probe_set, sep = "/"))
    original_volcano(results, stratum = stratum, fu = fu, probe_set = probe_set, ...)
  }, envir = .GlobalEnv)
  device_before <- dev.cur()
  generate_all_plots(list(analysis_change = results, metadata = list(shared_beta = TRUE)),
                     figures_dir = work_dir)
  expected_calls <- unlist(lapply(c("full", "filtered"), function(probe_set) {
    unlist(lapply(c("qq", "volcano"), function(kind) {
      c(paste(kind, "all", 1, probe_set, sep = "/"),
        paste(kind, "male", 2, probe_set, sep = "/"))
    }))
  }), use.names = FALSE)
  files <- list.files(work_dir, pattern = "\\.pdf$", full.names = TRUE)
  stopifnot(identical(calls, expected_calls), length(files) == 4L,
            all(file.info(files)$size > 0L), identical(dev.cur(), device_before))
  empty_results <- list(all = list(outcome_effects = effects[0L, ]),
                        male = NULL, female = NULL)
  empty_dir <- file.path(work_dir, "empty")
  quiet(generate_all_plots(empty_results, figures_dir = empty_dir))
  stopifnot(length(list.files(empty_dir)) == 0L, identical(dev.cur(), device_before))
  cat("PASS plot PDFs use actual stratum/FU cells and tolerate all-unavailable or empty results\n")

  unsupported <- results
  unsupported$all$outcome_effects$ARM <- "2"
  error <- tryCatch({ plot_volcano(unsupported); NULL }, error = conditionMessage)
  stopifnot(grepl("ARM tables are not supported", error),
            !"arm" %in% names(formals(plot_qq)),
            !"arm" %in% names(formals(plot_volcano)))
  cat("PASS plots use common analyte slopes and reject incompatible per-arm beta tables\n")

}

run_multiarm_report_tests()
