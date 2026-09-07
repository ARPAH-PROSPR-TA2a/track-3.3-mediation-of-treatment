.log_33 <- function(verbose, text) {
  if (isTRUE(verbose)) {
    message(paste0("[3.3] ", text))

    progress_log <- getOption("track33.progress_log")
    if (is.character(progress_log) && length(progress_log) == 1L &&
        !is.na(progress_log) && nzchar(progress_log)) {
      cat(
        "[3.3] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        " | ", text, "\n",
        file = progress_log, append = TRUE, sep = ""
      )
    }
  }

  invisible(NULL)
}


.format_duration_33 <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) {
    return("unknown")
  }
  if (seconds >= 3600) {
    return(sprintf("%.1f hr", seconds / 3600))
  }
  if (seconds >= 60) {
    return(sprintf("%.1f min", seconds / 60))
  }
  sprintf("%.1f sec", seconds)
}


.outcomewas_source_file <- sys.frame(1L)$ofile
if (!file.exists(.outcomewas_source_file) &&
    file.exists(basename(.outcomewas_source_file))) {
  .outcomewas_source_file <- basename(.outcomewas_source_file)
}
.outcomewas_dir <- dirname(normalizePath(.outcomewas_source_file, mustWork = TRUE))

source(file.path(.outcomewas_dir, "R", "validation_helpers.R"))
source(file.path(.outcomewas_dir, "R", "reporting_helpers.R"))
source(file.path(.outcomewas_dir, "R", "analysis_helpers.R"))
source(file.path(.outcomewas_dir, "R", "proteomics_translation_helpers.R"))

rm(.outcomewas_dir, .outcomewas_source_file)

FAST_outcome_WAS <- function(pheno,
                             omics,
                             omics_type = "Proteomics",
                             additional_covariates = NULL,
                             n_cores = NULL,
                             checkpoint_dir = NULL,
                             checkpoint_batch_size = 2000L,
                             verbose = FALSE) {

  analysis_started <- proc.time()[["elapsed"]]

  if (is.null(n_cores)) {
    detected_cores <- parallel::detectCores()
    if (length(detected_cores) != 1L || is.na(detected_cores) || detected_cores < 1L) {
      detected_cores <- future::availableCores()
    }
    n_cores <- max(1L, as.integer(detected_cores) - 1L)
  }

  # See TreatmentWAS note in the original scaffold. future snapshots the
  # current option state, and a NULL `exact` value can break restoration in
  # newer R versions.
  if ("exact" %in% names(.Options)) {
    .Options$exact <- NULL
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (n_cores > 1L) {
    future::plan(future::multicore(), workers = n_cores)
  } else {
    future::plan(future::sequential)
  }

  .log_33(verbose, "analysis validating inputs and preparing data (serial)")
  .validate_omics_type(omics_type)

  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(omics, pheno_list, verbose = verbose)

  filtered_probes <- NULL
  if (omics_type == "DNAm") {
    .log_33(verbose, "loading bundled DNAm probe lists (serial)")
    full_probes     <- readRDS("outcomewas/Data/FAST_epicv1_epicv2_probe_list.rds")
    filtered_probes <- readRDS("outcomewas/Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
    .validate_dnam_probe_coverage(full_probes, filtered_probes, omics_list$all$ANALYTE_NAME)
    omics_list <- .subset_omics_list(omics_list, full_probes)
  }

  backend <- if (n_cores > 1L) "multicore" else "sequential"
  checkpoint_status <- if (is.null(checkpoint_dir)) {
    "disabled"
  } else {
    paste0(checkpoint_dir, " (batch size=", checkpoint_batch_size, ")")
  }

  .log_33(
    verbose,
    paste0(
      "analysis starting: ", nrow(omics_list$all), " analytes, ",
      nrow(pheno_list$all), " phenotype rows; outcome=", pheno_list$outcome_type,
      "; backend=", backend, ", effective workers=", future::nbrOfWorkers(),
      "; checkpoints=", checkpoint_status
    )
  )

  analysis_change <- .run_stratified_analysis(
    pheno_list,
    omics_list,
    omics_type,
    pheno_list$outcome_type,
    additional_covariates,
    response_type = "change",
    filtered_probes = filtered_probes,
    checkpoint_dir = checkpoint_dir,
    checkpoint_batch_size = checkpoint_batch_size,
    verbose = verbose
  )

  analysis_level <- .run_stratified_analysis(
    pheno_list,
    omics_list,
    omics_type,
    pheno_list$outcome_type,
    additional_covariates,
    response_type = "level",
    filtered_probes = filtered_probes,
    checkpoint_dir = checkpoint_dir,
    checkpoint_batch_size = checkpoint_batch_size,
    verbose = verbose
  )

  .log_33(
    verbose,
    paste0(
      "analysis complete (",
      .format_duration_33(proc.time()[["elapsed"]] - analysis_started),
      ")"
    )
  )

  list(
    analysis_change = analysis_change,
    analysis_level = analysis_level
  )
}


FAST_outcome_WAS_reports <- function(pheno,
                                     omics,
                                     omics_type = "Proteomics",
                                     additional_covariates = NULL,
                                     verbose = FALSE) {

  reports_started <- proc.time()[["elapsed"]]
  .log_33(verbose, "reports validating inputs (serial)")

  .validate_omics_type(omics_type)

  pheno_list <- .validate_pheno(pheno, additional_covariates)
  omics_list <- .validate_omics(omics, pheno_list, verbose = verbose)

  if (omics_type == "DNAm") {
    .log_33(verbose, "loading bundled DNAm probe lists (serial)")
    full_probes     <- readRDS("outcomewas/Data/FAST_epicv1_epicv2_probe_list.rds")
    filtered_probes <- readRDS("outcomewas/Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds")
    .validate_dnam_probe_coverage(full_probes, filtered_probes, omics_list$all$ANALYTE_NAME)
    omics_list <- .subset_omics_list(omics_list, full_probes)
  }

  reports <- .generate_reports(
    pheno_list,
    omics_list,
    additional_covariates = additional_covariates,
    outcome_type = pheno_list$outcome_type,
    verbose = verbose
  )

  .log_33(
    verbose,
    paste0(
      "reports complete (",
      .format_duration_33(proc.time()[["elapsed"]] - reports_started),
      ")"
    )
  )

  reports
}
