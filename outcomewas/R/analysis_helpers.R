.build_group_key <- function(df, group_col) {
  if (length(group_col) == 1L) {
    return(as.character(df[[group_col]]))
  }

  do.call(
    paste,
    c(
      lapply(group_col, function(col) as.character(df[[col]])),
      sep = "\r"
    )
  )
}


.apply_multiple_testing_correction <- function(results_df, p_value_col = "P_VALUE",
                                               group_col = "COEFFICIENT") {

  if (is.null(results_df) || nrow(results_df) == 0) {
    return(results_df)
  }

  group_key <- .build_group_key(results_df, group_col)
  results_df$BH_P_VALUE <- NA_real_

  for (grp in unique(group_key)) {
    grp_idx <- which(group_key == grp)
    p_values_grp <- results_df[[p_value_col]][grp_idx]
    results_df$BH_P_VALUE[grp_idx] <- p.adjust(p_values_grp, method = "BH")
  }

  results_df
}


.prepare_analysis_data <- function(pheno_df, omics_df, pheno_baseline, omics_baseline) {

  omics_sample_ids <- colnames(omics_df)[-which(colnames(omics_df) == "ANALYTE_NAME")]
  shared_samples <- intersect(pheno_df$SAMPLE_ID, omics_sample_ids)
  pheno_merged <- pheno_df[pheno_df$SAMPLE_ID %in% shared_samples, ]

  baseline_subject_ids <- pheno_baseline$SUBJECT_ID
  omics_baseline_matrix <- as.matrix(omics_baseline)

  sample_subjects <- pheno_merged$SUBJECT_ID
  baseline_idx <- match(sample_subjects, baseline_subject_ids)
  baseline_col_idx <- match(pheno_baseline$SAMPLE_ID[baseline_idx], colnames(omics_baseline_matrix))

  list(
    pheno_merged = pheno_merged,
    shared_samples = shared_samples,
    baseline_col_idx = baseline_col_idx,
    omics_baseline_matrix = omics_baseline_matrix
  )
}


.drop_uninformative_covariates <- function(model_data, covariate_terms) {
  keep_terms <- character(0)

  for (term in covariate_terms) {
    term_values <- model_data[[term]]
    unique_values <- unique(term_values[!is.na(term_values)])

    if (length(unique_values) <= 1L) {
      if (term != "FEMALE") {
        warning("Dropping covariate '", term, "' because it has only one observed value in this stratum.")
      }
      next
    }

    keep_terms <- c(keep_terms, term)
  }

  keep_terms
}


.build_formula_string <- function(outcome_type, covariate_terms) {
  response <- if (outcome_type == "continuous") {
    "OUTCOME"
  } else {
    "survival::Surv(OUTCOME_TIME, OUTCOME_STATUS)"
  }

  formula_str <- paste(response, "~ analyte")
  if (length(covariate_terms) > 0) {
    formula_str <- paste(formula_str, paste(covariate_terms, collapse = " + "), sep = " + ")
  }

  formula_str
}


.perform_continuous_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline,
                                         fu_level,
                                         additional_covariates = NULL,
                                         response_type = c("change", "level"),
                                         checkpoint_dir = NULL,
                                         checkpoint_batch_size = 2000L,
                                         verbose = FALSE,
                                         progress_label = paste(
                                           response_type,
                                           paste0("FU", fu_level),
                                           sep = "/"
                                         )) {

  response_type <- match.arg(response_type)

  prep <- .prepare_analysis_data(pheno_df, omics_df, pheno_baseline, omics_baseline)
  pheno_merged <- prep$pheno_merged
  shared_samples <- prep$shared_samples
  baseline_col_idx <- prep$baseline_col_idx
  omics_baseline_matrix <- prep$omics_baseline_matrix

  model_data <- data.frame(pheno_merged)
  model_data$analyte <- NA_real_
  model_data$analyte_baseline <- NA_real_

  covariate_terms <- "analyte_baseline"
  adjustment_terms <- c("TREATMENT_GROUP", "FEMALE")
  if (!is.null(additional_covariates)) {
    adjustment_terms <- c(adjustment_terms, additional_covariates)
  }
  adjustment_terms <- .drop_uninformative_covariates(model_data, adjustment_terms)
  covariate_terms <- c(covariate_terms, adjustment_terms)
  formula_str <- .build_formula_string("continuous", covariate_terms)

  analyte_names <- omics_df$ANALYTE_NAME

  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }

  batches <- split(seq_along(analyte_names),
                   ceiling(seq_along(analyte_names) / checkpoint_batch_size))
  all_results <- vector("list", length(analyte_names))
  checkpoint_files <- if (!is.null(checkpoint_dir)) {
    file.path(checkpoint_dir, paste0("batch_", seq_along(batches), ".rds"))
  } else {
    character(0)
  }
  n_cached <- sum(file.exists(checkpoint_files))
  n_pending <- length(batches) - n_cached
  completed_pending <- 0L
  pending_elapsed <- 0
  new_failures <- 0L
  reported_worker_pids <- FALSE

  .log_33(
    verbose,
    paste0(
      progress_label, ": ", length(analyte_names), " analytes, ",
      length(batches), " batches (", n_cached, " cached, ",
      n_pending, " pending)"
    )
  )
  if (length(batches) > 0L && n_pending == 0L) {
    .log_33(verbose, paste0(progress_label, ": all batches cached; no workers launched"))
  }

  for (b in seq_along(batches)) {
    batch <- batches[[b]]
    batch_file <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, paste0("batch_", b, ".rds")) else NULL

    if (!is.null(batch_file) && file.exists(batch_file)) {
      all_results[batch] <- readRDS(batch_file)
      next
    }

    batch_started <- proc.time()[["elapsed"]]

    batch_items <- lapply(batch, function(i) list(
      analyte_name = analyte_names[i],
      fu_values = as.numeric(omics_df[i, shared_samples]),
      baseline_vals = omics_baseline_matrix[i, baseline_col_idx]
    ))

    n_workers <- min(future::nbrOfWorkers(), length(batch_items))
    .log_33(
      verbose,
      sprintf(
        "%s batch %d/%d: dispatching %d analytes across up to %d worker%s",
        progress_label, b, length(batches), length(batch), n_workers,
        if (n_workers == 1L) "" else "s"
      )
    )

    worker_results <- furrr::future_map(batch_items, function(item) {
      result <- tryCatch({
        md <- model_data
        if (response_type == "change") {
          md$analyte <- item$fu_values - item$baseline_vals
        } else {
          md$analyte <- item$fu_values
        }
        md$analyte_baseline <- item$baseline_vals

        fit <- lm(as.formula(formula_str), data = md)
        n_obs <- nrow(fit$model)
        coef_table <- summary(fit)$coefficients
        analyte_idx <- match("analyte", rownames(coef_table))
        if (is.na(analyte_idx)) {
          stop("analyte coefficient was not estimable")
        }

        coefs <- data.frame(
          ANALYTE_NAME = item$analyte_name,
          FU = fu_level,
          COEFFICIENT = rownames(coef_table),
          N_OBS = n_obs,
          EFFECT_SIZE = coef_table[, "Estimate"],
          SE = coef_table[, "Std. Error"],
          P_VALUE = coef_table[, "Pr(>|t|)"],
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        oe <- data.frame(
          ANALYTE_NAME = item$analyte_name,
          FU = fu_level,
          EFFECT_SIZE = coef_table[analyte_idx, "Estimate"],
          SE = coef_table[analyte_idx, "Std. Error"],
          P_VALUE = coef_table[analyte_idx, "Pr(>|t|)"],
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        list(coefficients = coefs, outcome_effects = oe)
      }, error = function(e) {
        warning("Error processing analyte '", item$analyte_name, "' at FU=", fu_level, ": ", e$message)
        NULL
      })

      list(
        result = result,
        worker_pid = Sys.getpid(),
        failed = is.null(result)
      )
    }, .options = furrr::furrr_options(seed = TRUE))

    worker_pids <- sort(unique(vapply(
      worker_results,
      function(result) result$worker_pid,
      integer(1)
    )))
    batch_failures <- sum(vapply(
      worker_results,
      function(result) result$failed,
      logical(1)
    ))
    batch_results <- lapply(worker_results, `[[`, "result")
    new_failures <- new_failures + batch_failures

    if (!reported_worker_pids) {
      .log_33(
        verbose,
        paste0(
          progress_label, " worker PIDs observed: ",
          paste(worker_pids, collapse = ", ")
        )
      )
      reported_worker_pids <- TRUE
    }

    if (!is.null(batch_file)) {
      saveRDS(batch_results, paste0(batch_file, ".tmp"))
      file.rename(paste0(batch_file, ".tmp"), batch_file)
    }

    all_results[batch] <- batch_results

    completed_pending <- completed_pending + 1L
    remaining_pending <- n_pending - completed_pending
    batch_elapsed <- proc.time()[["elapsed"]] - batch_started
    pending_elapsed <- pending_elapsed + batch_elapsed
    eta_text <- if (remaining_pending > 0L) {
      estimated_remaining <- pending_elapsed / completed_pending * remaining_pending
      paste0(
        "; estimated ", .format_duration_33(estimated_remaining),
        " remaining for ", remaining_pending, " pending batch",
        if (remaining_pending == 1L) "" else "es"
      )
    } else {
      ""
    }

    .log_33(
      verbose,
      paste0(
        progress_label, " batch ", b, "/", length(batches),
        " complete: ", length(worker_pids), " worker PID",
        if (length(worker_pids) == 1L) "" else "s",
        ", ", batch_failures, " failure",
        if (batch_failures == 1L) "" else "s",
        ", ", .format_duration_33(batch_elapsed),
        eta_text
      )
    )
  }

  if (new_failures > 0L) {
    .log_33(
      verbose,
      paste0(
        progress_label, ": ", new_failures, " newly attempted analyte",
        if (new_failures == 1L) "" else "s", " failed"
      )
    )
  }

  all_results <- Filter(Negate(is.null), all_results)
  if (length(all_results) == 0) {
    return(NULL)
  }

  coefficients <- do.call(rbind, lapply(all_results, `[[`, "coefficients"))
  outcome_effects <- do.call(rbind, lapply(all_results, `[[`, "outcome_effects"))
  row.names(coefficients) <- NULL
  row.names(outcome_effects) <- NULL

  list(
    coefficients = coefficients,
    outcome_effects = outcome_effects
  )
}


.perform_tte_analysis <- function(pheno_df, omics_df, pheno_baseline, omics_baseline,
                                  fu_level,
                                  additional_covariates = NULL,
                                  response_type = c("change", "level"),
                                  checkpoint_dir = NULL,
                                  checkpoint_batch_size = 2000L,
                                  verbose = FALSE,
                                  progress_label = paste(
                                    response_type,
                                    paste0("FU", fu_level),
                                    sep = "/"
                                  )) {

  response_type <- match.arg(response_type)

  prep <- .prepare_analysis_data(pheno_df, omics_df, pheno_baseline, omics_baseline)
  pheno_merged <- prep$pheno_merged
  shared_samples <- prep$shared_samples
  baseline_col_idx <- prep$baseline_col_idx
  omics_baseline_matrix <- prep$omics_baseline_matrix

  model_data <- data.frame(pheno_merged)
  model_data$analyte <- NA_real_
  model_data$analyte_baseline <- NA_real_

  if (sum(model_data$OUTCOME_STATUS == 1, na.rm = TRUE) == 0) {
    warning("No outcome events found in this stratum at FU=", fu_level, "; returning NULL for this follow-up.")
    return(NULL)
  }

  covariate_terms <- "analyte_baseline"
  adjustment_terms <- c("TREATMENT_GROUP", "FEMALE")
  if (!is.null(additional_covariates)) {
    adjustment_terms <- c(adjustment_terms, additional_covariates)
  }
  adjustment_terms <- .drop_uninformative_covariates(model_data, adjustment_terms)
  covariate_terms <- c(covariate_terms, adjustment_terms)
  formula_str <- .build_formula_string("tte", covariate_terms)

  analyte_names <- omics_df$ANALYTE_NAME

  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  }

  batches <- split(seq_along(analyte_names),
                   ceiling(seq_along(analyte_names) / checkpoint_batch_size))
  all_results <- vector("list", length(analyte_names))
  checkpoint_files <- if (!is.null(checkpoint_dir)) {
    file.path(checkpoint_dir, paste0("batch_", seq_along(batches), ".rds"))
  } else {
    character(0)
  }
  n_cached <- sum(file.exists(checkpoint_files))
  n_pending <- length(batches) - n_cached
  completed_pending <- 0L
  pending_elapsed <- 0
  new_failures <- 0L
  reported_worker_pids <- FALSE

  .log_33(
    verbose,
    paste0(
      progress_label, ": ", length(analyte_names), " analytes, ",
      length(batches), " batches (", n_cached, " cached, ",
      n_pending, " pending)"
    )
  )
  if (length(batches) > 0L && n_pending == 0L) {
    .log_33(verbose, paste0(progress_label, ": all batches cached; no workers launched"))
  }

  for (b in seq_along(batches)) {
    batch <- batches[[b]]
    batch_file <- if (!is.null(checkpoint_dir)) file.path(checkpoint_dir, paste0("batch_", b, ".rds")) else NULL

    if (!is.null(batch_file) && file.exists(batch_file)) {
      all_results[batch] <- readRDS(batch_file)
      next
    }

    batch_started <- proc.time()[["elapsed"]]

    batch_items <- lapply(batch, function(i) list(
      analyte_name = analyte_names[i],
      fu_values = as.numeric(omics_df[i, shared_samples]),
      baseline_vals = omics_baseline_matrix[i, baseline_col_idx]
    ))

    n_workers <- min(future::nbrOfWorkers(), length(batch_items))
    .log_33(
      verbose,
      sprintf(
        "%s batch %d/%d: dispatching %d analytes across up to %d worker%s",
        progress_label, b, length(batches), length(batch), n_workers,
        if (n_workers == 1L) "" else "s"
      )
    )

    worker_results <- furrr::future_map(batch_items, function(item) {
      result <- tryCatch({
        md <- model_data
        if (response_type == "change") {
          md$analyte <- item$fu_values - item$baseline_vals
        } else {
          md$analyte <- item$fu_values
        }
        md$analyte_baseline <- item$baseline_vals

        fit <- survival::coxph(as.formula(formula_str), data = md, ties = "efron")
        n_obs <- stats::nobs(fit)
        n_events <- fit$nevent
        coef_table <- summary(fit)$coefficients
        analyte_idx <- match("analyte", rownames(coef_table))
        if (is.na(analyte_idx)) {
          stop("analyte coefficient was not estimable")
        }

        coefs <- data.frame(
          ANALYTE_NAME = item$analyte_name,
          FU = fu_level,
          COEFFICIENT = rownames(coef_table),
          N_OBS = n_obs,
          N_EVENTS = n_events,
          EFFECT_SIZE = coef_table[, "coef"],
          HAZARD_RATIO = coef_table[, "exp(coef)"],
          SE = coef_table[, "se(coef)"],
          P_VALUE = coef_table[, "Pr(>|z|)"],
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        oe <- data.frame(
          ANALYTE_NAME = item$analyte_name,
          FU = fu_level,
          EFFECT_SIZE = coef_table[analyte_idx, "coef"],
          HAZARD_RATIO = coef_table[analyte_idx, "exp(coef)"],
          SE = coef_table[analyte_idx, "se(coef)"],
          P_VALUE = coef_table[analyte_idx, "Pr(>|z|)"],
          stringsAsFactors = FALSE,
          row.names = NULL
        )

        list(coefficients = coefs, outcome_effects = oe)
      }, error = function(e) {
        warning("Error processing analyte '", item$analyte_name, "' at FU=", fu_level, ": ", e$message)
        NULL
      })

      list(
        result = result,
        worker_pid = Sys.getpid(),
        failed = is.null(result)
      )
    }, .options = furrr::furrr_options(seed = TRUE, packages = "survival"))

    worker_pids <- sort(unique(vapply(
      worker_results,
      function(result) result$worker_pid,
      integer(1)
    )))
    batch_failures <- sum(vapply(
      worker_results,
      function(result) result$failed,
      logical(1)
    ))
    batch_results <- lapply(worker_results, `[[`, "result")
    new_failures <- new_failures + batch_failures

    if (!reported_worker_pids) {
      .log_33(
        verbose,
        paste0(
          progress_label, " worker PIDs observed: ",
          paste(worker_pids, collapse = ", ")
        )
      )
      reported_worker_pids <- TRUE
    }

    if (!is.null(batch_file)) {
      saveRDS(batch_results, paste0(batch_file, ".tmp"))
      file.rename(paste0(batch_file, ".tmp"), batch_file)
    }

    all_results[batch] <- batch_results

    completed_pending <- completed_pending + 1L
    remaining_pending <- n_pending - completed_pending
    batch_elapsed <- proc.time()[["elapsed"]] - batch_started
    pending_elapsed <- pending_elapsed + batch_elapsed
    eta_text <- if (remaining_pending > 0L) {
      estimated_remaining <- pending_elapsed / completed_pending * remaining_pending
      paste0(
        "; estimated ", .format_duration_33(estimated_remaining),
        " remaining for ", remaining_pending, " pending batch",
        if (remaining_pending == 1L) "" else "es"
      )
    } else {
      ""
    }

    .log_33(
      verbose,
      paste0(
        progress_label, " batch ", b, "/", length(batches),
        " complete: ", length(worker_pids), " worker PID",
        if (length(worker_pids) == 1L) "" else "s",
        ", ", batch_failures, " failure",
        if (batch_failures == 1L) "" else "s",
        ", ", .format_duration_33(batch_elapsed),
        eta_text
      )
    )
  }

  if (new_failures > 0L) {
    .log_33(
      verbose,
      paste0(
        progress_label, ": ", new_failures, " newly attempted analyte",
        if (new_failures == 1L) "" else "s", " failed"
      )
    )
  }

  all_results <- Filter(Negate(is.null), all_results)
  if (length(all_results) == 0) {
    return(NULL)
  }

  coefficients <- do.call(rbind, lapply(all_results, `[[`, "coefficients"))
  outcome_effects <- do.call(rbind, lapply(all_results, `[[`, "outcome_effects"))
  row.names(coefficients) <- NULL
  row.names(outcome_effects) <- NULL

  list(
    coefficients = coefficients,
    outcome_effects = outcome_effects
  )
}


.perform_analysis <- function(pheno_df, omics_df, omics_type, outcome_type,
                              additional_covariates = NULL,
                              response_type = c("change", "level"),
                              checkpoint_dir = NULL,
                              checkpoint_batch_size = 2000L,
                              verbose = FALSE,
                              progress_label = response_type) {

  response_type <- match.arg(response_type)

  pheno_baseline_all <- pheno_df[pheno_df$FU == 0, ]
  fu_levels <- sort(unique(as.integer(as.character(pheno_df$FU))))
  fu_levels <- fu_levels[fu_levels != 0]

  all_coefficients <- list()
  all_outcome_effects <- list()

  .log_33(
    verbose,
    paste0(
      progress_label, ": ", length(fu_levels), " follow-up level",
      if (length(fu_levels) == 1L) "" else "s", " to process"
    )
  )

  for (fu_level in fu_levels) {
    fu_started <- proc.time()[["elapsed"]]
    fu_progress_label <- paste0(progress_label, "/FU", fu_level)
    pheno_analysis <- pheno_df[as.integer(as.character(pheno_df$FU)) == fu_level, ]
    complete_subjects <- intersect(pheno_baseline_all$SUBJECT_ID, pheno_analysis$SUBJECT_ID)

    if (length(complete_subjects) == 0) {
      .log_33(verbose, paste0(fu_progress_label, ": skipped; no complete subjects"))
      next
    }

    .log_33(
      verbose,
      paste0(
        fu_progress_label, " starting: ", length(complete_subjects),
        " complete subjects; outcome=", outcome_type
      )
    )

    pheno_baseline <- pheno_baseline_all[pheno_baseline_all$SUBJECT_ID %in% complete_subjects, ]
    pheno_analysis <- pheno_analysis[pheno_analysis$SUBJECT_ID %in% complete_subjects, ]

    baseline_sample_ids <- pheno_baseline$SAMPLE_ID
    omics_baseline <- omics_df[, colnames(omics_df) %in% baseline_sample_ids, drop = FALSE]

    fu_checkpoint_dir <- if (!is.null(checkpoint_dir)) {
      file.path(checkpoint_dir, paste0("FU", fu_level))
    } else {
      NULL
    }

    fu_results <- if (outcome_type == "continuous") {
      .perform_continuous_analysis(
        pheno_analysis,
        omics_df,
        pheno_baseline,
        omics_baseline,
        fu_level = fu_level,
        additional_covariates = additional_covariates,
        response_type = response_type,
        checkpoint_dir = fu_checkpoint_dir,
        checkpoint_batch_size = checkpoint_batch_size,
        verbose = verbose,
        progress_label = fu_progress_label
      )
    } else {
      .perform_tte_analysis(
        pheno_analysis,
        omics_df,
        pheno_baseline,
        omics_baseline,
        fu_level = fu_level,
        additional_covariates = additional_covariates,
        response_type = response_type,
        checkpoint_dir = fu_checkpoint_dir,
        checkpoint_batch_size = checkpoint_batch_size,
        verbose = verbose,
        progress_label = fu_progress_label
      )
    }

    if (is.null(fu_results)) {
      .log_33(
        verbose,
        paste0(
          fu_progress_label, " complete: no estimable results (",
          .format_duration_33(proc.time()[["elapsed"]] - fu_started), ")"
        )
      )
      next
    }

    all_coefficients[[paste0("FU", fu_level)]] <- fu_results$coefficients
    all_outcome_effects[[paste0("FU", fu_level)]] <- fu_results$outcome_effects

    .log_33(
      verbose,
      paste0(
        fu_progress_label, " complete: ", nrow(fu_results$outcome_effects),
        " outcome-effect rows (",
        .format_duration_33(proc.time()[["elapsed"]] - fu_started), ")"
      )
    )
  }

  if (length(all_coefficients) == 0 || length(all_outcome_effects) == 0) {
    return(NULL)
  }

  coefficients <- do.call(rbind, all_coefficients)
  outcome_effects <- do.call(rbind, all_outcome_effects)
  row.names(coefficients) <- NULL
  row.names(outcome_effects) <- NULL

  .log_33(
    verbose,
    paste0(progress_label, ": model batches complete; applying BH correction (serial)")
  )

  coefficients <- .apply_multiple_testing_correction(
    coefficients,
    group_col = c("FU", "COEFFICIENT")
  )
  coefficients <- coefficients[
    order(coefficients$ANALYTE_NAME, coefficients$FU, coefficients$COEFFICIENT),
  ]

  outcome_effects <- .apply_multiple_testing_correction(
    outcome_effects,
    group_col = "FU"
  )
  outcome_effects <- outcome_effects[
    order(outcome_effects$ANALYTE_NAME, outcome_effects$FU),
  ]

  list(
    coefficients = coefficients,
    outcome_effects = outcome_effects
  )
}


.add_filtered_bh_column <- function(df, filtered_probes, group_col) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }

  group_key <- .build_group_key(df, group_col)
  df$BH_P_VALUE_FILTERED <- NA_real_

  for (grp in unique(group_key)) {
    idx <- which(group_key == grp & df$ANALYTE_NAME %in% filtered_probes)
    if (length(idx) > 0) {
      df$BH_P_VALUE_FILTERED[idx] <- p.adjust(df$P_VALUE[idx], method = "BH")
    }
  }

  df
}


.add_filtered_bh_correction <- function(outputs, filtered_probes) {
  for (stratum in c("all", "male", "female")) {
    if (is.null(outputs[[stratum]])) next

    outputs[[stratum]]$coefficients <- .add_filtered_bh_column(
      outputs[[stratum]]$coefficients,
      filtered_probes,
      group_col = c("FU", "COEFFICIENT")
    )

    outputs[[stratum]]$outcome_effects <- .add_filtered_bh_column(
      outputs[[stratum]]$outcome_effects,
      filtered_probes,
      group_col = "FU"
    )
  }

  outputs
}


.run_stratified_analysis <- function(pheno_list, omics_list, omics_type, outcome_type,
                                     additional_covariates,
                                     response_type = c("change", "level"),
                                     filtered_probes = NULL,
                                     checkpoint_dir = NULL,
                                     checkpoint_batch_size = 2000L,
                                     verbose = FALSE) {

  response_type <- match.arg(response_type)
  outputs <- list(all = NULL, male = NULL, female = NULL)

  for (dataset in c("all", "male", "female")) {
    if (is.null(pheno_list[[dataset]])) next
    stratum_started <- proc.time()[["elapsed"]]
    progress_label <- paste(response_type, dataset, sep = "/")

    .log_33(verbose, paste0(progress_label, " starting"))

    stratum_checkpoint_dir <- if (!is.null(checkpoint_dir)) {
      file.path(checkpoint_dir, response_type, dataset)
    } else {
      NULL
    }

    analysis_results <- .perform_analysis(
      pheno_list[[dataset]],
      omics_list[[dataset]],
      omics_type,
      outcome_type,
      additional_covariates,
      response_type,
      stratum_checkpoint_dir,
      checkpoint_batch_size,
      verbose = verbose,
      progress_label = progress_label
    )

    if (is.null(analysis_results)) {
      outputs[[dataset]] <- NULL
      .log_33(
        verbose,
        paste0(
          progress_label, " complete: no estimable results (",
          .format_duration_33(proc.time()[["elapsed"]] - stratum_started), ")"
        )
      )
      next
    }

    outputs[[dataset]] <- list(
      coefficients = analysis_results$coefficients,
      outcome_effects = analysis_results$outcome_effects
    )

    .log_33(
      verbose,
      paste0(
        progress_label, " complete (",
        .format_duration_33(proc.time()[["elapsed"]] - stratum_started), ")"
      )
    )
  }

  if (!is.null(filtered_probes)) {
    .log_33(
      verbose,
      paste0(response_type, ": applying filtered-probe BH correction (serial)")
    )
    outputs <- .add_filtered_bh_correction(outputs, filtered_probes)
  }

  outputs
}
