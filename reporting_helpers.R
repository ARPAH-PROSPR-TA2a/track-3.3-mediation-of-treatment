.create_pheno_data_report <- function(pheno_df) {

  groups <- unique(pheno_df[, c("FU", "FEMALE")])
  groups <- groups[order(groups$FU, groups$FEMALE), , drop = FALSE]

  report <- data.frame(
    FU = groups$FU,
    FEMALE = groups$FEMALE,
    N_SUBJECTS = NA_integer_,
    N_CONTROL = NA_integer_,
    N_TREATMENT = NA_integer_,
    N_SAMPLES = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(report))) {
    cell <- pheno_df[pheno_df$FU == report$FU[i] &
                       pheno_df$FEMALE == report$FEMALE[i], ]
    report$N_SUBJECTS[i] <- length(unique(cell$SUBJECT_ID))
    report$N_CONTROL[i] <- length(unique(cell$SUBJECT_ID[cell$TREATMENT_GROUP == 0]))
    report$N_TREATMENT[i] <- length(unique(cell$SUBJECT_ID[cell$TREATMENT_GROUP == 1]))
    report$N_SAMPLES[i] <- nrow(cell)
  }

  rownames(report) <- NULL
  report
}


.create_omics_data_report <- function(sample_ids, omics_df) {

  analyte_names <- omics_df$ANALYTE_NAME
  omics_numeric <- omics_df[, setdiff(names(omics_df), "ANALYTE_NAME"), drop = FALSE]
  omics_numeric <- omics_numeric[, colnames(omics_numeric) %in% sample_ids, drop = FALSE]

  report <- data.frame(
    ANALYTE_NAME = analyte_names,
    N_NONMISSING = NA_integer_,
    MEAN = NA_real_,
    MEDIAN = NA_real_,
    SD = NA_real_,
    MIN = NA_real_,
    MAX = NA_real_,
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


.create_addx_covariate_report <- function(pheno_df, covariate_names) {

  if (is.null(covariate_names) || length(covariate_names) == 0) {
    return(NULL)
  }

  results_list <- list(
    COVARIATE_NAME = character(),
    TYPE = character(),
    N_NA = integer(),
    SUMMARY = list()
  )

  for (i in seq_along(covariate_names)) {
    covar_name <- covariate_names[i]
    covar_data <- pheno_df[[covar_name]]
    n_na <- sum(is.na(covar_data))

    if (is.numeric(covar_data)) {
      covar_type <- "numeric"
      summary_stats <- list(
        mean = mean(covar_data, na.rm = TRUE),
        median = median(covar_data, na.rm = TRUE),
        sd = sd(covar_data, na.rm = TRUE),
        min = min(covar_data, na.rm = TRUE),
        max = max(covar_data, na.rm = TRUE)
      )
    } else if (is.factor(covar_data)) {
      covar_type <- "factor"
      level_counts <- table(covar_data, useNA = "no")
      summary_stats <- list(
        n_levels = nlevels(covar_data),
        level_names = levels(covar_data),
        counts = as.numeric(level_counts)
      )
    } else if (is.logical(covar_data)) {
      covar_type <- "logical"
      summary_stats <- list(
        n_true = sum(covar_data == TRUE, na.rm = TRUE),
        n_false = sum(covar_data == FALSE, na.rm = TRUE)
      )
    } else {
      covar_type <- "unknown"
      summary_stats <- list()
    }

    results_list$COVARIATE_NAME <- c(results_list$COVARIATE_NAME, covar_name)
    results_list$TYPE <- c(results_list$TYPE, covar_type)
    results_list$N_NA <- c(results_list$N_NA, n_na)
    results_list$SUMMARY[[i]] <- summary_stats
  }

  report <- data.frame(
    COVARIATE_NAME = results_list$COVARIATE_NAME,
    TYPE = results_list$TYPE,
    N_NA = results_list$N_NA,
    stringsAsFactors = FALSE
  )
  report$SUMMARY <- results_list$SUMMARY

  report
}


.subject_level_outcome_frame <- function(pheno_df) {
  fu_num <- as.integer(as.character(pheno_df$FU))
  ordered <- pheno_df[order(pheno_df$SUBJECT_ID, fu_num), ]
  ordered[!duplicated(ordered$SUBJECT_ID), ]
}


.create_analysis_sample_summary <- function(pheno_df) {
  fu_num <- as.integer(as.character(pheno_df$FU))
  baseline_subjects <- unique(pheno_df$SUBJECT_ID[fu_num == 0])
  fu_levels <- sort(unique(fu_num[fu_num != 0]))

  report <- data.frame(
    FU = fu_levels,
    N_SUBJECTS = NA_integer_,
    N_BASELINE_SAMPLES = NA_integer_,
    N_FOLLOWUP_SAMPLES = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(fu_levels)) {
    fu_level <- fu_levels[i]
    fu_subjects <- unique(pheno_df$SUBJECT_ID[fu_num == fu_level])
    complete_subjects <- intersect(baseline_subjects, fu_subjects)

    report$N_SUBJECTS[i] <- length(complete_subjects)
    report$N_BASELINE_SAMPLES[i] <- sum(fu_num == 0 & pheno_df$SUBJECT_ID %in% complete_subjects)
    report$N_FOLLOWUP_SAMPLES[i] <- sum(fu_num == fu_level & pheno_df$SUBJECT_ID %in% complete_subjects)
  }

  report
}


.create_continuous_outcome_report <- function(pheno_df) {

  pheno_subject <- .subject_level_outcome_frame(pheno_df)
  groups <- list(
    all = pheno_subject,
    control = pheno_subject[pheno_subject$TREATMENT_GROUP == 0, ],
    treatment = pheno_subject[pheno_subject$TREATMENT_GROUP == 1, ]
  )

  report <- data.frame(
    GROUP = names(groups),
    N_SUBJECTS = NA_integer_,
    MEAN = NA_real_,
    MEDIAN = NA_real_,
    SD = NA_real_,
    MIN = NA_real_,
    MAX = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(groups)) {
    values <- groups[[i]]$OUTCOME
    report$N_SUBJECTS[i] <- length(values)
    report$MEAN[i] <- mean(values, na.rm = TRUE)
    report$MEDIAN[i] <- median(values, na.rm = TRUE)
    report$SD[i] <- sd(values, na.rm = TRUE)
    report$MIN[i] <- min(values, na.rm = TRUE)
    report$MAX[i] <- max(values, na.rm = TRUE)
  }

  report
}


.create_tte_outcome_report <- function(pheno_df) {

  pheno_subject <- .subject_level_outcome_frame(pheno_df)
  groups <- list(
    all = pheno_subject,
    control = pheno_subject[pheno_subject$TREATMENT_GROUP == 0, ],
    treatment = pheno_subject[pheno_subject$TREATMENT_GROUP == 1, ]
  )

  report <- data.frame(
    GROUP = names(groups),
    N_SUBJECTS = NA_integer_,
    N_EVENTS = NA_integer_,
    N_CENSORED = NA_integer_,
    EVENT_RATE = NA_real_,
    TIME_MEAN = NA_real_,
    TIME_MEDIAN = NA_real_,
    TIME_MIN = NA_real_,
    TIME_MAX = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(groups)) {
    df <- groups[[i]]
    report$N_SUBJECTS[i] <- nrow(df)
    report$N_EVENTS[i] <- sum(df$OUTCOME_STATUS == 1, na.rm = TRUE)
    report$N_CENSORED[i] <- sum(df$OUTCOME_STATUS == 0, na.rm = TRUE)
    report$EVENT_RATE[i] <- mean(df$OUTCOME_STATUS == 1, na.rm = TRUE)
    report$TIME_MEAN[i] <- mean(df$OUTCOME_TIME, na.rm = TRUE)
    report$TIME_MEDIAN[i] <- median(df$OUTCOME_TIME, na.rm = TRUE)
    report$TIME_MIN[i] <- min(df$OUTCOME_TIME, na.rm = TRUE)
    report$TIME_MAX[i] <- max(df$OUTCOME_TIME, na.rm = TRUE)
  }

  report
}


.generate_reports <- function(pheno_list, omics_list, additional_covariates = NULL,
                              outcome_type = c("continuous", "tte")) {

  outcome_type <- match.arg(outcome_type)

  variable_summaries <- list(all = NULL, male = NULL, female = NULL)

  for (dataset in c("all", "male", "female")) {
    if (is.null(pheno_list[[dataset]])) next

    pheno_df <- pheno_list[[dataset]]
    omics_df <- omics_list[[dataset]]
    fu_levels <- as.character(sort(unique(as.integer(as.character(pheno_df$FU)))))
    tx_levels <- c("0", "1")
    cell_reports <- list()

    for (fu in fu_levels) {
      for (tx in tx_levels) {
        cell_key <- paste0("_FU", fu, "_Tx", tx)
        cell_mask <- as.character(pheno_df$FU) == fu &
          as.character(pheno_df$TREATMENT_GROUP) == tx
        cell_pheno <- pheno_df[cell_mask, ]

        if (nrow(cell_pheno) == 0) next

        cell_reports[[paste0("omics", cell_key)]] <-
          .create_omics_data_report(cell_pheno$SAMPLE_ID, omics_df)

        if (!is.null(additional_covariates)) {
          cell_reports[[paste0("covariates", cell_key)]] <-
            .create_addx_covariate_report(cell_pheno, additional_covariates)
        }
      }
    }

    variable_summaries[[dataset]] <- cell_reports
  }

  outcome_summary <- if (outcome_type == "continuous") {
    .create_continuous_outcome_report(pheno_list$all)
  } else {
    .create_tte_outcome_report(pheno_list$all)
  }

  outcome_reports <- list(
    outcome_type = outcome_type,
    analysis_sample_summary = .create_analysis_sample_summary(pheno_list$all),
    outcome_summary = outcome_summary
  )

  list(
    pheno_summary = .create_pheno_data_report(pheno_list$all),
    variable_summaries = variable_summaries,
    outcome_reports = outcome_reports
  )
}
