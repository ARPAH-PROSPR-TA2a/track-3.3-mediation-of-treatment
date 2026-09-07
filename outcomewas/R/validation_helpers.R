# Helper function for NZV detection
.is_near_zero_variance <- function(x) {
  if (!is.numeric(x)) return(FALSE)
  if (all(is.na(x))) return(FALSE)
  var(x, na.rm = TRUE) < 1e-8
}


.detect_outcome_type <- function(pheno) {
  has_continuous <- "OUTCOME" %in% names(pheno)
  has_tte_time   <- "OUTCOME_TIME" %in% names(pheno)
  has_tte_status <- "OUTCOME_STATUS" %in% names(pheno)
  has_tte        <- has_tte_time || has_tte_status

  if (has_continuous && has_tte) {
    stop("pheno must contain either OUTCOME or OUTCOME_TIME/OUTCOME_STATUS, not both.")
  }

  if (has_continuous) {
    return("continuous")
  }

  if (has_tte_time && has_tte_status) {
    return("tte")
  }

  if (has_tte_time || has_tte_status) {
    stop("Time-to-event outcomes require both OUTCOME_TIME and OUTCOME_STATUS.")
  }

  stop("pheno must contain either OUTCOME or OUTCOME_TIME and OUTCOME_STATUS.")
}


.subject_constant <- function(x, subject_ids) {
  split_values <- split(as.character(x), subject_ids)
  all(vapply(
    split_values,
    function(values) length(unique(values[!is.na(values)])) <= 1L,
    logical(1)
  ))
}


.validate_omics_type <- function(omics_type) {

  acceptable_types <- c("DNAm", "Proteomics", "Metabolomics", "other")

  if (!omics_type %in% acceptable_types) {
    stop(
      "Invalid omics_type '", omics_type, "'. ",
      "Must be one of: ", paste(acceptable_types, collapse = ", ")
    )
  }

  if (omics_type == "DNAm") {
    message("DNAm: beta values or M-values are used on the supplied scale; no transformation is applied.")
  } else if (omics_type == "Metabolomics") {
    message("Metabolomics: inputs should be log2-transformed prior to analysis.")
  } else if (omics_type == "Proteomics") {
    message("Proteomics: inputs should be log2-transformed prior to analysis.")
  }
}


.validate_pheno <- function(pheno, additional_covariates = NULL) {

  if (is.matrix(pheno)) {
    pheno <- as.data.frame(pheno)
  } else if (!is.data.frame(pheno)) {
    stop("pheno must be a data.frame or matrix")
  }

  if (!is.null(additional_covariates) && !is.character(additional_covariates)) {
    stop("additional_covariates must be NULL or a character vector")
  }

  outcome_type <- .detect_outcome_type(pheno)
  outcome_cols <- if (outcome_type == "continuous") {
    "OUTCOME"
  } else {
    c("OUTCOME_TIME", "OUTCOME_STATUS")
  }

  required_cols <- c("SAMPLE_ID", "FU", "SUBJECT_ID", "FEMALE", "TREATMENT_GROUP", outcome_cols)
  missing_cols <- setdiff(required_cols, names(pheno))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  fu_num <- suppressWarnings(as.integer(as.character(pheno$FU)))
  if (any(is.na(fu_num))) {
    stop("FU must be integer-valued (e.g., 0, 1, 2, ...).")
  }
  if (any(fu_num < 0)) {
    stop("FU must be non-negative (baseline is FU == 0).")
  }
  if (!any(fu_num == 0)) {
    stop("pheno must contain at least one baseline sample (FU == 0)")
  }
  if (!any(fu_num == 1)) {
    stop("pheno must contain at least one follow-up sample (FU == 1)")
  }

  unique_fu <- sort(unique(fu_num))
  expected_fu <- seq.int(0, max(unique_fu))
  if (!identical(unique_fu, expected_fu)) {
    stop(
      "FU must be encoded as consecutive integers starting at 0: ",
      paste(expected_fu, collapse = ", "),
      ". Found: ", paste(unique_fu, collapse = ", "),
      ". If your raw follow-up is in months (e.g., 0,3,6,12), recode to 0,1,2,3 before running."
    )
  }
  if (!is.factor(pheno$FU)) {
    warning("FU column is not a factor. Converting to factor.")
  }
  pheno$FU <- factor(fu_num)

  if (!all(pheno$FEMALE %in% 0:1)) {
    stop("FEMALE must contain only values 0 or 1")
  }
  if (!is.factor(pheno$FEMALE)) {
    warning("FEMALE column is not a factor. Converting to factor.")
  }
  pheno$FEMALE <- factor(pheno$FEMALE)

  if (!all(pheno$TREATMENT_GROUP %in% 0:1)) {
    stop("TREATMENT_GROUP must contain only values 0 or 1")
  }
  if (!any(pheno$TREATMENT_GROUP == 0) || !any(pheno$TREATMENT_GROUP == 1)) {
    stop("pheno must contain both control (TREATMENT_GROUP == 0) and treatment (TREATMENT_GROUP == 1)")
  }
  if (!is.factor(pheno$TREATMENT_GROUP)) {
    warning("TREATMENT_GROUP column is not a factor. Converting to factor.")
  }
  pheno$TREATMENT_GROUP <- factor(pheno$TREATMENT_GROUP)

  if (any(duplicated(pheno$SAMPLE_ID))) {
    stop("SAMPLE_ID contains duplicate values")
  }

  subject_fu_pairs <- paste(pheno$SUBJECT_ID, pheno$FU, sep = "_")
  if (any(duplicated(subject_fu_pairs))) {
    warning("Found duplicate SUBJECT_ID/FU pairs. Keeping first occurrence, discarding replicates.")
    pheno <- pheno[!duplicated(subject_fu_pairs), ]
  }

  if (!is.null(additional_covariates)) {
    missing_addl <- setdiff(additional_covariates, names(pheno))
    if (length(missing_addl) > 0) {
      stop("Additional covariates not found in pheno: ", paste(missing_addl, collapse = ", "))
    }
  }

  if (outcome_type == "continuous") {
    if (!is.numeric(pheno$OUTCOME)) {
      stop("OUTCOME must be numeric")
    }
  } else {
    if (!is.numeric(pheno$OUTCOME_TIME)) {
      stop("OUTCOME_TIME must be numeric")
    }
    if (any(pheno$OUTCOME_TIME < 0, na.rm = TRUE)) {
      stop("OUTCOME_TIME must be non-negative")
    }
    warning("OUTCOME_TIME is assumed to be measured in years. Convert source data before running if time is recorded in days, months, or another unit.")

    status_num <- suppressWarnings(as.integer(as.character(pheno$OUTCOME_STATUS)))
    non_missing_status <- !is.na(pheno$OUTCOME_STATUS)
    if (any(is.na(status_num[non_missing_status])) || !all(status_num[non_missing_status] %in% 0:1)) {
      stop("OUTCOME_STATUS must contain only values 0 or 1")
    }
    pheno$OUTCOME_STATUS <- status_num
  }

  for (covar in if (is.null(additional_covariates)) character(0) else additional_covariates) {
    col_data <- pheno[[covar]]

    if (!is.numeric(col_data) && !is.factor(col_data) && !is.logical(col_data)) {
      stop("Additional covariate '", covar, "' must be numeric, factor, or logical")
    }

    na_rows <- is.na(col_data)
    if (any(na_rows)) {
      message("Dropping ", sum(na_rows), " sample(s) with NA values in covariate '", covar, "'.")
      pheno <- pheno[!na_rows, ]
    }
  }

  missing_outcome_mask <- Reduce(`|`, lapply(outcome_cols, function(col) is.na(pheno[[col]])))
  if (any(missing_outcome_mask)) {
    message(
      "Dropping ", sum(missing_outcome_mask),
      " sample(s) with missing or incomplete outcome data; please check your data."
    )
    pheno <- pheno[!missing_outcome_mask, ]
  }

  constant_cols <- c("TREATMENT_GROUP", "FEMALE", outcome_cols)
  non_constant <- constant_cols[!vapply(
    constant_cols,
    function(col) .subject_constant(pheno[[col]], pheno$SUBJECT_ID),
    logical(1)
  )]
  if (length(non_constant) > 0) {
    stop(
      "The following columns must be constant within SUBJECT_ID: ",
      paste(non_constant, collapse = ", ")
    )
  }

  fu_num_current <- as.integer(as.character(pheno$FU))
  subjects_with_baseline <- unique(pheno$SUBJECT_ID[fu_num_current == 0])
  subjects_with_followup <- unique(pheno$SUBJECT_ID[fu_num_current > 0])
  complete_subjects <- intersect(subjects_with_baseline, subjects_with_followup)

  n_incomplete <- length(unique(pheno$SUBJECT_ID)) - length(complete_subjects)
  if (n_incomplete > 0) {
    message("Dropping ", n_incomplete, " subject(s) missing either a baseline or a follow-up sample.")
  }
  pheno <- pheno[pheno$SUBJECT_ID %in% complete_subjects, ]

  if (nrow(pheno) == 0) {
    stop("No subjects remain after requiring complete baseline/follow-up data, complete outcomes, and complete covariates.")
  }

  n_male <- sum(pheno$FEMALE == 0, na.rm = TRUE)
  n_female <- sum(pheno$FEMALE == 1, na.rm = TRUE)

  if (n_male == 0 || n_female == 0) {
    warning("Dataset contains only one gender. Male and female subsets will be NULL.")
    male_pheno <- NULL
    female_pheno <- NULL
  } else {
    male_pheno <- pheno[pheno$FEMALE == 0, ]
    female_pheno <- pheno[pheno$FEMALE == 1, ]
  }

  cols_to_keep <- c(required_cols, additional_covariates)
  cols_to_keep <- intersect(cols_to_keep, names(pheno))

  pheno <- pheno[, cols_to_keep, drop = FALSE]
  if (!is.null(male_pheno)) {
    male_pheno <- male_pheno[, cols_to_keep, drop = FALSE]
    female_pheno <- female_pheno[, cols_to_keep, drop = FALSE]
  }

  list(
    all = pheno,
    male = male_pheno,
    female = female_pheno,
    outcome_type = outcome_type,
    requires_mixed_effects = FALSE
  )
}


.validate_omics <- function(omics, pheno_list, verbose = FALSE,
                            progress_every = 100000L) {

  if (is.matrix(omics)) {
    omics <- as.data.frame(omics)
  } else if (!is.data.frame(omics)) {
    stop("omics must be a data.frame or matrix")
  }

  if (!"ANALYTE_NAME" %in% names(omics)) {
    stop("omics must contain ANALYTE_NAME column")
  }

  analyte_names <- omics$ANALYTE_NAME
  omics_numeric <- omics[, setdiff(names(omics), "ANALYTE_NAME"), drop = FALSE]

  if (!all(sapply(omics_numeric, is.numeric))) {
    stop("All columns in omics (except ANALYTE_NAME) must be numeric")
  }

  pheno_sample_ids <- pheno_list$all$SAMPLE_ID
  omics_sample_ids <- names(omics_numeric)
  shared_samples <- intersect(omics_sample_ids, pheno_sample_ids)

  if (length(shared_samples) == 0) {
    stop("No overlap between omics column names and pheno SAMPLE_IDs")
  }

  message("Found ", length(shared_samples), " samples shared between omics and pheno")
  message("  Omics only: ", length(setdiff(omics_sample_ids, pheno_sample_ids)))
  message("  Pheno only: ", length(setdiff(pheno_sample_ids, omics_sample_ids)))

  omics_numeric <- omics_numeric[, shared_samples, drop = FALSE]

  n_with_na <- 0
  n_with_nzv <- 0
  n_analytes <- length(analyte_names)
  quality_check_started <- proc.time()[["elapsed"]]

  .log_33(
    verbose,
    paste0("omics quality checks: 0/", n_analytes, " analytes (serial)")
  )

  for (i in seq_along(analyte_names)) {
    analyte_data <- as.numeric(omics_numeric[i, ])

    if (any(is.na(analyte_data))) {
      n_with_na <- n_with_na + 1
    }

    if (.is_near_zero_variance(analyte_data)) {
      n_with_nzv <- n_with_nzv + 1
    }

    if (i < n_analytes && i %% progress_every == 0L) {
      .log_33(
        verbose,
        sprintf(
          "omics quality checks: %d/%d analytes (%.1f%%; %s)",
          i,
          n_analytes,
          100 * i / n_analytes,
          .format_duration_33(proc.time()[["elapsed"]] - quality_check_started)
        )
      )
    }
  }

  .log_33(
    verbose,
    paste0(
      "omics quality checks: ", n_analytes, "/", n_analytes,
      " analytes complete (", n_with_na, " with NA, ",
      n_with_nzv, " near-zero variance; ",
      .format_duration_33(proc.time()[["elapsed"]] - quality_check_started),
      ")"
    )
  )

  if (n_with_na > 0) {
    warning(n_with_na, " analytes contain NA values")
  }

  if (n_with_nzv > 0) {
    warning(n_with_nzv, " analytes have near-zero variance")
  }

  omics_all <- cbind(ANALYTE_NAME = analyte_names, omics_numeric)

  omics_male <- NULL
  if (!is.null(pheno_list$male)) {
    male_sample_ids <- pheno_list$male$SAMPLE_ID
    male_cols <- intersect(colnames(omics_all), male_sample_ids)
    omics_male <- omics_all[, c("ANALYTE_NAME", male_cols), drop = FALSE]
  }

  omics_female <- NULL
  if (!is.null(pheno_list$female)) {
    female_sample_ids <- pheno_list$female$SAMPLE_ID
    female_cols <- intersect(colnames(omics_all), female_sample_ids)
    omics_female <- omics_all[, c("ANALYTE_NAME", female_cols), drop = FALSE]
  }

  list(
    all = omics_all,
    male = omics_male,
    female = omics_female
  )
}


.validate_dnam_probe_coverage <- function(full_probes, filtered_probes, available_probes) {
  full_present <- sum(full_probes %in% available_probes)
  full_missing <- length(full_probes) - full_present
  if (full_present == 0) {
    stop("No probes from full probe list found in data")
  }
  if (full_missing > 0) {
    warning(sprintf("%d of %d probes from full probe list not found in data",
                    full_missing, length(full_probes)))
  }

  filtered_present <- sum(filtered_probes %in% available_probes)
  filtered_missing <- length(filtered_probes) - filtered_present
  if (filtered_present == 0) {
    stop("No probes from filtered probe list found in data")
  }
  if (filtered_missing > 0) {
    warning(sprintf("%d of %d probes from filtered probe list not found in data",
                    filtered_missing, length(filtered_probes)))
  }
}


.subset_omics_list <- function(omics_list, analyte_subset) {
  if (is.null(analyte_subset)) {
    return(omics_list)
  }

  subsetted <- list()
  for (dataset in c("all", "male", "female")) {
    if (is.null(omics_list[[dataset]])) {
      subsetted[[dataset]] <- NULL
      next
    }

    omics_df <- omics_list[[dataset]]
    matching_analytes <- omics_df$ANALYTE_NAME %in% analyte_subset
    subsetted[[dataset]] <- omics_df[matching_analytes, ]
  }

  subsetted
}
