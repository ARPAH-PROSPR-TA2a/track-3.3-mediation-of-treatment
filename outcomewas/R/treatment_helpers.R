# Codes identify arms throughout the models; labels describe them in the outputs.
.prepare_treatments <- function(pheno, treatment_labels = NULL) {
  codes <- as.character(pheno$TREATMENT_GROUP)
  numbers <- suppressWarnings(as.integer(codes))
  if (anyNA(numbers) || any(numbers < 0L) || any(as.character(numbers) != codes)) {
    stop("TREATMENT_GROUP must contain non-missing non-negative integer codes: 0 = control; 1, 2, ... = active treatments.")
  }
  arms <- as.character(sort(unique(numbers)))
  if (!"0" %in% arms || length(arms) < 2L) {
    stop("pheno must contain control (TREATMENT_GROUP == 0) and at least one active treatment.")
  }

  # Validate before dropping duplicate visits or incomplete covariates.
  by_subject <- split(codes, pheno$SUBJECT_ID, drop = TRUE)
  if (any(vapply(by_subject, function(x) length(unique(x)) != 1L, logical(1)))) {
    stop("TREATMENT_GROUP must be constant within SUBJECT_ID.")
  }
  if (is.null(treatment_labels)) {
    if (length(arms) > 2L) {
      stop('Multiple active treatments require treatment_labels covering every code, including control 0; e.g., c(`0` = "Control", `1` = "Diet", `2` = "Exercise").')
    }
    treatment_labels <- setNames(c("Control", "Treatment"), arms)
  }
  if (!is.character(treatment_labels) || !is.null(dim(treatment_labels)) ||
      is.null(names(treatment_labels)) || anyNA(names(treatment_labels)) ||
      anyDuplicated(names(treatment_labels)) || !setequal(names(treatment_labels), arms)) {
    stop("treatment_labels must be a named character vector with unique code names matching every observed TREATMENT_GROUP code exactly, including 0.")
  }
  if (anyNA(treatment_labels) || any(!nzchar(trimws(treatment_labels))) ||
      anyDuplicated(treatment_labels)) {
    stop("treatment_labels must contain unique, non-missing, non-blank labels.")
  }
  treatment_labels <- treatment_labels[arms]
  active_arms <- setdiff(arms, "0")
  treatments <- data.frame(
    ARM = active_arms, TREATMENT = unname(treatment_labels[active_arms]),
    REFERENCE = unname(treatment_labels["0"]),
    CONTRAST = paste0(active_arms, "_vs_0"), stringsAsFactors = FALSE
  )
  pheno$TREATMENT_GROUP <- factor(codes, levels = arms)
  list(pheno = pheno, treatments = treatments, treatment_labels = treatment_labels)
}

# Keep the code/label contract above aligned with Track 1.1.1.
.outcome_metadata <- function(pheno_list, additional_covariates = NULL) {
  response <- if (pheno_list$outcome_type == "continuous") "OUTCOME" else {
    "survival::Surv(OUTCOME_TIME, OUTCOME_STATUS)"
  }
  predictors <- "analyte + analyte_baseline + TREATMENT_GROUP"
  list(
    schema_version = "2.0",
    treatment_labels = pheno_list$treatment_labels,
    treatments = pheno_list$treatments,
    outcome_type = pheno_list$outcome_type,
    additional_covariates = additional_covariates,
    analyte_slope = "common_across_treatment_arms",
    p_adjust = "BH within analysis, stratum and FU",
    model_formula = paste(response, "~", paste(c(predictors, "FEMALE", additional_covariates),
                                              collapse = " + ")),
    inference = if (pheno_list$outcome_type == "continuous") {
      "Model covariance; t tests with pooled model residual degrees of freedom"
    } else {
      "Model covariance; Wald z tests; Efron ties"
    },
    variance_assumption = if (pheno_list$outcome_type == "continuous") {
      "Common residual variance across treatment arms"
    } else {
      "Common Cox baseline hazard across treatment arms"
    }
  )
}
