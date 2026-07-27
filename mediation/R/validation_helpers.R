.supported_analyses <- c("analysis_change", "analysis_level")
.supported_strata <- c("all", "male", "female")


.validate_effect_table <- function(table, table_name, context, outcome_type = NULL) {
  if (!is.data.frame(table)) {
    stop(context, "$", table_name, " must be a data frame.")
  }

  required <- c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE")
  missing <- setdiff(required, names(table))
  if (length(missing) > 0L) {
    stop(
      context, "$", table_name, " is missing required columns: ",
      paste(missing, collapse = ", "), "."
    )
  }

  if (nrow(table) == 0L) {
    stop(context, "$", table_name, " must contain at least one row.")
  }
  if (anyNA(table$ANALYTE_NAME) || any(!nzchar(as.character(table$ANALYTE_NAME)))) {
    stop(context, "$", table_name, "$ANALYTE_NAME must be non-missing and non-empty.")
  }

  fu <- suppressWarnings(as.integer(as.character(table$FU)))
  if (anyNA(fu) || any(fu <= 0L)) {
    stop(context, "$", table_name, "$FU must contain positive integer follow-up levels.")
  }

  key <- paste(as.character(table$ANALYTE_NAME), fu, sep = "\r")
  if (anyDuplicated(key)) {
    stop(context, "$", table_name, " contains duplicate ANALYTE_NAME x FU keys.")
  }

  for (column in c("EFFECT_SIZE", "SE")) {
    value <- table[[column]]
    if (!is.numeric(value) || any(!is.finite(value))) {
      stop(context, "$", table_name, "$", column, " must contain finite numeric values.")
    }
  }
  if (any(table$SE <= 0)) {
    stop(context, "$", table_name, "$SE must be strictly positive.")
  }

  if (!is.null(outcome_type)) {
    has_hazard_ratio <- "HAZARD_RATIO" %in% names(table)
    if (outcome_type == "tte" && !has_hazard_ratio) {
      stop(context, "$", table_name, " must contain HAZARD_RATIO for outcome_type = 'tte'.")
    }
    if (outcome_type == "continuous" && has_hazard_ratio) {
      stop(context, "$", table_name, " contains HAZARD_RATIO but outcome_type = 'continuous'.")
    }
  }

  invisible(TRUE)
}


.validate_results_object <- function(results, role, outcome_type = NULL) {
  if (!is.list(results)) {
    stop(role, " results must be a list.")
  }

  missing_analyses <- setdiff(.supported_analyses, names(results))
  if (length(missing_analyses) > 0L) {
    stop(role, " results are missing: ", paste(missing_analyses, collapse = ", "), ".")
  }

  effect_name <- if (role == "treatment") "treatment_effects" else "outcome_effects"

  for (analysis in .supported_analyses) {
    analysis_results <- results[[analysis]]
    if (!is.list(analysis_results)) {
      stop(role, " results$", analysis, " must be a list.")
    }

    available_names <- names(analysis_results)
    if (is.null(available_names) || !any(available_names %in% .supported_strata)) {
      stop(
        role, " results$", analysis,
        " must name at least one supported stratum: ",
        paste(.supported_strata, collapse = ", "), "."
      )
    }
    unexpected_strata <- setdiff(available_names, .supported_strata)
    if (length(unexpected_strata) > 0L) {
      stop(
        role, " results$", analysis, " contains unsupported strata: ",
        paste(unexpected_strata, collapse = ", "), "."
      )
    }

    for (stratum in intersect(.supported_strata, available_names)) {
      entry <- analysis_results[[stratum]]
      if (is.null(entry)) next
      context <- paste0(role, " results$", analysis, "$", stratum)
      if (!is.list(entry) || is.null(entry[[effect_name]])) {
        stop(context, " must contain $", effect_name, ".")
      }
      .validate_effect_table(
        entry[[effect_name]],
        effect_name,
        context,
        outcome_type = if (role == "outcome") outcome_type else NULL
      )
    }
  }

  invisible(TRUE)
}
