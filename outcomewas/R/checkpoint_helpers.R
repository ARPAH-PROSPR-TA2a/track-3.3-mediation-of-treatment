# Each batch stores one result per analyte, before labels and BH correction.
.read_outcome_checkpoint <- function(batch_file, analyte_names, fu_level,
                                     outcome_type, checkpoint_state = NULL) {
  if (is.null(batch_file) || !file.exists(batch_file)) return(NULL)

  invalid <- function(reason) {
    stop("Cannot reuse checkpoint '", batch_file, "': ", reason,
         ". Use a fresh checkpoint_dir for this analysis.", call. = FALSE)
  }
  payload <- tryCatch(readRDS(batch_file), error = function(e) {
    invalid(paste("could not read the file -", conditionMessage(e)))
  })
  if (!is.list(payload) || is.data.frame(payload) || !is.null(names(payload)) ||
      length(payload) != length(analyte_names) || length(payload) == 0L) {
    invalid("unrecognized batch structure or incorrect number of analyte results")
  }

  coefficient_columns <- c("ANALYTE_NAME", "FU", "COEFFICIENT", "N_OBS",
                           "EFFECT_SIZE", "SE", "P_VALUE")
  legacy_effect_columns <- c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE", "P_VALUE")
  if (outcome_type == "tte") {
    coefficient_columns <- c(coefficient_columns, "N_EVENTS", "HAZARD_RATIO")
    legacy_effect_columns <- c(legacy_effect_columns, "HAZARD_RATIO")
  }
  effect_columns <- c(legacy_effect_columns, "STATUS")

  has_columns <- function(df, columns) {
    is.data.frame(df) && identical(sort(names(df)), sort(columns))
  }
  valid_table_types <- function(df, columns) {
    numeric_columns <- setdiff(columns, c("ANALYTE_NAME", "COEFFICIENT", "STATUS"))
    is.atomic(df$ANALYTE_NAME) && is.null(dim(df$ANALYTE_NAME)) &&
      !anyNA(df$ANALYTE_NAME) && all(vapply(df[numeric_columns], function(x) {
        is.numeric(x) && is.null(dim(x))
      }, logical(1)))
  }

  legacy_found <- FALSE
  for (i in seq_along(payload)) {
    result <- payload[[i]]
    # Historical fit failures occupied a NULL entry in the analyte list.
    if (is.null(result)) {
      legacy_found <- TRUE
      next
    }
    if (!is.list(result) || is.data.frame(result) || length(result) != 2L ||
        !setequal(names(result), c("coefficients", "outcome_effects"))) {
      invalid(paste("unrecognized result structure at analyte position", i))
    }

    coefficients <- result$coefficients
    effects <- result$outcome_effects
    legacy_effect <- has_columns(effects, legacy_effect_columns)
    if (!has_columns(coefficients, coefficient_columns) ||
        (!legacy_effect && !has_columns(effects, effect_columns)) ||
        nrow(effects) != 1L || (legacy_effect && nrow(coefficients) == 0L)) {
      invalid(paste("unrecognized coefficient or outcome-effect schema at analyte position", i))
    }
    if (!valid_table_types(coefficients, coefficient_columns) ||
        !valid_table_types(effects, legacy_effect_columns) ||
        !is.character(coefficients$COEFFICIENT) ||
        !is.null(dim(coefficients$COEFFICIENT)) || anyNA(coefficients$COEFFICIENT) ||
        anyDuplicated(coefficients$COEFFICIENT)) {
      invalid(paste("invalid column types or coefficient names at analyte position", i))
    }

    expected_analyte <- as.character(analyte_names[i])
    expected_fu <- as.numeric(as.character(fu_level))
    if (anyNA(coefficients$FU) || anyNA(effects$FU) ||
        !all(as.character(coefficients$ANALYTE_NAME) == expected_analyte) ||
        as.character(effects$ANALYTE_NAME) != expected_analyte ||
        !all(coefficients$FU == expected_fu) || effects$FU != expected_fu) {
      invalid(paste("analyte or follow-up does not match batch position", i))
    }

    if (legacy_effect) {
      legacy_found <- TRUE
    } else {
      if (!is.character(effects$STATUS) || !is.null(dim(effects$STATUS)) ||
          anyNA(effects$STATUS) ||
          !effects$STATUS %in% c("ok", "not_estimable", "fit_error")) {
        invalid(paste("unrecognized outcome-effect status at analyte position", i))
      }
      if (effects$STATUS == "ok" && !"analyte" %in% coefficients$COEFFICIENT) {
        invalid(paste("successful outcome effect has no analyte coefficient at position", i))
      }
    }
  }

  # Check every entry before accepting a legacy batch for recomputation.
  if (legacy_found) {
    if (is.null(checkpoint_state) || !isTRUE(checkpoint_state$notified)) {
      directory <- if (is.null(checkpoint_state)) dirname(batch_file) else checkpoint_state$directory
      message("Older-format checkpoints detected in '", directory,
              "'. FAST will recompute affected batches and replace them as they complete.")
      if (!is.null(checkpoint_state)) checkpoint_state$notified <- TRUE
    }
    return(NULL)
  }
  payload
}


.save_outcome_checkpoint <- function(batch_results, batch_file) {
  if (is.null(batch_file)) return(invisible(NULL))

  temporary_file <- paste0(batch_file, ".tmp")
  on.exit(unlink(temporary_file), add = TRUE)
  tryCatch(saveRDS(batch_results, temporary_file), error = function(e) {
    stop("Could not write checkpoint '", batch_file, "': ", conditionMessage(e),
         call. = FALSE)
  })
  if (!file.rename(temporary_file, batch_file)) {
    stop("Could not replace checkpoint '", batch_file,
         "'. The previous checkpoint was retained.", call. = FALSE)
  }
  invisible(NULL)
}
