.supported_analyses <- c("analysis_change", "analysis_level")
.supported_strata <- c("all", "male", "female")

.mediation_arm_codes <- function(codes, context, positive = FALSE) {
  codes <- as.character(codes)
  numbers <- suppressWarnings(as.integer(codes))
  minimum <- if (positive) 1L else 0L
  if (anyNA(numbers) || any(numbers < minimum) || any(as.character(numbers) != codes)) {
    stop(context, " must contain canonical ", if (positive) "positive" else "non-negative",
         " integer treatment codes.")
  }
  codes
}

.validate_mediation_labels <- function(labels, context) {
  if (anyNA(labels) || any(!nzchar(trimws(as.character(labels))))) {
    stop(context, " must contain non-missing, non-blank labels.")
  }
  as.character(labels)
}

.validate_effect_table <- function(table, table_name, context, outcome_type = NULL) {
  context <- paste0(context, "$", table_name)
  if (!is.data.frame(table)) stop(context, " must be a data frame.")
  missing <- setdiff(c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE"), names(table))
  if (length(missing)) stop(context, " is missing required columns: ", paste(missing, collapse = ", "), ".")
  if (nrow(table) == 0L) stop(context, " must contain at least one row.")
  if (anyNA(table$ANALYTE_NAME) || any(!nzchar(as.character(table$ANALYTE_NAME)))) {
    stop(context, "$ANALYTE_NAME must be non-missing and non-empty.")
  }
  table$ANALYTE_NAME <- as.character(table$ANALYTE_NAME)
  fu <- suppressWarnings(as.numeric(as.character(table$FU)))
  if (anyNA(fu) || any(!is.finite(fu)) || any(fu <= 0 | fu != floor(fu)) || any(fu > .Machine$integer.max)) {
    stop(context, "$FU must contain positive integer follow-up levels.")
  }
  table$FU <- as.integer(fu)
  key_columns <- c("ANALYTE_NAME", "FU")
  if (table_name == "treatment_effects") {
    supplied_labels <- intersect(c("TREATMENT", "REFERENCE"), names(table))
    if (!"CONTRAST" %in% names(table)) {
      # Historical Track 1.1.1 summaries supported only the 1-versus-0 study.
      table$CONTRAST <- "1_vs_0"
      if (!"TREATMENT" %in% names(table)) table$TREATMENT <- "Treatment"
      if (!"REFERENCE" %in% names(table)) table$REFERENCE <- "Control"
    } else if (!all(c("TREATMENT", "REFERENCE") %in% names(table))) {
      stop(context, " must supply TREATMENT and REFERENCE labels with CONTRAST.")
    }
    table$CONTRAST <- as.character(table$CONTRAST)
    if (anyNA(table$CONTRAST) || any(!grepl("^[1-9][0-9]*_vs_0$", table$CONTRAST))) {
      stop(context, "$CONTRAST must use stable positive-code-versus-control identifiers, e.g. 2_vs_0.")
    }
    .mediation_arm_codes(sub("_vs_0$", "", table$CONTRAST), paste0(context, "$CONTRAST"), TRUE)
    for (column in c("TREATMENT", "REFERENCE")) {
      table[[column]] <- .validate_mediation_labels(table[[column]], paste0(context, "$", column))
    }
    attr(table, "supplied_labels") <- supplied_labels
    key_columns <- c(key_columns, "CONTRAST")
  } else if ("ARM" %in% names(table)) {
    stop(context, " must contain one common analyte slope per ANALYTE_NAME x FU; ARM-specific beta tables are not supported.")
  }
  if (anyDuplicated(table[key_columns])) {
    stop(context, " contains duplicate ", paste(key_columns, collapse = " x "), " keys.")
  }
  if (!"STATUS" %in% names(table)) table$STATUS <- "ok"
  table$STATUS <- as.character(table$STATUS)
  if (anyNA(table$STATUS) || any(!table$STATUS %in% c("ok", "not_estimable", "fit_error"))) {
    stop(context, "$STATUS must be ok, not_estimable, or fit_error.")
  }
  available <- table$STATUS == "ok"
  for (column in c("EFFECT_SIZE", "SE")) {
    value <- table[[column]]
    if (!is.numeric(value) || any(!is.finite(value[available]))) {
      stop(context, "$", column, " must contain finite numeric values for STATUS = 'ok'.")
    }
    if (any(!is.na(value[!available]))) stop(context, "$", column, " must be NA for unavailable STATUS rows.")
  }
  if (any(table$SE[available] <= 0)) stop(context, "$SE must be strictly positive for STATUS = 'ok'.")
  if ("IN_FILTERED_SET" %in% names(table)) {
    if (!is.logical(table$IN_FILTERED_SET) || anyNA(table$IN_FILTERED_SET)) {
      stop(context, "$IN_FILTERED_SET must contain non-missing logical membership values.")
    }
  } else if ("BH_P_VALUE_FILTERED" %in% names(table)) {
    if (any(!available)) {
      stop(context, " has unavailable inference without IN_FILTERED_SET; regenerate this input with explicit filtered-probe membership.")
    }
    table$IN_FILTERED_SET <- !is.na(table$BH_P_VALUE_FILTERED)
  }
  if (!is.null(outcome_type)) {
    has_hazard_ratio <- "HAZARD_RATIO" %in% names(table)
    if (outcome_type == "tte" && !has_hazard_ratio) stop(context, " must contain HAZARD_RATIO for outcome_type = 'tte'.")
    if (outcome_type == "continuous" && has_hazard_ratio) stop(context, " contains HAZARD_RATIO but outcome_type = 'continuous'.")
  }
  table
}

.mediation_metadata_labels <- function(results, role) {
  labels <- results$metadata$treatment_labels
  if (is.null(labels)) return(NULL)
  if (!is.character(labels) || is.null(names(labels)) || anyDuplicated(names(labels))) {
    stop(role, " metadata$treatment_labels must be a named character vector with unique arm codes.")
  }
  codes <- .mediation_arm_codes(names(labels), paste0(role, " metadata$treatment_labels names"))
  .validate_mediation_labels(labels, paste0(role, " metadata$treatment_labels"))
  if (!"0" %in% codes || length(codes) < 2L || anyDuplicated(labels)) {
    stop(role, " metadata$treatment_labels must name control and active arms with distinct labels.")
  }
  labels
}

# Resolve visit labels separately from numeric joins and inference.
# Baseline may be metadata-only; effect rows remain positive FU.
.mediation_fu_labels <- function(results, role) {
  metadata <- results$metadata$fu_labels
  if (!is.null(metadata)) {
    codes <- names(metadata)
    numbers <- suppressWarnings(as.integer(codes))
    if (!is.character(metadata) || !is.null(dim(metadata)) || !length(metadata) ||
        is.null(codes) || anyNA(numbers) || any(numbers < 0L) ||
        any(as.character(numbers) != codes) || anyDuplicated(codes)) {
      stop(role, " metadata$fu_labels must be a named character vector with unique non-negative integer FU codes.")
    }
    .validate_mediation_labels(metadata, paste0(role, " metadata$fu_labels"))
  }

  effect_name <- if (role == "treatment") "treatment_effects" else "outcome_effects"
  observed <- character()
  rows <- list()
  for (analysis in .supported_analyses) {
    for (stratum in intersect(.supported_strata, names(results[[analysis]]))) {
      table <- results[[analysis]][[stratum]][[effect_name]]
      if (is.null(table)) next
      codes <- as.character(table$FU)
      observed <- union(observed, codes)
      if (!"FU_LABEL" %in% names(table)) next
      context <- paste(role, analysis, stratum, "FU_LABEL")
      if (!is.character(table$FU_LABEL)) stop(context, " must be character.")
      labels <- .validate_mediation_labels(table$FU_LABEL, context)
      rows[[length(rows) + 1L]] <- unique(data.frame(FU = codes, FU_LABEL = labels))
    }
  }
  if (!is.null(metadata) && any(!observed %in% names(metadata))) {
    stop(role, " metadata$fu_labels must cover every observed FU code in its effect tables.")
  }
  row_map <- NULL
  if (length(rows)) {
    rows <- unique(do.call(rbind, rows))
    if (anyDuplicated(rows$FU)) {
      stop(role, " FU_LABEL must map each FU consistently across analyses and strata.")
    }
    row_map <- setNames(rows$FU_LABEL, rows$FU)
  }
  if (!is.null(metadata) && !is.null(row_map) &&
      any(unname(metadata[names(row_map)]) != unname(row_map))) {
    stop(role, " FU_LABEL disagrees with metadata$fu_labels.")
  }
  labels <- if (!is.null(metadata)) metadata else row_map
  if (is.null(labels)) return(NULL)
  labels[order(as.integer(names(labels)))]
}

.resolve_mediation_fu_labels <- function(treatment_results, outcome_results) {
  alpha <- .mediation_fu_labels(treatment_results, "treatment")
  beta <- .mediation_fu_labels(outcome_results, "outcome")
  if (xor(is.null(alpha), is.null(beta))) {
    source <- if (is.null(alpha)) "beta (outcome)" else "alpha (treatment)"
    warning("Follow-up labels are supplied only by ", source,
            "; continuing with numeric FU matching and inheriting those labels. ",
            "Label agreement could not be checked.", call. = FALSE)
  }
  common <- intersect(names(alpha), names(beta))
  if (any(unname(alpha[common]) != unname(beta[common]))) {
    stop("Follow-up labels disagree between alpha and beta for FU: ",
         paste(common[unname(alpha[common]) != unname(beta[common])], collapse = ", "), ".")
  }
  labels <- c(alpha, beta[setdiff(names(beta), names(alpha))])
  if (is.null(labels)) return(NULL)
  labels[order(as.integer(names(labels)))]
}

.validate_results_object <- function(results, role, outcome_type = NULL) {
  if (!is.list(results)) stop(role, " results must be a list.")
  missing_analyses <- setdiff(.supported_analyses, names(results))
  if (length(missing_analyses)) stop(role, " results are missing: ", paste(missing_analyses, collapse = ", "), ".")
  effect_name <- if (role == "treatment") "treatment_effects" else "outcome_effects"
  if (role == "outcome") {
    # Guard against silently reusing results from the former interaction model.
    if ("arm_specific" %in% results$metadata$beta_model) {
      stop("Outcome metadata$beta_model = 'arm_specific' is not supported; supply one common analyte slope per ANALYTE_NAME x FU.")
    }
    slope <- results$metadata$analyte_slope
    if (!is.null(slope) && !identical(slope, "common_across_treatment_arms")) {
      stop("Outcome metadata$analyte_slope must be 'common_across_treatment_arms'.")
    }
  }
  label_map <- .mediation_metadata_labels(results, role)
  maps <- list()
  for (analysis in .supported_analyses) {
    analysis_results <- results[[analysis]]
    if (!is.list(analysis_results)) stop(role, " results$", analysis, " must be a list.")
    available_names <- names(analysis_results)
    if (is.null(available_names) || !any(available_names %in% .supported_strata)) {
      stop(role, " results$", analysis, " must name at least one supported stratum: ", paste(.supported_strata, collapse = ", "), ".")
    }
    unexpected <- setdiff(available_names, .supported_strata)
    if (length(unexpected)) stop(role, " results$", analysis, " contains unsupported strata: ", paste(unexpected, collapse = ", "), ".")
    for (stratum in intersect(.supported_strata, available_names)) {
      entry <- analysis_results[[stratum]]
      if (is.null(entry)) next
      context <- paste0(role, " results$", analysis, "$", stratum)
      if (!is.list(entry) || is.null(entry[[effect_name]])) stop(context, " must contain $", effect_name, ".")
      table <- .validate_effect_table(entry[[effect_name]], effect_name, context,
                                     if (role == "outcome") outcome_type else NULL)
      if (role == "treatment") {
        map <- unique(table[c("CONTRAST", "TREATMENT", "REFERENCE")])
        map$ARM <- sub("_vs_0$", "", map$CONTRAST)
        if (!is.null(label_map)) {
          if (any(!map$ARM %in% names(label_map))) stop(context, " contrasts disagree with treatment metadata arm codes.")
          supplied <- attr(table, "supplied_labels")
          if (("TREATMENT" %in% supplied && any(map$TREATMENT != unname(label_map[map$ARM]))) ||
              ("REFERENCE" %in% supplied && any(map$REFERENCE != unname(label_map["0"])))) {
            stop(context, " treatment labels disagree with metadata.")
          }
          table$TREATMENT <- unname(label_map[sub("_vs_0$", "", table$CONTRAST)])
          table$REFERENCE <- unname(label_map["0"])
          attr(table, "supplied_labels") <- c("TREATMENT", "REFERENCE")
          map <- unique(table[c("CONTRAST", "TREATMENT", "REFERENCE")])
          map$ARM <- sub("_vs_0$", "", map$CONTRAST)
        }
        maps[[length(maps) + 1L]] <- map
      }
      results[[analysis]][[stratum]][[effect_name]] <- table
    }
  }
  if (length(maps)) {
    map <- unique(do.call(rbind, maps))
    if (anyDuplicated(map$ARM) || anyDuplicated(map$TREATMENT) ||
        (role == "treatment" && (length(unique(map$REFERENCE)) != 1L || any(map$TREATMENT == map$REFERENCE)))) {
      stop(role, " treatment labels must map each arm consistently and distinctly across results.")
    }
    map <- map[order(as.integer(map$ARM)), , drop = FALSE]
    row.names(map) <- NULL
    attr(results, "treatment_map") <- map
  }
  results
}
