.protein_annotation_columns <- c(
  "SomaId",
  "TargetFullName",
  "Target",
  "UniProt",
  "EntrezGeneID",
  "EntrezGeneSymbol",
  "Uniprot_Unique",
  "Symbol_Unique"
)


.as_positive_integer_fu <- function(x, argument) {
  if (length(x) == 0L || !is.atomic(x)) {
    stop(argument, " must contain one or more positive integer values.")
  }

  numeric_x <- suppressWarnings(as.numeric(as.character(x)))
  invalid <- is.na(numeric_x) |
    !is.finite(numeric_x) |
    numeric_x <= 0 |
    numeric_x != floor(numeric_x)

  if (any(invalid)) {
    stop(argument, " must contain only positive integer values.")
  }

  as.integer(numeric_x)
}


.validate_fu_labels <- function(fu_labels) {
  if (is.null(fu_labels)) {
    return(NULL)
  }

  if (!is.character(fu_labels) || length(fu_labels) == 0L) {
    stop("fu_labels must be NULL or a non-empty named character vector.")
  }

  label_names <- names(fu_labels)
  if (is.null(label_names) ||
      any(is.na(label_names)) ||
      any(trimws(label_names) == "")) {
    stop("fu_labels must have a nonblank name for every label.")
  }

  label_fu <- .as_positive_integer_fu(label_names, "names(fu_labels)")
  normalized_names <- as.character(label_fu)
  if (anyDuplicated(normalized_names)) {
    stop("fu_labels must have unique FU names.")
  }

  if (any(is.na(fu_labels)) || any(trimws(fu_labels) == "")) {
    stop("fu_labels must contain only nonblank labels.")
  }

  names(fu_labels) <- normalized_names
  fu_labels
}


.validate_protein_annotation <- function(protein_annotation) {
  if (!is.data.frame(protein_annotation)) {
    stop("protein_annotation must be a data.frame.")
  }
  if (anyDuplicated(names(protein_annotation))) {
    stop("protein_annotation column names must be unique.")
  }

  required_columns <- c("AptName", .protein_annotation_columns)
  missing_columns <- setdiff(required_columns, names(protein_annotation))
  if (length(missing_columns) > 0L) {
    stop(
      "protein_annotation is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  apt_names <- as.character(protein_annotation$AptName)
  if (any(is.na(apt_names)) || any(trimws(apt_names) == "")) {
    stop("protein_annotation$AptName must not contain missing or blank values.")
  }
  if (anyDuplicated(apt_names)) {
    duplicated_name <- apt_names[duplicated(apt_names)][1L]
    stop(
      "protein_annotation$AptName must be unique; duplicate value: ",
      duplicated_name
    )
  }

  invisible(apt_names)
}


.annotate_result_table <- function(
    result_table,
    table_name,
    protein_annotation,
    annotation_keys,
    fu_labels = NULL
) {
  table_name <- match.arg(table_name, c("coefficients", "outcome_effects"))

  if (!is.data.frame(result_table)) {
    stop(table_name, " must be a data.frame.")
  }
  if (anyDuplicated(names(result_table))) {
    stop(table_name, " column names must be unique.")
  }

  reserved_columns <- c(.protein_annotation_columns, "FU_LABEL")
  conflicting_columns <- intersect(names(result_table), reserved_columns)
  if (length(conflicting_columns) > 0L) {
    stop(
      table_name, " contains reserved annotation columns: ",
      paste(conflicting_columns, collapse = ", ")
    )
  }

  required_columns <- if (table_name == "coefficients") {
    c(
      "ANALYTE_NAME", "FU", "COEFFICIENT", "N_OBS",
      "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE"
    )
  } else {
    c("ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")
  }
  missing_columns <- setdiff(required_columns, names(result_table))
  if (length(missing_columns) > 0L) {
    stop(
      table_name, " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  analyte_names <- as.character(result_table$ANALYTE_NAME)
  if (any(is.na(analyte_names)) || any(trimws(analyte_names) == "")) {
    stop(table_name, "$ANALYTE_NAME must not contain missing or blank values.")
  }

  result_fu <- .as_positive_integer_fu(
    result_table$FU,
    paste0(table_name, "$FU")
  )
  available_fu <- sort(unique(result_fu))

  if (!is.null(fu_labels)) {
    missing_labels <- setdiff(available_fu, as.integer(names(fu_labels)))
    if (length(missing_labels) > 0L) {
      stop(
        "fu_labels is missing labels for FU: ",
        paste(missing_labels, collapse = ", ")
      )
    }
  }

  key_columns <- c("ANALYTE_NAME", "FU")
  key_values <- data.frame(
    ANALYTE_NAME = analyte_names,
    FU = result_fu,
    stringsAsFactors = FALSE
  )
  if (table_name == "coefficients") {
    key_columns <- c(key_columns, "COEFFICIENT")
    key_values$COEFFICIENT <- result_table$COEFFICIENT
  }
  if (anyDuplicated(key_values)) {
    stop(table_name, " must have unique ", paste(key_columns, collapse = " x "), " rows.")
  }

  matched_rows <- match(analyte_names, annotation_keys)
  zero_coverage_fu <- available_fu[!vapply(
    available_fu,
    function(fu_value) any(!is.na(matched_rows[result_fu == fu_value])),
    logical(1)
  )]
  if (length(zero_coverage_fu) > 0L) {
    stop(
      "No unique ", table_name,
      " analytes matched protein_annotation$AptName for FU: ",
      paste(zero_coverage_fu, collapse = ", ")
    )
  }

  unmatched_analytes <- unique(analyte_names[is.na(matched_rows)])
  if (length(unmatched_analytes) > 0L) {
    warning(
      length(unmatched_analytes), " unique ", table_name,
      " analytes were not found in protein_annotation; ",
      "rows were retained with NA annotation values."
    )
  }

  annotation <- protein_annotation[
    matched_rows,
    .protein_annotation_columns,
    drop = FALSE
  ]
  remaining_columns <- setdiff(
    names(result_table),
    c("ANALYTE_NAME", "FU")
  )

  output_parts <- list(
    result_table[, "ANALYTE_NAME", drop = FALSE],
    annotation,
    result_table[, "FU", drop = FALSE]
  )
  if (!is.null(fu_labels)) {
    output_parts[[length(output_parts) + 1L]] <- data.frame(
      FU_LABEL = unname(fu_labels[as.character(result_fu)]),
      stringsAsFactors = FALSE
    )
  }
  output_parts[[length(output_parts) + 1L]] <- result_table[
    , remaining_columns,
    drop = FALSE
  ]

  annotated <- do.call(cbind, output_parts)
  row.names(annotated) <- NULL

  if (nrow(annotated) != nrow(result_table)) {
    stop("Internal error: protein annotation changed ", table_name, " row count.")
  }

  annotated
}


get_annotated_results <- function(
    results,
    protein_annotation,
    fu_labels = NULL
) {
  if (!is.list(results)) {
    stop("results must be a FAST_outcome_WAS result list.")
  }

  validated_labels <- .validate_fu_labels(fu_labels)
  annotation_keys <- .validate_protein_annotation(protein_annotation)
  annotated_results <- results
  annotated_table_count <- 0L

  for (analysis_name in c("analysis_change", "analysis_level")) {
    if (!analysis_name %in% names(results) || !is.list(results[[analysis_name]])) {
      stop("results is missing the required list: ", analysis_name)
    }

    for (group in c("all", "male", "female")) {
      if (!group %in% names(results[[analysis_name]])) {
        stop("results$", analysis_name, " is missing group: ", group)
      }

      result_group <- results[[analysis_name]][[group]]
      if (is.null(result_group)) {
        next
      }
      if (!is.list(result_group)) {
        stop("results$", analysis_name, "$", group, " must be a list or NULL.")
      }

      for (table_name in c("coefficients", "outcome_effects")) {
        if (!table_name %in% names(result_group)) {
          stop(
            "results$", analysis_name, "$", group,
            " is missing ", table_name, "."
          )
        }

        annotated_results[[analysis_name]][[group]][[table_name]] <-
          .annotate_result_table(
            result_table = result_group[[table_name]],
            table_name = table_name,
            protein_annotation = protein_annotation,
            annotation_keys = annotation_keys,
            fu_labels = validated_labels
          )
        annotated_table_count <- annotated_table_count + 1L
      }
    }
  }

  if (annotated_table_count == 0L) {
    stop("results does not contain any available result tables.")
  }

  annotated_results
}
