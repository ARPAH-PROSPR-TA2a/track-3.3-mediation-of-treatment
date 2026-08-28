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
  if (anyDuplicated(fu_labels)) {
    stop("fu_labels must contain unique labels.")
  }

  names(fu_labels) <- normalized_names
  fu_labels
}


.validate_protein_annotation <- function(protein_annotation) {
  if (!is.data.frame(protein_annotation)) {
    stop("protein_annotation must be a data.frame.")
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


.select_outcome_effects <- function(results, analysis_type, group) {
  analysis_name <- paste0("analysis_", analysis_type)

  if (!is.list(results)) {
    stop("results must be a FAST_outcome_WAS result list.")
  }
  if (!analysis_name %in% names(results) || !is.list(results[[analysis_name]])) {
    stop("results is missing the required list: ", analysis_name)
  }
  if (!group %in% names(results[[analysis_name]])) {
    stop("results$", analysis_name, " is missing group: ", group)
  }

  result_group <- results[[analysis_name]][[group]]
  if (is.null(result_group)) {
    return(NULL)
  }
  if (!is.list(result_group) || !"outcome_effects" %in% names(result_group)) {
    stop("results$", analysis_name, "$", group, " is missing outcome_effects.")
  }

  outcome_effects <- result_group$outcome_effects
  if (!is.data.frame(outcome_effects)) {
    stop(
      "results$", analysis_name, "$", group,
      "$outcome_effects must be a data.frame."
    )
  }

  outcome_effects
}


get_annotated_outcome_effects <- function(
    results,
    protein_annotation,
    analysis_type = c("change", "level"),
    group = c("all", "male", "female"),
    fu = NULL,
    fu_labels = NULL
) {
  analysis_type <- match.arg(analysis_type)
  group <- match.arg(group)
  outcome_effects <- .select_outcome_effects(results, analysis_type, group)
  if (is.null(outcome_effects)) {
    return(NULL)
  }

  required_effect_columns <- c(
    "ANALYTE_NAME", "FU", "EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE"
  )
  missing_effect_columns <- setdiff(required_effect_columns, names(outcome_effects))
  if (length(missing_effect_columns) > 0L) {
    stop(
      "outcome_effects is missing required columns: ",
      paste(missing_effect_columns, collapse = ", ")
    )
  }

  analyte_names <- as.character(outcome_effects$ANALYTE_NAME)
  if (any(is.na(analyte_names)) || any(trimws(analyte_names) == "")) {
    stop("outcome_effects$ANALYTE_NAME must not contain missing or blank values.")
  }

  outcome_fu <- .as_positive_integer_fu(outcome_effects$FU, "outcome_effects$FU")
  result_keys <- data.frame(
    ANALYTE_NAME = analyte_names,
    FU = outcome_fu,
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(result_keys)) {
    stop("outcome_effects must have unique ANALYTE_NAME x FU rows.")
  }

  if (!is.null(fu)) {
    requested_fu <- unique(.as_positive_integer_fu(fu, "fu"))
    absent_fu <- setdiff(requested_fu, unique(outcome_fu))
    if (length(absent_fu) > 0L) {
      stop(
        "Requested fu values are not present in outcome_effects: ",
        paste(absent_fu, collapse = ", ")
      )
    }

    keep <- outcome_fu %in% requested_fu
    outcome_effects <- outcome_effects[keep, , drop = FALSE]
    analyte_names <- analyte_names[keep]
    outcome_fu <- outcome_fu[keep]
  }

  validated_labels <- .validate_fu_labels(fu_labels)
  if (!is.null(validated_labels)) {
    missing_labels <- setdiff(unique(outcome_fu), as.integer(names(validated_labels)))
    if (length(missing_labels) > 0L) {
      stop(
        "fu_labels is missing labels for FU: ",
        paste(missing_labels, collapse = ", ")
      )
    }
  }

  annotation_keys <- .validate_protein_annotation(protein_annotation)
  matched_rows <- match(analyte_names, annotation_keys)
  unique_analytes <- unique(analyte_names)
  unique_matches <- match(unique_analytes, annotation_keys)

  if (length(unique_analytes) > 0L && all(is.na(unique_matches))) {
    stop("No unique outcome analytes matched protein_annotation$AptName.")
  }

  unmatched_analytes <- unique_analytes[is.na(unique_matches)]
  if (length(unmatched_analytes) > 0L) {
    warning(
      length(unmatched_analytes), " of ", length(unique_analytes),
      " unique outcome analytes were not found in protein_annotation; ",
      "rows were retained with NA annotation values."
    )
  }

  annotation <- protein_annotation[
    matched_rows,
    .protein_annotation_columns,
    drop = FALSE
  ]
  remaining_columns <- setdiff(
    names(outcome_effects),
    c("ANALYTE_NAME", "FU")
  )

  output_parts <- list(
    outcome_effects[, "ANALYTE_NAME", drop = FALSE],
    annotation,
    outcome_effects[, "FU", drop = FALSE]
  )
  if (!is.null(validated_labels)) {
    output_parts[[length(output_parts) + 1L]] <- data.frame(
      FU_LABEL = unname(validated_labels[as.character(outcome_fu)]),
      stringsAsFactors = FALSE
    )
  }
  output_parts[[length(output_parts) + 1L]] <- outcome_effects[
    , remaining_columns,
    drop = FALSE
  ]

  annotated <- do.call(cbind, output_parts)
  row.names(annotated) <- NULL

  if (nrow(annotated) != nrow(outcome_effects)) {
    stop("Internal error: protein annotation changed outcome_effects row count.")
  }

  annotated
}


write_annotated_outcome_effect_tables <- function(
    results,
    protein_annotation,
    output_rds,
    fu_labels = NULL
) {
  if (!is.character(output_rds) ||
      length(output_rds) != 1L ||
      is.na(output_rds) ||
      trimws(output_rds) == "") {
    stop("output_rds must be one nonblank file path.")
  }

  output_rds <- path.expand(output_rds)
  output_parent <- dirname(output_rds)
  if (!dir.exists(output_parent)) {
    stop("The output_rds parent directory does not exist: ", output_parent)
  }

  validated_labels <- .validate_fu_labels(fu_labels)
  tables <- list()

  for (analysis_type in c("change", "level")) {
    for (group in c("all", "male", "female")) {
      outcome_effects <- .select_outcome_effects(results, analysis_type, group)
      if (is.null(outcome_effects)) {
        next
      }
      if (!"FU" %in% names(outcome_effects)) {
        stop("outcome_effects is missing required columns: FU")
      }

      available_fu <- sort(unique(
        .as_positive_integer_fu(outcome_effects$FU, "outcome_effects$FU")
      ))
      for (fu_value in available_fu) {
        table <- get_annotated_outcome_effects(
          results = results,
          protein_annotation = protein_annotation,
          analysis_type = analysis_type,
          group = group,
          fu = fu_value,
          fu_labels = validated_labels
        )

        fu_suffix <- if (is.null(validated_labels)) {
          paste0("FU", fu_value)
        } else {
          validated_labels[[as.character(fu_value)]]
        }
        table_name <- paste(analysis_type, group, fu_suffix, sep = "_")
        tables[[table_name]] <- table
      }
    }
  }

  saveRDS(tables, output_rds)
  invisible(tables)
}
