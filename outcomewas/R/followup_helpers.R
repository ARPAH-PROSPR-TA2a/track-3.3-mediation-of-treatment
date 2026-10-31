# Visit codes identify models and checkpoints; labels describe the outputs.
.validate_outcome_fu_labels <- function(fu_labels, fu_codes) {
  if (is.null(fu_labels)) return(NULL)
  if (!is.character(fu_labels) || !is.null(dim(fu_labels)) ||
      is.null(names(fu_labels)) || anyNA(names(fu_labels)) ||
      anyDuplicated(names(fu_labels)) || !setequal(names(fu_labels), fu_codes)) {
    stop("fu_labels must be a named character vector with unique code names matching every observed FU code exactly, including baseline 0.")
  }
  if (anyNA(fu_labels) || any(!nzchar(trimws(fu_labels))) || anyDuplicated(fu_labels)) {
    stop("fu_labels must contain unique, non-missing, non-blank labels.")
  }
  fu_labels[fu_codes]
}

.attach_outcome_fu_labels <- function(table, fu_labels) {
  if (is.null(table) || is.null(fu_labels)) return(table)
  table$FU_LABEL <- unname(fu_labels[as.character(table$FU)])
  columns <- setdiff(names(table), "FU_LABEL")
  table[, append(columns, "FU_LABEL", after = match("FU", columns)), drop = FALSE]
}
