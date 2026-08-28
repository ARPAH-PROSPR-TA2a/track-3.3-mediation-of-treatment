translation_file <- "~/FAST/Data/CALERIE/Raw/Proteomics/CALERIE_cleaned_protein_translation_table.csv"
results_file <- "~/FAST/Outputs/3.3/Proteomics_3.3OWAS_INF/Proteomics_3.3OWAS_INF_results.rds"
reports_file <- "~/FAST/Outputs/3.3/Proteomics_3.3OWAS_INF/Proteomics_3.3OWAS_INF_reports.rds"
annotated_results_file <- "~/FAST/Outputs/3.3/Proteomics_3.3OWAS_INF/Proteomics_3.3OWAS_INF_results_annotated.rds"
annotated_reports_file <- "~/FAST/Outputs/3.3/Proteomics_3.3OWAS_INF/Proteomics_3.3OWAS_INF_reports_annotated.rds"

annotation_cols <- c(
  "AptName",
  # "SeqId",
  # "SomaId",
  # "TargetFullName",
  # "Target",
  # "UniProt",
  # "EntrezGeneID",
  # "EntrezGeneSymbol",
  "Uniprot_Unique",
  "Symbol_Unique"
)

read_annotation <- function(path) {
  protein_annotation <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  
  missing_cols <- setdiff(annotation_cols, names(protein_annotation))
  if (length(missing_cols) > 0) {
    stop("Translation table is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (anyDuplicated(protein_annotation$AptName)) {
    duplicated_apt_names <- unique(protein_annotation$AptName[duplicated(protein_annotation$AptName)])
    stop("AptName is not unique in translation table. Example duplicate: ", duplicated_apt_names[[1]])
  }
  
  protein_annotation[annotation_cols]
}

annotate_df <- function(df, protein_annotation) {
  if (!is.data.frame(df) || !"ANALYTE_NAME" %in% names(df)) {
    return(df)
  }
  
  matched_rows <- match(df$ANALYTE_NAME, protein_annotation$AptName)
  metadata <- protein_annotation[matched_rows, setdiff(annotation_cols, "AptName"), drop = FALSE]
  names(metadata) <- paste0("ANNOT_", names(metadata))
  
  insert_after <- match("ANALYTE_NAME", names(df))
  cbind(
    df[seq_len(insert_after)],
    metadata,
    df[setdiff(names(df), names(df)[seq_len(insert_after)])],
    stringsAsFactors = FALSE
  )
}

annotate_recursive <- function(x, protein_annotation) {
  if (is.data.frame(x)) {
    return(annotate_df(x, protein_annotation))
  }
  
  if (is.list(x)) {
    return(lapply(x, annotate_recursive, protein_annotation = protein_annotation))
  }
  
  x
}

summarize_annotations <- function(x) {
  counts <- list(tables = 0L, rows = 0L, unmatched = 0L)
  
  walk <- function(y) {
    if (is.data.frame(y) && "ANALYTE_NAME" %in% names(y)) {
      counts$tables <<- counts$tables + 1L
      counts$rows <<- counts$rows + nrow(y)
      if ("ANNOT_UniProt" %in% names(y)) {
        counts$unmatched <<- counts$unmatched + sum(is.na(y$ANNOT_UniProt) | y$ANNOT_UniProt == "")
      }
    } else if (is.list(y)) {
      invisible(lapply(y, walk))
    }
    invisible(NULL)
  }
  
  walk(x)
  counts
}

annotate_rds <- function(input_file, output_file, protein_annotation) {
  obj <- readRDS(input_file)
  annotated <- annotate_recursive(obj, protein_annotation)
  saveRDS(annotated, output_file)
  
  stats <- summarize_annotations(annotated)
  cat(
    sprintf(
      "%s -> %s: annotated %d tables, %d rows, %d unmatched rows\n",
      input_file, output_file, stats$tables, stats$rows, stats$unmatched
    )
  )
  
  invisible(annotated)
}

protein_annotation <- read_annotation(translation_file)
annotate_rds(results_file, annotated_results_file, protein_annotation)
annotate_rds(reports_file, annotated_reports_file, protein_annotation)
