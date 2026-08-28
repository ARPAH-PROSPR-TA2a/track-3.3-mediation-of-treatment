pipeline_repo <- path.expand("~/FAST/GitHub/track-3.3")
translation_file <- path.expand(
  "~/FAST/Data/CALERIE/Raw/Proteomics/CALERIE_cleaned_protein_translation_table.csv"
)
results_file <- path.expand(
  "~/FAST/Outputs/3.3/Proteomics_3.3OWAS_INF/Proteomics_3.3OWAS_INF_results.rds"
)
annotated_tables_file <- path.expand(
  "~/FAST/Outputs/3.3/Proteomics_3.3OWAS_INF/annotated_outcome_effect_tables.rds"
)
fu_labels <- c("1" = "3mo", "2" = "6mo", "3" = "12mo", "4" = "24mo")

source(file.path(pipeline_repo, "outcomewas", "main.R"), chdir = TRUE)

stopifnot(
  file.exists(translation_file),
  file.exists(results_file),
  exists("write_annotated_outcome_effect_tables")
)

results <- readRDS(results_file)
protein_annotation <- read.csv(
  translation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

write_annotated_outcome_effect_tables(
  results = results,
  protein_annotation = protein_annotation,
  output_rds = annotated_tables_file,
  fu_labels = fu_labels
)

cat("Annotated outcome-effect tables saved to: ", annotated_tables_file, "\n", sep = "")
