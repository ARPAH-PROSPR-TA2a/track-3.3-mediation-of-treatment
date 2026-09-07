Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

library(dplyr)
library(forcats)
library(readr)
library(readxl)
library(tibble)
library(tidyr)
print(getwd())
# -----------------------------
# Explicit paths and parameters
# -----------------------------
omics_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/Proteomics/CALERIE_cleaned_log2_soma_matrix.csv")
pheno_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/CALERIE_clinical_w_omics_crosswalks.csv")
genetic_pcs_path <- path.expand("~/FAST/Data/CALERIE/Raw/CALERIE_GeneticPCs_20200727.xlsx")
outcome_path <- path.expand("~/FAST/Data/CALERIE/Raw/MetS_CALERIE_24mo_Outcome.csv")
translation_path <- path.expand("~/FAST/Data/CALERIE/Raw/Proteomics/CALERIE_cleaned_protein_translation_table.csv")

pipeline_repo <- path.expand("~/FAST/GitHub/track-3.3")
out_dir <- path.expand("~/FAST/Outputs/3.3/Proteomics_3.3OWAS_MetS")
annotated_results_path <- file.path(out_dir, "Proteomics_3.3OWAS_MetS_results_annotated.rds")

omics_type <- "Proteomics"
n_cores <- 47 # 48-vCPU production host; leave one core for the OS
checkpoint_batch_size <- 2000L
fu_labels <- c("1" = "3mo", "2" = "6mo", "3" = "12mo", "4" = "24mo")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")
options(track33.progress_log = log_file)

log_status <- function(text) {
  line <- paste0(
    "[3.3] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ", text
  )
  message(line)
  cat(line, "\n", file = log_file, append = TRUE, sep = "")
}

run_started <- proc.time()[["elapsed"]]
log_status("START Proteomics run")
log_status(paste0(
  "parallel configuration: n_cores=", n_cores,
  "; BLAS/OpenMP threads per process=1"
))


# -----------------------------
# Read raw inputs
# -----------------------------
input_started <- proc.time()[["elapsed"]]
log_status("loading raw inputs (serial)")
omics_raw <- read.csv(
  omics_raw_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

pheno_raw <- read.csv(
  pheno_raw_path,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

PCs <- read_excel(genetic_pcs_path)

outcome <- read_csv(outcome_path, col_names = TRUE, show_col_types = FALSE) |>
  dplyr::select(-1) |>
  dplyr::rename(SUBJECT_ID = DEID, OUTCOME = 'MetS_Score')

log_status(paste0(
  "raw inputs loaded: ", nrow(omics_raw), " proteomics samples, ",
  nrow(pheno_raw), " phenotype rows, ", nrow(outcome), " outcome rows (",
  round(proc.time()[["elapsed"]] - input_started, 1), " sec)"
))
preparation_started <- proc.time()[["elapsed"]]
log_status("input QC and table preparation starting (serial)")

stopifnot(!any(is.na(omics_raw)))
stopifnot(!any(duplicated(colnames(omics_raw))))
stopifnot(!any(duplicated(omics_raw$SampleId)))

# -----------------------------
# Create pheno
# -----------------------------
pheno <- pheno_raw |>
  dplyr::select(DEID, CR, deidsite, agebl, female, bmistrat, fumo, proteomics_barcode) |>
  left_join(
    PCs |>
      dplyr::select(ID, PC1, PC2, PC3),
    by = c("DEID" = "ID")
  ) |>
  drop_na() |>
  group_by(DEID, fumo) |>
  dplyr::slice(1) |>
  ungroup() |>
  dplyr::rename(
    SAMPLE_ID = proteomics_barcode,
    SUBJECT_ID = DEID,
    FU = fumo,
    FEMALE = female,
    TREATMENT_GROUP = CR
  ) |>
  mutate(
    agebl = as.numeric(scale(agebl)),
    PC1 = as.numeric(scale(PC1)),
    PC2 = as.numeric(scale(PC2)),
    PC3 = as.numeric(scale(PC3))
  ) |>
  mutate(
    FEMALE = as_factor(FEMALE),
    TREATMENT_GROUP = as_factor(TREATMENT_GROUP),
    deidsite = as_factor(deidsite),
    bmistrat = as_factor(bmistrat)
  ) |>
  inner_join(outcome, by = "SUBJECT_ID") |>
  mutate(
    FU = case_when(
      FU == 0 ~ 0,
      FU == 3 ~ 1,
      FU == 6 ~ 2,
      FU == 12 ~ 3,
      FU == 24 ~ 4
    ),
    FU = factor(FU, levels = c(0, 1, 2, 3, 4))
  )

# -----------------------------
# Create omics
# -----------------------------
shared_sample <- intersect(omics_raw$SampleId, pheno$SAMPLE_ID)
stopifnot(length(shared_sample) > 0)

pheno <- pheno |>
  filter(SAMPLE_ID %in% shared_sample) |>
  arrange(match(SAMPLE_ID, omics_raw$SampleId))

omics <- omics_raw |>
  filter(SampleId %in% pheno$SAMPLE_ID) |>
  arrange(match(SampleId, pheno$SAMPLE_ID)) |>
  column_to_rownames("SampleId") |>
  t() |>
  as.data.frame(check.names = FALSE) |>
  rownames_to_column("ANALYTE_NAME")

# -----------------------------
# Define covariates and validate input
# -----------------------------
covariates <- c(
  c("PC1", "PC2", "PC3"),
  c("agebl", "deidsite", "bmistrat")
)

stopifnot(!any(is.na(pheno)))
stopifnot(all(covariates %in% colnames(pheno)))
stopifnot(identical(colnames(omics)[-1], pheno$SAMPLE_ID))

# -----------------------------
# Load pipeline functions
# -----------------------------



# Forcing the user to set WD is messing up the entire pipeline.


source(file.path(pipeline_repo, "outcomewas/main.R"), chdir = TRUE)

stopifnot(
  exists("FAST_outcome_WAS"),
  exists("FAST_outcome_WAS_reports"),
  exists("get_annotated_results"),
  exists("pheno"),
  exists("omics"),
  exists("covariates")
)

log_status(paste0(
  "inputs ready: ", nrow(omics), " analytes, ", nrow(pheno),
  " matched phenotype rows (",
  round(proc.time()[["elapsed"]] - preparation_started, 1), " sec)"
))

# -----------------------------
# Analysis
# -----------------------------
log_status(paste0(
  "analysis call starting: n_cores=", n_cores,
  ", checkpoint batch size=", checkpoint_batch_size,
  "; parallel work occurs inside uncached batches"
))

results <- FAST_outcome_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates,
  n_cores = n_cores,
  checkpoint_dir = file.path(out_dir, "Proteomics_3.3OWAS_MetS_checkpoints"),
  checkpoint_batch_size = checkpoint_batch_size,
  verbose = TRUE
)

saveRDS(
  results,
  file = file.path(out_dir, "Proteomics_3.3OWAS_MetS_results.rds")
)

log_status("analysis model fitting complete; raw results saved")

log_status("START protein annotation")

protein_annotation <- read.csv(
  translation_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

annotated_results <- get_annotated_results(
  results = results,
  protein_annotation = protein_annotation,
  fu_labels = fu_labels
)
saveRDS(annotated_results, annotated_results_path)

log_status(paste0(
  "DONE protein annotation; annotated results: ", annotated_results_path
))

log_status("DONE analysis")

# -----------------------------
# Reports
# -----------------------------
log_status("START report run")

reports <- FAST_outcome_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates,
  verbose = TRUE
)

saveRDS(
  reports,
  file = file.path(out_dir, "Proteomics_3.3OWAS_MetS_reports.rds")
)

log_status("DONE reports; report output saved")

# -----------------------------
# Plotting
# -----------------------------
source(file.path(pipeline_repo, "outcomewas", "R", "plotting_helpers.R"))

fig_change <- file.path(out_dir, "Figures", "change")
fig_level <- file.path(out_dir, "Figures", "level")

dir.create(fig_change, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_level, recursive = TRUE, showWarnings = FALSE)

log_status("change-result plotting starting (serial)")
generate_all_plots(
  results,
  figures_dir = fig_change,
  analysis = "analysis_change"
)
log_status("change-result plotting complete")

log_status("level-result plotting starting (serial)")
generate_all_plots(
  results,
  figures_dir = fig_level,
  analysis = "analysis_level"
)
log_status("level-result plotting complete")
log_status(paste0(
  "run complete (",
  round(proc.time()[["elapsed"]] - run_started, 1),
  " sec)"
))
