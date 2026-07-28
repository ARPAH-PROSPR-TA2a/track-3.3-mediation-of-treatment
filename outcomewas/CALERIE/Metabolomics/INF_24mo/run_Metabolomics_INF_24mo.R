library(dplyr)
library(forcats)
library(readr)
library(readxl)
library(tibble)
library(tidyr)

# -----------------------------
# Explicit paths and parameters
# -----------------------------
omics_raw_path <- path.expand("~/CALERIE/Metabolomics/data/OBT_log.csv")
pheno_raw_path <- path.expand("~/CALERIE/data/Proteomics/CALERIE_clinical_w_omics_crosswalks.csv")
genetic_pcs_path <- path.expand("~/CALERIE/data/Proteomics/CALERIE_GeneticPCs_20200727.xlsx")
outcome_path <- path.expand("~/CALERIE/Metabolomics/3.3/data/INF_CALERIE_24mo_Outcome.csv")

pipeline_repo <- path.expand("~/CALERIE/repos/3.3/track-3.3-mediation-of-treatment")
out_dir <- path.expand("~/CALERIE/Metabolomics/3.3/INF_24mo")

omics_type <- "Metabolomics"
n_cores <- 3

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")

cat("START Metabolomics run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Read raw inputs
# -----------------------------
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
  select(-1) |>
  rename(SUBJECT_ID = DEID, OUTCOME = INF_Score)

metabolite_cols <- names(omics_raw)[25:ncol(omics_raw)]

stopifnot(!any(is.na(omics_raw[, metabolite_cols])))
stopifnot(!any(duplicated(colnames(omics_raw))))
stopifnot(!any(duplicated(omics_raw$CLIENT_IDENTIFIER)))

# -----------------------------
# Create pheno
# -----------------------------
pheno <- pheno_raw |>
  select(DEID, CR, deidsite, agebl, female, bmistrat, fumo, metabolomics_barcode) |>
  left_join(
    PCs |>
      select(ID, PC1, PC2, PC3),
    by = c("DEID" = "ID")
  ) |>
  drop_na() |>
  group_by(DEID, fumo) |>
  slice(1) |>
  ungroup() |>
  rename(
    SAMPLE_ID = metabolomics_barcode,
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
  filter(FU %in% c(0, 24)) |>
  inner_join(outcome, by = "SUBJECT_ID") |>
  mutate(
    FU = case_when(
      FU == 0 ~ 0,
      FU == 24 ~ 1
    ),
    FU = factor(FU, levels = c(0, 1))
  )

stopifnot(!any(is.na(pheno)))
stopifnot(!any(duplicated(pheno$SAMPLE_ID)))

# -----------------------------
# Create omics
# -----------------------------
shared_sample <- intersect(omics_raw$CLIENT_IDENTIFIER, pheno$SAMPLE_ID)
stopifnot(length(shared_sample) > 0)

pheno <- pheno |>
  filter(SAMPLE_ID %in% shared_sample) |>
  arrange(match(SAMPLE_ID, omics_raw$CLIENT_IDENTIFIER))

omics <- omics_raw |>
  filter(CLIENT_IDENTIFIER %in% pheno$SAMPLE_ID) |>
  arrange(match(CLIENT_IDENTIFIER, pheno$SAMPLE_ID)) |>
  select(CLIENT_IDENTIFIER, all_of(metabolite_cols)) |>
  column_to_rownames("CLIENT_IDENTIFIER") |>
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

stopifnot(all(covariates %in% colnames(pheno)))
stopifnot(identical(colnames(omics)[-1], pheno$SAMPLE_ID))

# -----------------------------
# Load pipeline functions
# -----------------------------
source(file.path(pipeline_repo, "main.R"), chdir = TRUE)

stopifnot(
  exists("FAST_outcome_WAS"),
  exists("FAST_outcome_WAS_reports"),
  exists("pheno"),
  exists("omics"),
  exists("covariates")
)

# -----------------------------
# Analysis
# -----------------------------
results <- FAST_outcome_WAS(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates,
  n_cores = n_cores,
  checkpoint_dir = file.path(out_dir, "checkpoints")
)

saveRDS(
  results,
  file = file.path(out_dir, "results.rds")
)

cat("DONE analysis: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Reports
# -----------------------------
cat("START report run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

reports <- FAST_outcome_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates
)

saveRDS(
  reports,
  file = file.path(out_dir, "reports.rds")
)

cat("DONE reports: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

# -----------------------------
# Plotting
# -----------------------------
source(file.path(pipeline_repo, "plotting_helpers.R"), chdir = TRUE)

fig_change <- file.path(out_dir, "Figures", "change")
fig_level <- file.path(out_dir, "Figures", "level")

dir.create(fig_change, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_level, recursive = TRUE, showWarnings = FALSE)

generate_all_plots(
  results,
  figures_dir = fig_change,
  analysis = "analysis_change"
)

generate_all_plots(
  results,
  figures_dir = fig_level,
  analysis = "analysis_level"
)

