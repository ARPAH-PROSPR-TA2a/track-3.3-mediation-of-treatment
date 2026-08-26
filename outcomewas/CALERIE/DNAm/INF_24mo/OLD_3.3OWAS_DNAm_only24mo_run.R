library(dplyr)
library(forcats)
library(haven)
library(readr)
library(tidyr)

# -----------------------------
# Explicit paths and parameters
# -----------------------------
betas_raw_path <- path.expand("~/CALERIE/DNAm/data/omics/GRSet_fully_filtered_bmiq_chunk_mvals.rds")
pheno_raw_path <- path.expand("~/CALERIE/DNAm/data/pheno/CALERIE_CPR_processed_pheno.rds")
control_pc_path <- path.expand("~/CALERIE/DNAm/data/covariate/CALERIE_control_pcs_rgset_goodsamples.csv")
cell_pcs_path <- path.expand("~/CALERIE/DNAm/data/covariate/Cell_PCs.csv")
outcome_path <- path.expand("~/CALERIE/DNAm/3.3/data/INF_CALERIE_24mo_Outcome.csv")

pipeline_repo <- path.expand("~/CALERIE/repos/3.3/track-3.3-mediation-of-treatment")
out_dir <- path.expand("~/CALERIE/DNAm/3.3/OutcomeWAS/INF_24mo")

omics_type <- "DNAm"
n_cores <- 3

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")

# -----------------------------
# Read raw inputs
# -----------------------------
betas_raw <- readRDS(betas_raw_path)
pheno_raw <- readRDS(pheno_raw_path)
control_pc <- read.csv(control_pc_path, stringsAsFactors = FALSE)
cell_PCs <- read_csv(cell_pcs_path, show_col_types = FALSE)

outcome <- read_csv(outcome_path, col_names = TRUE, show_col_types = FALSE) |>
  select(-1) |>
  rename(SUBJECT_ID = DEID, OUTCOME = INF_Score)

stopifnot(!any(duplicated(rownames(betas_raw))))
stopifnot(!any(duplicated(colnames(betas_raw))))
stopifnot(!any(is.na(betas_raw)))

stopifnot(!any(duplicated(pheno_raw$Barcode)))
stopifnot(!any(duplicated(control_pc$filenames)))
stopifnot(!any(duplicated(cell_PCs$SAMPLE_ID)))

# -----------------------------
# Create pheno
# -----------------------------
pheno <- pheno_raw |>
  select(
    Participant_ID, Time_Point, Barcode, fu, CR, deidsite, agebl, female,
    bmistrat, snppc1.x, snppc2.x, snppc3.x
  ) |>
  left_join(
    control_pc |> select(filenames, paste0("PC", 1:20)),
    by = c("Barcode" = "filenames")
  ) |>
  group_by(Participant_ID, Time_Point) |>
  slice(1) |>
  ungroup() |>
  mutate(
    agebl = as.numeric(scale(agebl)),
    snppc1.x = as.numeric(scale(snppc1.x)),
    snppc2.x = as.numeric(scale(snppc2.x)),
    snppc3.x = as.numeric(scale(snppc3.x))
  ) |>
  mutate(
    fu = factor(zap_labels(fu), levels = c(0, 1, 2)),
    female = as_factor(female),
    CR = as_factor(CR),
    deidsite = as_factor(deidsite),
    bmistrat = as_factor(bmistrat)
  ) |>
  rename(
    SAMPLE_ID = Barcode,
    SUBJECT_ID = Participant_ID,
    FU = fu,
    FEMALE = female,
    CONTROL_STATUS = CR
  ) |>
  select(-any_of(paste0("PC", 4:20))) |>
  left_join(
    cell_PCs |> select(SAMPLE_ID, any_of(paste0("cell_PC", 1:4))),
    by = "SAMPLE_ID"
  ) |>
  rename(TREATMENT_GROUP = CONTROL_STATUS) |>
  inner_join(outcome, by = "SUBJECT_ID") |>
  mutate(
    FU = case_when(
      Time_Point == "base" ~ 0,
      Time_Point == "12" ~ 1,
      Time_Point == "24" ~ 2
    ),
    FU = factor(FU, levels = c(0, 1, 2))
  ) |>
  drop_na()

# -----------------------------
# Create omics
# -----------------------------
shared_sample <- intersect(colnames(betas_raw), pheno$SAMPLE_ID)
stopifnot(length(shared_sample) > 0)

pheno <- pheno |>
  filter(SAMPLE_ID %in% shared_sample) |>
  arrange(match(SAMPLE_ID, colnames(betas_raw)))

omics <- betas_raw[, pheno$SAMPLE_ID, drop = FALSE] |>
  data.frame(check.names = FALSE)

omics <- data.frame(
  ANALYTE_NAME = rownames(omics),
  omics,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# -----------------------------
# Define covariates and validate input
# -----------------------------
covariates <- c(
  c("snppc1.x", "snppc2.x", "snppc3.x"),
  c("agebl", "deidsite", "bmistrat"),
  paste0("cell_PC", 1:4),
  paste0("PC", 1:3)
)

stopifnot(all(covariates %in% colnames(pheno)))
stopifnot(identical(colnames(omics)[-1], pheno$SAMPLE_ID))

# -----------------------------
# Load pipeline functions
# -----------------------------
setwd(pipeline_repo)
source(
  file.path("outcomewas", "main.R")
)

stopifnot(
  exists("FAST_outcome_WAS"),
  exists("FAST_outcome_WAS_reports")
)

# -----------------------------
# Analysis
# -----------------------------
cat("START DNAm run: ", as.character(Sys.time()), "\n",
    file = log_file, append = TRUE, sep = "")

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
source(file.path(pipeline_repo, "outcomewas", "R", "plotting_helpers.R"))

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