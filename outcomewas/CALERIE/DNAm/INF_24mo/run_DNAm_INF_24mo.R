Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

library(dplyr)
library(forcats)
library(readr)

# -----------------------------
# Paths and parameters
# -----------------------------
omics_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/GRSet_fully_filtered_bmiq_chunk.rds")
pheno_raw_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_CPR_processed_pheno.rds")
control_pc_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/CALERIE_control_pcs_rgset_goodsamples.csv")
cell_pcs_path <- path.expand("~/FAST/Data/CALERIE/Raw/DNAm/Cell_PCs.csv")
outcome_path <- path.expand("~/FAST/Data/CALERIE/Raw/INF_CALERIE_24mo_Outcome.csv")

pipeline_repo <- path.expand("~/FAST/GitHub/track-3.3")
out_dir <- path.expand("~/FAST/Outputs/3.3/DNAm_betas_3.3OWAS_INF")
output_prefix <- "DNAm_betas_3.3OWAS_INF"

# DNAm probe lists are resolved relative to the repository root.
setwd(pipeline_repo)

omics_type <- "DNAm"
n_cores <- 47 # 48-vCPU production host; leave one core for the OS
checkpoint_batch_size <- 2000L

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
log_status("START DNAm beta run")
log_status(paste0(
  "parallel configuration: n_cores=", n_cores,
  "; BLAS/OpenMP threads per process=1"
))

# -----------------------------
# Read inputs
# -----------------------------
input_started <- proc.time()[["elapsed"]]
log_status("loading raw inputs (serial)")
omics_raw <- readRDS(omics_raw_path)
pheno_raw <- readRDS(pheno_raw_path)
control_pc <- read.csv(control_pc_path, stringsAsFactors = FALSE)
cell_PCs <- read_csv(cell_pcs_path, show_col_types = FALSE)
log_status(paste0(
  "raw inputs loaded: ", nrow(omics_raw), " CpGs, ",
  ncol(omics_raw), " methylation samples, ", nrow(pheno_raw),
  " phenotype rows (", round(proc.time()[["elapsed"]] - input_started, 1),
  " sec)"
))
preparation_started <- proc.time()[["elapsed"]]
log_status("input QC and table preparation starting (serial)")

# -----------------------------
# Input QC
# -----------------------------
stopifnot(!any(duplicated(rownames(omics_raw))))
stopifnot(!any(duplicated(colnames(omics_raw))))
stopifnot(!anyNA(omics_raw))
# Reject measurements outside the finite beta range.
beta_range <- range(as.matrix(omics_raw))
stopifnot(all(is.finite(beta_range)), beta_range[1] >= 0, beta_range[2] <= 1)
log_status(paste0("methylation scale: beta; range=", paste(beta_range, collapse = " to ")))

stopifnot(!any(duplicated(pheno_raw$Barcode)))
stopifnot(!any(duplicated(control_pc$filenames)))
stopifnot(!any(duplicated(cell_PCs$SAMPLE_ID)))

# -----------------------------
# Build phenotype table
# -----------------------------
pheno <- pheno_raw |>
  select(
    Participant_ID, Time_Point, Barcode, fu, CR, deidsite, agebl, female,
    bmistrat, snppc1.x, snppc2.x, snppc3.x
  ) |>
  left_join(
    control_pc |> select(filenames, paste0("PC", 1:20)),
    by = c("Barcode" = "filenames")
  )

pheno <- pheno |>
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
    fu = factor(haven::zap_labels(fu), levels = c(0, 1, 2)),
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
    TREATMENT_GROUP = CR
  ) |>
  select(
    -any_of(paste0("PC", 4:20))
  ) |>
  left_join(
    cell_PCs |> select(SAMPLE_ID, any_of(paste0("cell_PC", 1:4))),
    by = "SAMPLE_ID"
  )

outcome <- read_csv(outcome_path, col_names = TRUE, show_col_types = FALSE) |>
  select(SUBJECT_ID = DEID, OUTCOME = INF_Score)

stopifnot(!anyNA(outcome$SUBJECT_ID), !any(duplicated(outcome$SUBJECT_ID)))

# Preserve the 24-month-only analysis: FU=1 here is FU=2 in Track 1.1.1.
# Align visit keys explicitly before using these results for mediation.
pheno <- pheno |>
  filter(Time_Point %in% c("base", "24")) |>
  inner_join(outcome, by = "SUBJECT_ID") |>
  mutate(
    FU = case_when(
      Time_Point == "base" ~ 0,
      Time_Point == "24" ~ 1
    ),
    FU = factor(FU, levels = c(0, 1))
  )

# -----------------------------
# Define covariates and retain matched, complete samples (as in Track 1.1.1)
# -----------------------------
covariates <- c(
  c("snppc1.x", "snppc2.x", "snppc3.x"),
  c("agebl", "deidsite", "bmistrat"),
  paste0("cell_PC", 1:4),
  paste0("PC", 1:3)
)

required_pheno_cols <- c(
  "SAMPLE_ID", "SUBJECT_ID", "FU", "FEMALE", "TREATMENT_GROUP",
  "OUTCOME", covariates
)
stopifnot(all(required_pheno_cols %in% colnames(pheno)))

pheno <- pheno |>
  filter(SAMPLE_ID %in% colnames(omics_raw)) |>
  filter(if_all(all_of(required_pheno_cols), ~ !is.na(.x))) |>
  select(all_of(required_pheno_cols))

stopifnot(nrow(pheno) > 0)
stopifnot(!any(duplicated(pheno$SAMPLE_ID)))
stopifnot(all(pheno$SAMPLE_ID %in% colnames(omics_raw)))

# -----------------------------
# Build omics table
# -----------------------------
omics <- data.frame(
  ANALYTE_NAME = rownames(omics_raw),
  data.frame(omics_raw[, pheno$SAMPLE_ID, drop = FALSE], check.names = FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
stopifnot(identical(colnames(omics)[-1], pheno$SAMPLE_ID))

# -----------------------------
# Load pipeline functions
# -----------------------------
source(file.path(pipeline_repo, "outcomewas", "main.R"), chdir = TRUE)

stopifnot(
  exists("FAST_outcome_WAS"),
  exists("FAST_outcome_WAS_reports")
)

# -----------------------------
# Analysis
# -----------------------------
log_status(paste0(
  "inputs ready: ", nrow(omics), " CpGs, ", nrow(pheno),
  " matched phenotype rows (",
  round(proc.time()[["elapsed"]] - preparation_started, 1), " sec)"
))
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
  checkpoint_dir = file.path(out_dir, paste0(output_prefix, "_checkpoints")),
  checkpoint_batch_size = checkpoint_batch_size,
  verbose = TRUE
)

saveRDS(
  results,
  file = file.path(out_dir, paste0(output_prefix, "_results.rds"))
)

log_status("analysis complete; results saved")

# -----------------------------
# Reports
# -----------------------------
log_status("reports starting (serial)")

reports <- FAST_outcome_WAS_reports(
  pheno = pheno,
  omics = omics,
  omics_type = omics_type,
  additional_covariates = covariates,
  verbose = TRUE
)

saveRDS(
  reports,
  file = file.path(out_dir, paste0(output_prefix, "_reports.rds"))
)

log_status("reports complete; report output saved")

# -----------------------------
# Plotting
# -----------------------------
source(file.path(pipeline_repo, "outcomewas", "R", "plotting_helpers.R"), chdir = TRUE)

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
  "run complete (", round(proc.time()[["elapsed"]] - run_started, 1), " sec)"
))
