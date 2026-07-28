source(file.path("outcomewas", "main.R"))


run_checks <- function(checks) {
  all_pass <- TRUE
  for (name in names(checks)) {
    status <- if (isTRUE(checks[[name]])) "PASS" else "FAIL"
    cat(status, " ", name, "\n", sep = "")
    if (!isTRUE(checks[[name]])) {
      all_pass <- FALSE
      if (is.character(checks[[name]])) {
        cat("       ", checks[[name]], "\n", sep = "")
      }
    }
  }
  all_pass
}


normalize_frame <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }

  order_cols <- intersect(c("ANALYTE_NAME", "FU", "COEFFICIENT"), names(df))
  if (length(order_cols) > 0) {
    df <- df[do.call(order, df[order_cols]), , drop = FALSE]
  }

  row.names(df) <- NULL
  df
}


compare_results <- function(r1, r2, tol = 1e-10) {
  for (response in c("analysis_change", "analysis_level")) {
    for (stratum in c("all", "male", "female")) {
      obj1 <- r1[[response]][[stratum]]
      obj2 <- r2[[response]][[stratum]]

      if (is.null(obj1) || is.null(obj2)) {
        if (!identical(is.null(obj1), is.null(obj2))) {
          return(paste0(response, "$", stratum, ": NULL mismatch"))
        }
        next
      }

      for (slot in c("coefficients", "outcome_effects")) {
        d1 <- normalize_frame(obj1[[slot]])
        d2 <- normalize_frame(obj2[[slot]])
        eq <- all.equal(d1, d2, tolerance = tol, check.attributes = FALSE)
        if (!isTRUE(eq)) {
          return(paste0(response, "$", stratum, "$", slot, ": ", eq[1]))
        }
      }
    }
  }

  TRUE
}


count_tmp_files <- function(checkpoint_dir) {
  length(list.files(
    checkpoint_dir,
    pattern = "\\.tmp$",
    recursive = TRUE,
    full.names = TRUE
  ))
}


expected_batches <- function(n_analytes, batch_size) {
  ceiling(n_analytes / batch_size)
}


restrict_to_complete_subjects <- function(pheno_df) {
  fu_num <- as.integer(as.character(pheno_df$FU))
  subjects_with_baseline <- unique(pheno_df$SUBJECT_ID[fu_num == 0])
  subjects_with_followup <- unique(pheno_df$SUBJECT_ID[fu_num > 0])
  complete_subjects <- intersect(subjects_with_baseline, subjects_with_followup)
  pheno_df[pheno_df$SUBJECT_ID %in% complete_subjects, ]
}


create_subject_level_outcomes <- function(pheno_df, omics_df, preferred_fu = 1L) {
  fu_num <- as.integer(as.character(pheno_df$FU))
  pheno_baseline <- pheno_df[fu_num == 0, ]
  pheno_followup <- pheno_df[fu_num != 0, ]

  followup_priority <- ifelse(
    as.integer(as.character(pheno_followup$FU)) == preferred_fu,
    -1L,
    as.integer(as.character(pheno_followup$FU))
  )
  pheno_followup <- pheno_followup[order(pheno_followup$SUBJECT_ID, followup_priority), ]
  pheno_followup <- pheno_followup[!duplicated(pheno_followup$SUBJECT_ID), ]

  subjects <- intersect(pheno_baseline$SUBJECT_ID, pheno_followup$SUBJECT_ID)
  pheno_baseline <- pheno_baseline[match(subjects, pheno_baseline$SUBJECT_ID), ]
  pheno_followup <- pheno_followup[match(subjects, pheno_followup$SUBJECT_ID), ]

  signal_baseline <- as.numeric(omics_df[1, pheno_baseline$SAMPLE_ID])
  signal_followup <- as.numeric(omics_df[1, pheno_followup$SAMPLE_ID])
  signal_change <- signal_followup - signal_baseline

  treatment <- as.numeric(as.character(pheno_followup$TREATMENT_GROUP))
  bmi_scaled <- as.numeric(scale(pheno_followup$mbmi))
  signal_scaled <- as.numeric(scale(signal_change))

  set.seed(20260518)
  continuous_outcome <- 0.9 * signal_scaled + 0.4 * treatment + 0.25 * bmi_scaled +
    rnorm(length(signal_scaled), sd = 0.7)

  lp <- 0.7 * signal_scaled + 0.35 * treatment + 0.2 * bmi_scaled
  event_rate <- 0.08 * exp(lp)
  event_time <- rexp(length(lp), rate = pmax(event_rate, 1e-4))
  censor_time <- rexp(length(lp), rate = 0.06)
  outcome_status <- as.integer(event_time <= censor_time)
  outcome_time <- pmin(event_time, censor_time)

  data.frame(
    SUBJECT_ID = pheno_followup$SUBJECT_ID,
    OUTCOME = continuous_outcome,
    OUTCOME_TIME = outcome_time,
    OUTCOME_STATUS = outcome_status,
    stringsAsFactors = FALSE
  )
}


attach_subject_level_outcome <- function(pheno_df, subject_outcomes, type = c("continuous", "tte")) {
  type <- match.arg(type)
  out <- merge(pheno_df, subject_outcomes, by = "SUBJECT_ID", all.x = TRUE, sort = FALSE)
  out <- out[match(pheno_df$SAMPLE_ID, out$SAMPLE_ID), ]

  if (type == "continuous") {
    out$OUTCOME_TIME <- NULL
    out$OUTCOME_STATUS <- NULL
  } else {
    out$OUTCOME <- NULL
  }

  out
}


cat("OutcomeWAS Parallelization and Checkpointing Test Suite\n")
cat("=======================================================\n\n")

pheno_raw <- readRDS("outcomewas/Examples/ExampleData/pheno_example.rds")
omics_raw <- readRDS("outcomewas/Examples/ExampleData/proteomics_log2.rds")

if (any(sapply(pheno_raw, function(x) inherits(x, "haven_labelled")))) {
  for (col in names(pheno_raw)) {
    if (inherits(pheno_raw[[col]], "haven_labelled")) {
      raw_values <- as.vector(pheno_raw[[col]])
      numeric_attempt <- suppressWarnings(as.numeric(raw_values))
      if (!all(is.na(numeric_attempt))) {
        pheno_raw[[col]] <- numeric_attempt
      } else {
        pheno_raw[[col]] <- as.character(raw_values)
      }
    }
  }
}

pheno_raw$SUBJECT_ID <- as.character(pheno_raw$SUBJECT_ID)
pheno_raw$FU <- factor(as.integer(as.character(pheno_raw$FU)))
pheno_raw$FEMALE <- factor(as.integer(as.character(pheno_raw$FEMALE)))
pheno_raw$TREATMENT_GROUP <- factor(as.integer(as.character(pheno_raw$TREATMENT_GROUP)))
if ("CONTROL_STATUS" %in% names(pheno_raw)) {
  pheno_raw$CONTROL_STATUS <- NULL
}
pheno_raw[intersect(c("OUTCOME", "OUTCOME_TIME", "OUTCOME_STATUS"), names(pheno_raw))] <- NULL

required_base_cols <- c(
  "SAMPLE_ID", "FU", "SUBJECT_ID", "FEMALE", "TREATMENT_GROUP",
  "agebl", "agevis", "ethnic", "race3", "mbmi"
)
pheno <- pheno_raw[complete.cases(pheno_raw[, required_base_cols]), ]
pheno <- pheno[!duplicated(paste(pheno$SUBJECT_ID, pheno$FU)), ]

omics_full <- as.data.frame(omics_raw, check.names = FALSE)
if (!"ANALYTE_NAME" %in% names(omics_full)) {
  omics_full <- cbind(
    ANALYTE_NAME = rownames(omics_full),
    omics_full,
    stringsAsFactors = FALSE
  )
}
sample_names <- setdiff(colnames(omics_full), "ANALYTE_NAME")

shared_samples <- intersect(pheno$SAMPLE_ID, sample_names)
pheno <- pheno[pheno$SAMPLE_ID %in% shared_samples, ]
pheno <- restrict_to_complete_subjects(pheno)
omics_full <- omics_full[, c("ANALYTE_NAME", shared_samples), drop = FALSE]

omics_complete_idx <- rowSums(is.na(omics_full[, shared_samples, drop = FALSE])) == 0
omics_small <- omics_full[omics_complete_idx, , drop = FALSE]
omics_small <- omics_small[seq_len(min(30L, nrow(omics_small))), , drop = FALSE]

if (nrow(omics_small) < 30L) {
  stop("Need at least 30 complete analytes in the tracked example data.")
}

pheno_single_fu <- pheno[pheno$FU %in% c(0, 1), ]
pheno_single_fu <- restrict_to_complete_subjects(pheno_single_fu)
pheno_multi_fu <- pheno

subject_outcomes <- create_subject_level_outcomes(pheno_multi_fu, omics_small, preferred_fu = 1L)
pheno_single_continuous <- attach_subject_level_outcome(pheno_single_fu, subject_outcomes, type = "continuous")
pheno_multi_continuous <- attach_subject_level_outcome(pheno_multi_fu, subject_outcomes, type = "continuous")
pheno_multi_tte <- attach_subject_level_outcome(pheno_multi_fu, subject_outcomes, type = "tte")

additional_covariates <- c("agebl", "agevis", "ethnic", "race3", "mbmi")
fu_levels <- sort(unique(as.integer(as.character(pheno_multi_fu$FU))))
fu_levels <- fu_levels[fu_levels != 0]
parallel_workers <- 2L

cat("Data prepared\n")
cat("  Single-FU subjects: ", length(unique(pheno_single_fu$SUBJECT_ID)), "\n")
cat("  Multi-FU subjects:  ", length(unique(pheno_multi_fu$SUBJECT_ID)), "\n")
cat("  Multi-FU levels:    ", paste(fu_levels, collapse = ", "), "\n")
cat("  Small analytes:     ", nrow(omics_small), "\n")
cat("  TTE events:         ", sum(subject_outcomes$OUTCOME_STATUS), "\n\n")


# =============================================================================
# PARALLELIZATION TESTS
# =============================================================================

cat("PARALLELIZATION TESTS\n")
cat("=====================\n\n")

cat("P1: Serial vs parallel - single-FU continuous\n")

results_single_serial <- FAST_outcome_WAS(
  pheno = pheno_single_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1
)
results_single_parallel <- FAST_outcome_WAS(
  pheno = pheno_single_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers
)

p1_pass <- run_checks(list(
  "Results are numerically identical" = compare_results(results_single_serial, results_single_parallel)
))
cat("\n")

cat("P2: Serial vs parallel - multi-FU continuous\n")

results_multi_cont_serial <- FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1
)
results_multi_cont_parallel <- FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers
)

p2_pass <- run_checks(list(
  "Results are numerically identical" = compare_results(results_multi_cont_serial, results_multi_cont_parallel)
))
cat("\n")

cat("P3: Serial vs parallel - multi-FU TTE\n")

results_multi_tte_serial <- suppressWarnings(FAST_outcome_WAS(
  pheno = pheno_multi_tte,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1
))
results_multi_tte_parallel <- suppressWarnings(FAST_outcome_WAS(
  pheno = pheno_multi_tte,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers
))

p3_pass <- run_checks(list(
  "Results are numerically identical" = compare_results(results_multi_tte_serial, results_multi_tte_parallel)
))
cat("\n")

cat("P4: n_cores = NULL auto-detection\n")

p4_error <- tryCatch({
  FAST_outcome_WAS(
    pheno = pheno_single_continuous,
    omics = omics_small,
    omics_type = "Proteomics",
    additional_covariates = additional_covariates,
    n_cores = NULL
  )
  NULL
}, error = function(e) e$message)

p4_pass <- run_checks(list(
  "Runs without error" = is.null(p4_error)
))
cat("\n")

cat("P5: future::plan() restoration\n")

old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::sequential)

invisible(FAST_outcome_WAS(
  pheno = pheno_single_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers
))

p5_pass <- run_checks(list(
  "Plan restored to sequential" = inherits(future::plan(), "sequential")
))
cat("\n")

cat("P6: Large omics does not exceed future globals size limit\n")

n_large_probes <- 5000L
n_samples_single <- length(unique(pheno_single_continuous$SAMPLE_ID))

set.seed(999)
omics_large_mat <- matrix(
  runif(n_large_probes * n_samples_single),
  nrow = n_large_probes,
  ncol = n_samples_single,
  dimnames = list(NULL, unique(pheno_single_continuous$SAMPLE_ID))
)
omics_large <- cbind(
  ANALYTE_NAME = paste0("probe_", seq_len(n_large_probes)),
  as.data.frame(omics_large_mat),
  stringsAsFactors = FALSE
)

single_matrix_bytes <- n_large_probes * n_samples_single * 8
test_max_size <- floor(single_matrix_bytes * 0.6)
old_max_size <- getOption("future.globals.maxSize")
options(future.globals.maxSize = test_max_size)
on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)

p6_error <- tryCatch({
  FAST_outcome_WAS(
    pheno = pheno_single_continuous,
    omics = omics_large,
    omics_type = "Proteomics",
    n_cores = 2
  )
  NULL
}, error = function(e) e$message)

options(future.globals.maxSize = old_max_size)

p6_pass <- run_checks(list(
  "Large omics runs within tight globals size limit" = is.null(p6_error)
))
if (!is.null(p6_error)) {
  cat("       Error: ", p6_error, "\n", sep = "")
}
cat("\n")


# =============================================================================
# CHECKPOINTING TESTS
# =============================================================================

cat("CHECKPOINTING TESTS\n")
cat("===================\n\n")

BATCH_SIZE <- 10L
N_ANALYTES <- nrow(omics_small)
N_BATCHES <- expected_batches(N_ANALYTES, BATCH_SIZE)
STRATA <- c("all", "male", "female")
RESPONSES <- c("change", "level")

cat("C1/C2: Directory structure and batch file counts\n")

ckpt_dir_c1 <- file.path(tempdir(), "owas_ckpt_c1")
unlink(ckpt_dir_c1, recursive = TRUE, force = TRUE)
on.exit(unlink(ckpt_dir_c1, recursive = TRUE, force = TRUE), add = TRUE)

invisible(FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers,
  checkpoint_dir = ckpt_dir_c1,
  checkpoint_batch_size = BATCH_SIZE
))

structure_checks <- list()
for (resp in RESPONSES) {
  for (stratum in STRATA) {
    for (fu_level in fu_levels) {
      subdir <- file.path(ckpt_dir_c1, resp, stratum, paste0("FU", fu_level))
      n_files <- length(list.files(subdir, pattern = "^batch_[0-9]+\\.rds$"))
      structure_checks[[paste0(resp, "/", stratum, "/FU", fu_level, " directory exists")]] <-
        dir.exists(subdir)
      structure_checks[[paste0(resp, "/", stratum, "/FU", fu_level, " has ", N_BATCHES, " batch files")]] <-
        n_files == N_BATCHES
    }
  }
}

c1c2_pass <- run_checks(structure_checks)
cat("\n")

cat("C3: No .tmp files after clean run\n")

c3_pass <- run_checks(list(
  "No .tmp files remaining" = count_tmp_files(ckpt_dir_c1) == 0
))
cat("\n")

cat("C4: Checkpointed results match non-checkpointed\n")

results_from_cache <- FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = 1,
  checkpoint_dir = ckpt_dir_c1,
  checkpoint_batch_size = BATCH_SIZE
)

c4_pass <- run_checks(list(
  "Checkpointed results match reference" = compare_results(results_multi_cont_serial, results_from_cache)
))
cat("\n")

cat("C5/C6: Resume - partial crash simulation\n")

ckpt_dir_c5 <- file.path(tempdir(), "owas_ckpt_c5")
unlink(ckpt_dir_c5, recursive = TRUE, force = TRUE)
on.exit(unlink(ckpt_dir_c5, recursive = TRUE, force = TRUE), add = TRUE)

results_before_crash <- FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers,
  checkpoint_dir = ckpt_dir_c5,
  checkpoint_batch_size = BATCH_SIZE
)

crash_subdir <- file.path(ckpt_dir_c5, "change", "all", paste0("FU", max(fu_levels)))
all_batch_files <- sort(list.files(
  crash_subdir,
  pattern = "^batch_[0-9]+\\.rds$",
  full.names = TRUE
))
files_to_delete <- tail(all_batch_files, 2L)
kept_files <- setdiff(all_batch_files, files_to_delete)
mtimes_before <- file.mtime(kept_files)

invisible(file.remove(files_to_delete))
Sys.sleep(1)

results_after_resume <- FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers,
  checkpoint_dir = ckpt_dir_c5,
  checkpoint_batch_size = BATCH_SIZE
)

mtimes_after <- file.mtime(kept_files)

c5c6_pass <- run_checks(list(
  "Deleted batch files were recreated" = all(file.exists(files_to_delete)),
  "Resumed results match pre-crash results" = compare_results(results_before_crash, results_after_resume),
  "Untouched batches were not recomputed" = all(mtimes_before == mtimes_after)
))
cat("\n")

cat("C7: checkpoint_batch_size larger than n_analytes\n")

ckpt_dir_c7 <- file.path(tempdir(), "owas_ckpt_c7")
unlink(ckpt_dir_c7, recursive = TRUE, force = TRUE)
on.exit(unlink(ckpt_dir_c7, recursive = TRUE, force = TRUE), add = TRUE)

invisible(FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers,
  checkpoint_dir = ckpt_dir_c7,
  checkpoint_batch_size = N_ANALYTES + 100L
))

c7_checks <- list()
for (resp in RESPONSES) {
  for (stratum in STRATA) {
    for (fu_level in fu_levels) {
      subdir <- file.path(ckpt_dir_c7, resp, stratum, paste0("FU", fu_level))
      n_files <- length(list.files(subdir, pattern = "^batch_[0-9]+\\.rds$"))
      c7_checks[[paste0(resp, "/", stratum, "/FU", fu_level, " has exactly 1 batch file")]] <-
        n_files == 1L
    }
  }
}

c7_pass <- run_checks(c7_checks)
cat("\n")

cat("C8: checkpoint_batch_size = 1\n")

ckpt_dir_c8 <- file.path(tempdir(), "owas_ckpt_c8")
unlink(ckpt_dir_c8, recursive = TRUE, force = TRUE)
on.exit(unlink(ckpt_dir_c8, recursive = TRUE, force = TRUE), add = TRUE)

results_batch1 <- FAST_outcome_WAS(
  pheno = pheno_multi_continuous,
  omics = omics_small,
  omics_type = "Proteomics",
  additional_covariates = additional_covariates,
  n_cores = parallel_workers,
  checkpoint_dir = ckpt_dir_c8,
  checkpoint_batch_size = 1L
)

c8_checks <- list(
  "Results match reference" = compare_results(results_multi_cont_serial, results_batch1)
)
for (resp in RESPONSES) {
  for (stratum in STRATA) {
    for (fu_level in fu_levels) {
      subdir <- file.path(ckpt_dir_c8, resp, stratum, paste0("FU", fu_level))
      n_files <- length(list.files(subdir, pattern = "^batch_[0-9]+\\.rds$"))
      c8_checks[[paste0(resp, "/", stratum, "/FU", fu_level, " has ", N_ANALYTES, " batch files")]] <-
        n_files == N_ANALYTES
    }
  }
}

c8_pass <- run_checks(c8_checks)
cat("\n")


# =============================================================================
# SUMMARY
# =============================================================================

cat("FINAL SUMMARY\n")
cat("=============\n")
cat("P1 Serial vs parallel single-FU continuous: ", if (p1_pass) "PASS" else "FAIL", "\n", sep = "")
cat("P2 Serial vs parallel multi-FU continuous:  ", if (p2_pass) "PASS" else "FAIL", "\n", sep = "")
cat("P3 Serial vs parallel multi-FU TTE:         ", if (p3_pass) "PASS" else "FAIL", "\n", sep = "")
cat("P4 n_cores = NULL auto-detect:              ", if (p4_pass) "PASS" else "FAIL", "\n", sep = "")
cat("P5 future::plan() restoration:              ", if (p5_pass) "PASS" else "FAIL", "\n", sep = "")
cat("P6 Tight globals size regression:           ", if (p6_pass) "PASS" else "FAIL", "\n", sep = "")
cat("C1/C2 Directory structure and counts:       ", if (c1c2_pass) "PASS" else "FAIL", "\n", sep = "")
cat("C3 No .tmp files:                           ", if (c3_pass) "PASS" else "FAIL", "\n", sep = "")
cat("C4 Checkpointed results match reference:    ", if (c4_pass) "PASS" else "FAIL", "\n", sep = "")
cat("C5/C6 Resume behavior:                      ", if (c5c6_pass) "PASS" else "FAIL", "\n", sep = "")
cat("C7 Oversized batch size:                    ", if (c7_pass) "PASS" else "FAIL", "\n", sep = "")
cat("C8 checkpoint_batch_size = 1:               ", if (c8_pass) "PASS" else "FAIL", "\n", sep = "")

all_pass <- p1_pass && p2_pass && p3_pass && p4_pass && p5_pass && p6_pass &&
  c1c2_pass && c3_pass && c4_pass && c5c6_pass && c7_pass && c8_pass

cat("\n")
if (all_pass) {
  cat("ALL TESTS PASSED\n")
} else {
  cat("SOME TESTS FAILED\n")
  stop("One or more OutcomeWAS parallel/checkpoint tests failed.")
}
