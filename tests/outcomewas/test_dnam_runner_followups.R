# Run the production preparation, models, reports, and plots on synthetic raw files.
# Only paths and worker/batch counts are replaced in temporary runner copies.
repo <- normalizePath(".", mustWork = TRUE)
fixture_dir <- tempfile("track33-dnam-followups-")
dir.create(fixture_dir)

run_checks <- function() {
  set.seed(20260909)
  n <- 120L
  person <- rep(seq_len(n), each = 3L)
  visit <- rep(c("base", "12", "24"), n)
  barcodes <- paste(sprintf("synthetic-%03d", person), visit, sep = "-")
  raw_pheno <- data.frame(
    Participant_ID = sprintf("synthetic-%03d", person), Time_Point = visit,
    Barcode = barcodes, fu = rep(0:2, n), CR = (person - 1L) %% 2L,
    female = ((person - 1L) %/% 2L) %% 2L,
    deidsite = rep(sample(1:3, n, replace = TRUE), each = 3L),
    agebl = rep(runif(n, 25, 50), each = 3L),
    bmistrat = rep(sample(0:1, n, replace = TRUE), each = 3L),
    snppc1.x = rep(rnorm(n), each = 3L), snppc2.x = rep(rnorm(n), each = 3L),
    snppc3.x = rep(rnorm(n), each = 3L)
  )
  controls <- data.frame(filenames = barcodes)
  for (i in 1:20) controls[[paste0("PC", i)]] <- rnorm(length(barcodes))
  cells <- data.frame(SAMPLE_ID = barcodes)
  for (i in 1:4) cells[[paste0("cell_PC", i)]] <- rnorm(length(barcodes))
  cells$cell_PC1[cells$SAMPLE_ID == "synthetic-004-base"] <- NA_real_
  probes <- head(readRDS("outcomewas/Data/FAST_epicv1_epicv2_sugden_TruD_probe_list.rds"), 4L)
  # Subject 5 has only FU1; subject 7 has only FU2. Both must contribute.
  beta_samples <- c(setdiff(barcodes, c("synthetic-005-24", "synthetic-007-12")),
                    "synthetic-omics-only")
  betas <- matrix(runif(4L * length(beta_samples), 0.05, 0.95), nrow = 4L,
                  dimnames = list(probes, beta_samples))
  outcomes <- data.frame(DEID = sprintf("synthetic-%03d", seq_len(n)),
                          INF_Score = rnorm(n), MetS_Score = rnorm(n))
  outcomes[6L, c("INF_Score", "MetS_Score")] <- NA_real_
  saveRDS(raw_pheno, file.path(fixture_dir, "pheno.rds"))
  saveRDS(betas, file.path(fixture_dir, "betas.rds"))
  write.csv(controls, file.path(fixture_dir, "controls.csv"), row.names = FALSE)
  write.csv(cells, file.path(fixture_dir, "cells.csv"), row.names = FALSE)
  write.csv(outcomes, file.path(fixture_dir, "outcomes.csv"), row.names = FALSE)
  saveRDS(raw_pheno[raw_pheno$Time_Point != "12", ], file.path(fixture_dir, "no12.rds"))

  run_runner <- function(runner, label, pheno_file = "pheno.rds") {
    overrides <- list(
      omics_raw_path = file.path(fixture_dir, "betas.rds"),
      pheno_raw_path = file.path(fixture_dir, pheno_file),
      control_pc_path = file.path(fixture_dir, "controls.csv"),
      cell_pcs_path = file.path(fixture_dir, "cells.csv"),
      outcome_path = file.path(fixture_dir, "outcomes.csv"), pipeline_repo = repo,
      out_dir = file.path(fixture_dir, label), n_cores = 2L, checkpoint_batch_size = 2L
    )
    code <- parse(file.path(repo, runner))
    replaced <- character()
    for (i in seq_along(code)) {
      expr <- code[[i]]
      if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
          is.symbol(expr[[2L]]) && as.character(expr[[2L]]) %in% names(overrides)) {
        name <- as.character(expr[[2L]])
        code[[i]][[3L]] <- overrides[[name]]
        replaced <- c(replaced, name)
      }
    }
    stopifnot(setequal(replaced, names(overrides)))
    snapshot_path <- file.path(fixture_dir, paste0(label, "-prepared.rds"))
    snapshot <- substitute(saveRDS(list(pheno = pheno, omics = omics,
      output_prefix = output_prefix), TARGET), list(TARGET = snapshot_path))
    script_path <- file.path(fixture_dir, paste0(label, "-runner.R"))
    writeLines(c(unlist(lapply(code, deparse)), deparse(snapshot)), script_path)
    log_path <- file.path(fixture_dir, paste0(label, "-console.log"))
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script_path),
                      stdout = log_path, stderr = log_path)
    list(status = status, log = readLines(log_path), out_dir = overrides$out_dir,
         prepared = if (file.exists(snapshot_path)) readRDS(snapshot_path) else NULL)
  }

  # Independently assemble covariates from raw files; do not use pipeline helpers.
  manual <- raw_pheno
  for (column in c("agebl", "snppc1.x", "snppc2.x", "snppc3.x")) {
    manual[[column]] <- as.numeric(scale(raw_pheno[[column]]))
  }
  for (column in c("deidsite", "bmistrat")) manual[[column]] <- factor(manual[[column]])
  manual$TREATMENT_GROUP <- factor(raw_pheno$CR)
  manual$FEMALE <- factor(raw_pheno$female)
  for (column in paste0("PC", 1:3)) {
    manual[[column]] <- controls[[column]][match(manual$Barcode, controls$filenames)]
  }
  for (column in paste0("cell_PC", 1:4)) {
    manual[[column]] <- cells[[column]][match(manual$Barcode, cells$SAMPLE_ID)]
  }
  covariates <- c("snppc1.x", "snppc2.x", "snppc3.x", "agebl", "deidsite", "bmistrat",
                  paste0("cell_PC", 1:4), paste0("PC", 1:3))
  expected_samples <- setdiff(barcodes, c("synthetic-004-base", "synthetic-005-24",
    "synthetic-007-12", paste0("synthetic-006-", c("base", "12", "24"))))
  runners <- c(INF = "outcomewas/CALERIE/DNAm/INF_24mo/run_DNAm_INF_24mo.R",
               MetS = "outcomewas/CALERIE/DNAm/MetS_24mo/run_DNAm_Mets_24mo.R")

  for (label in names(runners)) {
    run <- run_runner(runners[[label]], label)
    if (run$status != 0L) stop(paste(run$log, collapse = "\n"))
    prepared <- run$prepared
    pheno <- prepared$pheno
    raw_rows <- match(pheno$SAMPLE_ID, raw_pheno$Barcode)
    score <- outcomes[[paste0(label, "_Score")]]
    stopifnot(setequal(pheno$SAMPLE_ID, expected_samples), all(complete.cases(pheno)),
      identical(colnames(prepared$omics)[-1L], pheno$SAMPLE_ID),
      identical(sort(unique(as.integer(as.character(pheno$FU)))), 0:2),
      all(as.integer(as.character(pheno$FU)) == raw_pheno$fu[raw_rows]),
      isTRUE(all.equal(pheno$OUTCOME, score[match(pheno$SUBJECT_ID, outcomes$DEID)],
                       tolerance = 1e-12)),
      all(vapply(split(pheno$OUTCOME, pheno$SUBJECT_ID), function(x) length(unique(x)) == 1L,
                 logical(1))))
    prefix <- paste0("DNAm_betas_3.3OWAS_", label)
    stopifnot(identical(prepared$output_prefix, prefix))
    results <- readRDS(file.path(run$out_dir, paste0(prefix, "_results.rds")))
    reports <- readRDS(file.path(run$out_dir, paste0(prefix, "_reports.rds")))
    stopifnot(dir.exists(file.path(run$out_dir, paste0(prefix, "_checkpoints"))),
      file.exists(file.path(run$out_dir, "run.log")),
      length(list.files(file.path(run$out_dir, "Figures"), pattern = "\\.pdf$", recursive = TRUE)) == 8L,
      identical(sort(unique(reports$outcome_reports$analysis_sample_summary$FU)), 1:2),
      all(paste0("omics_FU", 0:2, "_Tx0") %in% names(reports$variable_summaries$all)))

    model_data <- manual
    model_data$OUTCOME <- score[match(model_data$Participant_ID, outcomes$DEID)]
    model_data <- model_data[model_data$Barcode %in% colnames(betas) &
      complete.cases(model_data[, c("OUTCOME", "FEMALE", "TREATMENT_GROUP", covariates)]), ]
    baseline <- model_data[model_data$fu == 0L, ]
    for (fu in 1:2) {
      followup <- model_data[model_data$fu == fu &
        model_data$Participant_ID %in% baseline$Participant_ID, ]
      stopifnot(nrow(followup) == 117L,
        ("synthetic-005" %in% followup$Participant_ID) == (fu == 1L),
        ("synthetic-007" %in% followup$Participant_ID) == (fu == 2L))
      baseline_samples <- baseline$Barcode[match(followup$Participant_ID, baseline$Participant_ID)]
      followup$analyte_baseline <- as.numeric(betas[1L, baseline_samples])
      for (analysis in c("analysis_change", "analysis_level")) {
        followup$analyte <- as.numeric(betas[1L, followup$Barcode])
        if (analysis == "analysis_change") followup$analyte <- followup$analyte - followup$analyte_baseline
        for (stratum in c("all", "male", "female")) {
          md <- followup
          if (stratum != "all") {
            sex <- if (stratum == "female") "1" else "0"
            md <- md[as.character(md$FEMALE) == sex, ]
          }
          adjustment <- c("analyte_baseline", "TREATMENT_GROUP",
                           if (stratum == "all") "FEMALE", covariates)
          fit <- lm(reformulate(c("analyte", adjustment), response = "OUTCOME"), data = md)
          expected <- coef(summary(fit))["analyte", c("Estimate", "Std. Error")]
          result <- results[[analysis]][[stratum]]$outcome_effects
          coefficients <- results[[analysis]][[stratum]]$coefficients
          stopifnot(nrow(result) == 8L, identical(sort(unique(result$FU)), 1:2),
                    identical(sort(unique(coefficients$FU)), 1:2),
                    "N_OBS" %in% names(coefficients), all(is.finite(result$EFFECT_SIZE)),
                    all(coefficients$N_OBS[coefficients$FU == fu] == nrow(md)))
          actual <- result[result$FU == fu & result$ANALYTE_NAME == probes[1L], c("EFFECT_SIZE", "SE")]
          stopifnot(nrow(actual) == 1L,
            isTRUE(all.equal(as.numeric(actual[1L, ]), unname(expected), tolerance = 1e-10)))
        }
      }
    }
    cat("PASS ", label, ": raw FU0/1/2, fixed 24-month outcome, independent visit pairing, manual beta fits, reports and plots\n", sep = "")

    no12 <- run_runner(runners[[label]], paste0(label, "-no12"), pheno_file = "no12.rds")
    stopifnot(no12$status != 0L, is.null(no12$prepared),
      any(grepl("Expected baseline-paired methylation at both 12 months (FU=1) and 24 months (FU=2).",
                no12$log, fixed = TRUE)),
      !any(grepl("analysis call starting", no12$log, fixed = TRUE)),
      length(list.files(no12$out_dir, pattern = "(_results\\.rds$|batch.*\\.rds$)", recursive = TRUE)) == 0L)
    cat("PASS ", label, ": absent 12-month data rejected before model fitting\n", sep = "")
  }
}

tryCatch(run_checks(), finally = unlink(fixture_dir, recursive = TRUE))
