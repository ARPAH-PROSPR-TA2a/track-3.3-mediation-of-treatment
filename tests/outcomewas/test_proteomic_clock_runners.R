# Exercise production preparation, fits, reports, and plots on shuffled long-format
# clock inputs. Temporary runner copies replace only paths, workers, and the XLSX reader.
repo <- normalizePath(".", mustWork = TRUE)
fixture_dir <- tempfile("track33-proteomic-clocks-")
dir.create(fixture_dir)

run_checks <- function() {
  set.seed(20261022)
  n <- 96L
  months <- c(0L, 3L, 6L, 12L, 24L)
  person <- rep(seq_len(n), each = length(months))
  visit <- rep(months, n)
  barcodes <- paste(sprintf("synthetic-%03d", person), visit, sep = "-")
  raw_pheno <- data.frame(
    DEID = sprintf("synthetic-%03d", person), fumo = visit,
    proteomics_barcode = barcodes, CR = (person - 1L) %% 2L,
    female = ((person - 1L) %/% 2L) %% 2L,
    deidsite = rep(sample(1:3, n, replace = TRUE), each = length(months)),
    agebl = rep(runif(n, 25, 50), each = length(months)),
    bmistrat = rep(sample(0:1, n, replace = TRUE), each = length(months))
  )
  pcs <- data.frame(ID = sprintf("synthetic-%03d", seq_len(n)),
                    PC1 = rnorm(n), PC2 = rnorm(n), PC3 = rnorm(n))
  outcomes <- data.frame(index = seq_len(n), DEID = pcs$ID,
                         INF_Score = rnorm(n), MetS_Score = rnorm(n))
  raw_clocks <- data.frame(
    Organ = rep(c("Brain", "WholeBody"), each = length(barcodes)),
    proteomics_barcode = rep(barcodes, 2L),
    Adjusted_Predicted_Age = rnorm(2L * length(barcodes))
  )
  # Each FU must use its own available baseline pairs, not complete five-visit cases.
  missing_samples <- c("synthetic-005-3", "synthetic-007-12", "synthetic-009-0")
  raw_clocks <- raw_clocks[!raw_clocks$proteomics_barcode %in% missing_samples, ]
  raw_clocks <- raw_clocks[sample(nrow(raw_clocks)), ]
  write.csv(raw_clocks, file.path(fixture_dir, "clocks.csv"), row.names = FALSE)
  write.csv(rbind(raw_clocks, raw_clocks[1L, ]),
            file.path(fixture_dir, "duplicate-clocks.csv"), row.names = FALSE)
  write.csv(raw_pheno[sample(nrow(raw_pheno)), ],
            file.path(fixture_dir, "pheno.csv"), row.names = FALSE)
  write.csv(pcs[sample(nrow(pcs)), ], file.path(fixture_dir, "pcs.csv"), row.names = FALSE)
  write.csv(outcomes, file.path(fixture_dir, "outcomes.csv"), row.names = FALSE)

  run_runner <- function(runner, label, duplicate = FALSE) {
    code <- parse(file.path(repo, runner))
    overrides <- list(
      omics_raw_path = file.path(fixture_dir, if (duplicate) "duplicate-clocks.csv" else "clocks.csv"),
      pheno_raw_path = file.path(fixture_dir, "pheno.csv"),
      genetic_pcs_path = file.path(fixture_dir, "pcs.csv"),
      outcome_path = file.path(fixture_dir, "outcomes.csv"), pipeline_repo = repo,
      n_cores = 1L
    )
    # Preserve the production leaf so a legacy checkpoint is present alongside it.
    output_assignment <- Filter(function(x) is.call(x) &&
      identical(x[[1L]], as.name("<-")) && identical(x[[2L]], as.name("out_dir")), code)
    output_leaf <- basename(eval(output_assignment[[1L]][[3L]]))
    stopifnot(identical(output_leaf, paste0(label, "_24mo_allFU_v1")))
    overrides$out_dir <- file.path(fixture_dir, if (duplicate) "duplicates" else "outputs", output_leaf)
    replaced <- character()
    reader_replaced <- FALSE
    for (i in seq_along(code)) {
      expr <- code[[i]]
      if (!is.call(expr) || !identical(expr[[1L]], as.name("<-")) || !is.symbol(expr[[2L]])) next
      name <- as.character(expr[[2L]])
      if (name %in% names(overrides)) {
        code[[i]][[3L]] <- overrides[[name]]
        replaced <- c(replaced, name)
      } else if (identical(name, "PCs")) {
        stopifnot(identical(expr[[3L]][[1L]], as.name("read_excel")))
        code[[i]][[3L]] <- quote(read.csv(genetic_pcs_path, stringsAsFactors = FALSE))
        reader_replaced <- TRUE
      }
    }
    stopifnot(setequal(replaced, names(overrides)), reader_replaced)
    name <- paste0(label, if (duplicate) "-duplicates" else "")
    snapshot_path <- file.path(fixture_dir, paste0(name, "-prepared.rds"))
    snapshot <- substitute(saveRDS(list(pheno = pheno, omics = omics,
      omics_type = omics_type), TARGET), list(TARGET = snapshot_path))
    script_path <- file.path(fixture_dir, paste0(name, "-runner.R"))
    writeLines(c(unlist(lapply(code, deparse)), deparse(snapshot)), script_path)
    log_path <- file.path(fixture_dir, paste0(name, "-console.log"))
    status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script_path),
                      stdout = log_path, stderr = log_path)
    list(status = status, log = readLines(log_path), out_dir = overrides$out_dir,
         prepared = if (file.exists(snapshot_path)) readRDS(snapshot_path) else NULL)
  }

  # Independently reconstruct the adjusted model from the raw files.
  manual <- raw_pheno
  manual$agebl <- as.numeric(scale(manual$agebl))
  for (column in c("PC1", "PC2", "PC3")) {
    manual[[column]] <- as.numeric(scale(pcs[[column]][match(manual$DEID, pcs$ID)]))
  }
  for (column in c("CR", "female", "deidsite", "bmistrat")) {
    manual[[column]] <- factor(manual[[column]])
  }
  manual$FU <- match(manual$fumo, months) - 1L
  manual <- manual[!manual$proteomics_barcode %in% missing_samples, ]
  brain <- raw_clocks[raw_clocks$Organ == "Brain", ]
  manual$clock <- brain$Adjusted_Predicted_Age[match(manual$proteomics_barcode, brain$proteomics_barcode)]
  runners <- c(
    INF = "outcomewas/CALERIE/Proteomics/Proteomic_Clocks/INF_24mo/run_Proteomic_Clocks_INF_24mo.R",
    MetS = "outcomewas/CALERIE/Proteomics/Proteomic_Clocks/MetS_24mo/run_Proteomic_Clock_MetS_24mo.R"
  )

  for (label in names(runners)) {
    legacy_file <- file.path(fixture_dir, "outputs", paste0(label, "_24mo"),
                             "checkpoints", "batch_1.rds")
    dir.create(dirname(legacy_file), recursive = TRUE)
    saveRDS("legacy checkpoint must remain untouched", legacy_file)
    legacy_hash <- tools::md5sum(legacy_file)
    run <- run_runner(runners[[label]], label)
    if (run$status != 0L) stop(paste(run$log, collapse = "\n"))
    pheno <- run$prepared$pheno
    omics <- run$prepared$omics
    score <- outcomes[[paste0(label, "_Score")]]
    stopifnot(identical(run$prepared$omics_type, "other"),
      setequal(pheno$SAMPLE_ID, setdiff(barcodes, missing_samples)),
      !anyDuplicated(pheno$SAMPLE_ID), identical(colnames(omics)[-1L], pheno$SAMPLE_ID),
      identical(sort(unique(as.integer(as.character(pheno$FU)))), 0:4),
      all(as.integer(as.character(pheno$FU)) ==
        match(raw_pheno$fumo[match(pheno$SAMPLE_ID, raw_pheno$proteomics_barcode)], months) - 1L),
      isTRUE(all.equal(pheno$OUTCOME, score[match(pheno$SUBJECT_ID, outcomes$DEID)], tolerance = 1e-12)),
      all(vapply(split(pheno$OUTCOME, pheno$SUBJECT_ID), function(x) length(unique(x)) == 1L, logical(1))),
      all(vapply(omics[-1L], is.numeric, logical(1))),
      identical(tools::md5sum(legacy_file), legacy_hash))
    for (organ in unique(raw_clocks$Organ)) {
      rows <- raw_clocks[raw_clocks$Organ == organ, ]
      expected <- rows$Adjusted_Predicted_Age[match(pheno$SAMPLE_ID, rows$proteomics_barcode)]
      stopifnot(isTRUE(all.equal(as.numeric(omics[omics$ANALYTE_NAME == organ, -1L]),
                                expected, tolerance = 1e-12)))
    }
    results <- readRDS(file.path(run$out_dir, "results.rds"))
    reports <- readRDS(file.path(run$out_dir, "reports.rds"))
    sample_summary <- reports$outcome_reports$analysis_sample_summary
    stopifnot(dir.exists(file.path(run$out_dir, "checkpoints")),
      file.exists(file.path(run$out_dir, "run.log")),
      length(list.files(file.path(run$out_dir, "Figures"), pattern = "\\.pdf$", recursive = TRUE)) == 4L,
      identical(sort(sample_summary$FU), 1:4),
      all(paste0("omics_FU", 0:4, "_Tx0") %in% names(reports$variable_summaries$all)))

    model_data <- manual
    model_data$OUTCOME <- score[match(model_data$DEID, outcomes$DEID)]
    baseline <- model_data[model_data$FU == 0L, ]
    for (fu in 1:4) {
      followup <- model_data[model_data$FU == fu & model_data$DEID %in% baseline$DEID, ]
      followup$analyte_baseline <- baseline$clock[match(followup$DEID, baseline$DEID)]
      stopifnot(nrow(followup) == if (fu %in% c(1L, 3L)) 94L else 95L,
        !("synthetic-009" %in% followup$DEID),
        ("synthetic-005" %in% followup$DEID) == (fu != 1L),
        ("synthetic-007" %in% followup$DEID) == (fu != 3L),
        sample_summary$N_SUBJECTS[sample_summary$FU == fu] == nrow(followup))
      for (analysis in c("analysis_change", "analysis_level")) {
        followup$analyte <- followup$clock
        if (analysis == "analysis_change") followup$analyte <- followup$analyte - followup$analyte_baseline
        for (stratum in c("all", "male", "female")) {
          md <- followup
          if (stratum != "all") md <- md[as.character(md$female) == if (stratum == "female") "1" else "0", ]
          adjustment <- c("analyte_baseline", "CR", if (stratum == "all") "female",
                           "PC1", "PC2", "PC3", "agebl", "deidsite", "bmistrat")
          fit <- lm(reformulate(c("analyte", adjustment), response = "OUTCOME"), data = md)
          expected <- coef(summary(fit))["analyte", c("Estimate", "Std. Error")]
          effects <- results[[analysis]][[stratum]]$outcome_effects
          coefficients <- results[[analysis]][[stratum]]$coefficients
          stopifnot(nrow(effects) == 8L, identical(sort(unique(effects$FU)), 1:4),
            is.numeric(effects$EFFECT_SIZE), all(is.finite(effects$EFFECT_SIZE)),
            all(effects$STATUS == "ok"), all(coefficients$N_OBS[coefficients$FU == fu] == nrow(md)))
          actual <- effects[effects$FU == fu & effects$ANALYTE_NAME == "Brain", c("EFFECT_SIZE", "SE")]
          stopifnot(nrow(actual) == 1L,
            isTRUE(all.equal(as.numeric(actual[1L, ]), unname(expected), tolerance = 1e-10)))
        }
      }
    }
    cat("PASS ", label, ": shuffled clock inputs, FU1/2/3/4, fixed outcome, independent baseline pairs, manual fits, reports and plots\n", sep = "")

    duplicate <- run_runner(runners[[label]], label, duplicate = TRUE)
    stopifnot(duplicate$status != 0L, is.null(duplicate$prepared),
      any(grepl("duplicat", duplicate$log, ignore.case = TRUE)),
      !file.exists(file.path(duplicate$out_dir, "results.rds")),
      !dir.exists(file.path(duplicate$out_dir, "checkpoints")))
    cat("PASS ", label, ": duplicate Organ/barcode rejected before fitting\n", sep = "")
  }
}

tryCatch(run_checks(), finally = unlink(fixture_dir, recursive = TRUE))
