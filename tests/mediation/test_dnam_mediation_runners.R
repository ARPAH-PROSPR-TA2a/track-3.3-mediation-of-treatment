source(file.path("mediation", "main.R"))


run_runner_tests <- function(outcome_name) {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
  runner_name <- paste0("run_DNAm_", outcome_name, "_24mo.R")
  runner <- file.path(
    repo_root, "mediation", "CALERIE", "DNAm", paste0(outcome_name, "_24mo"), runner_name
  )
  if (!file.exists(runner)) stop("DNAm ", outcome_name, " mediation runner not found: ", runner)

  temp_root <- tempfile(paste0("dnam-", tolower(outcome_name), "-mediation-"))
  dir.create(temp_root)
  on.exit(unlink(temp_root, recursive = TRUE), add = TRUE)
  analyses <- c("analysis_change", "analysis_level")
  strata <- c("all", "male", "female")
  result_name <- paste0("DNAm_betas_3.3Med_", outcome_name, "_results.rds")

  expect_true <- function(value, label) {
    if (!isTRUE(value)) stop("FAIL: ", label)
    cat("PASS ", label, "\n", sep = "")
  }

  # Check production defaults before fixture paths replace the editable values.
  runner_code <- parse(runner, keep.source = FALSE)
  assignment_index <- function(name) {
    index <- which(vapply(runner_code, function(expression) {
      is.call(expression) && identical(expression[[1]], as.name("<-")) &&
        identical(expression[[2]], as.name(name))
    }, logical(1)))
    if (length(index) != 1L) stop("Expected one editable assignment for ", name)
    index
  }
  default_paths <- list(
    repo = "~/FAST/GitHub/track-3.3",
    treatment_path = "~/FAST/Outputs/1.1.1/DNAm_Bvals_1.1.1/DNAm_betas_1.1.1.rds",
    outcome_path = paste0("~/FAST/Outputs/3.3/DNAm_betas_3.3OWAS_", outcome_name,
                          "/DNAm_betas_3.3OWAS_", outcome_name, "_results.rds"),
    out_dir = paste0("~/FAST/Outputs/3.3/DNAm_betas_3.3Med_", outcome_name)
  )
  cat("DNAm beta / 24-month ", outcome_name, " mediation runner tests\n\n", sep = "")
  for (name in names(default_paths)) {
    expect_true(
      identical(runner_code[[assignment_index(name)]][[3]], default_paths[[name]]),
      paste("Default", name, "matches the server layout for", outcome_name)
    )
  }

  alpha <- data.frame(
    ANALYTE_NAME = rep(sprintf("cg%08d", 1:4), times = 2L),
    FU = rep(1:2, each = 4L),
    EFFECT_SIZE = c(0.04, -0.015, 0.001, 0, 0.06, -0.02, 0.005, 0),
    SE = c(0.003, 0.003, 0.01, 0.01, 0.005, 0.004, 0.004, 0.01),
    stringsAsFactors = FALSE
  )
  beta <- alpha
  beta$EFFECT_SIZE <- c(2, 0.8, -0.05, 0, 1.5, 1.1, -0.1, 0)
  beta$SE <- c(0.2, 0.25, 0.2, 0.2, 0.15, 0.2, 0.1, 0.2)
  for (table_name in c("alpha", "beta")) {
    table <- get(table_name)
    table$P_VALUE <- 2 * pnorm(-abs(table$EFFECT_SIZE / table$SE))
    table$BH_P_VALUE <- .adjust_within_fu(table$P_VALUE, table$FU)
    table$BH_P_VALUE_FILTERED <- ifelse(
      table$ANALYTE_NAME %in% sprintf("cg%08d", c(1, 3)),
      table$BH_P_VALUE,
      NA_real_
    )
    assign(table_name, table)
  }

  make_results <- function(table, effect_name) {
    result <- setNames(vector("list", length(analyses)), analyses)
    for (analysis in analyses) {
      result[[analysis]] <- setNames(vector("list", length(strata)), strata)
      for (stratum in strata) {
        effects <- table
        multiplier <- match(analysis, analyses) / match(stratum, strata)
        effects$EFFECT_SIZE <- effects$EFFECT_SIZE * multiplier
        result[[analysis]][[stratum]] <- setNames(list(effects), effect_name)
        # Production inputs also contain large tables that mediation does not use.
        result[[analysis]][[stratum]]$coefficients <- data.frame(
          TERM = c("(Intercept)", "AGE"), ESTIMATE = c(1, 0.1)
        )
      }
    }
    result
  }

  treatment <- make_results(alpha, "treatment_effects")
  outcome <- make_results(beta[c(8, 2, 5, 3, 7, 1, 6, 4), ], "outcome_effects")

  launch <- function(label, treatment_input = treatment, outcome_input = outcome,
                     preexisting = FALSE) {
    case_dir <- file.path(temp_root, label)
    dir.create(case_dir)
    treatment_path <- file.path(case_dir, "treatment.rds")
    outcome_path <- file.path(case_dir, "outcome.rds")
    output_dir <- file.path(case_dir, "output")
    saveRDS(treatment_input, treatment_path)
    saveRDS(outcome_input, outcome_path)
    input_paths <- c(treatment_path, outcome_path)
    before <- tools::md5sum(input_paths)

    existing_path <- file.path(output_dir, result_name)
    if (preexisting) {
      dir.create(output_dir)
      saveRDS(list(existing_result = "preserve me"), existing_path)
      existing_hash <- tools::md5sum(existing_path)
    }

    # Edit only the four user-configurable path assignments in a temporary copy.
    fixture_paths <- list(
      repo = repo_root, treatment_path = treatment_path,
      outcome_path = outcome_path, out_dir = output_dir
    )
    code <- runner_code
    for (name in names(fixture_paths)) {
      code[[assignment_index(name)]][[3]] <- fixture_paths[[name]]
    }
    fixture_runner <- file.path(case_dir, runner_name)
    writeLines(unlist(lapply(code, deparse)), fixture_runner)
    workspace_path <- file.path(case_dir, "workspace.rds")
    expression <- bquote({
      setwd(.(case_dir))
      source(.(runner_name))
      stopifnot(identical(getwd(), .(repo_root)))
      saveRDS(mget(c("treatment_results", "outcome_results", "mediation_results",
                     "mediation_summary"), envir = .GlobalEnv), .(workspace_path))
    })
    output <- suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      args = c("-e", shQuote(paste(deparse(expression), collapse = "\n"))),
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    if (!identical(before, tools::md5sum(input_paths))) {
      stop("FAIL: runner modified its input files for ", label)
    }
    if (preexisting && !identical(existing_hash, tools::md5sum(existing_path))) {
      stop("FAIL: runner replaced the pre-existing result")
    }
    list(status = status, output = output, directory = output_dir,
         workspace = workspace_path)
  }

  successful <- launch("success")
  if (successful$status != 0L) {
    stop("FAIL: valid runner invocation\n", paste(successful$output, collapse = "\n"))
  }
  artifacts <- c(result_name, "summary.tsv", "provenance.rds", "run.log")
  expect_true(
    all(file.exists(file.path(successful$directory, artifacts))),
    "Runner saves results, summary, provenance, and run.log"
  )
  log <- readLines(file.path(successful$directory, "run.log"))
  other_outcome <- if (outcome_name == "INF") "MetS" else "INF"
  expect_true(
    any(grepl(paste0("START DNAm beta / ", outcome_name, " mediation"), log, fixed = TRUE)) &&
      any(grepl(paste0("fixed 24mo ", outcome_name, " outcome"), log, fixed = TRUE)) &&
      !any(grepl(paste0("START DNAm beta / ", other_outcome, " mediation"), log, fixed = TRUE)) &&
      !any(grepl(paste0("fixed 24mo ", other_outcome, " outcome"), log, fixed = TRUE)),
    paste("Run log correctly identifies", outcome_name, "as the 24-month outcome")
  )
  expect_true(
    all(vapply(file.path(dirname(successful$directory), c("treatment.rds", "outcome.rds")),
               function(path) any(grepl(path, log, fixed = TRUE)), logical(1))),
    "Run log records both input paths"
  )
  expect_true(
    any(grepl("Calculating indirect effects", log, fixed = TRUE)) &&
      any(grepl(file.path(successful$directory, result_name), log, fixed = TRUE)) &&
      any(grepl("DONE", log, fixed = TRUE)),
    "Run log records calculation, saved result path, and completion"
  )
  expected_preflight <- expand.grid(role = c("treatment", "outcome"), analysis = analyses,
                                    stratum = strata, fu = 1:2, stringsAsFactors = FALSE)
  expected_preflight <- with(expected_preflight, sprintf(
    "preflight %s/%s/%s/FU%d: 4 rows; 2 filtered", role, analysis, stratum, fu
  ))
  expect_true(
    all(vapply(expected_preflight, function(line) any(grepl(line, log, fixed = TRUE)),
               logical(1))),
    "Run log records counts for every treatment/outcome analysis, stratum, and FU"
  )
  expect_true(
    sum(grepl("8 matched; 0 treatment-only; 0 outcome-only", log, fixed = TRUE)) == 6L &&
      any(grepl("BH_LT_0_05", log, fixed = TRUE)),
    "Run log records all complete joins and the significance summary"
  )

  actual <- readRDS(file.path(successful$directory, result_name))
  expected <- FAST_mediation(treatment, outcome, outcome_type = "continuous")
  for (analysis in analyses) {
    for (stratum in strata) {
      expect_true(
        isTRUE(all.equal(actual[[analysis]][[stratum]], expected[[analysis]][[stratum]],
                         tolerance = 1e-12)),
        paste("Exact key alignment matches FAST_mediation for", analysis, stratum)
      )
    }
  }
  expect_true(
    identical(attr(actual, "outcome_type"), "continuous") &&
      identical(attr(actual, "inference"), "sobel"),
    "Runner preserves continuous Sobel metadata"
  )

  workspace <- readRDS(successful$workspace)
  expect_true(
    identical(names(workspace), c("treatment_results", "outcome_results",
                                  "mediation_results", "mediation_summary")) &&
      isTRUE(all.equal(workspace$mediation_results, expected, tolerance = 1e-12)),
    "Sourcing leaves the inputs, mediation results, and summary available for inspection"
  )

  summary <- read.delim(file.path(successful$directory, "summary.tsv"),
                        stringsAsFactors = FALSE, check.names = FALSE)
  columns <- c("ANALYSIS", "STRATUM", "FU", "CONTRAST", "TREATMENT", "REFERENCE",
               "METHYLATION_MONTH", "OUTCOME_MONTH", "N_REQUESTED",
               "N_TESTED", "N_UNAVAILABLE", "N_FILTERED", "P_LT_0_05", "BH_LT_0_05",
               "FILTERED_BH_LT_0_05")
  expect_true(identical(names(summary), columns) && nrow(summary) == 12L,
              "Summary has the specified columns and all 12 analysis/stratum/FU rows")
  expect_true(isTRUE(all.equal(summary, workspace$mediation_summary,
                               check.attributes = FALSE)),
              "Saved summary matches the inspectable mediation_summary object")
  for (analysis in analyses) {
    for (stratum in strata) {
      effects <- expected[[analysis]][[stratum]]$mediation_effects
      for (fu in 1:2) {
        observed <- summary[summary$ANALYSIS == analysis &
                              summary$STRATUM == stratum & summary$FU == fu, ]
        if (nrow(observed) != 1L) stop("FAIL: missing or duplicate summary key")
        subset <- effects[effects$FU == fu, ]
        values <- c(
          METHYLATION_MONTH = fu * 12L,
          OUTCOME_MONTH = 24L,
          N_TESTED = nrow(subset),
          N_FILTERED = sum(!is.na(subset$BH_P_VALUE_FILTERED)),
          P_LT_0_05 = sum(subset$P_VALUE < 0.05),
          BH_LT_0_05 = sum(subset$BH_P_VALUE < 0.05),
          FILTERED_BH_LT_0_05 = sum(subset$BH_P_VALUE_FILTERED < 0.05, na.rm = TRUE)
        )
        if (!isTRUE(all.equal(as.numeric(observed[names(values)]),
                              as.numeric(values)))) {
          stop("FAIL: summary timing or counts differ for ", analysis, "/", stratum, "/", fu)
        }
      }
    }
  }
  cat("PASS Summary timing and significance counts match saved mediation effects\n")

  # Actual runner integration with nonconsecutive arms and unavailable filtered
  # probes, reusing one common outcome slope for both treatment contrasts.
  labels <- c(`0` = "Control", `2` = "Diet", `7` = "Exercise")
  multi_alpha <- do.call(rbind, lapply(c("2", "7"), function(arm) {
    x <- alpha
    x$CONTRAST <- paste0(arm, "_vs_0")
    x$TREATMENT <- unname(labels[arm])
    x$REFERENCE <- "Control"
    x$STATUS <- "ok"
    x$IN_FILTERED_SET <- x$ANALYTE_NAME %in% sprintf("cg%08d", c(1, 3))
    x
  }))
  multi_beta <- beta
  multi_beta$STATUS <- "ok"
  multi_beta$IN_FILTERED_SET <- multi_beta$ANALYTE_NAME %in% sprintf("cg%08d", c(1, 3))
  missing_alpha <- with(multi_alpha, FU == 1 & CONTRAST == "2_vs_0" & ANALYTE_NAME == "cg00000001")
  multi_alpha[missing_alpha, c("EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE", "BH_P_VALUE_FILTERED")] <- NA_real_
  multi_alpha$STATUS[missing_alpha] <- "fit_error"
  missing_beta <- with(multi_beta, FU == 2 & ANALYTE_NAME == "cg00000003")
  multi_beta[missing_beta, c("EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE", "BH_P_VALUE_FILTERED")] <- NA_real_
  multi_beta$STATUS[missing_beta] <- "not_estimable"
  multi_treatment <- make_results(multi_alpha, "treatment_effects")
  multi_outcome <- make_results(multi_beta, "outcome_effects")
  multi_treatment$metadata <- list(treatment_labels = labels)
  multi_outcome$metadata <- list(analyte_slope = "common_across_treatment_arms", treatment_labels = labels)
  multi <- launch("multiarm-unavailable", multi_treatment, multi_outcome)
  if (multi$status != 0L) stop("FAIL: multi-arm runner\n", paste(multi$output, collapse = "\n"))
  multi_summary <- read.delim(file.path(multi$directory, "summary.tsv"))
  expect_true(nrow(multi_summary) == 24L &&
                setequal(multi_summary$CONTRAST, c("2_vs_0", "7_vs_0")) &&
                all(multi_summary$N_REQUESTED == 4L) &&
                all(multi_summary$N_FILTERED == 2L) &&
                sum(multi_summary$N_UNAVAILABLE) == 18L &&
                all(multi_summary$N_TESTED + multi_summary$N_UNAVAILABLE == multi_summary$N_REQUESTED),
              "Multi-arm runner separates contrasts and retains unavailable filtered members")
  multi_results <- readRDS(file.path(multi$directory, result_name))
  multi_provenance <- readRDS(file.path(multi$directory, "provenance.rds"))
  expect_true(identical(multi_provenance$outcome_metadata, multi_outcome$metadata) &&
                identical(multi_provenance$mediation_metadata$metadata,
                          attr(multi_results, "metadata")),
              "Runner provenance records the common analyte slope and mediation estimand")

  # A contrast may be requested at only one follow-up. Summaries must describe
  # observed combinations rather than filling the FU-by-contrast grid.
  sparse_treatment <- multi_treatment
  for (analysis in analyses) {
    for (stratum in strata) {
      effects <- sparse_treatment[[analysis]][[stratum]]$treatment_effects
      sparse_treatment[[analysis]][[stratum]]$treatment_effects <-
        effects[!(effects$FU == 2L & effects$CONTRAST == "7_vs_0"), ]
    }
  }
  sparse <- launch("multiarm-sparse-followup", sparse_treatment, multi_outcome)
  if (sparse$status != 0L) stop("FAIL: sparse multi-arm runner\n", paste(sparse$output, collapse = "\n"))
  sparse_summary <- read.delim(file.path(sparse$directory, "summary.tsv"))
  expected_sparse <- multi_summary[!(multi_summary$FU == 2L & multi_summary$CONTRAST == "7_vs_0"), ]
  expect_true(nrow(sparse_summary) == 18L &&
                !any(sparse_summary$FU == 2L & sparse_summary$CONTRAST == "7_vs_0") &&
                all(sparse_summary$N_REQUESTED > 0L),
              "Sparse follow-up summaries contain only observed contrast combinations")
  expect_true(!anyNA(sparse_summary[c("TREATMENT", "REFERENCE")]) &&
                isTRUE(all.equal(sparse_summary, expected_sparse, check.attributes = FALSE)),
              "Sparse summaries preserve follow-up order, labels, counts and significance values")

  expect_rejection <- function(label, treatment_input = treatment, outcome_input = outcome,
                              preexisting = FALSE) {
    result <- launch(label, treatment_input, outcome_input, preexisting)
    expect_true(result$status != 0L, paste("Rejects", label))
    patterns <- c(
      "missing-FU2" = "FU1.*FU2",
      "missing-stratum" = "female.*must contain",
      "missing-analysis" = "analysis_level.*must contain",
      "partial-keys" = "Partial mediation join",
      "filtered-membership-mismatch" = "Filtered-probe membership disagrees",
      "missing-filtered-column" = "must contain.*BH_P_VALUE_FILTERED",
      "nonfinite-SE" = "SE.*finite numeric",
      "nonpositive-SE" = "SE.*strictly positive",
      "nonfinite-effect" = "EFFECT_SIZE.*finite numeric",
      "fractional-FU" = "FU1.*FU2",
      "duplicate-keys" = "duplicate ANALYTE_NAME x FU",
      "overflow-in-derived-estimates" = "Non-finite mediation output",
      "empty-filtered-FU" = "[Ff]iltered.*FU|FU.*[Ff]iltered",
      "pre-existing-result" = "Output already exists"
    )
    if (!any(grepl(patterns[[label]], result$output))) {
      stop("FAIL: unexpected rejection for ", label, "\n",
           paste(result$output, collapse = "\n"))
    }
    if (!preexisting) {
      expect_true(!file.exists(file.path(result$directory, result_name)),
                  paste("Does not publish a result for", label))
    }
  }

  bad <- treatment
  table <- bad$analysis_change$all$treatment_effects
  bad$analysis_change$all$treatment_effects <- table[table$FU == 1L, ]
  expect_rejection("missing-FU2", treatment_input = bad)

  bad <- outcome
  bad$analysis_level$female <- NULL
  expect_rejection("missing-stratum", outcome_input = bad)

  bad <- outcome
  bad$analysis_level <- NULL
  expect_rejection("missing-analysis", outcome_input = bad)

  bad <- outcome
  bad$analysis_change$all$outcome_effects <- bad$analysis_change$all$outcome_effects[-1, ]
  expect_rejection("partial-keys", outcome_input = bad)

  bad <- outcome
  table <- bad$analysis_change$all$outcome_effects
  index <- which(!is.na(table$BH_P_VALUE_FILTERED))[1]
  bad$analysis_change$all$outcome_effects$BH_P_VALUE_FILTERED[index] <- NA_real_
  expect_rejection("filtered-membership-mismatch", outcome_input = bad)

  bad <- treatment
  bad$analysis_change$all$treatment_effects$BH_P_VALUE_FILTERED <- NULL
  expect_rejection("missing-filtered-column", treatment_input = bad)

  bad <- treatment
  bad$analysis_change$all$treatment_effects$SE[1] <- Inf
  expect_rejection("nonfinite-SE", treatment_input = bad)

  bad <- outcome
  bad$analysis_change$all$outcome_effects$SE[1] <- 0
  expect_rejection("nonpositive-SE", outcome_input = bad)

  bad <- outcome
  bad$analysis_level$female$outcome_effects$EFFECT_SIZE[1] <- NA_real_
  expect_rejection("nonfinite-effect", outcome_input = bad)

  bad <- treatment
  bad$analysis_change$all$treatment_effects$FU[1] <- 1.5
  expect_rejection("fractional-FU", treatment_input = bad)

  bad <- treatment
  table <- bad$analysis_change$all$treatment_effects
  bad$analysis_change$all$treatment_effects <- rbind(table, table[1, ])
  expect_rejection("duplicate-keys", treatment_input = bad)

  bad <- treatment
  bad$analysis_change$all$treatment_effects$EFFECT_SIZE[1] <- .Machine$double.xmax
  expect_rejection("overflow-in-derived-estimates", treatment_input = bad)

  bad_treatment <- treatment
  bad_outcome <- outcome
  alpha_fu2 <- bad_treatment$analysis_change$all$treatment_effects$FU == 2L
  beta_fu2 <- bad_outcome$analysis_change$all$outcome_effects$FU == 2L
  bad_treatment$analysis_change$all$treatment_effects$BH_P_VALUE_FILTERED[alpha_fu2] <- NA_real_
  bad_outcome$analysis_change$all$outcome_effects$BH_P_VALUE_FILTERED[beta_fu2] <- NA_real_
  expect_rejection("empty-filtered-FU", treatment_input = bad_treatment,
                   outcome_input = bad_outcome)

  expect_rejection("pre-existing-result", preexisting = TRUE)
  cat("PASS Existing results and every input RDS remain byte-identical\n")
  cat("\nAll DNAm ", outcome_name, " mediation runner tests passed.\n", sep = "")
}

for (outcome_name in c("INF", "MetS")) run_runner_tests(outcome_name)
