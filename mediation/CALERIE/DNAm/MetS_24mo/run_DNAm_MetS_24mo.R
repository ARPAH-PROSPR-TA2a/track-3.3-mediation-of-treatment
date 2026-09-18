# -----------------------------
# Paths: edit these for your machine
# -----------------------------
repo <- "~/FAST/GitHub/track-3.3"
treatment_path <- "~/FAST/Outputs/1.1.1/DNAm_Bvals_1.1.1/DNAm_betas_1.1.1.rds"
outcome_path <- "~/FAST/Outputs/3.3/DNAm_betas_3.3OWAS_MetS/DNAm_betas_3.3OWAS_MetS_results.rds"
out_dir <- "~/FAST/Outputs/3.3/DNAm_betas_3.3Med_MetS"

setwd(repo)
source(file.path(repo, "mediation", "main.R"))

# Mediation uses saved effect summaries and runs in one R process.
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
)
started <- Sys.time()
analyses <- c("analysis_change", "analysis_level")
strata <- c("all", "male", "female")
result_path <- file.path(out_dir, "DNAm_betas_3.3Med_MetS_results.rds")
if (any(file.exists(c(result_path, file.path(out_dir, c("summary.tsv", "provenance.rds")))))) {
  stop("Output already exists; choose a new out_dir.")
}
stopifnot(file.exists(treatment_path), file.exists(outcome_path))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "run.log")
log_status <- function(text) {
  line <- paste0("[3.3 mediation] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", text)
  message(line)
  cat(line, "\n", file = log_file, append = TRUE, sep = "")
}
log_status("START DNAm beta / MetS mediation; one R process; BLAS/OpenMP threads=1")
log_status(paste("Pipeline checkout:", normalizePath(repo)))
log_status("FU1=12mo and FU2=24mo DNAm; fixed 24mo MetS outcome")

# -----------------------------
# Read inputs
# -----------------------------
load_started <- Sys.time()
log_status(sprintf("Loading treatment input: %s (%.2f GiB on disk)",
                   treatment_path, file.info(treatment_path)$size / 1024^3))
treatment_results <- readRDS(treatment_path)
log_status(sprintf("Treatment loaded in %.1f seconds; full object %.2f GiB",
                   as.numeric(difftime(Sys.time(), load_started, units = "secs")),
                   as.numeric(object.size(treatment_results)) / 1024^3))
# Coefficient tables are large and unused by mediation. Release them before
# loading the next file; keep the treatment_effects tables for inspection.
for (analysis in analyses) {
  for (stratum in strata) {
    treatment_results[[analysis]][[stratum]]$coefficients <- NULL
  }
}
invisible(gc())
log_status(sprintf("Treatment coefficient tables released; retained object %.2f GiB",
                   as.numeric(object.size(treatment_results)) / 1024^3))

load_started <- Sys.time()
log_status(sprintf("Loading outcome input: %s (%.2f GiB on disk)",
                   outcome_path, file.info(outcome_path)$size / 1024^3))
outcome_results <- readRDS(outcome_path)
log_status(sprintf("Outcome loaded in %.1f seconds; full object %.2f GiB",
                   as.numeric(difftime(Sys.time(), load_started, units = "secs")),
                   as.numeric(object.size(outcome_results)) / 1024^3))
for (analysis in analyses) {
  for (stratum in strata) {
    outcome_results[[analysis]][[stratum]]$coefficients <- NULL
  }
}
invisible(gc())
log_status(sprintf("Outcome coefficient tables released; retained object %.2f GiB",
                   as.numeric(object.size(outcome_results)) / 1024^3))

# -----------------------------
# Check the intended analysis
# -----------------------------
# FU1 = 12-month DNAm; FU2 = 24-month DNAm; MetS is measured at 24 months.
# Both inputs must use beta units and compatible coding/covariates. MetS requires
# an observed outcome, so its cohort may differ from the treatment analysis.
log_status("Checking both analyses, all three strata, FU1/FU2 and filtered probe counts")
for (analysis in analyses) {
  for (stratum in strata) {
    effect_tables <- list(treatment = treatment_results[[analysis]][[stratum]]$treatment_effects,
                          outcome = outcome_results[[analysis]][[stratum]]$outcome_effects)
    for (role in names(effect_tables)) {
      effects <- effect_tables[[role]]
      if (!is.data.frame(effects) ||
          !all(c("FU", "BH_P_VALUE_FILTERED") %in% names(effects))) {
        stop(analysis, "/", stratum, " must contain both effect tables with FU and BH_P_VALUE_FILTERED.")
      }
      fu <- suppressWarnings(as.numeric(as.character(effects$FU)))
      if (anyNA(fu) || !setequal(fu, c(1, 2))) {
        stop(analysis, "/", stratum, " must contain exactly FU1=12mo and FU2=24mo.")
      }
      for (visit in 1:2) {
        if (!any(fu == visit & !is.na(effects$BH_P_VALUE_FILTERED))) {
          stop("No filtered probes for ", analysis, "/", stratum, "/FU", visit, ".")
        }
        log_status(sprintf("preflight %s/%s/%s/FU%d: %d rows; %d filtered",
                           role, analysis, stratum, visit, sum(fu == visit),
                           sum(fu == visit & !is.na(effects$BH_P_VALUE_FILTERED))))
      }
    }
  }
}
rm(effects, effect_tables)

# -----------------------------
# Mediation: alpha x beta, Sobel inference, BH within each FU/stratum/analysis
# -----------------------------
mediation_started <- Sys.time()
log_status("Calculating indirect effects: product of coefficients, Sobel inference, within-FU BH")
# The core validates keys, estimates, SEs and filtered membership. Treat its
# partial-join/stratum warnings as errors so omitted results get reviewed.
mediation_results <- withCallingHandlers(
  FAST_mediation(
    treatment_results = treatment_results,
    outcome_results = outcome_results,
    outcome_type = "continuous"
  ),
  warning = function(w) stop(conditionMessage(w), call. = FALSE),
  error = function(e) log_status(paste("FAILED:", conditionMessage(e)))
)
log_status(sprintf("Mediation complete in %.1f seconds; result object %.2f GiB",
                   as.numeric(difftime(Sys.time(), mediation_started, units = "secs")),
                   as.numeric(object.size(mediation_results)) / 1024^3))

# -----------------------------
# Summarize and save
# -----------------------------
summary_rows <- list()
for (analysis in analyses) {
  for (stratum in strata) {
    effects <- mediation_results[[analysis]][[stratum]]$mediation_effects
    join <- mediation_results[[analysis]][[stratum]]$join_diagnostics
    log_status(sprintf("join %s/%s: %d matched; %d treatment-only; %d outcome-only",
                       analysis, stratum, join$matched_rows,
                       join$treatment_only_rows, join$outcome_only_rows))
    numeric_columns <- c("INDIRECT_EFFECT", "INDIRECT_SE", "Z_VALUE", "P_VALUE",
                         "CI_LOWER", "CI_UPPER", "BH_P_VALUE")
    if (any(vapply(effects[numeric_columns], function(x) any(!is.finite(x)), logical(1)))) {
      stop("Non-finite mediation output for ", analysis, "/", stratum, ".")
    }
    for (fu in 1:2) {
      x <- effects[effects$FU == fu, , drop = FALSE]
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        ANALYSIS = analysis, STRATUM = stratum, FU = fu,
        METHYLATION_MONTH = 12L * fu, OUTCOME_MONTH = 24L,
        N_TESTED = nrow(x), N_FILTERED = sum(!is.na(x$BH_P_VALUE_FILTERED)),
        P_LT_0_05 = sum(x$P_VALUE < 0.05),
        BH_LT_0_05 = sum(x$BH_P_VALUE < 0.05),
        FILTERED_BH_LT_0_05 = sum(x$BH_P_VALUE_FILTERED < 0.05, na.rm = TRUE)
      )
    }
  }
}
mediation_summary <- do.call(rbind, summary_rows)
rm(effects, x, summary_rows)

save_started <- Sys.time()
log_status(paste("Saving mediation results:", result_path))
saveRDS(mediation_results, result_path)
log_status(sprintf("Results RDS saved in %.1f seconds (%.2f GiB on disk)",
                   as.numeric(difftime(Sys.time(), save_started, units = "secs")),
                   file.info(result_path)$size / 1024^3))
write.table(mediation_summary, file.path(out_dir, "summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(
  list(treatment_path = normalizePath(treatment_path),
       outcome_path = normalizePath(outcome_path), repo = normalizePath(repo),
       started = started, finished = Sys.time(),
       methylation_months = c(`1` = 12L, `2` = 24L), outcome_month = 24L,
       method = "product_of_coefficients", inference = "sobel",
       session_info = sessionInfo()),
  file.path(out_dir, "provenance.rds")
)
log_status(paste(capture.output(print(mediation_summary, row.names = FALSE)), collapse = "\n"))
log_status(paste("Summary saved to", file.path(out_dir, "summary.tsv")))
log_status(paste("Provenance saved to", file.path(out_dir, "provenance.rds")))
log_status(sprintf("DONE: results saved to %s; total %.1f seconds", result_path,
                   as.numeric(difftime(Sys.time(), started, units = "secs"))))
