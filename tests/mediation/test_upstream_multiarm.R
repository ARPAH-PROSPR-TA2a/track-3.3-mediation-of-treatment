# Optional integration with the actual local (including uncommitted) Track 1.1.1.
track111_dir <- Sys.getenv("TRACK111_DIR", unset = file.path(dirname(getwd()), "Track1.1.1"))
if (!file.exists(file.path(track111_dir, "main.R"))) {
  cat("SKIP live multi-arm integration: set TRACK111_DIR to a Track 1.1.1 checkout.\n")
  quit(status = 0L)
}
source(file.path("outcomewas", "main.R"))
source(file.path("tests", "outcomewas", "synthetic_multiarm.R"))
fixture <- .make_multiarm_outcome_fixture(n_analytes = 2L)
upstream <- new.env(parent = globalenv())
source(file.path(normalizePath(track111_dir), "main.R"), local = upstream)
treatment <- suppressWarnings(suppressMessages(upstream$FAST_omics_WAS(
  fixture$pheno, fixture$omics, omics_type = "other",
  additional_covariates = "AGE", n_cores = 1L,
  treatment_labels = fixture$treatment_labels
)))
outcome <- suppressWarnings(suppressMessages(FAST_outcome_WAS(
  fixture$pheno, fixture$omics, omics_type = "other",
  additional_covariates = "AGE", n_cores = 1L,
  treatment_labels = fixture$treatment_labels
)))
source(file.path("mediation", "main.R"))
results <- FAST_mediation(treatment, outcome)
for (analysis in c("analysis_change", "analysis_level")) {
  for (stratum in c("all", "male", "female")) {
    alpha <- treatment[[analysis]][[stratum]]$treatment_effects
    beta <- outcome[[analysis]][[stratum]]$outcome_effects
    indirect <- results[[analysis]][[stratum]]$mediation_effects
    alpha_key <- paste(alpha$ANALYTE_NAME, alpha$FU, alpha$CONTRAST)
    indirect_key <- paste(indirect$ANALYTE_NAME, indirect$FU, indirect$CONTRAST)
    ai <- match(indirect_key, alpha_key)
    bi <- match(paste(indirect$ANALYTE_NAME, indirect$FU),
                paste(beta$ANALYTE_NAME, beta$FU))
    stopifnot(nrow(indirect) == 8L, all(indirect$STATUS == "ok"),
              setequal(indirect$CONTRAST, c("2_vs_0", "7_vs_0")),
              all(indirect$TREATMENT == alpha$TREATMENT[ai]),
              isTRUE(all.equal(indirect$INDIRECT_EFFECT, alpha$EFFECT_SIZE[ai] * beta$EFFECT_SIZE[bi])),
              results[[analysis]][[stratum]]$join_diagnostics$outcome_only_rows == 0L,
              results[[analysis]][[stratum]]$join_diagnostics$matched_outcome_rows == 4L)
  }
}
cat("PASS local Track 1.1.1 joint mixed models feed common-slope OutcomeWAS and separate mediation contrasts across both FUs/analyses/all strata\n")
