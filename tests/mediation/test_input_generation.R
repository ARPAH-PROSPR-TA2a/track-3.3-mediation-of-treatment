repo_root <- normalizePath(getwd(), mustWork = TRUE)
default_track111 <- file.path(dirname(repo_root), "Track1.1.1")
track111_dir <- Sys.getenv("TRACK111_DIR", unset = default_track111)

if (!dir.exists(track111_dir)) {
  cat(
    "SKIP live Track 1.1.1 integration: checkout not found at ",
    track111_dir,
    ". Set TRACK111_DIR to enable it.\n",
    sep = ""
  )
  quit(status = 0L)
}

temp_output <- tempfile("mediation-inputs-")
dir.create(temp_output, recursive = TRUE)
on.exit(unlink(temp_output, recursive = TRUE), add = TRUE)

rscript <- file.path(R.home("bin"), "Rscript")
generator <- file.path(
  repo_root,
  "mediation",
  "Examples",
  "generate_example_inputs.R"
)
status <- system2(
  rscript,
  c(generator, temp_output),
  env = paste0("TRACK111_DIR=", normalizePath(track111_dir, mustWork = TRUE))
)
if (status != 0L) stop("Live upstream input generation failed.")

generated_treatment <- readRDS(file.path(temp_output, "proteomics_treatment_results.rds"))
generated_outcome <- readRDS(file.path(temp_output, "proteomics_outcome_results.rds"))
expected_dir <- file.path(
  repo_root,
  "mediation",
  "Examples",
  "ExampleInputs"
)
expected_treatment <- readRDS(file.path(expected_dir, "proteomics_treatment_results.rds"))
expected_outcome <- readRDS(file.path(expected_dir, "proteomics_outcome_results.rds"))

compare_effect_tables <- function(generated, expected, table_name, label) {
  for (analysis in c("analysis_change", "analysis_level")) {
    for (stratum in c("all", "male", "female")) {
      actual <- generated[[analysis]][[stratum]][[table_name]]
      reference <- expected[[analysis]][[stratum]][[table_name]]
      key_columns <- c("ANALYTE_NAME", "FU")
      actual <- actual[order(actual$ANALYTE_NAME, actual$FU), ]
      reference <- reference[order(reference$ANALYTE_NAME, reference$FU), ]
      if (!identical(actual[key_columns], reference[key_columns])) {
        stop(label, " keys changed for ", analysis, "/", stratum, ".")
      }
      for (column in c("EFFECT_SIZE", "SE", "P_VALUE", "BH_P_VALUE")) {
        if (!isTRUE(all.equal(actual[[column]], reference[[column]], tolerance = 1e-6))) {
          stop(label, "$", column, " changed for ", analysis, "/", stratum, ".")
        }
      }
    }
  }
}

compare_effect_tables(
  generated_treatment,
  expected_treatment,
  "treatment_effects",
  "Track 1.1.1 treatment effects"
)
cat("PASS Track 1.1.1 regenerated the expected alpha inputs\n")

compare_effect_tables(
  generated_outcome,
  expected_outcome,
  "outcome_effects",
  "OutcomeWAS outcome effects"
)
cat("PASS OutcomeWAS regenerated the expected beta inputs\n")

source(file.path(repo_root, "mediation", "main.R"))
mediation <- FAST_mediation(
  generated_treatment,
  generated_outcome,
  outcome_type = "continuous"
)

for (analysis in c("analysis_change", "analysis_level")) {
  for (stratum in c("all", "male", "female")) {
    diagnostics <- mediation[[analysis]][[stratum]]$join_diagnostics
    if (diagnostics$matched_rows != 120L ||
        diagnostics$treatment_only_rows != 0L ||
        diagnostics$outcome_only_rows != 0L) {
      stop("Incomplete upstream join for ", analysis, "/", stratum, ".")
    }
  }
}
cat("PASS All upstream ANALYTE_NAME x FU keys join completely\n")

effects <- mediation$analysis_change$all$mediation_effects
planted <- effects[effects$ANALYTE_NAME == "PROT_00001" & effects$FU == 1L, ]
if (nrow(planted) != 1L || planted$BH_P_VALUE >= 0.05) {
  stop("The deterministic PROT_00001/FU1 mediation check failed.")
}
cat("PASS Deterministic PROT_00001/FU1 mediation result\n")
cat("\nMediation live input-generation test passed.\n")
