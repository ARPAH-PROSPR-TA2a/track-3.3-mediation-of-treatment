test_scripts <- c(
  file.path("tests", "mediation", "test_input_generation.R"),
  file.path("tests", "mediation", "test_comprehensive.R"),
  file.path("tests", "outcomewas", "test_protein_annotation.R"),
  file.path("tests", "outcomewas", "test_comprehensive.R"),
  file.path("tests", "outcomewas", "test_dnam_betas.R"),
  file.path("tests", "outcomewas", "test_parallel_checkpoint.R"),
  file.path("tests", "outcomewas", "test_verbose_progress.R")
)

rscript <- file.path(R.home("bin"), "Rscript")
statuses <- vapply(test_scripts, function(script) {
  cat("\nRunning ", script, "\n", sep = "")
  system2(rscript, script)
}, integer(1))

if (any(statuses != 0L)) {
  failed <- test_scripts[statuses != 0L]
  stop("Test scripts failed: ", paste(failed, collapse = ", "))
}

cat("\nAll Track 3.3 test scripts passed.\n")
