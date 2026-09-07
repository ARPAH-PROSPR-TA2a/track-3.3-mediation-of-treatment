source(file.path("outcomewas", "main.R"))

# Keep the tracked examples unchanged: construct a small beta-scale fixture in memory.
pheno <- readRDS("outcomewas/Examples/ExampleData/pheno_example.rds")
for (column in names(pheno)) {
  if (inherits(pheno[[column]], "haven_labelled")) {
    pheno[[column]] <- as.vector(pheno[[column]])
  }
}
pheno <- pheno[, c("SAMPLE_ID", "SUBJECT_ID", "FU", "FEMALE", "TREATMENT_GROUP")]
pheno$SUBJECT_ID <- as.character(pheno$SUBJECT_ID)
for (column in c("FU", "FEMALE", "TREATMENT_GROUP")) {
  pheno[[column]] <- factor(as.integer(as.character(pheno[[column]])))
}
pheno <- pheno[complete.cases(pheno) & pheno$FU %in% c(0, 1), ]
pheno <- pheno[!duplicated(paste(pheno$SUBJECT_ID, pheno$FU)), ]

omics <- as.data.frame(readRDS("outcomewas/Examples/ExampleData/dnam_mvalues.rds"),
                       check.names = FALSE)
if (!"ANALYTE_NAME" %in% names(omics)) {
  omics <- cbind(ANALYTE_NAME = rownames(omics), omics, stringsAsFactors = FALSE)
}
pheno <- pheno[pheno$SAMPLE_ID %in% names(omics), ]
paired_subjects <- intersect(pheno$SUBJECT_ID[pheno$FU == 0],
                             pheno$SUBJECT_ID[pheno$FU == 1])
pheno <- pheno[pheno$SUBJECT_ID %in% paired_subjects, ]
samples <- as.character(pheno$SAMPLE_ID)
valid_probes <- readRDS("outcomewas/Data/FAST_epicv1_epicv2_probe_list.rds")
keep <- omics$ANALYTE_NAME %in% valid_probes &
  rowSums(!is.finite(as.matrix(omics[, samples, drop = FALSE]))) == 0L
omics <- omics[head(which(keep), 4L), c("ANALYTE_NAME", samples), drop = FALSE]
stopifnot(nrow(omics) == 4L)
omics[, samples] <- plogis(as.matrix(omics[, samples, drop = FALSE]) * log(2))
stopifnot(all(as.matrix(omics[, samples]) >= 0), all(as.matrix(omics[, samples]) <= 1))

baseline <- pheno[pheno$FU == 0, ]
followup <- pheno[pheno$FU == 1, ]
baseline <- baseline[match(followup$SUBJECT_ID, baseline$SUBJECT_ID), ]
beta_baseline <- as.numeric(omics[1, as.character(baseline$SAMPLE_ID)])
beta_followup <- as.numeric(omics[1, as.character(followup$SAMPLE_ID)])
set.seed(20260907)
outcome <- 3 * (beta_followup - beta_baseline) + beta_baseline +
  0.3 * as.numeric(as.character(followup$TREATMENT_GROUP)) +
  rnorm(nrow(followup), sd = 0.2)
pheno$OUTCOME <- outcome[match(pheno$SUBJECT_ID, followup$SUBJECT_ID)]

serial <- FAST_outcome_WAS(pheno, omics, omics_type = "DNAm", n_cores = 1L)
reports <- FAST_outcome_WAS_reports(pheno, omics, omics_type = "DNAm")

# Independent fits use the original beta fixture, not pipeline preparation helpers.
manual_data <- followup
manual_data$OUTCOME <- outcome
manual_data$analyte_baseline <- beta_baseline
for (response in c("change", "level")) {
  manual_data$analyte <- if (response == "change") {
    beta_followup - beta_baseline
  } else {
    beta_followup
  }
  fit <- lm(OUTCOME ~ analyte + analyte_baseline + TREATMENT_GROUP + FEMALE,
            data = manual_data)
  expected <- coef(summary(fit))["analyte", c("Estimate", "Std. Error")]
  effects <- serial[[paste0("analysis_", response)]]$all$outcome_effects
  actual <- effects[effects$ANALYTE_NAME == omics$ANALYTE_NAME[1] & effects$FU == 1,
                    c("EFFECT_SIZE", "SE")]
  stopifnot(nrow(actual) == 1L,
            isTRUE(all.equal(as.numeric(actual[1, ]), unname(expected), tolerance = 1e-10)))
  cat("PASS beta-scale ", response, " effect and SE match manual lm\n", sep = "")
}

cell <- pheno$FU == 0 & pheno$TREATMENT_GROUP == 0
reported <- reports$variable_summaries$all$omics_FU0_Tx0
expected_means <- rowMeans(as.matrix(omics[, as.character(pheno$SAMPLE_ID[cell]), drop = FALSE]))
stopifnot(isTRUE(all.equal(reported$MEAN[match(omics$ANALYTE_NAME, reported$ANALYTE_NAME)],
                           unname(expected_means), tolerance = 1e-12)))
for (stratum in reports$variable_summaries) {
  for (summary in stratum[startsWith(names(stratum), "omics_")]) {
    stopifnot(nrow(summary) == 4L, all(summary$MEAN >= 0 & summary$MEAN <= 1))
  }
}
cat("PASS reports retain beta-scale summaries\n")

checkpoint_dir <- tempfile("outcomewas-beta-checkpoint-")
tryCatch({
  parallel_result <- FAST_outcome_WAS(pheno, omics, omics_type = "DNAm", n_cores = 2L,
                                      checkpoint_dir = checkpoint_dir, checkpoint_batch_size = 2L)
  resumed <- FAST_outcome_WAS(pheno, omics, omics_type = "DNAm", n_cores = 1L,
                              checkpoint_dir = checkpoint_dir, checkpoint_batch_size = 2L)
  for (response in c("analysis_change", "analysis_level")) {
    for (stratum in c("all", "male", "female")) {
      for (slot in c("coefficients", "outcome_effects")) {
        expected <- serial[[response]][[stratum]][[slot]]
        stopifnot(isTRUE(all.equal(expected, parallel_result[[response]][[stratum]][[slot]],
                                  tolerance = 1e-10)),
                  isTRUE(all.equal(expected, resumed[[response]][[stratum]][[slot]],
                                  tolerance = 1e-10)))
      }
    }
  }
}, finally = unlink(checkpoint_dir, recursive = TRUE))
cat("PASS beta-scale serial, two-worker, and checkpoint-resumed results agree\n")
