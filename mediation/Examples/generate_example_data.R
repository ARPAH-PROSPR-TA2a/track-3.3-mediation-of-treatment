repo_root <- normalizePath(getwd(), mustWork = TRUE)
data_dir <- file.path(repo_root, "mediation", "Examples", "ExampleData")

pheno <- readRDS(file.path(data_dir, "pheno_base.rds"))
omics <- readRDS(file.path(data_dir, "proteomics_log2.rds"))

fu <- as.integer(as.character(pheno$FU))
baseline <- pheno[fu == 0L, , drop = FALSE]
followup_1 <- pheno[fu == 1L, , drop = FALSE]
subjects <- intersect(baseline$SUBJECT_ID, followup_1$SUBJECT_ID)
baseline <- baseline[match(subjects, baseline$SUBJECT_ID), , drop = FALSE]
followup_1 <- followup_1[match(subjects, followup_1$SUBJECT_ID), , drop = FALSE]

signal_analyte <- "PROT_00001"
signal_row <- match(signal_analyte, omics$ANALYTE_NAME)
if (is.na(signal_row)) stop("Missing planted signal analyte: ", signal_analyte)

signal_change <-
  as.numeric(omics[signal_row, followup_1$SAMPLE_ID]) -
  as.numeric(omics[signal_row, baseline$SAMPLE_ID])
treatment <- as.numeric(as.character(followup_1$TREATMENT_GROUP))
bmi_scaled <- as.numeric(scale(followup_1$mbmi))
signal_scaled <- as.numeric(scale(signal_change))

set.seed(20260518)
outcome <-
  0.9 * signal_scaled +
  0.4 * treatment +
  0.25 * bmi_scaled +
  rnorm(length(signal_scaled), sd = 0.7)

outcome_by_subject <- setNames(outcome, as.character(subjects))
pheno$OUTCOME <- unname(outcome_by_subject[as.character(pheno$SUBJECT_ID)])
if (anyNA(pheno$OUTCOME)) stop("Generated outcome is missing for one or more subjects.")

saveRDS(pheno, file.path(data_dir, "pheno_example.rds"))

cat("Mediation example phenotype generated\n")
cat("  Subjects:       ", length(subjects), "\n", sep = "")
cat("  Samples:        ", nrow(pheno), "\n", sep = "")
cat("  Signal analyte: ", signal_analyte, "\n", sep = "")
cat("  Seed:           20260518\n")
